/*

COPYRIGHT (c)2016 THE OHIO STATE UNIVERSITY
ALL RIGHTS RESERVED

PERMISSION IS GRANTED TO USE, COPY, CREATE DERIVATIVE WORKS AND
REDISTRIBUTE THIS SOFTWARE AND SUCH DERIVATIVE WORKS FOR NONCOMMERCIAL
EDUCATION AND RESEARCH PURPOSES, SO LONG AS NO FEE IS CHARGED, AND SO
LONG AS THE COPYRIGHT NOTICE ABOVE, THIS GRANT OF PERMISSION, AND THE
DISCLAIMER BELOW APPEAR IN ALL COPIES MADE; AND SO LONG AS THE NAME OF
THE OHIO STATE UNIVERSITY IS NOT USED IN ANY ADVERTISING OR PUBLICITY
PERTAINING TO THE USE OR DISTRIBUTION OF THIS SOFTWARE WITHOUT
SPECIFIC, WRITTEN PRIOR AUTHORIZATION.

THIS SOFTWARE IS PROVIDED AS IS, WITHOUT REPRESENTATION FROM THE OHIO
STATE UNIVERSITY AS TO ITS FITNESS FOR ANY PURPOSE, AND WITHOUT
WARRANTY BY THE OHIO STATE UNIVERSITY OF ANY KIND, EITHER EXPRESS OR
IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE OHIO STATE
UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES, INCLUDING SPECIAL,
INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, WITH RESPECT TO ANY
CLAIM ARISING OUT OF OR IN CONNECTION WITH THE USE OF THE SOFTWARE,
EVEN IF IT HAS BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.
*/

/* Homopolymer repeat calling */

/* implementation of site description */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "msdesc.h"
#include "mssequence.h"

/* helper functions for XML interpretation */
static int ParsePatterns(PMSDESCRIPTION pDescription, xmlNode *pNode)
{
  int i;                   /* loop variable                        */
  char repeatedNucleotide; /* nucleotide in homopolymer repeat     */
  char *pattern;           /* all parts of the pattern together    */
  xmlChar *text;           /* text of a node                       */
  xmlChar *flank5Text;     /* text with 5' flanking sequence       */
  xmlChar *flank3Text;     /* text with 3' flanking sequence       */
  xmlNode *pMyNode;        /* current node                         */

  /* determine total length */
  repeatedNucleotide=' ';
  flank5Text=NULL;
  flank3Text=NULL;
  pDescription->Length=0;
  for(pMyNode=pNode->children;pMyNode;pMyNode=pMyNode->next) {
    if (pNode->type==XML_ELEMENT_NODE &&
	(!strcmp((const char*)pMyNode->name,"FivePrimeFlanking") ||
	 !strcmp((const char*)pMyNode->name,"HomopolymerRun") ||
	 !strcmp((const char*)pMyNode->name,"ThreePrimeFlanking"))) {
      text=xmlNodeGetContent(pMyNode);
      if (!text) {
	fprintf(stderr,"Problem reading site pattern!\n");
	pDescription->Length=0;
        if (flank5Text) xmlFree(flank5Text);
        if (flank3Text) xmlFree(flank3Text);
	return(1);
      }
      if (!strcmp((const char*)pMyNode->name,"HomopolymerRun")) {
	for(i=1;i<strlen((const char*)text);i++) {
	  if (text[0]!=text[i]) {
	    fprintf(stderr,"All nucleotides in repeat must be equal!\n");
	    pDescription->Length=0;
	    xmlFree(text);
	    if (flank5Text) xmlFree(flank5Text);
	    if (flank3Text) xmlFree(flank3Text);
	    return(1);
	  }
        }
        repeatedNucleotide=(char)text[0];
	pDescription->RepeatLength=strlen((const char*)text);
        xmlFree(text);
        pDescription->Length++;
      } else {
        pDescription->Length+=strlen((const char*)text);
        if (!strcmp((const char*)pMyNode->name,"FivePrimeFlanking")) {
          pDescription->RepeatBasePosition=strlen((const char*)text);
          flank5Text=text;
        } else
          flank3Text=text;
      }
    } // end of: one of the known pattern nodes?
  } // end of loop over nodesin XML document
  if (!flank3Text || !flank5Text || repeatedNucleotide==' ') {
    fprintf(stderr,"A site pattern has to contain both flanks and a homopolymer run.\n");
    pDescription->Length=0;
    if (flank5Text) xmlFree(flank5Text);
    if (flank3Text) xmlFree(flank3Text);
    return(1);
  }

  /* acquire memory */
  pattern=(char *)malloc((pDescription->Length+1)*sizeof(char));
  if (!pattern) {
    fprintf(stderr,"Could not get memory for repeat sequence!\n");
    pDescription->Length=0;
    xmlFree(flank5Text);
    xmlFree(flank3Text);
    return(1);
  }

  /* save pattern sequence */
  strcpy(pattern, (const char*)flank5Text);
  pattern[pDescription->RepeatBasePosition]=repeatedNucleotide;
  pattern[pDescription->RepeatBasePosition+1]='\0';
  strcat(pattern, (const char*)flank3Text);
  pDescription->Sequence=CreateInternalSequence(pattern);
  xmlFree(flank5Text);
  xmlFree(flank3Text);
  free((void*)pattern);
  if (!pDescription->Sequence) {
    fprintf(stderr,"Site pattern illegal!\n");
    pDescription->Length=0;
    return(1);
  }

  return(0);
}

typedef struct tagPARAMETERDESCRIPTION {
  char *name;     /* name of parameter                  */
  int found;      /* non-zero if parameter was found    */
  int *variable;  /* pointer to variable to be filled   */
  int mandatory;  /* non-zero if parameter is mandatory */
} PARAMETERDESCRIPTION, *PPARAMETERDESCRIPTION;

static int ParseScoring(PMSDESCRIPTION pDescription, xmlNode *pNode)
{
  int i;           /* loop variable         */
  xmlChar *text;   /* text of a node        */
  PARAMETERDESCRIPTION Parameters[]=
    {
      {"Minimum5Flank",0,&(pDescription->Minimal5Flank),1},
      {"Minimum3Flank",0,&(pDescription->Minimal3Flank),1},
      {"MaximumMismatches",0,&(pDescription->MaximalMismatch),1},
      {"MaximumMismatches5Required",0,&(pDescription->MaxMismatch5Req),1},
      {"MaximumMismatches3Required",0,&(pDescription->MaxMismatch3Req),1},
    };
  for(pNode=pNode->children;pNode;pNode=pNode->next) {
    if (pNode->type != XML_ELEMENT_NODE) continue;
    for(i=0;i<sizeof(Parameters)/sizeof(Parameters[0]);i++) {
      if (!strcmp((const char*)pNode->name,Parameters[i].name)) {
	if (Parameters[i].found) {
	  fprintf(stderr,"Two instances of scoring parameter \"%s\"!\n",
		  Parameters[i].name);
	  return(1);
	}
	Parameters[i].found=1;
	text=xmlNodeGetContent(pNode);
	if (!text) {
	  fprintf(stderr,
		  "Problem reading value of scoring parameter \"%s\"!\n",
		  Parameters[i].name);
	  return(1);
	}
	*(Parameters[i].variable)=atoi((const char*)text);
	xmlFree(text);
	break;
      }
    }
  }

  for(i=0;i<sizeof(Parameters)/sizeof(Parameters[0]);i++) {
    if (Parameters[i].mandatory && !Parameters[i].found) {
      fprintf(stderr,"Scoring parameter \"%s\" is missing!\n",
	      Parameters[i].name);
      return(1);
    }
  }

  return(0);
}

static int ParseQuality(PMSDESCRIPTION pDescription, xmlNode *pNode)
{
  int i;           /* loop variable         */
  xmlChar *text;   /* text of a node        */
  PARAMETERDESCRIPTION Parameters[]=
    {
      {"NCutoff",0,&(pDescription->CutoffForNs),1},
      {"MaximumNs",0,&(pDescription->MaximumNs),1},
    };
  for(pNode=pNode->children;pNode;pNode=pNode->next) {
    if (pNode->type != XML_ELEMENT_NODE) continue;
    for(i=0;i<sizeof(Parameters)/sizeof(Parameters[0]);i++) {
      if (!strcmp((const char*)pNode->name,Parameters[i].name)) {
	if (Parameters[i].found) {
	  fprintf(stderr,"Two instances of scoring parameter \"%s\"!\n",
		  Parameters[i].name);
	  return(1);
	}
	Parameters[i].found=1;
	text=xmlNodeGetContent(pNode);
	if (!text) {
	  fprintf(stderr,
		  "Problem reading value of scoring parameter \"%s\"!\n",
		  Parameters[i].name);
	  return(1);
	}
	*(Parameters[i].variable)=atoi((const char*)text);
	xmlFree(text);
	break;
      }
    }
  }

  for(i=0;i<sizeof(Parameters)/sizeof(Parameters[0]);i++) {
    if (Parameters[i].mandatory && !Parameters[i].found) {
      fprintf(stderr,"Scoring parameter \"%s\" is missing!\n",
	      Parameters[i].name);
    }
  }

  return(0);
}

/* load a description from a file */
PMSDESCRIPTION LoadSiteDescription(const char *Filename)
{
  int error;                     /* non-zero if error occurred          */
  int foundPattern;              /* non-zero if pattern found           */
  int foundScoring;              /* non-zero if scoring parameters found*/
  int foundQuality;              /* non-zero if quality filter found    */
  PMSDESCRIPTION pDescription;   /* site description to be loaded       */
  xmlDocPtr pXMLDocument;        /* the XML document in the file        */
  xmlNode *pNode;                /* current node in XML document        */

  /* parse XML file */
  xmlInitParser();
  pXMLDocument=xmlParseFile(Filename);
  if (!pXMLDocument) {
    fprintf(stderr,"Could not parse site description file!\n");
    return(NULL);
  }

  /* get memory for site description */
  pDescription=(PMSDESCRIPTION)malloc(sizeof(MSDESCRIPTION));
  if (!pDescription) {
    fprintf(stderr,"Could not allocate memory for site description!\n");
    xmlFreeDoc(pXMLDocument);
    xmlCleanupParser();
    return(NULL);
  }
  pDescription->Length=0;
  pDescription->Sequence=NULL;
  pDescription->RepeatBasePosition=-1;
  pDescription->Minimal5Flank=5;
  pDescription->Minimal3Flank=5;
  pDescription->MaximalMismatch=10;
  pDescription->MaxMismatch5Req=1;
  pDescription->MaxMismatch3Req=1;
  pDescription->CutoffForNs=20;
  pDescription->MaximumNs=30;

  error=0;
  do { /* region of controlled exit */

    for(pNode=xmlDocGetRootElement(pXMLDocument);pNode;pNode=pNode->next)
      if (pNode->type==XML_ELEMENT_NODE &&
	  !strcmp((const char*)pNode->name,"SiteDescription"))
	break;
    if (!pNode) {
      fprintf(stderr,"Could not find <SiteDescription> in description file!\n");
      error=1;
      break;
    }
    foundPattern=0;
    foundScoring=0;
    foundQuality=0;
    for(pNode=pNode->children;pNode;pNode=pNode->next) {
      if (pNode->type != XML_ELEMENT_NODE) continue;
      if (!strcmp((const char*)pNode->name,"Pattern")) {
	if (foundPattern) {
	  fprintf(stderr,"Site description contains two patterns!\n");
	  error=1;
	  break;
	}
	foundPattern=1;
	error=ParsePatterns(pDescription,pNode);
      }
      if (!strcmp((const char*)pNode->name,"ScoringRules")) {
	if (foundScoring) {
	  fprintf(stderr,"Site description contains two scoring parameter "
		  "sets!\n");
	  error=1;
	  break;
	}
	foundScoring=1;
	error=ParseScoring(pDescription,pNode);
      }
      if (!strcmp((const char*)pNode->name,"QualityFilter")) {
	if (foundQuality) {
	  fprintf(stderr,"Site description contains two quality filter "
		  "parameter sets!\n");
	  error=1;
	  break;
	}
	foundQuality=1;
	error=ParseQuality(pDescription,pNode);
      }
      if (error) break;
    }

    if (!foundPattern) {
      fprintf(stderr,"Site description does not contain a pattern!\n");
      error=1;
      break;
    }

    if (pDescription->RepeatBasePosition<pDescription->Minimal5Flank) {
      fprintf(stderr,"5' flank is shorter than requested minimal 5' flank "
	      "length!\n");
      error=1;
      break;
    }
    
    if (pDescription->Length-pDescription->RepeatBasePosition-1
	<pDescription->Minimal3Flank) {
      fprintf(stderr,"3' flank is shorter than requested minimal 3' flank "
	      "length!\n");
      error=1;
      break;
    }

  } while(0); /* end of region of controlled exit */
  
  xmlFreeDoc(pXMLDocument);
  xmlCleanupParser();

  if (error) DestroySiteDescription(&pDescription);
  
  return(pDescription);
}

int DestroySiteDescription(PMSDESCRIPTION *ppDescription)
{
  if (!ppDescription || !*ppDescription) return(0);
  if ((*ppDescription)->Sequence) free((void*)(*ppDescription)->Sequence);

  free((void *)*ppDescription);
  *ppDescription=NULL;

  return(1);
}
