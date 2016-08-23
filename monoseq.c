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

/* implementation of main program */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msdesc.h"
#include "mssequence.h"
#include "mshmm.h"
#include "mscall.h"

int DealWithOneRead(PMSDESCRIPTION pDescription, int *Workspace,
		    const int *Sequence, int *pLength)
{
  int score;        /* score of sequence                         */

  if (pLength) *pLength=0;

  /* check if read contains site */
  score=ReadContainsSite(Sequence, pDescription, Workspace);
  if (!score) return(0);

  /* find repeat length in read */
  score=FindRepeatLength(Sequence, pDescription, pLength);
  if (!score) return(0);
  
  return(0);
}

int DealWithOneReadData(PMSDESCRIPTION pDescription, int *Workspace,
		        char *sequence, char *quality, char lowestQuality,
		        int *Histogram, int LengthOfHistogram)
{
  int i;            /* loop variable                             */
  int Length;       /* length of homopolymer run which was found */
  int *fsequence;   /* forward sequence                          */
  int *rsequence;   /* reverse complement sequence               */
  char *p;

  /* set low quality bases and dotted bases to N */
  for(i=0,p=sequence;*p;p++,i++) {
    if (*p=='.') *p='N';
    if (!quality[i]) {
      fprintf(stderr,"Missing quality scores in line: %s\n", sequence);
      return(1);
    }
    if (((int)(quality[i]-lowestQuality))<pDescription->CutoffForNs) *p='N';
  }

  /* eliminate trailing N's */
  for(--p;*p=='N' && p>=sequence;p--)
    *p='\0';

  /* count Ns and exit if too big */
  for(p=sequence,i=0;*p;p++)
    if (*p=='N') i++;
  if (i>pDescription->MaximumNs) {
    return(0);
  }
  
  /* get forward version of sequence */
  fsequence=CreateInternalSequence(sequence);
  if (!fsequence) {
    fprintf(stderr,"Sequence %s could not be interpreted!\n",sequence);
    return(1);
  }
  DealWithOneRead(pDescription, Workspace, fsequence, &Length);
  if (Length && Histogram && Length<LengthOfHistogram)
    Histogram[Length]++;

  rsequence=CreateReverseComplement(fsequence);
  free((void *)fsequence);

  /* get reverse complement version of sequence */
  if (!rsequence) {
    fprintf(stderr,
	    "Reverse complement of sequence %s could not be interpreted!\n",
	    sequence);
    return(1);
  }
  DealWithOneRead(pDescription, Workspace, rsequence, &Length);
  if (Length && Histogram && Length<LengthOfHistogram)
    Histogram[Length]++;
  free((void *)rsequence);
  
  return(0);
}


int DealWithOneQseqLine(PMSDESCRIPTION pDescription, int *Workspace,
			char *Line, int *Histogram, int LengthOfHistogram)
{
  int i;            /* loop variable                             */
  char *sequence;   /* beginning of actual read sequence         */
  char *quality;    /* beginning of quality score                */
 
  /* find sequence and quality score in qseq line */
  for(i=0,sequence=Line;*sequence && i<8;sequence++) {
    if (*sequence=='\t') i++;
  }
  if (!*sequence) {
    fprintf(stderr,"Could not find read sequence in line: %s\n", Line);
    return(1);
  }
  for(quality=sequence;*quality && *quality!='\t';quality++);
  if (!*quality) {
    fprintf(stderr,"Could not find quality scores in line: %s\n", Line);
    return(1);
  }
  *(quality++)='\0';
  return(DealWithOneReadData(pDescription, Workspace, sequence, quality, '@',
			     Histogram, LengthOfHistogram));
}

void usage(void)
{
  fprintf(stderr,"monoseq usage: monoseq [options] <site_des> [<fastq_file>]\n");
  fprintf(stderr,"          -q : use qseq file\n");
}

int main(int argc, char **argv)
{
  int i;                        /* loop variable                     */
  int coverage;                 /* total coverage                    */
  int maxlen;                   /* maximal homopolymerlength         */
  int k0,k1;                    /* optimal true homopolymer lengths  */
  double f;                     /* optimal variant frequency         */
  int FirstArg;                 /* first real command line argument  */
  int ReadStdin;                /* non-zero if input from stdin      */
  int WantFastqFile;            /* non-zero if fastq format          */
  int *Workspace;               /* workspace for dynamic programming */
  int Histogram[2000];          /* homopolymer length histogram      */
  char Line[10000];             /* one line from input file          */
  char Sequence[10000];         /* sequence for current record       */
  PMSDESCRIPTION pDescription;  /* description of site               */
  FILE *fp;                     /* qseq or fastq file                */

  /* interpret command line arguments */
  WantFastqFile=1;
  for(FirstArg=1;FirstArg<argc && argv[FirstArg][0]=='-' && argv[FirstArg][1];
      FirstArg++) {
    switch(argv[FirstArg][1]) {
    case 'q':
      WantFastqFile=0;
      break;
    default:
      usage();
      return(1);
    }
  }
  if (argc<FirstArg+1 || argc>FirstArg+2) {
    usage();
    return(1);
  }
  
  /* open sequence file */
  if (argc>FirstArg+1) {
    ReadStdin=0;
    fp=fopen(argv[FirstArg+1],"r");
    if (!fp) {
      fprintf(stderr,"Could not open input file %s!\n",argv[FirstArg+1]);
      return(2);
    }
  } else {
    ReadStdin=1;
    fp=stdin;
  }

  /* read site description */
  pDescription=LoadSiteDescription(argv[FirstArg]);
  if (!pDescription) {
    fprintf(stderr,"Could not load site description from file %s\n",
	    argv[FirstArg]);
    if (!ReadStdin) fclose(fp);
    return(3);
  }

  /* get memory for dynamics programming workspace */
  Workspace=(int *)malloc((pDescription->Length+1)*2*sizeof(int));
  if (!Workspace) {
    fprintf(stderr,"Could not allocate memory for workspace!\n");
    if (!ReadStdin) fclose(fp);
    DestroySiteDescription(&pDescription);
    return(4);
  }

  /* initialize histogram */
  for(i=0;i<sizeof(Histogram)/sizeof(Histogram[0]);i++)
    Histogram[i]=0;

  i=0;
  while(!feof(fp) && fgets(Line,sizeof(Line),fp)) {
    if (WantFastqFile) {
      /* fastq file */
      if (i%2==1) {
	Line[strlen(Line)-1]='\0';
      }
      if (++i==2) {
	strcpy(Sequence,Line);
      } else if (i==4) {
	if (DealWithOneReadData(pDescription, Workspace, Sequence, Line,
				'!', Histogram,
				sizeof(Histogram)/sizeof(Histogram[0]))) {
	  fprintf(stderr,"Error interpreting fastq file line %s!\n",Line);
	  if (!ReadStdin) fclose(fp);
	  DestroySiteDescription(&pDescription);
	  free((void*)Workspace);
	  return(5);
	}
	i=0;
      }
    } else {
      /* qseq file */
      if (DealWithOneQseqLine(pDescription, Workspace, Line, Histogram,
			      sizeof(Histogram)/sizeof(Histogram[0]))) {
	fprintf(stderr,"Error interpreting qseq file line %s!\n",Line);
	if (!ReadStdin) fclose(fp);
	DestroySiteDescription(&pDescription);
	free((void*)Workspace);
	return(5);
      }
    }
  }

  /* find the best two alleles and their relative frequency */
  f=FindVariantFrequency(Histogram, &k0, &k1);

  /* find total number of reads and maximum homopolymer length */
  maxlen=0;
  coverage=0;
  for(i=0;i<sizeof(Histogram)/sizeof(Histogram[0]);i++) {
    if (Histogram[i]) {
      maxlen=i;
      coverage+=Histogram[i];
    }
  }
  if (k0>maxlen) maxlen=k0;
  if (k1>maxlen) maxlen=k1;
  
  /* output results */
  printf("# coverage max_homopol wt_homopol variant_frequency counts raw_frequencies called_frequencies\n");
  printf("%1d %1d %1d ",coverage,maxlen,pDescription->RepeatLength);
  if (pDescription->RepeatLength==k0) {
    printf("%5.3f",f);
  } else if (pDescription->RepeatLength==k1) {
    printf("%5.3f",1-f);
  } else {
    printf("%5.3f",1.0);
  }
  for(i=0;i<=maxlen;i++) {
    printf(" %4d",Histogram[i]);
  }
  for(i=0;i<=maxlen;i++) {
    printf(" %5.3f",(double)Histogram[i]/(double)coverage);
  }
  for(i=0;i<=maxlen;i++) {
    if (i==k0) {
      printf(" %5.3f",1-f);
    } else if (i==k1) {
      printf(" %5.3f",f);
    } else {
      printf(" %5.3f",0.0);
    }
  }
  printf("\n");
  
  /* release all resources */
  if (!ReadStdin) fclose(fp);
  DestroySiteDescription(&pDescription);
  free((void*)Workspace);
  
  return(0);
}


