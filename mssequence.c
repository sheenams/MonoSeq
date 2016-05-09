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

/* implementation of sequence functions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mssequence.h"

int IUPACToInternal(char c)
{
  int r;  /* result */
  
  switch(toupper(c)) {
  case 'A':
    r=MONOSEQ_BASE_A;
    break;
  case 'C':
    r=MONOSEQ_BASE_C;
    break;
  case 'G':
    r=MONOSEQ_BASE_G;
    break;
  case 'T':
  case 'U':
    r=MONOSEQ_BASE_T;
    break;
  case 'R':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_G;
    break;
  case 'Y':
    r=MONOSEQ_BASE_C | MONOSEQ_BASE_T;
    break;
  case 'S':
    r=MONOSEQ_BASE_C | MONOSEQ_BASE_G;
    break;
  case 'W':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_T;
    break;
  case 'K':
    r=MONOSEQ_BASE_T | MONOSEQ_BASE_G;
    break;
  case 'M':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_C;
    break;
  case 'B':
    r=MONOSEQ_BASE_C | MONOSEQ_BASE_G | MONOSEQ_BASE_T;
    break;
  case 'D':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_G | MONOSEQ_BASE_T;
    break;
  case 'H':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_C | MONOSEQ_BASE_T;
    break;
  case 'V':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_C | MONOSEQ_BASE_G;
    break;
  case 'N':
    r=MONOSEQ_BASE_A | MONOSEQ_BASE_C | MONOSEQ_BASE_G | MONOSEQ_BASE_T;
    break;
  default:
    r=-1;
    break;
  }
  return(r);
}

/* calculate complementary nucleotide */
int ComplementaryNucleotide(int n)
{
  int r;  /* result */

  r=0;
  if (n & MONOSEQ_BASE_A) r|=MONOSEQ_BASE_T;
  if (n & MONOSEQ_BASE_T) r|=MONOSEQ_BASE_A;
  if (n & MONOSEQ_BASE_G) r|=MONOSEQ_BASE_C;
  if (n & MONOSEQ_BASE_C) r|=MONOSEQ_BASE_G;

  return(r);
}

int *CreateInternalSequence(const char *IUPACString)
{
  const char *p;     /* pointer to current IUPAC character */
  int *Sequence;     /* sequence to be constructed         */
  int *q;            /* pointer to current result          */

  Sequence=(int*)malloc((strlen(IUPACString)+1)*sizeof(int));
  if (!Sequence) return(NULL);

  for(p=IUPACString,q=Sequence;*p;p++,q++) {
    *q=IUPACToInternal(*p);
    if (*q<0) {
      free((void*)Sequence);
      return(NULL);
    }
  }
  *q=0;

  return(Sequence);
}

/* create reverse complement */
int *CreateReverseComplement(const int *forward)
{
  int i;             /* loop variable                      */
  const int *p;      /* pointer to current forward base    */
  int *Sequence;     /* sequence to be constructed         */
  int *q;            /* pointer to current result          */

  /* count length */
  for(p=forward,i=0;*p;p++,i++);

  /* get memory */
  Sequence=(int*)malloc((i+1)*sizeof(int));
  if (!Sequence) return(NULL);

  for(p=forward,q=Sequence+i-1;*p;p++,q--) {
    *q=ComplementaryNucleotide(*p);
    if (*q<0) {
      free((void*)Sequence);
      return(NULL);
    }
  }
  Sequence[i]=0;

  return(Sequence);
}

/* print sequence */
void PrintSequence(const int *sequence, FILE *fp)
{
  const int *p;                      /* pointer to current nucleotide */
  char bases[17]=" ATWGRKDCMYHSVBN"; /* IUPAC codes                   */

  for(p=sequence;*p;p++) {
    if (*p>0 && *p<16) fprintf(fp,"%c",bases[*p]);
  }
}
