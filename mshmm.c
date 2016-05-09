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

/* implementation of business logic */

#include <stdio.h>
#include <stdlib.h>
#include "mshmm.h"

/* very large negative number */
#define NEGINFINITY -1000000000

/* check if read contains site */
int ReadContainsSite(const int *Read, PMSDESCRIPTION pDescription,
		     int *Workspace)
{
  /* Workspace has to be 2*(length+1) */
  int i,j;          /* loop variable                            */
  int *current;     /* current scores                           */
  int *new;         /* new scores                               */
  int *swap;        /* score swap variable                      */
  const int *pCRead;/* current letter in read                   */
  int *pPattern;    /* current letter in pattern                */
  int *pLast;       /* score at last read base                  */
  int *pNew;        /* score at new base                        */
  int score;        /* score for next position                  */
  int maxscore;     /* maximal score at read position           */
  int newscore;     /* trial score                              */
  int rlen;         /* read length                              */

  /* determine read length */
  for(rlen=0,pCRead=Read;*pCRead;pCRead++,rlen++);
  
  /* initialization */
  current=Workspace;
  new=Workspace+pDescription->Length+1;
  for(i=0;i<=pDescription->RepeatBasePosition-pDescription->Minimal5Flank;i++)
    current[i]=0;
  for(;i<=pDescription->Length;i++)
    current[i]=NEGINFINITY;
    
  /* actual recursion */
  maxscore=score=0;
  for(i=1,pCRead=Read;*pCRead;i++,pCRead++) {
    new[0]=i;
    for(j=1,pPattern=pDescription->Sequence,pLast=current,pNew=new+1;
	*pPattern;pPattern++,j++,pLast++,pNew++) {
      score=*pLast;
      if (j==pDescription->RepeatBasePosition+1) {
	newscore=current[j];
	if (newscore>score) score=newscore;
      }
      if ((*pPattern) & (*pCRead)) {
	score++;
      }
      *pNew=score;
    } /* end of loop over positions in site */
    score+=rlen-i;
    if (score>maxscore) maxscore=score;
    /* switch working arrays */
    swap=new;
    new=current;
    current=swap;
  } /* end of loop over positions in read */

  /* check alignments that run into the end of the read */
  for(i=pDescription->RepeatBasePosition+pDescription->Minimal3Flank+1;
      i<pDescription->Length;i++) {
    if (current[i]>maxscore) maxscore=current[i];
  }
  
  if (score<rlen-pDescription->MaximalMismatch) score=0;

  if (pDescription->Length-maxscore>2*pDescription->MaximalMismatch)
     maxscore=0;

  /* 0 if site not found, actual number of matches otherwise */
  return(maxscore);
}

/* find repeat length in read */
int FindRepeatLength(const int *Read, PMSDESCRIPTION pDescription,
		     int *pRepeatLength)
{
  int i,j;          /* loop variable                            */
  int rpos;         /* position of maximum in read              */
  int ppos;         /* position of maximum in pattern           */
  int mmcount5;     /* mismatch count in 5' flank               */
  int mmcounttot;   /* total mismatch count in                  */
  int mmcount3;     /* mismatch count in 3' flank               */
  int **w;          /* score matrix                             */
  int **o;          /* origins for all the recursion positions  */
  const int *pCRead;/* current letter in read                   */
  int *pPattern;    /* current letter in pattern                */
  int *pLast;       /* score at last read base                  */
  int *pNew;        /* score at new base                        */
  int score;        /* score for next position                  */
  int newscore;     /* trial score                              */
  int rlen;         /* read length                              */

  /* initialize return values */
  if (pRepeatLength) *pRepeatLength=0;
  
  /* determine read length */
  for(rlen=0,pCRead=Read;*pCRead;pCRead++,rlen++);

  /* get memory for scoring array */
  for(pCRead=Read,rlen=0;*pCRead;pCRead++,rlen++);
  w=(int **)malloc(2*(rlen+1)*sizeof(int *));
  if (!w) return(0);
  w[0]=(int *)malloc(2*(rlen+1)*(pDescription->Length+1)*sizeof(int));
  if (!w[0]) {
    free((void *)w);
    return(0);
  }
  for(i=1;i<=rlen;i++)
    w[i]=w[0]+i*(pDescription->Length+1);
  o=w+rlen+1;
  for(i=0;i<=rlen;i++)
    o[i]=w[0]+(i+rlen+1)*(pDescription->Length+1);
  
  /* initialization */
  for(i=0;i<=pDescription->RepeatBasePosition-pDescription->Minimal5Flank;i++) {
    w[0][i]=0;
    o[0][i]=-1;
  }
  for(;i<=pDescription->Length;i++) {
    w[0][i]=NEGINFINITY;
    o[0][i]=-1;
  }
    
  /* actual recursion */
  for(i=1,pCRead=Read;*pCRead;i++,pCRead++) {
    w[i][0]=i;
    o[i][0]=-1;
    for(j=1,pPattern=pDescription->Sequence,pLast=w[i-1],pNew=w[i]+1;
	*pPattern;pPattern++,j++,pLast++,pNew++) {
      score=*pLast;
      o[i][j]=j-1;
      if (j==pDescription->RepeatBasePosition+1) {
	newscore=w[i-1][j];
	if (newscore>score) {
	  score=newscore;
	  o[i][j]=j;
	}
      }
      if ((*pPattern) & (*pCRead)) {
	score++;
      }
      *pNew=score;
    } /* end of loop over positions in site */
  } /* end of loop over positions in read */

  /* find high score */
  score=NEGINFINITY;
  rpos=-1;
  ppos=-1;
  for(i=0;i<=rlen;i++) {
    if (w[i][pDescription->Length]+rlen-i>score) {
      score=w[i][pDescription->Length]+rlen-i;
      rpos=i;
      ppos=pDescription->Length;
    }
  }
  for(i=pDescription->RepeatBasePosition+pDescription->Minimal3Flank+1;
      i<pDescription->Length;i++) {
    if (w[rlen][i]>score) {
      score=w[rlen][i];
      rpos=rlen;
      ppos=i;
    }
  }

  if (score<rlen-pDescription->MaximalMismatch)
    score=0;
  
  if (score) {
    mmcount3=0;
    mmcount5=0;
    mmcounttot=0;
    rlen=0;
    for(i=ppos,pCRead=Read+rpos-1;
	rpos>=0 && i>=0;i=o[rpos--][i],pCRead--) {
      if (i==pDescription->RepeatBasePosition+1) rlen++;
      if (!(pDescription->Sequence[i-1] & *pCRead)) mmcounttot++;
      if (!(pDescription->Sequence[i-1] & *pCRead) || *pCRead==15) {
	if (i>pDescription->RepeatBasePosition+1 &&
	    i<=pDescription->RepeatBasePosition+
	    pDescription->Minimal3Flank+1) {
	  mmcount3++;
	} else if (i<=pDescription->RepeatBasePosition &&
		   i>pDescription->RepeatBasePosition-
		   pDescription->Minimal3Flank) {
	  mmcount5++;
	}
      }
    }
    
    /* check if maximal number of mismatches in required flank regions OK */
    if (mmcounttot>pDescription->MaximalMismatch) score=0;
    if (mmcount5>pDescription->MaxMismatch5Req) score=0;
    if (mmcount3>pDescription->MaxMismatch3Req) score=0;

    if (score && pRepeatLength) *pRepeatLength=rlen;
  }

  /* get rid of all the memory */
  free((void *)w[0]);
  free((void *)w);

  /* 0 if site not found, actual number of matches otherwise */
  return(score);
}
