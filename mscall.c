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

/* implementation of actual variant frequency calling */

#include <math.h>
#include "mscall.h"

/* observed distributions for every length */
double ObservedProbability[MAXTRUEREPEATLENGTH][MAXOBSERVEDREPEATLENGTH]={
  {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 1 */
  {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 2 */
  {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 3 */
  {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 4 */
  {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 5 */
  {0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},  /* 6 */
  {0.000000,0.000000,0.000029,0.000000,0.000206,0.017441,0.974029,0.007912,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},  /* 7 */
  {0.000000,0.000083,0.000000,0.000000,0.000000,0.001583,0.043583,0.939083,0.015250,0.000250,0.000000,0.000000,0.0,0.0,0.0},  /* 8 */
  {0.000000,0.000000,0.000000,0.000000,0.000200,0.000200,0.003800,0.078200,0.876200,0.038800,0.002200,0.000000,0.000000,0.000000,0.000000},  /* 9 */
  {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.001333,0.016000,0.174000,0.725333,0.075667,0.007000,0.000667,0.000000,0.000000}, /* 10 */
  {0.000000,0.000000,0.000000,0.000000,0.000200,0.000000,0.000000,0.002600,0.032800,0.202600,0.642600,0.104600,0.012600,0.001200,0.000000}, /* 11 */
};

/* non-zero if observed probabilities have already been rescaled */
int ObservedProbabilitiesRescaled=0;

/* calculate the log-likelihood given two variants k0 and k1 and
   frequency f of k1 */
double CalculateLogLikelihood(int *N, int k0, int k1, double f)
{
  int i;        /* loop variable                */
  double s;     /* summation variable           */
  double d;     /* intermediate result variable */

  s=0.0;
  for(i=0;i<MAXOBSERVEDREPEATLENGTH;i++) {
    if (N[i]>0) {
      d=(1.0-f)*ObservedProbability[k0][i]+f*ObservedProbability[k1][i];
      s+=((double)N[i])*log(d);
    }
  }
  return(s);
}

/* calculate the derivative of the log-likelihood given two variants
   k0 and k1 and frequency f of k1 */
double CalculateDerivative(int *N, int k0, int k1, double f)
{
  int i;        /* loop variable                */
  double s;     /* summation variable           */

  s=0.0;
  for(i=0;i<MAXOBSERVEDREPEATLENGTH;i++) {
    if (N[i]>0) {
      if (ObservedProbability[k0][i]==0.0) {
        s+=((double)N[i])/f;
      } else {
        s+=((double)N[i])*
          (ObservedProbability[k1][i]-ObservedProbability[k0][i])/
          ((1.0-f)*ObservedProbability[k0][i]+f*ObservedProbability[k1][i]);
      }
    }
  }
  return(s);
}

/* maximize the log-likelihood if k0 and k1 are the two real alleles */
double FindOptimalMixture(int *N, int k0, int k1)
{
  double der;          /* derivative               */
  double fmin,fmax;    /* current limits           */
  double f;            /* current fraction         */

  /* initialize limits */
  fmin=0.00001;
  f=fmax=0.99999;

  /* check that initial derivative is positive */
  der=CalculateDerivative(N, k0, k1, fmin);
  if (der>0.0) {
    /* check that final derivative is negative */
    der=CalculateDerivative(N, k0, k1, fmax);
    if (der<0.0) {
      while(fmax-fmin>0.00001) {
        f=(fmax+fmin)/2.0;
        der=CalculateDerivative(N, k0, k1, f);
        if (der==0.0) break;
        if (der>0.0)
          fmin=f;
        else
          fmax=f;
      } /* end of loop over refinements */
      return(f);
    }
  }
  /* at this point the derivative did not change sign - the maximum must
     be at a boundary */
  der=CalculateLogLikelihood(N, k0, k1, fmin);
  f=CalculateLogLikelihood(N, k0, k1, fmax);
  return(der>f?fmin:fmax);
}

/* find the best pair of alleles and their variant frequency */
double FindVariantFrequency(int *N, int *pk0, int *pk1)
{
  int i;           /* loop variable                             */
  int k0,k1;       /* current putative true homopolymer lengths */
  int k0opt,k1opt; /* so far best true homopolymer lengths      */
  double f,fopt;   /* current and so far optimal frequency      */
  double L,Lopt;   /* current and so far optimal log-likelihood */

  /* if the observed probabilities were never rescaled we have to do
     this now, since zero probabilities result in numerical problems */
  if (!ObservedProbabilitiesRescaled) {
    for(k0=0;k0<MAXTRUEREPEATLENGTH;k0++) {
      f=0.0;
      for(i=0;i<MAXOBSERVEDREPEATLENGTH;i++) {
	ObservedProbability[k0][i]+=0.0001;
	f+=ObservedProbability[k0][i];
      }
      for(i=0;i<MAXOBSERVEDREPEATLENGTH;i++) {
	ObservedProbability[k0][i]/=f;
      }
    }
    ObservedProbabilitiesRescaled=1;
  }
  
  Lopt=-1e38;
  k0opt=k1opt=0;
  fopt=0.0;

  for(k0=4;k0<MAXTRUEREPEATLENGTH-1;k0++) {
    for(k1=k0+1;k1<MAXTRUEREPEATLENGTH;k1++) {
      f=FindOptimalMixture(N,k0,k1);
      L=CalculateLogLikelihood(N,k0,k1,f);
      if (L>Lopt) {
        Lopt=L;
        k0opt=k0;
        k1opt=k1;
        fopt=f;
      }
    }
  }

  if (pk0) *pk0=k0opt;
  if (pk1) *pk1=k1opt;

  return(fopt);
}

/*
int main(int argc, char **argv)
{

  for(i=0;i<MAXTRUEREPEATLENGTH;i++) {
    if (i==k0)
      printf("%f ",1-f);
    else if (i==k1)
      printf("%f ",f);
    else
      printf("0.0 ");
  }
  printf("\n");

  return(0);
}

*/
