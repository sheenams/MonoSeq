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

/* declaration of STR description */

#ifndef __MSDESC_H
#define __MSDESC_H

/* description of an str */

typedef struct tagMSDESCRIPTION {
  /* actual pattern */
  int Length;      /* total length of sequence                              */
  int *Sequence;   /* actual sequence from beginning to end                 */
  int RepeatBasePosition; /* position of base repeated in homopolymer run   */
  int RepeatLength;/* length of repeated base in reference                  */

  /* matching conditions */
  int Minimal5Flank;/* minimal number of bases in 5' flanking sequence      */
  int Minimal3Flank;/* minimal number of bases in 3' flanking sequence      */
  int MaximalMismatch;/* maximal number of allowed mismatches               */
  int MaxMismatch5Req;/* maximal number of misatches in 5' required region  */
  int MaxMismatch3Req;/* maximal number of misatches in 3' required region  */
  
  /* quality filter parameters */
  int CutoffForNs;  /* minimal quality score to not turn to N               */
  int MaximumNs;    /* maximal number of Ns in a read                       */
  
} MSDESCRIPTION, *PMSDESCRIPTION;

/* load a description from a file */
PMSDESCRIPTION LoadSiteDescription(const char *Filename);

/* get rid of description */
int DestroySiteDescription(PMSDESCRIPTION *pDescription);

#endif /* ndef __MSDESC_H */
