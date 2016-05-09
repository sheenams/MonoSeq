#
# COPYRIGHT (c)2016 THE OHIO STATE UNIVERSITY
# ALL RIGHTS RESERVED
#
# PERMISSION IS GRANTED TO USE, COPY, CREATE DERIVATIVE WORKS AND
# REDISTRIBUTE THIS SOFTWARE AND SUCH DERIVATIVE WORKS FOR NONCOMMERCIAL
# EDUCATION AND RESEARCH PURPOSES, SO LONG AS NO FEE IS CHARGED, AND SO
# LONG AS THE COPYRIGHT NOTICE ABOVE, THIS GRANT OF PERMISSION, AND THE
# DISCLAIMER BELOW APPEAR IN ALL COPIES MADE; AND SO LONG AS THE NAME OF
# THE OHIO STATE UNIVERSITY IS NOT USED IN ANY ADVERTISING OR PUBLICITY
# PERTAINING TO THE USE OR DISTRIBUTION OF THIS SOFTWARE WITHOUT
# SPECIFIC, WRITTEN PRIOR AUTHORIZATION.
#
# THIS SOFTWARE IS PROVIDED AS IS, WITHOUT REPRESENTATION FROM THE OHIO
# STATE UNIVERSITY AS TO ITS FITNESS FOR ANY PURPOSE, AND WITHOUT
# WARRANTY BY THE OHIO STATE UNIVERSITY OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE OHIO STATE
# UNIVERSITY SHALL NOT BE LIABLE FOR ANY DAMAGES, INCLUDING SPECIAL,
# INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, WITH RESPECT TO ANY
# CLAIM ARISING OUT OF OR IN CONNECTION WITH THE USE OF THE SOFTWARE,
# EVEN IF IT HAS BEEN OR IS HEREAFTER ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.
#

DEBUG =
#DEBUG = -g

# path to header files for libxml2 - adjust to fit
XMLHEADERPATH= /usr/include/libxml2

CCFLAGS = -Wall -O3 ${DEBUG}
LNFLAGS = -Wall -O3 ${DEBUG}

FILES = monoseq.o mscall.o mshmm.o msdesc.o mssequence.o

default:	monoseq

clean:		
		/bin/rm -rf *.o monoseq

distclean:	clean
		/bin/rm -rf *~

monoseq:	${FILES}
		gcc ${LNFLAGS} -o monoseq ${FILES} -lm -lxml2

monoseq.o:	monoseq.c msdesc.h mshmm.h mssequence.h
		gcc -c ${CCFLAGS} monoseq.c

msdesc.o:	msdesc.c msdesc.h mssequence.h
		gcc -c ${CCFLAGS} -I ${XMLHEADERPATH} msdesc.c

mshmm.o:	mshmm.c mshmm.h msdesc.h
		gcc -c ${CCFLAGS} mshmm.c

mscall.o:	mscall.c mscall.h
		gcc -c ${CCFLAGS} mscall.c

mssequence.o:	mssequence.c mssequence.h
		gcc -c ${CCFLAGS} mssequence.c


