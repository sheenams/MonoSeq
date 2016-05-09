# MonoSeq
A variant caller for homopolymer runs

MonoSeq is a tool that determines variations in homopolymer run length in Illumina sequencing data.  It is provided with raw reads and a definition of the homopolymer run of interest and reports a variant frequency that is adjusted for the variability of homopolymer run length intrinsic to the Illumina platform.

This tool was developed at The Ohio State University by Ralf Bundschuh, Chris Walker, and Paul Goodfellow.

### Installation ###

Download the source files to a directory on a unix based system (including Mac).  The only dependency of the software is the XML library libxml2 which should be installed on most distributions. Typing

    make

should compile the source code.  You may have to adjust the line with XMLHEADERPATH if your XML library headers are in a different directory or replace gcc with your compiler's name if the gnu C compiler is not installed on your system.

