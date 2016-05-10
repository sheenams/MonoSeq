# MonoSeq
A variant caller for homopolymer runs

MonoSeq is a tool that determines variations in homopolymer run length in Illumina sequencing data.  It is provided with raw reads and a definition of the homopolymer run of interest and reports a variant frequency that is adjusted for the variability of homopolymer run length intrinsic to the Illumina platform.

This tool was developed at The Ohio State University by Ralf Bundschuh, Chris Walker, and Paul Goodfellow.

### Installation ###

Download the source files to a directory on a unix based system (including Mac).  The only dependency of the software is the XML library libxml2.  The library itself should be installed on most distributions but you may have to install its header files in a package typically called libxml2-dev or something similar. Typing

    make

should compile the source code.  You may have to adjust the line with XMLHEADERPATH if your XML library headers are in a different directory or the lines with CC and LN if the gnu C compiler is not installed on your system.  Compilation creates the executable monoseq.

### Usage ###

Monoseq analyzes one homopolymer repeat of one sample at a time.  Monoseq is called as

    monoseq [-q] <site_description_file> [<fastq_file>]
    
Here, the site description file is an XML file that defines the homopolymer to analyze.  The file chr10_125780753.xml is included with the source code as an example and contains comments explaining the meaning of every entry in the file.  Essentially, a homopolymer is specified by its flanking sequences and the homopolymer repeat itself.  In addition, the homopolymer definition file contains a few details how reads should be selected.  If in doubt, these parameters should be left at their default values in the sample file.

The second (optional) argument is a file in FASTQ (or - if the "-q" option is given - qseq) format that contains the reads that should be analyzed.  If no second argument is given, monoseq will read this information from standard input.  Monoseq will identify reads containing the homopolymer and a meaningful amount of flanking sequence itself, thus, in principle full sequencing files could be submitted.  However, to decrease the likelihood of wrongly identified reads and more importantly in order to significantly increase speed, it is recommended to extract only the relevant reads that have been aligned to a region surrounding the homopolymer repeat. This can be achieved using samtools with a command line such as

    samtools view -b <bam_file> chr10:125780653-125780853 | samtools fastq - | monoseq chr10_125780753.xml
    
The output of monoseq consists of two lines with the following format:

    # coverage max_homopol wt_homopol variant_frequency counts raw_frequencies called_frequencies
    13  10   9 0.096    0    0    0    0    0    0    0    0    2    9    2 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.154 0.692 0.154 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.904 0.09

that sumarize the homopolymer call at the requested location.  These output values have the following meaning:

1. The first number (here 13) is the total number of reads that are informative for the length of this homopolymer repeat.
2. The second number (here 10) is the number of repeated nucleotides in the longest version of the homopolymer repeat found in the submitted reads.
3. The third number (here 9) is the number of repeated nucleotides in the wild type homopolymer repeat (i.e., the number of nucleotides in the homopolymer description file)
4. The forth number (here 0.096) is the adjusted variant frequency of the homopolymer repeat, i.e., the variant frequency after taking into account sequencing artifacts.
5. The following numbers are the number of reads supporting each homopolymer length.  Since in this example the longest observed length is 10 (as indicated by the 10 under 2.), there are 11 of these numbers corresponding to homopolymer lengths 0 through 10.  In this example, 2 reads supported a homopolymer length of 8, 9 reads supported a homopolymer length of 9, and 2 reads supported a homopolymer length of 10.
6. The following (in this example with a 10 as the maximal under 2. again 11) numbers are the unadjusted frequencies, i.e., the set of the numbers under 5. divided by the total number of reads from 1.
7. The last set of  (in this example with a 10 as the maximal under 2. again 11) numbers are the adjusted frequencies.  These are calculated by taking into account the artifacts of an Illumina sequencer when sequencing homopolymer repeats under the assumption that only two alleles have non-zero true frequency.  The adjusted variant frequency reported in 4. is one minus the adjusted frequency reported for the wild type homopolymer length.
