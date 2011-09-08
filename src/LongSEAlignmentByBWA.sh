#!/bin/bash
# 2011-9-7 bwa alignment through unix pipe for long-read single-end data
noOfAlnThreads=4
if test $# -lt 3
then
	echo "Usage: $0 refFastaFname fastqF1 outputBamFnamePrefix"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "	2. This script uses $noOfAlnThreads threads in aln."
	echo "	3. Unmapped reads will be removed."
	echo "	4. Final bam output is sorted. sorting uses 2G mem."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta individual_sequence/1_VRC_ref_454/GINC65G04.fastq.gz 1_VRC_ref_454/GINC65G04"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
refFastaFname=$1
fastqF1=$2
outputBamFnamePrefix=$3

$bwaPath bwasw  -t $noOfAlnThreads $refFastaFname $fastqF1 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix
