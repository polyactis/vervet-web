#!/bin/bash
# 2011-9-7 bwa alignment through unix pipe for paired-end short-read
noOfAlnThreads=4
if test $# -lt 4
then
	echo "Usage: $0 refFastaFname fastqF1 fastqF2 outputBamFnamePrefix"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "	2. This script uses $noOfAlnThreads threads in both aln."
	echo "	3. Unmapped reads will be removed."
	echo "	4. Final bam output is sorted. sorting uses 2G mem."
	echo "	5. shell has to be bash. (...) is used."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta individual_sequence/435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT_1.fastq.gz individual_sequence/435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT_2.fastq.gz 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
refFastaFname=$1
fastqF1=$2
fastqF2=$3
outputBamFnamePrefix=$4

$bwaPath sampe -P $refFastaFname <($bwaPath aln -t $noOfAlnThreads $refFastaFname $fastqF1) <($bwaPath aln -t $noOfAlnThreads $refFastaFname $fastqF2) $fastqF1 $fastqF2 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix
