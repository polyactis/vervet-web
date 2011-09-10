#!/bin/bash
# 2011-9-7 bwa alignment through unix pipe for short-read single-end data
noOfAlnThreads=4
if test $# -lt 4
then
	echo "Usage: $0 refFastaFname saiF1 fastqF1 outputBamFnamePrefix"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "	2. Unmapped reads will be removed."
	echo "	3. Final bam output is sorted. sorting uses 2G mem."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta individual_sequence/435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT_1.fastq.gz 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
refFastaFname=$1
saiF1=$2
fastqF1=$3
outputBamFnamePrefix=$4

$bwaPath samse $refFastaFname $saiF1 $fastqF1 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix
