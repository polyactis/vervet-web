#!/bin/bash
# 2011-11-20 run samtools calmd on input bam
# http://sourceforge.net/mailarchive/forum.php?thread_name=201010141156.33921.ashish%40strandls.com&forum_name=samtools-announce

if test $# -lt 3
then
	echo "Usage: $0 inputBam refFastaFname outputBam"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "Example:"
	echo "	$0 NoBQTag.bam 120_480K_supercontigs.fasta withBQTag.bam"
exit
fi
#	echo "	4. shell has to be bash. (...) is used."

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
picardPath=~/script/picard/dist/
inputBam=$1
refFastaFname=$2
outputBam=$3


$samtoolsPath calmd -br $inputBam $refFastaFname >$outputBam

