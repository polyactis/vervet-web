#!/bin/bash
# 2011-9-7 bwasw alignment piped to gzipped sam for long-read single-end data
noOfAlnThreads=4
if test $# -lt 3
then
	echo "Usage: $0 refFastaFname fastqF1 outputSamGzipFname [noOfAlnThreads]"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "	2. This script uses $noOfAlnThreads threads in bwasw by default."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta individual_sequence/1_VRC_ref_454/GINC65G04.fastq.gz 1_VRC_ref_454/GINC65G04.sam.gz"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
refFastaFname=$1
fastqF1=$2
outputSamGzipFname=$3

if test -n "$4"
then
	noOfAlnThreads=$4
fi

#$bwaPath bwasw  -t $noOfAlnThreads $refFastaFname $fastqF1 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix
$bwaPath bwasw  -t $noOfAlnThreads $refFastaFname $fastqF1 | gzip >$outputSamGzipFname

exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
echo exit code: $exitCode, $exitCode2
if test "$exitCode" = "0" && test "$exitCode2" = "0"
then
	#normal exit
	exit
else
	if test "$exitCode" = "0"
	then
		finalExitCode=$exitCode2
	else
		finalExitCode=$exitCode
	fi
	exit $finalExitCode
fi
