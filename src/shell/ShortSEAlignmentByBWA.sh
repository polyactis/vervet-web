#!/bin/bash
# 2011-9-7 bwa sampe piped to gzipped sam for single-end short-read
noOfAlnThreads=4
if test $# -lt 4
then
	echo "Usage: $0 refFastaFname saiF1 fastqF1 outputSamGzipFname"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta individual_sequence/435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT_1.fastq.gz 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT.sam.gz"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
refFastaFname=$1
saiF1=$2
fastqF1=$3
outputSamGzipFname=$4

#$bwaPath samse $refFastaFname $saiF1 $fastqF1 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix
$bwaPath samse $refFastaFname $saiF1 $fastqF1 | gzip >$outputSamGzipFname
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
