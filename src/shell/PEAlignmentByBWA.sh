#!/bin/bash
# 2011-9-7 bwa sampe piped to gzipped sam for paired-end short-read
noOfAlnThreads=3
if test $# -lt 6
then
	echo "Usage: $0 refFastaFname saiF1 saiF2 fastqF1 fastqF2 outputSamGzipFname"
	echo
	echo "Note:"
	echo "	1. refFastaFname must have been indexed by bwa."
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta gerald_81GPBABXX_8_CGATGT_1.sai gerald_81GPBABXX_8_CGATGT_1.fastq.gz gerald_81GPBABXX_8_CGATGT_2.sai gerald_81GPBABXX_8_CGATGT_2.fastq.gz 435_6136_2004030_GA_5/gerald_81GPBABXX_8_CGATGT.sam.gz"
exit
fi
#	echo "	4. shell has to be bash. (...) is used."

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
picardPath=~/script/picard/dist/
refFastaFname=$1
saiF1=$2
saiF2=$3
fastqF1=$4
fastqF2=$5
outputSamGzipFname=$6

#$bwaPath sampe -P $refFastaFname $saiF1 $saiF2 $fastqF1 $fastqF2 | $samtoolsPath view  -F 4 -bSh - | $samtoolsPath sort -m 2000000000 - $outputBamFnamePrefix

$bwaPath sampe -P $refFastaFname $saiF1 $saiF2 $fastqF1 $fastqF2 | gzip >$outputSamGzipFname

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
#$bwaPath sampe -P $refFastaFname $saiF1 $saiF2 $fastqF1 $fastqF2 | java -jar $picardPath/ViewSam.jar I=/dev/stdin ALIGNMENT_STATUS= Aligned | gzip >$outputSamGzipFname

#java -jar $picardPath/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT=/dev/stdout | java -jar $picardPath/SortSam.jar SORT_ORDER=coordinate INPUT=/dev/sdin OUTPUT=$outputBamFnamePrefix.bam
