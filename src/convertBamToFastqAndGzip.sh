#!/bin/sh

if test $# -lt 2
then
	echo "Usage: $0 INPUTFNAME OUTPUTFnamePREFIX"
	echo
	echo "	Given a sam or bam (*.?am) file which contains paired-end reads,"
	echo "	  1. call picard to convert them to two fastq files."
	echo "	  2. gzip the fastq files"
	echo "Note:
		"
	echo
	echo "Examples:	~/script//shell/convertBamToFastqAndGzip.sh ./fastq/gerald_81GPJABXX_5_TTAGGC.bam db/individual_sequence/1_3"
exit
fi

inputFname=$1
outputFnamePrefix=$2

parameter=''
#while test -n "$6"
#do
#parameter=$parameter' '$6
#shift
#done
picard_tool_path=~/script/picard/dist/
java -jar $picard_tool_path/SamToFastq.jar INPUT=$inputFname F=$outputFnamePrefix\_1.fastq F2=$outputFnamePrefix\_2.fastq

gzip $outputFnamePrefix\_1.fastq
gzip $outputFnamePrefix\_2.fastq


#for i in `ls $inputDir/*.?am`; do echo $i;
#	fname_prefix=`echo $i|awk -F . '{print $1}'`
#	if test -z $fname_prefix	# . is the first letter in the directory name
#	then
#		fname_prefix=.`echo $i|awk -F . '{print $2}'`
#	fi
#	input_fname=$i
#	output_prefix=$fname_prefix
#	java -jar ~/script/vervet/bin/picard-tools/SamToFastq.jar INPUT=$input_fname F=$output_prefix\_1.fastq F2=$output_prefix\_2.fastq
#done
#
#for i in `ls $inputDir/*.fastq`; do echo $i;
#	input_fname=$i
#	gzip $input_fname
#done
