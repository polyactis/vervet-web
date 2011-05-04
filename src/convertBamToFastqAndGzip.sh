#!/bin/sh
#$ -pe mpich 5

if test $# -lt 1
then
	echo "Usage: $0 INPUTDIR"
	echo
	echo "	Given sam or bam (*.?am) files in INPUTDIR which contain paired-end reads,"
	echo "	  1. call picard to convert them to fastq."
	echo "	  2. gzip all fastq files"
	echo "Note:
	1. '.' could only be included once in the INPUTDIR, i.e. in the beginning.
	   because it is used as a separator in awk to extract filename prefix.
		"
	echo
	echo "Examples:	~/script//shell/convertBamToFastqAndGzip.sh ./fastq/"
exit
fi

inputDir=$1
#outputDir=$2

parameter=''
#while test -n "$6"
#do
#parameter=$parameter' '$6
#shift
#done

for i in `ls $inputDir/*.?am`; do echo $i;
	fname_prefix=`echo $i|awk -F . '{print $1}'`
	if test -z $fname_prefix	# . is the first letter in the directory name
	then
		fname_prefix=.`echo $i|awk -F . '{print $2}'`
	fi
	input_fname=$i
	output_prefix=$fname_prefix
	java -jar ~/script/vervet-web/bin/picard-tools/SamToFastq.jar INPUT=$input_fname F=$output_prefix\_1.fastq F2=$output_prefix\_2.fastq
done

for i in `ls $inputDir/*.fastq`; do echo $i;
	input_fname=$i
	gzip $input_fname
done
