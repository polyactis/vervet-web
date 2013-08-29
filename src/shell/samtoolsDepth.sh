#!/bin/bash
# 2012.5.6
minMapQDefault=0
minBaseQDefault=0
if test $# -lt 3
then
	echo "Usage: $0 samtoolsPath inputBam outputDepthFname [minMapQ] [minBaseQ] [moreSAMtoolsArguments]"
	echo
	echo "Note:"
	echo "	1. run 'samtools depth -q minMapQ -Q minBaseQ' on the input bam."
	echo "		The depth output is chr,pos,depth, tab-delimited. could be very large."
	echo "	2. minMapQ is the minimum mapping quality for a read to be considered. default is $minMapQDefault."
	echo "	3. minBaseQ is the minimum base quality for a base to be considered. default is $minBaseQDefault."
	echo "	4. outputDepthFname could be gzipped or not. It detects it and output accordingly."
	echo
	echo "Example:"
	echo "	$0 ~/bin/samtools input.bam output.depth.tsv.gz "
	echo "	$0 ~/bin/samtools input.bam output.depth.tsv "
exit
fi

samtoolsPath=$1
inputBam=$2
outputDepthFname=$3
minMapQ=$4
minBaseQ=$5
shift
shift
shift
shift
shift
#after 5 shift, all arguments left
arguments=$*

if [ -z $minMapQ ]; then
	minMapQ=$minMapQDefault
fi

if [ -z $minBaseQ ]; then
	minBaseQ=$minBaseQDefault
fi
echo minMapQ: $minMapQ
echo minBaseQ: $minBaseQ

commandline="$samtoolsPath depth -q $minMapQ -Q $minBaseQ $arguments $inputBam"
#commandline="$samtoolsPath depth  $arguments $inputBam"
outputFilenameLength=`expr length $outputDepthFname`
gzSuffixStartPosition=`echo $outputFilenameLength-3+1|bc`
gzSuffix=`expr substr $outputDepthFname $gzSuffixStartPosition 3`
echo gzSuffix is $gzSuffix.

if test "$gzSuffix" = ".gz"; then
	$commandline | gzip > $outputDepthFname
	exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"	#must be together in one line. PIPESTATUS[1] in subsequent lines has different meaning.
	exitCode1=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
	exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
	
	echo "exit codes: $exitCode1, $exitCode2"
	
	if test "$exitCode1" = "0" && test "$exitCode2" = "0"
	then
		exit 0
	else
		exit 3
	fi
else
	$commandline > $outputDepthFname
fi
