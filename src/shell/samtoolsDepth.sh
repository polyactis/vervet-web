#!/bin/bash
# 2012.5.6
minMapQDefault=20
minBaseQDefault=30
if test $# -lt 3
then
	echo "Usage: $0 samtoolsPath inputBam outputDepthFname [minMapQ] [minBaseQ] [moreSAMtoolsArguments]"
	echo
	echo "Note:"
	echo "	1. run 'samtools depth -q minMapQ -Q minBaseQ' on the input bam."
	echo "		The depth output is chr,pos,depth, tab-delimited. could be very large."
	echo "	2. minMapQ is the minimum mapping quality for a read to be considered. default is $minMapQDefault."
	echo "	3. minBaseQ is the minimum base quality for a base to be considered. default is $minBaseQDefault."
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
#after 5 shift, all arguments left are for trioCaller
arguments=$*

if [ -z $minMapQ ]; then
	minMapQ=$minMapQDefault
fi

if [ -z $minBaseQ ]; then
	minBaseQ=$minBaseQDefault
fi
echo minMapQ: $minMapQ
echo minBaseQ: $minBaseQ

#$samtoolsPath depth -q $minMapQ -Q $minBaseQ $arguments $inputBam | gzip > $outputDepthFname
$samtoolsPath depth -q $minMapQ -Q $minBaseQ $arguments $inputBam > $outputDepthFname
