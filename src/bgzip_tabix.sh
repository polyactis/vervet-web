#!/bin/sh
# 2011-9-14 
if test $# -lt 2
then
	echo "Usage: $0 inputVCF outputVCF"
	echo
	echo "Note:"
	echo
	echo "Example:"
	echo "	$0 Contig0_1_2000000.vcf Contig0_1_2000000.vcf.gz"
exit
fi
#	echo "	4. shell has to be bash. (...) is used."

bgzipPath=~/bin/bgzip
tabixPath=~/bin/tabix


inputVCF=$1
outputVCF=$2

$bgzipPath $inputVCF
exitCode=$?
if test "$exitCode" = "0"
then
	$tabixPath -p vcf $outputVCF
else
	exit $exitCode
fi
