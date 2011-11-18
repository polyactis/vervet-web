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

if test -r $inputVCF
then
	$bgzipPath $inputVCF
	exitCode=$?
	if test "$exitCode" = "0"
	then
		$tabixPath -p vcf $outputVCF
	else
		exit $exitCode
	fi
else
	#2011-11-12 fake the output for the pipeline to go through (temporary fix)
	touch $outputVCF
	touch $outputVCF.tbi
fi
