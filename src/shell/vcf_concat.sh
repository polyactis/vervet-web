#!/bin/bash
# 2011-9-14 
if test $# -lt 2
then
	echo "Usage: $0 outputVCF inputVCF1 inputVCF2 ..."
	echo
	echo "Note:"
	echo "	1. every input vcf file is gzipped and tabix indexed."
	echo "	2. outputFile will be bgzipped and tabix-indexed."
	echo
	echo "Example:"
	echo "	$0 Contig0_1_3000000.concat.vcf.gz Contig0_1_2000000.gatk.vcf.gz Contig0_2000000_3000000.vcf.gz"
exit
fi

bgzipPath=~/bin/bgzip
tabixPath=~/bin/tabix
#export below will actually fail the perl library search
#export PERL5LIB=~/script/vcftools/lib/
concatPath=~/bin/vcftools/vcf-concat

outputVCF=$1
shift
no_of_inputs=$#
$concatPath $*|$bgzipPath -c >$outputVCF

exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
echo exit code: $exitCode, $exitCode2
if test "$exitCode" = "0" && test "$exitCode2" = "0"
then
	$tabixPath -p vcf $outputVCF
else
	if test "$exitCode" = "0"
	then
		exit $exitCode2
	else
		exit $exitCode
	fi
fi
