#!/bin/bash
# 2011-9-14 
if test $# -lt 2
then
	echo "Usage: $0 outputVCF inputVCF1 inputVCF2 ..."
	echo
	echo "a script taking intersection (AND operator) of >=2 vcf files. using vcftools/vcf-isec"
	echo "Note:"
	echo "	1. every input vcf file is gzipped and tabix indexed."
	echo "	2. outputFile will be bgzipped and tabix-indexed."
	echo
	echo "Example:"
	echo "	$0 Contig0_1_2000000.isec.vcf.gz Contig0_1_2000000.gatk.vcf.gz Contig0_1_2000000.vcf.gz"
exit
fi

bgzipPath=~/bin/bgzip
tabixPath=~/bin/tabix
#export below will actually fail the perl library search
#export PERL5LIB="~/script/vcftools/lib/"
isecPath=~/bin/vcftools/vcf-isec

outputVCF=$1
shift
no_of_inputs=$#
#echo $no_of_inputs
echo $isecPath -n +$no_of_inputs $*
$isecPath -n +$no_of_inputs $*|$bgzipPath -c >$outputVCF
exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`

echo exit code: $exitCode, $exitCode2
#exitCode=255 or 2 of vcf-isec is a warning
#Warning: Read 0 lines from samtools/Contig687_1_1258568.vcf.gz, the tabix index may be broken.
#at /home/crocea/bin/vcftools/vcf-isec line 21
#main::error('Warning: Read 0 lines from samtools/Contig687_1_1258568.vcf.g...') called at /home/crocea/bin/vcftools/vcf-isec line 87
#main::__ANON__('Warning: Read 0 lines from samtools/Contig687_1_1258568.vcf.g...') called at /home/crocea/bin/vcftools/vcf-isec line 272
#main::vcf_isec('HASH(0x2339900)') called at /home/crocea/bin/vcftools/vcf-isec line 12

if ( test "$exitCode" = "0" || test "$exitCode" = "255"  || test "$exitCode" = "2" ) && test "$exitCode2" = "0"
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
