#!/bin/bash
if test $# -lt 3
then
	echo "Usage: $0 vcftoolsPath vcftoolsArgument1 vcftoolsArgument12 ..."
	echo
	echo "Note:"
	echo "	1. 2011-11-12 wrapper of vcftools"
	echo "	2. it's fail-safe wrapper around vcftools, in case the input vcf is empty (no loci, regardless of header)"
	echo "		it'll do following things with the purpose of making downstream job go through.:"
	echo "		1. it'll copy the input vcf if --recode option is on."
	echo "		2. it'll generate an empty .TiTv output."
	echo "		3. it'll generate an empty .windowed.pi file."
	echo
	echo "Example:"
	echo "	$0 ~/bin/vcftools/vcftools --gzvcf Contig0.vcf.gz --positions /tmp/snps_to_be_included.tsv --mac 2 --geno 1 --recode --out /tmp/gatkContig0 --freq --counts --TsTv 100000 --window-pi 100000 --depth --het --hardy"
exit
fi

vcftoolsPath=$1
shift
#after shift, all arguments left are for vcftools
arguments=$*

shellDir=`dirname $0`
source $shellDir/commonWrapper.sh

vcfInputFname=`findValueGivenAnOptionName --vcf`
echo vcfInputFname: $vcfInputFname

gzvcfInputFname=`findValueGivenAnOptionName --gzvcf`
echo gzvcfInputFname: $gzvcfInputFname
if [ -n "$vcfInputFname" ]
then
	isVCFEmpty=`checkVCFFileIfEmpty $vcfInputFname`
else
	isVCFEmpty=`checkVCFFileIfEmpty $gzvcfInputFname`
fi
echo isVCFEmpty: $isVCFEmpty
TsTvWindowSize=`findValueGivenAnOptionName --TsTv`
echo TsTvWindowSize: $TsTvWindowSize
PiWindowSize=`findValueGivenAnOptionName --window-pi`
echo PiWindowSize: $PiWindowSize
snpDensityWindowSize=`findValueGivenAnOptionName --SNPdensity`
echo snpDensityWindowSize: $snpDensityWindowSize
outputNamePrefix=`findValueGivenAnOptionName --out`
echo outputNamePrefix: $outputNamePrefix
outputVCFFname=$outputNamePrefix.recode.vcf

recodeStringPosition=`echo $arguments|awk -F ' ' '{i=1; while (i<=NF){if ($i=="--recode") {print i}; ++i}}'`
echo recodeStringPosition: $recodeStringPosition


if test $isVCFEmpty -eq 0
then
	$vcftoolsPath $arguments
	exitCode=$?
	if test $exitCode -ne 0; then
		echo "error in running vcftools, exiting ..."
		exit $exitCode;
	fi
	if test -n "$recodeStringPosition" && test -n "$outputNamePrefix"
	then
		recodeVCF=$outputNamePrefix.recode.vcf
		if test -r $recodeVCF
		then
			echo "recode vcf output is available"
		else
			outputEmptyVCFWithInputHeader
		fi
	fi
else
	echo "vcftools is not run due to empty/inexistent VCF. need to touch some files so downstream jobs wouldn't die.";
	touch $outputNamePrefix.log
	if test -n "$TsTvWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.TsTv
		touch $outputNamePrefix.TsTv.summary
	fi
	if test -n "$PiWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.windowed.pi
	fi
	if test -n "$recodeStringPosition" && test -n "$outputNamePrefix"
	then
		outputEmptyVCFWithInputHeader
	fi
	if test -n "$snpDensityWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.snpden
	fi
fi
