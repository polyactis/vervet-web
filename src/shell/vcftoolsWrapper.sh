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
vcftoolsArguments=$*

findValueGivenAnOptionName () {
	if [ -z "$1" ]
	then
		echo "Option Name is not provided."
		echo ;
	else
		optionNamePosition=`echo $vcftoolsArguments|awk -F ' ' '{i=1; while (i<=NF){if ($i=="'$1'") {print i}; ++i}}'`
		#echo $1 position: $optionNamePosition
		if [ -z "$optionNamePosition" ]
		then
			echo;
		else
			optionValuePosition=`echo $optionNamePosition+1|bc`
			#echo optionValuePosition $optionValuePosition
			optionValue=`echo $vcftoolsArguments|awk -F ' ' '{ if ('$optionValuePosition'<=NF) {print $'$optionValuePosition'} else print }'`
			echo $optionValue
		fi
	fi
}

checkVCFFileIfEmpty () {
	fname=$1
	if test -r $fname
	then
		suffix=`echo $fname|awk -F '.' '{print $NF}'`
		if test "$suffix" = "gz"
		then
			numberOfLoci=`gunzip -c $fname|grep -v "^#"|wc -l|awk -F ' ' '{print $1}'`
		else
			numberOfLoci=`grep -v "^#" $fname|wc -l|awk -F ' ' '{print $1}'`
		fi
		if test $numberOfLoci -gt 0
		then
			echo 0;
		else
			echo 1;
		fi
	else
		echo 1;
	fi
}

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

recodeStringPosition=`echo $vcftoolsArguments|awk -F ' ' '{i=1; while (i<=NF){if ($i=="--recode") {print i}; ++i}}'`
echo recodeStringPosition: $recodeStringPosition

if test $isVCFEmpty -eq 0
then
	$vcftoolsPath $vcftoolsArguments
else
	echo "vcftools is not run due to empty/inexistent VCF but touch some files.";
	if test -n "$TsTvWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.TsTv
	fi
	if test -n "$PiWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.windowed.pi
	fi
	if test -n "$recodeStringPosition" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.recode.vcf
	fi
	if test -n "$snpDensityWindowSize" && test -n "$outputNamePrefix"
	then
		touch $outputNamePrefix.snpden
	fi
fi
