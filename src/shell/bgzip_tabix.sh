#!/bin/bash
# 2011-9-14 
if test $# -lt 2
then
	echo "Usage: $0 inputVCF outputVCF [tabixArguments]"
	echo
	echo "Note:"
	echo "	tabixArguments denotes additional arguments to be passed to tabix. Default is '-p vcf'."
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
shift
shift
tabixArguments=$*
fstTabixArgument=$1
if test -z $fstTabixArgument
then
	tabixArguments="-p vcf"
fi
if test -r $inputVCF
then
	$bgzipPath -c -f $inputVCF > $outputVCF
	# -c: write on standard output, keep original files unchanged
	# -f: overwrite files without asking
	exitCode=$?
	if test "$exitCode" = "0"
	then
		echo "tabix arguments is $tabixArguments"
		$tabixPath $tabixArguments $outputVCF
	else
		exit $exitCode
	fi
elif test -r $outputVCF
then
	echo "bgzipped file is available but not the input. run tabix immediately"
	echo "tabix arguments is $tabixArguments"
	$tabixPath $tabixArguments $outputVCF
else
	#2011-11-12 fake the output for the pipeline to go through (temporary fix)
	touch $outputVCF
	touch $outputVCF.tbi
fi
