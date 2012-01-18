#!/bin/bash
if test $# -lt 3
then
	echo "Usage: $0 trioCallerPath vcftoolsArgument1 vcftoolsArgument12 ..."
	echo
	echo "Note:"
	echo "	1. 2011-12-13 wrapper of TrioCaller"
	echo "	2. it's fail-safe wrapper around trioCaller, in case the input vcf is empty (no loci, header may exist)"
	echo "		it'll do following things with the purpose of making downstream job go through.:"
	echo "		1. it'll copy the input vcf into the output (to preserve the header)."
	echo
	echo "Example:"
	echo "	$0 ~/bin/trioCaller/TrioCaller --shotgun test.vcf --pedfile test.ped --states 50 --randomPhase --rounds 30 --prefix result"
exit
fi

trioCallerPath=$1
shift
#after shift, all arguments left are for trioCaller
arguments=$*
shellDir=`dirname $0`
source $shellDir/commonWrapper.sh

vcfInputFname=`findValueGivenAnOptionName --shotgun`
echo vcfInputFname: $vcfInputFname

isVCFEmpty=`checkVCFFileIfEmpty $vcfInputFname`
echo isVCFEmpty: $isVCFEmpty

outputNamePrefix=`findValueGivenAnOptionName --prefix`
echo outputNamePrefix: $outputNamePrefix
outputVCFFname=$outputNamePrefix.vcf


if test $isVCFEmpty -eq 0
then
	$trioCallerPath $arguments
	exitCode=$?
	echo "exit code is $exitCode"
	if test -r $outputVCFFname
	then
		echo "vcf output is available"
		exit $exitCode
	else
		outputEmptyVCFWithInputHeader
		exit $exitCode
	fi
else
	echo "trioCaller is not run due to empty/inexistent VCF. need to touch some files so downstream jobs wouldn't die.";
	outputEmptyVCFWithInputHeader
fi
