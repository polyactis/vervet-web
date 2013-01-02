#!/bin/bash
# 2011-9-14 
if test $# -lt 4
then
	echo "Usage: $0 tabixPath inputFile outputFile region [tabixArguments]"
	echo
	echo "A shell wrapper around tabix to hide the unix pipe, required for tabix output."
	echo "Note:"
	echo "	1. tabixArguments denotes additional arguments to be passed to tabix. Default is '-h'."
	echo "	2. If more than one region, supply a BED file in place of 'region' argument and set tabixArguments to '-h -B'."
	echo
	echo "Example:"
	echo "	$0 ~/bin/tabix Contig0.vcf.gz Contig0_1_2000000.vcf.gz Contig0:1-2000000"
exit
fi


tabixPath=$1
inputFile=$2
outputFile=$3
region=$4
shift
shift
shift
shift
tabixArguments=$*
if test -z $tabixArguments
then
	tabixArguments="-h"	#add header in retrieval
fi
if test -r $inputFile
then
	echo "tabix arguments is $tabixArguments"
	$tabixPath $tabixArguments $inputFile $region > $outputFile
else
	echo "Error: $inputFile is not readable"
	exit 4
fi
