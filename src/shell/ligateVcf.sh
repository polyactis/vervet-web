#!/bin/bash
# 2012.5.9
if test $# -lt 3
then
	echo "Usage: $0 pathToligateVcf.pl outputVCF input1VCF input2VCF ..."
	echo
	echo "This program calls umake's ligateVcf.pl to concatenate multiple phased VCFs through the overlapping SNPs between two adjacent input VCFs. The input VCFs are in the order of the chromosome/contig."
	echo "Note:"
	echo "	1. input vcf file could be gzipped or not. tbi file (index file of the gzipped vcf) is not used."
	echo "	2. outputFile is plain VCF, no bgzip/tabix-index."
	echo
	echo "Example:"
	echo "	$0  ~/bin/umake/scripts/ligateVcf.pl Contig0.vcf Contig0_1_20000.vcf Contig0_19000_40000.vcf"
exit
fi


executablePath=$1
outputVCF=$2
shift
shift
inputVCFList=$*

echo $executablePath $inputVCFList
echo outputVCF $outputVCF
$executablePath $inputVCFList >$outputVCF
