#!/bin/sh

if test $# -lt 1 ; then
	echo "$0 INPUT_VCF_Folder"
	echo
	echo "This script counts the number of variants in each VCF file in the VCF folder and sums them all up."
	echo
	exit 1
fi
inputFolder=$1
n=0
for i in `ls $inputFolder/*vcf.gz`; do
	noOfVariants=`zcat $i |grep -v ^#|wc -l|awk -F ' ' '{print $1}'`
	n=`echo $n+$noOfVariants|bc`
done
echo $n
