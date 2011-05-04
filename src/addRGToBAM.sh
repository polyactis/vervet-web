#!/bin/bash

if test $# -lt 3
then
	echo "Usage: $0 RG_ID RG_PL InputFname"
	echo
	echo $# parameters given.
	echo
	echo "Given a bam file without read group (RG) information, this script add RGs to the file."
	echo "	  1. output RG info into /tmp/rg.txt."
	echo "	  2. call samtools view and awk to add RG and RG_ID to the bam file."
	
	echo "Note:
		Output goes to stdout.
		If InputFname is -, it i stdin"
	echo
	echo "Examples:	"
	echo "	addRGToBAM.sh 454_vs_BAC LS454 454_1MbBAC.bam >454_1MbBAC.RG.bam"
	echo "	cat aethiops_vs_1MBAC.bam | addRGToBAM.sh aethiops_vs_BAC ILLUMINA - >aethiops_vs_1MBAC.RG.bam"
exit
fi


rg_id=$1
rg_sample=$1
platform=$2
if test -n $3
then
	inputFname=$3
else
	inputFname=-
fi
rgFname=/tmp/rg.txt
echo -e "@RG\tID:$rg_id\tSM:$rg_sample\tLB:$platform\tPL:$platform" > $rgFname

samtools view -h $inputFname | cat $rgFname - | awk '{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:'$rg_id'\n",$0; }' | samtools view -bhS -

rm $rgFname
