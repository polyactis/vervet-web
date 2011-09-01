#!/bin/bash

if test $# -lt 5
then
	echo "Usage: $0 RG_ID RG_PL InputFname TmpRGFname OutputFname"
	echo
	echo $# parameters given.
	echo
	echo "Given a bam file without read group (RG) information, this script add RGs to the file."
	echo "	  1. output RG info into TmpRGFname"
	echo "	  2. call samtools view and awk to add RG info and append RG_ID to the end of each line."
	echo "	  3. It doesn't replace the old RG if InputFname contains some RG already."
	
	echo "Note:
		Output goes to OutputFname.
		If InputFname is -, it i stdin.
		TmpRGFname will be created and deleted."
	echo
	echo "Examples:	"
	echo "	addRGToBAM.sh 454_vs_BAC LS454 454_1MbBAC.bam /tmp/rg.txt 454_1MbBAC.RG.bam"
	echo "	cat aethiops_vs_1MBAC.bam | addRGToBAM.sh aethiops_vs_BAC ILLUMINA - /tmp/rg.txt aethiops_vs_1MBAC.RG.bam"
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
TmpRGFname=$4
OutputFname=$5

echo -e "@RG\tID:$rg_id\tSM:$rg_sample\tLB:$platform\tPL:$platform" > $TmpRGFname

samtools view -h $inputFname | cat $TmpRGFname - | awk '{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:'$rg_id'\n",$0; }' | samtools view -bhS -o $OutputFname -

rm $TmpRGFname
