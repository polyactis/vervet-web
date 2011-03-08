#!/bin/bash

if test $# -lt 3
then
	echo "Usage: $0 RG_ID RG_PL InputFnamePrefix"
	echo
	echo "Given a bam file without read group (RG) information, this script add RGs to the file."
	echo "	  1. output RG info into /tmp/rg.txt."
	echo "	  2. call samtools view and awk to add RG and RG_ID to the bam file."
	
	echo "Note:
		InputFnamePrefix is part of the filename without .bam."
	echo
	echo "Examples:	addRGToBAM.sh 454_vs_BAC LS454 input_prefix"
	echo "	addRGToBAM.sh aethiops_vs_BAC ILLUMINA input_prefix"
exit
fi


rg_id=$1
rg_sample=$1
platform=$2
inputFnamePrefix=$3
rgFname=/tmp/rg.txt
echo -e "@RG\tID:$rg_id\tSM:$rg_sample\tLB:$platform\tPL:$platform" > $rgFname

samtools view -h $inputFnamePrefix\.bam | cat $rgFname - | awk '{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:'$rg_id'\n",$0; }' | samtools view -bhS - > $inputFnamePrefix\.RG.bam
