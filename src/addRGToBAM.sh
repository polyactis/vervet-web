#!/bin/sh
#
# $0 454_vs_BAC LS454 input_prefix
# $0 aethiops_vs_BAC ILLUMINA input_prefix
#

rg_id=$1
rg_sample=$1
platform=$2
inputFnamePrefix=$3
rgFname=/tmp/rg.txt
echo -e "@RG\tID:$rg_id\tSM:$rg_sample\tLB:$platform\tPL:$platform" > $rgFname

samtools view -h $inputFnamePrefix\.bam | cat $rgFname - | awk '{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:'$rg_id'\n",$0; }' | samtools view -bhS - > $inputFnamePrefix\.RG.bam
