#!/bin/sh
#PBS -q cmb -j oe -S /bin/bash
#PBS -l walltime=23:55:00,mem=2G
#PBS -d /home/cmb-03/mn/yuhuang/qjob_output
#PBS -k eo
#PBS -l nodes=1:myri:ppn=1
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_NAME.joblog.$JOB_ID
#$ -j y
#$ -l time=23:55:00
#$ -l h_data=4G
#$ -r y
#$ -V
source ~/.bash_profile

inputFname=$1
inputFnamePrefix=`echo $inputFname|awk -F . '{print $1}'`
DupRemovedOutputF=$inputFnamePrefix\_dupRemoved.bam
javaPath=~/bin/jdk/bin/java
samtoolsPath=~/bin/samtools


BuildBamIndexFilesJar=~/script/picard/dist/BuildBamIndex.jar
bai_output=$DupRemovedOutputF.bai
$samtoolsPath rmdup $inputFname $DupRemovedOutputF
$javaPath -Xms128m -Xmx2500m -jar $BuildBamIndexFilesJar VALIDATION_STRINGENCY=LENIENT INPUT=$DupRemovedOutputF OUTPUT=$bai_output
