#!/bin/sh
#PBS -q cmb -j oe -S /bin/bash
#PBS -l walltime=100:00:00,mem=5G
#PBS -d /home/cmb-01/yuhuang/qjob_output
#PBS -k eo
#PBS -l nodes=1:myri:ppn=1

#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_NAME.joblog.$JOB_ID
#$ -j y
#$ -l h_rt=80:00:00
#$ -l highp
#$ -r y
#$ -pe shared* 8
#$ -V
source ~/.bash_profile
# shared* 8 is used to lock down the whole machine to prevent memory abuse by other jobs.
# also MarkDuplicatesJar use multi-threads.
# remove "-l h_data=5G" because no suitable queues were found.
inputFname=$1
inputFnamePrefix=`echo $inputFname|awk -F . '{print $1}'`
MarkDupOutputF=$inputFnamePrefix\_dupMarked.bam
MarkDupOutputMetricF=$inputFnamePrefix\_dupMarked.metric
MarkDuplicatesJar=~/script/picard/dist/MarkDuplicates.jar
javaPath=~/bin/jdk/bin/java
tmpDir=/u/scratch/polyacti/
$javaPath -Xms128m -Xmx5000m -jar $MarkDuplicatesJar MAX_FILE_HANDLES=200 VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true INPUT=$inputFname OUTPUT=$MarkDupOutputF M=$MarkDupOutputMetricF MAX_RECORDS_IN_RAM=500000 TMP_DIR=$tmpDir

BuildBamIndexFilesJar=~/script/picard/dist/BuildBamIndex.jar
bai_output=$MarkDupOutputF.bai
$javaPath -Xms128m -Xmx2500m -jar $BuildBamIndexFilesJar VALIDATION_STRINGENCY=LENIENT INPUT=$MarkDupOutputF OUTPUT=$bai_output
