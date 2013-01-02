#!/bin/bash
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

javaPath=~/bin/jdk/bin/java
javaPath=$1
XmxMemoryInMegabyte=$2
jarPath=$3
inputFastq=$4
outputFnamePrefix=$5
maxNoOfReads=$6
logFilename=$7


$javaPath -Xms128m -Xmx$XmxMemoryInMegabyte\m -jar $jarPath VALIDATION_STRINGENCY=LENIENT I=$inputFastq O=$outputFnamePrefix M=$maxNoOfReads >& $logFilename
