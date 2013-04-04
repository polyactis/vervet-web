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
XmxMemoryInMegabyte=$2m
jarPath=$3
inputFastq=$4
outputFnamePrefix=$5
maxNoOfReads=$6
logFilename=$7


commandline="$javaPath -Xms1280m -Xmx$XmxMemoryInMegabyte -jar $jarPath VALIDATION_STRINGENCY=LENIENT I=$inputFastq O=$outputFnamePrefix M=$maxNoOfReads"
date
echo commandline is $commandline
echo "stdout & stderr are redirected to $logFilename."
$commandline 2>&1 | tee $logFilename

exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"	#must be together in one line. PIPESTATUS[1] in subsequent lines has different meaning.
exitCode1=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`

echo "exit codes: $exitCode1, $exitCode2"

if test "$exitCode1" = "0" && test "$exitCode2" = "0"
then
	exit 0
else
	echo "Non-zero exit after running picard's SamToFastq.jar."
	exit 3
fi
date
