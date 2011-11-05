#!/bin/bash
# 2011-9-14 call variants through samtools
if test $# -lt 5
then
	echo "Usage: $0 refFastaFname interval outputVCF siteType inputBAM1 inputBAM2 ..."
	echo
	echo "Note:"
	echo "	1. bamListFile contains a list of bam filenames, one file per line."
	echo "	2. refFastaFname must be faidx-indexed."
	echo "	3. shell has to be bash. PIPESTATUS is used."
	echo "	4. siteType. 1: all sites; 2: variants only"
	echo
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta Contig0:1-2000000 Contig0_1_2000000.vcf 1 1_vs_120.bam 2_vs_120.bam"
exit
fi

bwaPath=~/bin/bwa
samtoolsPath=~/bin/samtools
bcftoolsPath=~/bin/bcftools
vcfutilsPath=~/bin/vcfutils.pl

picardPath=~/script/picard/dist/

refFastaFname=$1
interval=$2
outputVCF=$3
siteType=$4
shift
shift
shift
shift
#after shifting for 3 times, all arguments left are bam files
bamFiles=$*

#-I:	do not perform indel calling
#-D:	instructs mpileup to output per-sample read depth.
#-S:	Output per-sample Phred-scaled strand bias P-value
#-q 20:	Minimum mapping quality for an alignment to be used
#-Q 20: Minimum base quality for a base to be considered
#-C 50:	Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [0]
#-d 250:	At a position, read maximally INT reads per input BAM.
#-u in -ug:	Similar to -g except that the output is uncompressed BCF, which is preferred for piping.
#-g in -ug:	Compute genotype likelihoods and output them in the binary call format (BCF)
#-r:	

#-c in "-vcg":	Call variants using Bayesian inference. This option automatically invokes option -e.
#-v in "-vcg":	Output variant sites only (force -c)
#-g in "-vcg":	Call per-sample genotypes at variant sites (force -c)
#-D100: maximum depth=100

if test "$siteType" = "1"
then
	bcftoolsArguments="-bcg"
else
	bcftoolsArguments="-bvcg"
fi

$samtoolsPath mpileup -S -D -q 30 -Q 20 -ug -r $interval -f $refFastaFname $bamFiles | $bcftoolsPath view $bcftoolsArguments - > $outputVCF.bcf

exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
echo exit code: $exitCode, $exitCode2
if test "$exitCode" = "0" && test "$exitCode2" = "0"
then
	$bcftoolsPath view $outputVCF.bcf | $vcfutilsPath varFilter -w 10 -d 3 -D5000 > $outputVCF
	exitCode=$?
	rm $outputVCF.bcf
	exit $exitCode
else
	if test "$exitCode" = "0"
	then
		exit $exitCode2
	else
		exit $exitCode
	fi
fi
