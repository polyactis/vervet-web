#!/bin/bash
# 2011-9-14 call variants through samtools
maxDPDefault=5000
if test $# -lt 5
then
	echo "Usage: $0 refFastaFname mpileupInterval bcftoolsInterval outputVCF siteType maxDP inputBAM1 inputBAM2 ..."
	echo
	echo "Note:"
	echo "	1. interval could be a string or BED file restricting samtools."
	echo "	2. refFastaFname must be faidx-indexed."
	echo "	3. shell has to be bash. PIPESTATUS is used."
	echo "	4. siteType. 1: all sites; 2: variants only"
	echo "	5. maxDP is the maximum depth a site could be of, should be set according to alignment depth."
	echo
	echo "Example:"
	echo "	$0 120_480K_supercontigs.fasta Contig0:1-2000000 Contig0_1_2000000.vcf 1 2000 1_vs_120.bam 2_vs_120.bam"
exit
fi
source $HOME/.bash_profile
bwaPath=$HOME/bin/bwa
samtoolsPath=$HOME/bin/samtools
bcftoolsPath=$HOME/bin/bcftools
vcfutilsPath=$HOME/bin/vcfutils.pl

picardPath=$HOME/script/picard/dist/

refFastaFname=$1
mpileupInterval=$2
bcftoolsInterval=$3
outputVCF=$4
siteType=$5
maxDP=$6	#2012.8.6
shift
shift
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
###### aruments for $vcfutilsPath
#-D100: maximum depth=100
#-w 10:	SNP within 10 bp around a gap to be filtered.
#-d 3:	minimum read depth is 3.
#-e 0:	min P-value for HWE (plus F<0). setting it to zero removes this filter.

if test "$siteType" = "1"
then
	bcftoolsArguments="-bcg"
else
	bcftoolsArguments="-bvcg"
fi


if test -r "$mpileupInterval"
then
	#it's a file, BED format.
	echo "$mpileupInterval is treated as a BED file."
	mpileupIntervalArgument="-l $mpileupInterval"
else
	echo "$mpileupInterval" is treated as a region-defining string.
	mpileupIntervalArgument="-r $mpileupInterval"
fi
if test -r "$bcftoolsInterval"
then
	#it's a file, BED format.
	echo "$bcftoolsInterval is treated as a BED file."
	bcftoolsIntervalArgument="-l $bcftoolsInterval"
else
	echo "$bcftoolsInterval is treated as a region-defining string and ignored."
	bcftoolsIntervalArgument=""
fi

#2011-11-04 from Vasily. lower threshold for sixth column (QUAL) in VCF.
low_quality_thresh=5

# 2012.5.4 "-S -D" from mpileup increases running time significantly
$samtoolsPath mpileup -S -D -q 30 -Q 20 -ug $mpileupIntervalArgument -f $refFastaFname $bamFiles | $bcftoolsPath view $bcftoolsArguments $bcftoolsIntervalArgument - > $outputVCF.bcf
exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
echo exit code: $exitCode, $exitCode2

indelSNPVCF=$outputVCF.indel_snp.vcf
indelVCF=$outputVCF.indel.vcf

if test "$exitCode" = "0" && test "$exitCode2" = "0"
then
	$bcftoolsPath view $outputVCF.bcf | $vcfutilsPath varFilter -w 10 -d 3 -D $maxDP -e 0 -2 0 > $indelSNPVCF
	exitCode=$?
	rm $outputVCF.bcf
	if test "$exitCode" = "0"
	then
		#split into INDEL and SNP-only VCF
		egrep "^#" $indelSNPVCF 1>$outputVCF
		egrep -v "^#" $indelSNPVCF | egrep -v INDEL |awk "{if (\$6>=$low_quality_thresh) print}">>$outputVCF
		exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]} ${PIPESTATUS[2]}"
		exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
		exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
		exitCode3=`echo $exitCodeAll|awk -F ' ' '{print $3}'`
		echo exit code: $exitCodeAll
		#if test "$exitCode" != "0" || test "$exitCode2" != "0" || test "$exitCode3" != "0"
		#then
		#	#non-zero exit if any non-zero exit happens along the pipe
		#	#sometimes, the $indelSNPVCF is empty (only vcf headers) and egrep would exit non-zero because no line matched the pattern.	#so don't do this for now
		#	exit 3;
		#fi
		
		egrep "^#" $indelSNPVCF 1>$indelVCF
		egrep -v "^#" $indelSNPVCF | egrep INDEL  1>>$indelVCF
		exitCodeAll="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
		exitCode=`echo $exitCodeAll|awk -F ' ' '{print $1}'`
		exitCode2=`echo $exitCodeAll|awk -F ' ' '{print $2}'`
		echo exit code: $exitCodeAll
		#if test "$exitCode" != "0" || test "$exitCode2" != "0"
		#then
		#	#non-zero exit if any non-zero exit happens along the pipe
		#	exit 0;	#should be non-zero, but not care about indels right now
		#fi
		
		rm $indelSNPVCF
	else
		exit $exitCode
	fi
else
	if test "$exitCode" = "0"
	then
		exit $exitCode2
	else
		exit $exitCode
	fi
fi
