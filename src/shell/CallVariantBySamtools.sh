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
#################
#crocea@crocea:~$ samtools mpileup
#
#Usage: samtools mpileup [options] in1.bam [in2.bam [...]]
#
#Input options:
#
#       -6           assume the quality is in the Illumina-1.3+ encoding
#       -A           count anomalous read pairs
#       -B           disable BAQ computation
#       -b FILE      list of input BAM files [null]
#       -C INT       parameter for adjusting mapQ; 0 to disable [0]
#       -d INT       max per-BAM depth to avoid excessive memory usage [250]
#       -f FILE      faidx indexed reference sequence file [null]
#       -G FILE      exclude read groups listed in FILE [null]
#       -l FILE      list of positions (chr pos) or regions (BED) [null]
#       -M INT       cap mapping quality at INT [60]
#       -r STR       region in which pileup is generated [null]
#       -R           ignore RG tags
#       -q INT       skip alignments with mapQ smaller than INT [0]
#       -Q INT       skip bases with baseQ/BAQ smaller than INT [13]
#
#Output options:
#
#       -D           output per-sample DP in BCF (require -g/-u)
#       -g           generate BCF output (genotype likelihoods)
#       -O           output base positions on reads (disabled by -g/-u)
#       -s           output mapping quality (disabled by -g/-u)
#       -S           output per-sample strand bias P-value in BCF (require -g/-u)
#       -u           generate uncompress BCF output
#
#SNP/INDEL genotype likelihoods options (effective with `-g' or `-u'):
#
#       -e INT       Phred-scaled gap extension seq error probability [20]
#       -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
#       -h INT       coefficient for homopolymer errors [100]
#       -I           do not perform indel calling
#       -L INT       max per-sample depth for INDEL calling [250]
#       -m INT       minimum gapped reads for indel candidates [1]
#       -o INT       Phred-scaled gap open sequencing error probability [40]
#       -P STR       comma separated list of platforms for indels [all]
#
#Notes: Assuming diploid individuals.


###### aruments for $vcfutilsPath varFilter
#Usage:   vcfutils.pl varFilter [options] <in.vcf>
#
#Options: -Q INT    minimum RMS mapping quality for SNPs [10]
#         -d INT    minimum read depth [2]
#         -D INT    maximum read depth [10000000]
#         -a INT    minimum number of alternate bases [2]
#         -w INT    SNP within INT bp around a gap to be filtered [3]
#         -W INT    window size for filtering adjacent gaps [10]
#         -1 FLOAT  min P-value for strand bias (given PV4) [0.0001]
#         -2 FLOAT  min P-value for baseQ bias [1e-100]
#         -3 FLOAT  min P-value for mapQ bias [0]
#         -4 FLOAT  min P-value for end distance bias [0.0001]
#		 -e FLOAT  min P-value for HWE (plus F<0) [0.0001]
#         -p        print filtered variants
#
#Note: Some of the filters rely on annotations generated by SAMtools/BCFtools.
#
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
