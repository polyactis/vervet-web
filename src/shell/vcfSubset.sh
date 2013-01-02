#!/bin/bash
# 2012.5.9
extraArgumentsDefault="-e -f"
if test $# -lt 4
then
	echo "Usage: $0 vcfSubsetPath sampleIDFile inputVCF outputVCF [extraArguments]"
	echo
	echo "Note:"
	echo "	0. sampleIDFile is a file containing individual ucla_id in each row. one column with header UCLAID."
	echo "	1. input vcf file could be gzipped or not. tbi file (index file of the gzipped vcf) is not used."
	echo "	2. outputFile is plain VCF, no bgzip/tabix-index."
	echo "	3. Default for extraArguments is $extraArgumentsDefault."
	echo
	echo "Example:"
	echo "	$0  ~/bin/vcftools/vcf-subset 10SampleID.tsv Contig0.vcf.gz Contig0.10Sample.vcf"
exit
fi

# vcf-subset Usage: vcf-subset [OPTIONS] in.vcf.gz > out.vcf
#Options:
#   -c, --columns <string>           File or comma-separated list of columns to keep in the vcf file. If file, one column per row
#   -a, --trim-alt-alleles           Remove alternate alleles if not found in the subset
#   -e, --exclude-ref                Exclude rows not containing variants.
#   -f, --force                      Proceed anyway even if VCF does not contain some of the samples.
#   -p, --private                    Print only rows where only the subset columns carry an alternate allele.
#   -r, --replace-with-ref           Replace the excluded types with reference allele instead of dot.
#   -t, --type <list>                Comma-separated list of variant types to include: SNPs,indels.
#   -u, --keep-uncalled              Do not exclude rows without calls.
#   -h, -?, --help                   This help message.
#Examples:
#   cat in.vcf | vcf-subset -r -t indels -e -c SAMPLE1 > out.vcf

executablePath=$1
sampleIDFile=$2
inputVCF=$3
outputVCF=$4
shift
shift
shift
shift
extraArguments=$*
if [ -z $extraArguments ]; then
	extraArguments=$extraArgumentsDefault
fi
echo "extraArguments is $extraArguments."

$executablePath $extraArguments -c $sampleIDFile $inputVCF >$outputVCF

