#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i AlignmentToTrioCallPipeline_VRC_Aln559_600_Trio620_632_648_top2Contigs.2011.12.14T1432/trioCaller/Contig1.vcf.gz 
		-O /tmp/Contig1 -o /tmp/Contig1_out.txt

Description:
	2012.2.16

	There are 4 output files, appended after outputFnamePrefix.
		.geno
			The genotype file contains 1 line per SNP.
				Each line contains 1 character per individual:
				0 means zero copies of reference allele.
				1 means one copy of reference allele.
				2 means two copies of reference allele.
				9 means missing data.
		
		.ind
			The indiv file contains 1 line per individual.  There are 3 columns:
				1st column is sample ID
				2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
				3rd column is a label which might refer to Case or Control status, or
					might be a population group label.  If this entry is set to "Ignore",
					then that individual and all genotype data from that individual will be
					removed from the data set in all convertf output.
		.snp
			The snp file contains 1 line per SNP.  There are 6 columns (last 2 optional):
				1st column is SNP name
				2nd column is chromosome.  X chromosome is encoded as 23.
					Also, Y is encoded as 24, mtDNA is encoded as 90, and XY is encoded as 91.
					Note: SNPs with illegal chromosome values, such as 0, will be removed
				3rd column is genetic position (in Morgans).  If unknown, ok to set to 0.0.
				4th column is physical position (in bases)
					Optional 5th and 6th columns are reference and variant alleles.
					For monomorphic SNPs, the variant allele can be encoded as X.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile
from AbstractVCFMapper import AbstractVCFMapper

class ConvertVCF2EigenStrat(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.update({
							('missingCallAsRefBase', 0, int):[0, '', 0, 'toggle to regard all missing calls as homozygous reference'],\
						})
	def __init__(self, **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.contig_id_pattern = re.compile(r'Contig(\d+).*')
		
		#2012.9.11 temporary to assign all missing calls to homozygous reference.
		#self.missingCallAsRefBase = 1

	def convertVCF2EigenStrat(self, inputFname, outputFnamePrefix, chromosome=None, chrLength=None, minDepth=1):
		"""
		2012.2.16
			
		"""
		sys.stderr.write("Converting %s from VCF to EigenStrat ...\n"%(inputFname))
		from pymodule.VCFFile import VCFFile
		vcfFile = VCFFile(inputFname=inputFname, minDepth=minDepth)
		
		genotype_f = open('%s.geno'%outputFnamePrefix, 'w')
		
		ind_writer = csv.writer(open('%s.ind'%outputFnamePrefix, 'w'), delimiter='\t')
		snp_writer = csv.writer(open('%s.snp'%outputFnamePrefix, 'w'), delimiter='\t')
		
		no_of_samples = 0
		for sample_id in vcfFile.sample_id_ls:
			if sample_id=='ref':
				continue
			no_of_samples += 1
			ind_writer.writerow([sample_id, 'U', 'Case'])
			
		no_of_snps = 0
		for vcfRecord in vcfFile:
			chr = vcfRecord.chr
			contig_id_pattern_result = self.contig_id_pattern.search(vcfRecord.chr)
			if contig_id_pattern_result:
				chr = contig_id_pattern_result.group(1)
			pos = vcfRecord.pos
			pos = int(pos)
			refBase = vcfRecord.refBase
			
			snp_id = '%s_%s'%(chr, pos)
			snp_writer.writerow([snp_id, chr, 0.0, pos, vcfRecord.alleleLs[0], vcfRecord.alleleLs[1]])
			no_of_snps += 1
			geno_line = ''
			for callData in vcfRecord.data_row[1:]:	#[0] is the ref base
				if callData:
					callForThisSample = callData['GT']
				else:
					if self.missingCallAsRefBase:
						callForThisSample = "NA"
					else:
						callForThisSample = '%s%s'%(refBase, refBase)
				if not callForThisSample or callForThisSample=='NA':
					#missing
					geno_line += '9'
				elif callForThisSample[0]==refBase and callForThisSample[1]==refBase:
					#homozygous reference allele
					geno_line += '2'
				elif callForThisSample[0]==callForThisSample[1] and callForThisSample[0]!=refBase:
					#homozygous alternative allele
					geno_line += '0'
				elif callForThisSample[0]!=callForThisSample[1]:
					#heterozygous
					geno_line += '1'
				else:	#missing
					geno_line += '9'
			geno_line += '\n'
			genotype_f.write(geno_line)
			
		del genotype_f, ind_writer, snp_writer
		sys.stderr.write("%s samples X %s SNPs.\n"%(no_of_samples, no_of_snps))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.convertVCF2EigenStrat(self.inputFname, self.outputFnamePrefix, chromosome=None, chrLength=None, minDepth=self.minDepth)

if __name__ == '__main__':
	main_class = ConvertVCF2EigenStrat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()