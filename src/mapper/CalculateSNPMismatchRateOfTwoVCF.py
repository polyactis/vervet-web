#!/usr/bin/env python
"""
Examples:
	%s 
	
	#2011-11-11 tolerate zero mismatch
	%s -i gatk/Contig799.vcf.gz -j samtools/Contig799.vcf.gz -m 0 -o /tmp/gatk_vs_samtools_contig0.tsv

Description:
	2011-11-7
		program look at the intersection SNP sites of two VCF
			calculate the mismatch rate between two VCF
			output sites whose mismatch rate is >=minMismatchRate
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile


class CalculateSNPMismatchRateOfTwoVCF(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						('jnputFname', 1, ): ['', 'j', 1, '2nd VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("maxMismatchRate", 0, float): [0, 'm', 1, 'maximum mismatch rate'],\
						('outputFname', 1, ): [None, 'o', 1, 'output file, tab-delimited: chr, position of SNPs whose mimatch rate is <=maxMismatchRate'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		vcfFile1 = VCFFile(inputFname=self.inputFname)
		vcfFile1.parseFile()
		vcfFile2 = VCFFile(inputFname=self.jnputFname)
		vcfFile2.parseFile()
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['chromosome', 'position']
		
		
		no_of_sites_of_input1 = len(vcfFile1.locus_id_ls)
		no_of_sites_of_input2 = len(vcfFile2.locus_id_ls)
		overlapping_sites_set = set(vcfFile1.locus_id_ls)&set(vcfFile2.locus_id_ls)
		no_of_overlapping_sites = len(overlapping_sites_set)
		no_of_total_sites = no_of_sites_of_input1+no_of_sites_of_input2-no_of_overlapping_sites
		
		no_of_samples = len(vcfFile1.sample_id2index)
		no_of_samples_in_vcf2 = len(vcfFile2.sample_id2index)
		if no_of_samples!=no_of_samples_in_vcf2:
			sys.stderr.write("Warning: sample size in %s is %s, in %s is %s. not matching.\n"%\
							(self.inputFname, no_of_samples, self.jnputFname, no_of_samples_in_vcf2))
		
		no_of_samples_to_compare = min(no_of_samples, no_of_samples_in_vcf2)
		
		writer.writerow(header)
		
		locus_id2mismatchData = {}
		for locus_id in overlapping_sites_set:
			row_index1 = vcfFile1.locus_id2row_index[locus_id]
			row_index2 = vcfFile2.locus_id2row_index[locus_id]
			no_of_mismatches = 0
			no_of_non_NA_pairs = 0.0
			for j in xrange(no_of_samples_to_compare):
				call1 = vcfFile1.genotype_call_matrix[row_index1][j]
				call2 = vcfFile2.genotype_call_matrix[row_index2][j]
				if call1!='NA' and call2!='NA':
					no_of_non_NA_pairs += 1
					if call1!=call2:
						no_of_mismatches += 1
					else:
						#do nothing
						pass
			if no_of_non_NA_pairs>0:
				mismatchRate = no_of_mismatches/float(no_of_non_NA_pairs)
			else:
				mismatchRate = -1
			locus_id2mismatchData[locus_id] = [mismatchRate, no_of_mismatches, no_of_non_NA_pairs]
		
		counter = 0
		locus_id_ls = locus_id2mismatchData.keys()
		locus_id_ls.sort()
		for locus_id in locus_id_ls:
			mismatchData = locus_id2mismatchData.get(locus_id)
			mismatchRate = mismatchData[0]
			if mismatchRate<=self.maxMismatchRate:
				counter += 1
				chr, pos = locus_id[:2]
				writer.writerow([chr, pos, mismatchRate])
		sys.stderr.write("%s loci passed the maxMismatchRate out of %s overlapped loci.\n"%(counter, len(overlapping_sites_set)))

if __name__ == '__main__':
	main_class = CalculateSNPMismatchRateOfTwoVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()