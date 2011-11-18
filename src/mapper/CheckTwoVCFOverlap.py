#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i gatk/Contig799.vcf.gz -j samtools/Contig799.vcf.gz -l 1000000 -c Contig799 -o /tmp/output

Description:
	2011-11-7
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


class CheckTwoVCFOverlap(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						('jnputFname', 1, ): ['', 'j', 1, '2nd VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("chromosome", 1, ): [None, 'c', 1, 'chromosome name for these two VCF.'],\
						("chrLength", 1, int): [0, 'l', 1, 'length of the reference used for the input VCF file.'],\
						('outputFname', 1, ): [None, 'o', 1, 'output file'],\
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
		header = ['chromosome', 'length', '#sitesInInput1', '#sitesInInput2', '#overlapping', 'overlappingOverTotal', \
				'overlappingOverInput1', 'overlappingOverInput2', '#segregatingSitesNormalized', ]
		
		
		no_of_sites_of_input1 = len(vcfFile1.locus_id_ls)
		no_of_sites_of_input2 = len(vcfFile2.locus_id_ls)
		overlapping_sites_set = set(vcfFile1.locus_id_ls)&set(vcfFile2.locus_id_ls)
		no_of_overlapping_sites = len(overlapping_sites_set)
		no_of_total_sites = no_of_sites_of_input1+no_of_sites_of_input2-no_of_overlapping_sites
		if no_of_total_sites>0:
			overlapping_fraction = no_of_overlapping_sites/float(no_of_total_sites)
		else:
			overlapping_fraction = -1
		
		if no_of_sites_of_input1>0:
			overlappingOverInput1 = no_of_overlapping_sites/float(no_of_sites_of_input1)
		else:
			overlappingOverInput1 = -1
		
		if no_of_sites_of_input2>0:
			overlappingOverInput2 = no_of_overlapping_sites/float(no_of_sites_of_input2)
		else:
			overlappingOverInput2 = -1
		
		no_of_samples = len(vcfFile1.sample_id2index)
		no_of_samples_in_vcf2 = len(vcfFile2.sample_id2index)
		if no_of_samples!=no_of_samples_in_vcf2:
			sys.stderr.write("Warning: sample size in %s is %s, in %s is %s. not matching.\n"%\
							(self.inputFname, no_of_samples, self.jnputFname, no_of_samples_in_vcf2))
		
		#exclude the ref sample in the 1st column
		if no_of_samples>1:
			normalizingConstant = float(utils.sumOfReciprocals(no_of_samples*2-1))
		else:
			normalizingConstant = 1
		noOfSegregatesSitesNormalized = no_of_overlapping_sites/(normalizingConstant*self.chrLength)
		
		sys.stderr.write("Finding matches for each sample at overlapping sites ...")
		no_of_samples_to_compare = min(no_of_samples, no_of_samples_in_vcf2)
		
		header_ls_for_no_of_matches = []
		header_ls_for_no_of_non_NA_pairs = []
		header_ls_for_matchFraction = []
		for j in xrange(no_of_samples_to_compare):
			sample_id = vcfFile1.sample_id_ls[j]
			header_ls_for_no_of_matches.append('no_of_matches_for_%s'%(sample_id))
			header_ls_for_no_of_non_NA_pairs.append('no_of_non_NA_pairs_for_%s'%(sample_id))
			header_ls_for_matchFraction.append('matchFraction_for_%s'%(sample_id))
		writer.writerow(header + header_ls_for_no_of_matches + header_ls_for_no_of_non_NA_pairs + header_ls_for_matchFraction)
		
		no_of_matches_per_sample_ls = [0]*no_of_samples_to_compare
		no_of_non_NA_pairs_per_sample_ls = [0]*no_of_samples_to_compare
		
		for locus_id in overlapping_sites_set:
			row_index1 = vcfFile1.locus_id2row_index[locus_id]
			row_index2 = vcfFile2.locus_id2row_index[locus_id]
			for j in xrange(no_of_samples_to_compare):
				call1 = vcfFile1.genotype_call_matrix[row_index1][j]
				call2 = vcfFile2.genotype_call_matrix[row_index2][j]
				if call1!='NA' and call2!='NA':
					no_of_non_NA_pairs_per_sample_ls[j] += 1
					if call1==call2:
						no_of_matches_per_sample_ls[j] += 1
					else:
						#do nothing
						pass
		matchFractionLs = [-1]*no_of_samples_to_compare
		for j in xrange(no_of_samples_to_compare):
			if no_of_non_NA_pairs_per_sample_ls[j]>0:
				matchFractionLs[j] = no_of_matches_per_sample_ls[j]/float(no_of_non_NA_pairs_per_sample_ls[j])
		sys.stderr.write("Done.\n")
		
		"""
		#reformat for output
		no_of_matches_per_sample_ls = map(repr, no_of_matches_per_sample_ls)
		no_of_non_NA_pairs_per_sample_ls = map(repr, no_of_non_NA_pairs_per_sample_ls)
		matchFractionLs = map(repr, matchFractionLs)
		"""
		writer.writerow([self.chromosome, self.chrLength, no_of_sites_of_input1, no_of_sites_of_input2, no_of_overlapping_sites, \
						overlapping_fraction, overlappingOverInput1, overlappingOverInput2, \
						noOfSegregatesSitesNormalized] + no_of_matches_per_sample_ls + \
						no_of_non_NA_pairs_per_sample_ls + matchFractionLs)

if __name__ == '__main__':
	main_class = CheckTwoVCFOverlap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()