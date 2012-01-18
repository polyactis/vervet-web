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
from pymodule import SNP
from AbstractVCFMapper import AbstractVCFMapper

class CheckTwoVCFOverlap(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.update({
						('jnputFname', 1, ): ['', 'j', 1, '2nd VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
	
	def outputOverlapSites(self, overlapping_sites_set, outputFname):
		"""
		2011-12.9
			overlapping_sites_set is a set of (chr, pos) tuples.
			output is tab-delimited, 3-column. Last column is always 0 to mimic output of CalculateSNPMismatchRateOfTwoVCF.py
				chromosome	position	0
		"""
		sys.stderr.write("Outputting overlap %s sites ..."%(len(overlapping_sites_set)))
		header = ['chromosome', 'position', 'random']
		overlapping_sites_list = list(overlapping_sites_set)
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(header)
		overlapping_sites_list.sort()
		for chr, pos in overlapping_sites_list:
			writer.writerow([chr, pos, 0])
		sys.stderr.write("%s sites.\n"%(len(overlapping_sites_list)))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		vcfFile1 = VCFFile(inputFname=self.inputFname, minDepth=self.minDepth)
		vcfFile1.parseFile()
		vcfFile2 = VCFFile(inputFname=self.jnputFname, minDepth=self.minDepth)
		vcfFile2.parseFile()
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['#chromosome', 'length', '#sitesInInput1', '#sitesInInput2', '#overlapping', 'overlappingOverTotal', \
				'overlappingOverInput1', 'overlappingOverInput2', '#segregatingSitesNormalized', ]
		
		
		no_of_sites_of_input1 = len(vcfFile1.locus_id_ls)
		no_of_sites_of_input2 = len(vcfFile2.locus_id_ls)
		overlapping_sites_set = set(vcfFile1.locus_id_ls)&set(vcfFile2.locus_id_ls)
		if self.outputFnamePrefix:
			outputFname = "%s_overlapSitePos.tsv"%(self.outputFnamePrefix)
			self.outputOverlapSites(overlapping_sites_set, outputFname)
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
		overlapping_sample_id_set = set(vcfFile1.sample_id2index.keys()) & set(vcfFile2.sample_id2index.keys())
		
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
		no_of_samples_to_compare = len(overlapping_sample_id_set)
		
		header_ls_for_no_of_matches = []
		header_ls_for_no_of_non_NA_pairs = []
		header_ls_for_matchFraction = []
		overlapping_sample_id_list = list(overlapping_sample_id_set)
		overlapping_sample_id_list.sort()
		for sample_id in overlapping_sample_id_list:
			header_ls_for_no_of_matches.append('no_of_matches_for_%s'%(sample_id))
			header_ls_for_no_of_non_NA_pairs.append('no_of_non_NA_pairs_for_%s'%(sample_id))
			header_ls_for_matchFraction.append('matchFraction_for_%s'%(sample_id))
		writer.writerow(header + header_ls_for_no_of_matches + header_ls_for_no_of_non_NA_pairs + header_ls_for_matchFraction)
		
		no_of_matches_per_sample_ls = [0]*no_of_samples_to_compare
		no_of_non_NA_pairs_per_sample_ls = [0]*no_of_samples_to_compare
		
		for locus_id in overlapping_sites_set:
			row_index1 = vcfFile1.locus_id2row_index[locus_id]
			row_index2 = vcfFile2.locus_id2row_index[locus_id]
			for j in xrange(len(overlapping_sample_id_list)):
				sample_id = overlapping_sample_id_list[j]
				col_index1 = vcfFile1.sample_id2index.get(sample_id)
				col_index2 = vcfFile2.sample_id2index.get(sample_id)
				#2012.1.17 bugfix below. so that 'AG' and 'GA' are same.
				call1 = SNP.nt2number[vcfFile1.genotype_call_matrix[row_index1][col_index1]]
				call2 = SNP.nt2number[vcfFile2.genotype_call_matrix[row_index2][col_index2]]
				if call1!=0 and call2!=0:
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