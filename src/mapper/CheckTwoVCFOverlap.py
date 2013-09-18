#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i gatk/Contig799.vcf.gz -j samtools/Contig799.vcf.gz -l 1000000 -c Contig799 -o /tmp/overlapSummary.tsv.gz
		--overlappingSitesOutputFname /tmp/overlapSite.tsv.gz
		--perSampleConcordanceOutputFname /tmp/Contig799_perSampleConcordance.tsv.gz

Description:
	2013.09.10 all three output files need to be explicitly specified on the commandline.
	2012.7.29
	2011-11-7
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, PassingData, utils
from pymodule import VCFFile, MatrixFile
from pymodule.yhio.SNP import nt2number
from AbstractVCFMapper import AbstractVCFMapper

class CheckTwoVCFOverlap(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.update({
				('jnputFname', 1, ): ['', 'j', 1, '2nd VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
				('overlappingSitesOutputFname', 0, ): ['', '', 1, 'output filename that would contain locations of all overlap sites between two VCF. no output if not given.'],\
				('perSampleConcordanceOutputFname', 0, ): ['', '', 1, 'output filename for the per-sample matching result.\
			Per-sample concordance will not be run if this is not given.', ],\
				})

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
	
	def outputOverlapSites(self, overlapping_sites_set=None, outputFname=None):
		"""
		2011-12.9
			overlapping_sites_set is a set of (chromosome, pos) tuples.
			output is tab-delimited, 3-column. Last column is always 0 to mimic output of CalculateSNPMismatchRateOfTwoVCF.py
				chromosome	position	0
		"""
		sys.stderr.write("Outputting overlap %s sites ..."%(len(overlapping_sites_set)))
		header = ['chromosome', 'position', 'random']
		overlapping_sites_list = list(overlapping_sites_set)
		writer = MatrixFile(outputFname, openMode='w', delimiter='\t')
		writer.writerow(header)
		overlapping_sites_list.sort()
		for chromosome, pos in overlapping_sites_list:
			writer.writerow([chromosome, pos, 0])
		sys.stderr.write("%s sites.\n"%(len(overlapping_sites_list)))
	
	def calculateOverlappingSites(self, vcfFile1=None, vcfFile2=None, outputFname=None, overlappingSitesOutputFname=None,\
						chromosome=None, chrLength=None):
		"""
		2013.09.10
			added argument overlappingSitesOutputFname
		2013.07.17 vcf files are no longer pre-loaded. read in locus ids first. 
		2012.8.16
		"""
		writer = MatrixFile(outputFname, openMode='w', delimiter='\t')
		header = ['#chromosome', 'length', '#sitesInInput1', '#sitesInInput2', '#overlapping', 'overlappingOverTotal', \
				'overlappingOverInput1', 'overlappingOverInput2', '#segregatingSitesNormalized', ]
		
		vcf1_locus_id_list = []
		for row in vcfFile1.reader:
			vcf1_locus_id_list.append((row[0], row[1]))
		vcf2_locus_id_list = []
		for row in vcfFile2.reader:
			vcf2_locus_id_list.append((row[0], row[1]))
		
		no_of_sites_of_input1 = len(vcf1_locus_id_list)
		no_of_sites_of_input2 = len(vcf2_locus_id_list)
		overlapping_sites_set = set(vcf1_locus_id_list)&set(vcf2_locus_id_list)
		if overlappingSitesOutputFname:
			#outputFname = "%s_overlapSitePos.tsv"%(outputFnamePrefix)
			self.outputOverlapSites(overlapping_sites_set=overlapping_sites_set, outputFname=overlappingSitesOutputFname)
		
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
							(vcfFile1.inputFname, no_of_samples, vcfFile2.inputFname, no_of_samples_in_vcf2))
		
		#exclude the ref sample in the 1st column
		if no_of_samples>1:
			normalizingConstant = float(utils.sumOfReciprocals(no_of_samples*2-1))
		else:
			normalizingConstant = 1
		noOfSegregatesSitesNormalized = no_of_overlapping_sites/(normalizingConstant*chrLength)
		
		writer.writerow(header)
		"""
		#reformat for output
		no_of_matches_per_sample_ls = map(repr, no_of_matches_per_sample_ls)
		no_of_non_NA_pairs_per_sample_ls = map(repr, no_of_non_NA_pairs_per_sample_ls)
		matchFractionLs = map(repr, matchFractionLs)
		"""
		writer.writerow([chromosome, chrLength, no_of_sites_of_input1, no_of_sites_of_input2, no_of_overlapping_sites, \
						overlapping_fraction, overlappingOverInput1, overlappingOverInput2, \
						noOfSegregatesSitesNormalized])
		del writer
		return PassingData(overlapping_sample_id_set=overlapping_sample_id_set,overlapping_sites_set=overlapping_sites_set) 
	
	
	def calculatePerSampleMismatchFraction(self, vcfFile1=None, vcfFile2=None, outputFname=None, overlapping_sample_id_set=None,\
										NA_call_encoding_set = set(['.', 'NA'])):
		"""
		2013.08.13 bugfix, derive overlapping_sites_set by itself, rather than use calculateOverlappingSites()
		2013.07.17 vcf files are no longer pre-loaded.
		2012.8.16
		"""
		sys.stderr.write("Finding matches for each sample at overlapping sites ...")
		writer = MatrixFile(outputFname, openMode='w', delimiter='\t')
		header = ['sample_id', 'no_of_matches', 'no_of_non_NA_pairs', 'matchFraction']
		no_of_samples_to_compare = len(overlapping_sample_id_set)
		
		vcfFile1._resetInput()
		vcfFile1.parseFile()
		vcfFile2._resetInput()
		vcfFile2.parseFile()
		
		overlapping_sites_set = set(vcfFile1.locus_id_ls) & set(vcfFile2.locus_id_ls)
		sys.stderr.write(" %s overlapping loci, "%(len(overlapping_sites_set)))
		
		header_ls_for_no_of_matches = []
		header_ls_for_no_of_non_NA_pairs = []
		header_ls_for_matchFraction = []
		overlapping_sample_id_list = list(overlapping_sample_id_set)
		overlapping_sample_id_list.sort()
		"""
		for sample_id in overlapping_sample_id_list:
			header_ls_for_no_of_matches.append('no_of_matches_for_%s'%(sample_id))
			header_ls_for_no_of_non_NA_pairs.append('no_of_non_NA_pairs_for_%s'%(sample_id))
			header_ls_for_matchFraction.append('matchFraction_for_%s'%(sample_id))
		
		#header = header + header_ls_for_no_of_matches + header_ls_for_no_of_non_NA_pairs + header_ls_for_matchFraction
		"""
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
				call1 = vcfFile1.genotype_call_matrix[row_index1][col_index1]
				call2 = vcfFile2.genotype_call_matrix[row_index2][col_index2]
				if call1 not in NA_call_encoding_set and call2 not in NA_call_encoding_set:
					no_of_non_NA_pairs_per_sample_ls[j] += 1
					if nt2number[call1]==nt2number[call2]:	#2013.07.03 bugfix, 'AT' and 'TA' should be same. no phase
						no_of_matches_per_sample_ls[j] += 1
					else:
						#do nothing
						pass
		matchFractionLs = [-1]*no_of_samples_to_compare
		for j in xrange(no_of_samples_to_compare):
			if no_of_non_NA_pairs_per_sample_ls[j]>0:
				matchFractionLs[j] = no_of_matches_per_sample_ls[j]/float(no_of_non_NA_pairs_per_sample_ls[j])
		
		writer.writerow(header)
		for i in xrange(no_of_samples_to_compare):
			data_row = [overlapping_sample_id_list[i], no_of_matches_per_sample_ls[i], no_of_non_NA_pairs_per_sample_ls[i],\
					matchFractionLs[i]]
			writer.writerow(data_row)
		del writer
		sys.stderr.write("%s samples.\n"%(no_of_samples_to_compare))
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		vcfFile1 = VCFFile(inputFname=self.inputFname, minDepth=self.minDepth)
		vcfFile2 = VCFFile(inputFname=self.jnputFname, minDepth=self.minDepth)
		"""
		if self.outputFnamePrefix:
			outputFnamePrefix = self.outputFnamePrefix
		elif self.outputFname:
			outputFnamePrefix = os.path.splitext(self.outputFname)[0]	#2012.8.20 bugfix, was using os.path.split()
		else:
			sys.stderr.write("could not get outputFnamePrefix from self.outputFnamePrefix %s or self.outputFname %s.\n"%\
							(self.outputFnamePrefix, self.outputFname))
			sys.exit(1)
		"""
		#overallOverlapOutputFname = '%s.tsv'%(outputFnamePrefix)
		#perSampleConcordanceOutputFname = '%s_perSample.tsv'%(outputFnamePrefix)
		
		pdata = self.calculateOverlappingSites(vcfFile1=vcfFile1, vcfFile2=vcfFile2, outputFname=self.outputFname,
							overlappingSitesOutputFname=self.overlappingSitesOutputFname, \
							chromosome=self.chromosome, chrLength=self.chrLength)
		if self.perSampleConcordanceOutputFname:
			self.calculatePerSampleMismatchFraction(vcfFile1=vcfFile1, vcfFile2=vcfFile2, \
												outputFname=self.perSampleConcordanceOutputFname,\
												overlapping_sample_id_set=pdata.overlapping_sample_id_set)
		
		

if __name__ == '__main__':
	main_class = CheckTwoVCFOverlap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()