#!/usr/bin/env python
"""
Examples:
	%s -i input.vcf -o -o selectedStKitts.vcf --country_id_ls 144 --tax_id_ls 60711
	
	%s -i input.vcf -o -o selectedNevis.vcf --country_id_ls 148 --tax_id_ls 60711

	%s -i input.vcf.gz -o selected.vcf --country_id_ls 135,136,144,148,151 --tax_id_ls 60711
		

Description:
	2012.10.5 program that extracts samples from a VCF and form a new VCF.
		need to re-calculate the AC/AF values of each variant.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions
from pymodule import VCFFile
from pymodule.utils import returnZeroFunc
#used in getattr(individual_site_id_set, '__len__', returnZeroFunc)()
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper

class ExtractSamplesFromVCF(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict
	option_default_dict.update({
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						("country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						("tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						("min_coverage", 0, int): [None, '', 1, 'minimum coverage for a sample to be extracted'],\
						("max_coverage", 0, int): [None, '', 1, 'maximum coverage for a sample to be extracted'],\
						("is_contaminated", 0, int): [None, '', 1, 'if not None, require individual_alignment.individual_sequence.is_contaminated be same value'],\
						("outputFormat", 1, int): [1, '', 1, 'output format. 1: a subset VCF file; \n\
	2: file with a list of sample IDs (one per line), with header; \n\
	3: file with a list of sample IDs (one per line), without header'],\
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
		listArgumentName_data_type_ls = [("site_id_ls", int), ('country_id_ls', int), ('tax_id_ls', int)]
		listArgumentName2hasContent = ProcessOptions.processListArguments(listArgumentName_data_type_ls,\
												emptyContent=[], class_to_have_attr=self)
	
	
	def extractSamples(self, db_vervet=None, inputFname=None, outputFname=None, \
					tax_id_set=None, site_id_set=None, country_id_set=None, \
					min_coverage=None, max_coverage=None, outputFormat=1, is_contaminated=None,\
					**keywords):
		"""
		2013.07.03 added argument is_contaminated (whether to fetch contaminated samples or not)
		2013.04.30 added argument min_coverage, max_coverage
		2012.10.10
			added argument outputFormat. 
		2012.10.5
			
		"""
		sys.stderr.write("Extracting samples from %s, %s sites & %s countries & %s taxonomies, min_coverage=%s, max_coverage=%s, outputFormat=%s, is_contaminated=%s ...\n"%\
							(inputFname,\
							getattr(site_id_set, '__len__', returnZeroFunc)(),\
							getattr(country_id_set, '__len__', returnZeroFunc)(),\
							getattr(tax_id_set, '__len__', returnZeroFunc)(), min_coverage, max_coverage,\
							outputFormat, is_contaminated ))
		vcfFile = VCFFile(inputFname=inputFname)
		
		oldHeader = vcfFile.header
		oldHeaderLength = len(oldHeader)
		newHeader = oldHeader[:vcfFile.sampleStartingColumn]	#anything before the samples are same
		no_of_samples = 0
		col_index2sampleID = {}	#this structure stores the selected samples and their column index 
		for col_index, individual_name in vcfFile.get_col_index_individual_name_ls():
			individualAlignment = db_vervet.parseAlignmentReadGroup(individual_name).individualAlignment
			if individualAlignment is not None:
				filteredAlignmentList = db_vervet.filterAlignments(alignmentLs=[individualAlignment], min_coverage=min_coverage, \
						max_coverage=max_coverage, individual_site_id=None, \
						sequence_filtered=None, individual_site_id_set=site_id_set, \
						mask_genotype_method_id=None, parent_individual_alignment_id=None,\
						country_id_set=country_id_set, tax_id_set=tax_id_set, excludeContaminant=False, \
						is_contaminated=is_contaminated, excludeTissueIDSet=None,\
						local_realigned=None, reduce_reads=None, report=False)
				if filteredAlignmentList:	#non-empty, passed the filter
					newHeader.append(individual_name)
					no_of_samples += 1
					col_index2sampleID[col_index] = individual_name
			else:
				sys.stderr.write("Warning: no individualAlignment for sample %s.\n"%(individual_name))
				sys.exit(3)
		
		no_of_snps = 0
		if outputFormat==1:
			outVCFFile = VCFFile(outputFname=outputFname)
			outVCFFile.metaInfoLs = vcfFile.metaInfoLs
			outVCFFile.header = newHeader
			outVCFFile.writeMetaAndHeader()
			
			newHeaderLength = len(newHeader)
			for vcfRecord in vcfFile:
				data_row =vcfRecord.row[:vcfFile.sampleStartingColumn]
				for i in xrange(vcfFile.sampleStartingColumn, oldHeaderLength):
					if i in col_index2sampleID:
						data_row.append(vcfRecord.row[i])
				outVCFFile.writer.writerow(data_row)
				no_of_snps += 1
			outVCFFile.close()
		elif outputFormat in [2,3]:
			outf = open(outputFname, 'w')
			if outputFormat==2:
				outf.write("sampleID\n")
			for col_index, sampleID in col_index2sampleID.iteritems():
				outf.write("%s\n"%(sampleID))
			outf.close()
		vcfFile.close()
		sys.stderr.write("%s samples X %s SNPs.\n"%(no_of_samples, no_of_snps))
		
	
	def run(self):
		"""
		2012.10.5
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.extractSamples(db_vervet=self.db_vervet, inputFname=self.inputFname, outputFname=self.outputFname, \
					tax_id_set=set(self.tax_id_ls), site_id_set=set(self.site_id_ls), country_id_set=set(self.country_id_ls),\
					min_coverage=self.min_coverage, max_coverage=self.max_coverage, outputFormat=self.outputFormat,
					is_contaminated=self.is_contaminated)
		

if __name__ == '__main__':
	main_class = ExtractSamplesFromVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()