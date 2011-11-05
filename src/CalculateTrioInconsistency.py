#!/usr/bin/env python
"""
Examples:
	# all sites
	%s -i AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/gatk/Contig119.vcf.gz
		-o /tmp/Contig119.trio.75_17_86.het.homo.inconsistency -t 75,17,86
	
	# homo ony
	%s -i AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/gatk/Contig119.vcf.gz
		-o /tmp/Contig119.trio.75_17_86.homo.only.inconsistency -t 75,17,86 -m
	

Description:
	2011-9-27
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

import subprocess, cStringIO
import VervetDB, csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData
from pymodule.VCFFile import VCFFile


class CalculateTrioInconsistency(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						('trio_isq_id_ls', 1, ): ['', 't', 1, 'a comma-separated list of fa_isq_id,mo_isq_id,child_isq_id', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("refSize", 0, int): [0, '', 1, 'size of the reference used for the input VCF file. NOT used now.'],\
						("windowSize", 1, int): [200000, 'w', 1, 'calculate inconsistency within each window'],\
						('outputFnamePrefix', 1, ): [None, 'o', 1, '%s.window.%s.tsv is window-based inconsistency rate. \
							%s.summary.tsv is inconsistency rate over whole input file.'],\
						('homoOnly', 0, int):[0, 'm', 0, 'toggle to look at homozygous-only sites'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		#record it as string
		self.trio_set_str = self.trio_isq_id_ls
		if self.trio_isq_id_ls:
			self.trio_isq_id_ls = getListOutOfStr(self.trio_isq_id_ls, data_type=int)
		

	def findTrioIndex(self, sample_id2index, trio_isq_id_ls, ):
		"""
		2011-9-27
			sample_id looks like 556_16_1985088_GA_vs_524 (aln-id_isq-id_ind-code_platform_vs_ref-isq-id)
		"""
		#default to -1 (non-existent)
		father_index = -1
		mother_index = -1
		child_index = -1
		
		isq_id2list_index = {}
		for isq_id in trio_isq_id_ls:
			isq_id2list_index[isq_id] = len(isq_id2list_index)
		for sample_id, col_index in sample_id2index.iteritems():
			multi_id_ls = sample_id.split('_')
			if len(multi_id_ls)>=2:
				isq_id = int(multi_id_ls[1])
				isq_list_index = isq_id2list_index.get(isq_id)
				if isq_list_index==0:
					father_index = col_index
				elif isq_list_index==1:
					mother_index = col_index
				elif isq_list_index==2:
					child_index = col_index
		return PassingData(father_index=father_index, mother_index=mother_index, child_index=child_index)
		
		
	def outputWindowKey2Data(self, windowOutputFname, trio_set_str=None, windowKey2data=None, header=None, windowSize=200000):
		"""
		2011-9-28
		"""
		windowOutputWriter = csv.writer(open(windowOutputFname, 'w'), delimiter='\t')
		if header is None:
			header = ['trio_set', 'chromosome', 'start', 'stop', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
		
		windowOutputWriter.writerow(header)
		
		windowKeyLs = windowKey2data.keys()
		windowKeyLs.sort()
		for windowKey in windowKeyLs:
			data = windowKey2data.get(windowKey)
			chr, windowNo = windowKey[:2]
			inconsistency = data[0]/data[1]
			windowStartPos = windowNo*windowSize + 1
			windowStopPos = (windowNo+1)*windowSize
			data = [trio_set_str, chr, windowStartPos, windowStopPos, data[0], data[1], inconsistency]
			windowOutputWriter.writerow(data)
		del windowOutputWriter
		
	def outputFrequencyKey2Data(self, frequencyOutputFname, trio_set_str=None, frequencyKey2data=None):
		"""
		2011-9-28
		"""
		frequencyOutputWriter = csv.writer(open(frequencyOutputFname, 'w'), delimiter='\t')
		header = ['trio_set', 'chromosome', 'startFrequency', 'stopFrequency', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
		frequencyOutputWriter.writerow(header)
		frequencyKeyLs = frequencyKey2data.keys()
		frequencyKeyLs.sort()
		for frequencyKey in frequencyKeyLs:
			data = frequencyKey2data.get(frequencyKey)
			chr, frequencyRounded = frequencyKey[:2]
			inconsistency = data[0]/data[1]
			startFrequency = frequencyRounded
			stopFrequency = frequencyRounded + 0.1
			data = [trio_set_str, chr, startFrequency, stopFrequency, data[0], data[1], inconsistency]
			frequencyOutputWriter.writerow(data)
		del frequencyOutputWriter
		
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		vcfFile = VCFFile(inputFname=self.inputFname)
		trio_col_index_data = self.findTrioIndex(vcfFile.sample_id2index, self.trio_isq_id_ls)
		father_index = trio_col_index_data.father_index
		mother_index = trio_col_index_data.mother_index
		child_index = trio_col_index_data.child_index
		if father_index==-1 or mother_index==-1 or child_index==-1:
			sys.stderr.write("no complete trio (%s,%s,%s) found in this vcf file.\n"%(father_index, mother_index, child_index))
			sys.exit(3)
		
		windowOutputFname = "%s.window.%s.tsv"%(self.outputFnamePrefix, self.windowSize)
		frequencyOutputFname = "%s.frequency.tsv"%(self.outputFnamePrefix)
		summaryOutputFname = "%s.summary.tsv"%(self.outputFnamePrefix)
		depthOutputFname = "%s.vs.depth.tsv"%(self.outputFnamePrefix)
		
		depthOutputWriter = csv.writer(open(depthOutputFname, 'w'), delimiter='\t')
		depthOutputWriter.writerow(['depthOfFather','depthOfMother', 'depthOfChild', 'isInconsistent'])
		
		
		summaryOutputWriter = csv.writer(open(summaryOutputFname, 'w'), delimiter='\t')
		header = ['trio_set', 'chromosome', 'start', 'stop', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
		
		summaryOutputWriter.writerow(header)
		
		windowKey2data = {}	#(chr,window No.) as key, [no_of_inconsistent, no_of_total, ratio] as value
		frequencyKey2data = {}	#(chr, frequency) as key, [no_of_inconsistent, no_of_total, ratio] as value
		depth2data = {}	#depth is key, [no_of_inconsistent, no_of_total, ratio] as value
		#for summary statistic
		no_of_inconsistent = 0.
		no_of_total = 0.
		chr_set = set()
		minStart = None
		for vcfRecord in vcfFile.parseIter():
			locus = vcfRecord.locus
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			if pos>self.refSize:
				self.refSize = pos
			if minStart is None or pos<minStart:
				minStart = pos
			
			chr_set.add(chr)
			
			fa_genotypeData = vcfRecord.data_row[father_index]
			mo_genotypeData = vcfRecord.data_row[mother_index]
			child_genotypeData = vcfRecord.data_row[child_index]
			if fa_genotypeData and mo_genotypeData and child_genotypeData:
				fa_genotype = fa_genotypeData['GT']
				fa_depth = fa_genotypeData['DP']
				mo_genotype = mo_genotypeData['GT']
				mo_depth = mo_genotypeData['DP']
				child_genotype = child_genotypeData['GT']
				child_depth = child_genotypeData['DP']
				
				if self.homoOnly and ( (len(fa_genotype)>1 and fa_genotype[0]!=fa_genotype[1]) or \
						(len(mo_genotype)>1 and mo_genotype[0]!=mo_genotype[1]) or \
						(len(child_genotype)>1 and child_genotype[0]!=child_genotype[1]) ):
					#2011-9-28 ignore loci with het call in one of the trio
					continue
				potential_child_genotype_set = set()
				for fg in fa_genotype:
					for mg in mo_genotype:
						potential_child_genotype_set.add('%s%s'%(fg,mg))
				
				windowNo = int(pos/self.windowSize)
				windowKey = (chr, windowNo)
				if windowKey not in windowKey2data:
					windowKey2data[windowKey] = [0., 0., 0.] #[no_of_inconsistent, no_of_total, ratio]
				windowKey2data[windowKey][1] += 1
				
				frequency = vcfRecord.info_tag2value.get("AF", vcfRecord.info_tag2value.get("AF1", None))
				if frequency is not None:
					frequency = float(frequency)
				else:
					frequency = -0.1
				frequencyRounded = int(frequency*10)/10.0	#round it by 0.1 interval
				frequencyKey = (chr, frequencyRounded)
				if frequencyKey not in frequencyKey2data:
					frequencyKey2data[frequencyKey] = [0., 0., 0.] #[no_of_inconsistent, no_of_total, ratio]
				frequencyKey2data[frequencyKey][1] += 1
				
				no_of_total += 1
				reverse_child_genotype = '%s%s'%(child_genotype[1], child_genotype[0])	#reversing nucleotide is still same call. no phasing
				isInconsistent = 0
				if child_genotype not in potential_child_genotype_set and reverse_child_genotype not in potential_child_genotype_set:
					no_of_inconsistent += 1
					windowKey2data[windowKey][0] += 1
					frequencyKey2data[frequencyKey][0] += 1
					isInconsistent = 1
				depthOutputWriter.writerow([fa_depth, mo_depth, child_depth, isInconsistent])
		
		del depthOutputWriter
		
		inconsistency = no_of_inconsistent/no_of_total
		sys.stderr.write("inconsistent rate: %s/%s=%s \n"%(no_of_inconsistent, no_of_total, inconsistency))
		chr_ls = list(chr_set)
		summaryOutputWriter.writerow([self.trio_set_str, ','.join(chr_ls), minStart, self.refSize, \
									no_of_inconsistent, no_of_total, inconsistency])
		del summaryOutputWriter
		
		self.outputWindowKey2Data(windowOutputFname, trio_set_str=self.trio_set_str, windowKey2data=windowKey2data, \
								header=header, windowSize=self.windowSize)
		
		self.outputFrequencyKey2Data(frequencyOutputFname, trio_set_str=self.trio_set_str, frequencyKey2data=frequencyKey2data)
		

if __name__ == '__main__':
	main_class = CalculateTrioInconsistency
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
