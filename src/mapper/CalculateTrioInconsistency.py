#!/usr/bin/env python
"""
Examples:
	# all sites
	%s -i AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/gatk/Contig119.vcf.gz
		--outputFnamePrefix /tmp/Contig119.trio.75_17_86.het.homo.inconsistency -t 552_75,554_17,555_86
	
	# homo ony
	%s -i AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/gatk/Contig119.vcf.gz
		--outputFnamePrefix /tmp/Contig119.trio.75_17_86.homo.only.inconsistency -t 552_75,554,555_86 -m
	

Description:
	2011-9-27
	two output files based on outputFnamePrefix:
		1. $outputFnamePrefix.window.$windowSize.tsv is window-based inconsistency rate. \
		2. $outputFnamePrefix.summary.tsv is inconsistency rate over whole input file.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, re, copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, NextGenSeq
from pymodule import VCFFile
from pymodule import AbstractVCFMapper

class CalculateTrioInconsistency(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVCFMapper.option_default_dict)
	option_default_dict.update({
						('trio_id', 1, ): ['', 't', 1, 'a comma-separated list of fa_id,mo_id,child_id. Each id is in the format of check IndividualAlignment.getCompositeID().\
							Only the first digit, alignment.id, is required. \
							Use 0 for missing parent.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("refSize", 0, int): [0, '', 1, 'size of the reference used for the input VCF file. NOT used now.'],\
						("windowSize", 1, int): [200000, 'w', 1, 'calculate inconsistency within each window'],\
						('homoOnly', 0, int):[0, 'm', 0, 'toggle to look at homozygous-only sites'],\
						('minDepth', 1, float): [1, '', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		#record it as string
		self.trio_set_str = self.trio_id
		if self.trio_id:
			self.trio_id = getListOutOfStr(self.trio_id, data_type=str)
		
		self.allNumberPattern = re.compile(r'^\d+')

	def findTrioIndex(self, sample_id2index, trio_id):
		"""
		2012.5.2
			bugfix. sample_id2index could contain "ref" (col_index=0) which doesn't have alignment_id embedded at all.
		2012.1.25
			trio_id is a coma-separated list of IndividualAlignment.getCompositeID().
			Use alignment.id as primary member identifier.
		2011-9-27
			sample_id in sample_id2index is something like 556_16_1985088_GA_vs_524 (aln.id_isq.id_ind.code_platform_vs_ref-isq.id)
				Its format is determined by IndividualAlignment.getReadGroup().
			
		"""
		#default to -1 (non-existent)
		father_index = -1
		mother_index = -1
		child_index = -1
		
		aln_id2list_index = {}
		for member_id in trio_id:
			aln_id = member_id.split('_')[0]
			aln_id = int(aln_id)
			aln_id2list_index[aln_id] = len(aln_id2list_index)
		
		for sample_id, col_index in sample_id2index.iteritems():
			multi_id_ls = sample_id.split('_')
			if len(multi_id_ls)>=1:
				aln_id_in_str = multi_id_ls[0]
				if self.allNumberPattern.search(aln_id_in_str):
					aln_id = int(aln_id_in_str)
				else:
					aln_id = aln_id_in_str
				trio_member_index = aln_id2list_index.get(aln_id)
				if trio_member_index==0:
					father_index = col_index
				elif trio_member_index==1:
					mother_index = col_index
				elif trio_member_index==2:
					child_index = col_index
		return PassingData(father_index=father_index, mother_index=mother_index, child_index=child_index)
		
		
	def outputWindowKey2Data(self, windowOutputFname, trio_set_str=None, windowKey2data=None, header=None, windowSize=200000):
		"""
		2011-9-28
		"""
		windowOutputWriter = csv.writer(open(windowOutputFname, 'w'), delimiter='\t')
		if header is None:
			header = ['#trio_set', 'chromosome', 'start', 'stop', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
		
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
		header = ['#trio_set', 'chromosome', 'startFrequency', 'stopFrequency', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
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
	
	def openOutputFiles(self, outputFnamePrefix, windowSize=None):
		"""
		2011-12-13
		"""
		windowOutputFname = "%s.window.%s.tsv"%(outputFnamePrefix, windowSize)
		
		frequencyOutputFname = "%s.frequency.tsv"%(outputFnamePrefix)
		
		depthOutputFname = "%s.vs.depth.tsv"%(outputFnamePrefix)
		depthOutputWriter = csv.writer(open(depthOutputFname, 'w'), delimiter='\t')
		depthOutputWriter.writerow(['#trio_set', 'chromosome', 'pos', 'depthOfFather','depthOfMother', 'depthOfChild', 'isInconsistent'])
		#2011.12.16 .filename is not supported for csv.writer
		#depthOutputWriter.filename = depthOutputFname
		
		summaryOutputFname = "%s.summary.tsv"%(outputFnamePrefix)
		summaryOutputWriter = csv.writer(open(summaryOutputFname, 'w'), delimiter='\t')
		header = ['#trio_set', 'chromosome', 'start', 'stop', 'no_of_inconsistent', 'no_of_total', 'inconsistency']
		summaryOutputWriter.writerow(header)
		#2011.12.16 .filename is not supported for csv.writer
		#summaryOutputWriter.filename = summaryOutputFname
		
		return PassingData(windowOutputFname=windowOutputFname, frequencyOutputFname=frequencyOutputFname, depthOutputWriter=depthOutputWriter, \
						summaryOutputWriter=summaryOutputWriter)
	
	def isGenotypeHet(self, genotype):
		"""
		2011.12.13
			genotype is of two letters (diploid). i.e. AG, AT,...
		"""
		if genotype and len(genotype)>1 and genotype[0]!=genotype[1]:
			return True
		else:
			return False
	
	
	def _calculateForDuo(self, vcfFile, outputDStruc=None, trio_col_index_data=None):
		"""
		2012.5.4
			make outputDStruc & trio_col_index_data an argument, rather than generate it within itself
		2011-12.15
		"""
		father_index = trio_col_index_data.father_index
		mother_index = trio_col_index_data.mother_index
		child_index = trio_col_index_data.child_index
		if father_index==-1:
			parent_index = mother_index
		else:
			parent_index =  father_index
		if parent_index ==-1 or child_index==-1:
			sys.stderr.write("eith parent or child index (%s,%s) is not found in this vcf file.\n"%(parent_index, child_index))
			sys.exit(3)
		
		depthOutputWriter = outputDStruc.depthOutputWriter
		summaryOutputWriter = outputDStruc.summaryOutputWriter
		
		windowKey2data = {}	#(chr,window No.) as key, [no_of_inconsistent, no_of_total, ratio] as value
		frequencyKey2data = {}	#(chr, frequency) as key, [no_of_inconsistent, no_of_total, ratio] as value
		depth2data = {}	#depth is key, [no_of_inconsistent, no_of_total, ratio] as value
		#for summary statistic
		no_of_inconsistent = 0.
		no_of_total = 0.
		chr_set = set()
		minStart = None
		for vcfRecord in vcfFile.parseIter():
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			if pos>self.refSize:
				self.refSize = pos
			if minStart is None or pos<minStart:
				minStart = pos
			
			chr_set.add(chr)
			
			pa_genotypeData = vcfRecord.data_row[parent_index]
			child_genotypeData = vcfRecord.data_row[child_index]
			if pa_genotypeData  and child_genotypeData:
				pa_genotype = pa_genotypeData['GT']
				pa_depth = pa_genotypeData.get('DP',None)
				child_genotype = child_genotypeData['GT']
				child_depth = child_genotypeData.get('DP',None)
				
				if self.homoOnly and ( self.isGenotypeHet(pa_genotype) or \
										self.isGenotypeHet(child_genotype) ):
					#2011-9-28 ignore loci with het call in one of the trio if self.homoOnly
					continue
				potential_child_genotype_set = set()
				for parent_allele in pa_genotype:
					for anotherAllele in 'ATCG':	#the other allele could be anything even when pa_genotype is het call. 
						potential_child_genotype_set.add('%s%s'%(parent_allele, anotherAllele))
				
				windowNo = int(pos/self.windowSize)
				windowKey = (chr, windowNo)
				if windowKey not in windowKey2data:
					windowKey2data[windowKey] = [0., 0., 0.] #[no_of_inconsistent, no_of_total, ratio]
				windowKey2data[windowKey][1] += 1
				
				frequency = vcfRecord.getAAF()
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
				depthOutputWriter.writerow([self.trio_set_str, chr, pos, pa_depth, pa_depth, child_depth, isInconsistent])
		
		del depthOutputWriter
		
		inconsistency = no_of_inconsistent/no_of_total
		sys.stderr.write("inconsistent rate: %s/%s=%s \n"%(no_of_inconsistent, no_of_total, inconsistency))
		chr_ls = list(chr_set)
		summaryOutputWriter.writerow([self.trio_set_str, ','.join(chr_ls), minStart, self.refSize, \
									no_of_inconsistent, no_of_total, inconsistency])
		del summaryOutputWriter
		
		self.outputWindowKey2Data(outputDStruc.windowOutputFname, trio_set_str=self.trio_set_str, windowKey2data=windowKey2data, \
								header=None, windowSize=self.windowSize)
		
		self.outputFrequencyKey2Data(outputDStruc.frequencyOutputFname, trio_set_str=self.trio_set_str, frequencyKey2data=frequencyKey2data)
		
	
	def _calculateForTrio(self, vcfFile, outputDStruc=None, trio_col_index_data=None):
		"""
		2012.5.4
			make outputDStruc & trio_col_index_data an argument, rather than generate it within itself
		2011.12.15
			from run()
		"""
		father_index = trio_col_index_data.father_index
		mother_index = trio_col_index_data.mother_index
		child_index = trio_col_index_data.child_index
		
		if father_index==-1 or mother_index==-1 or child_index==-1:
			sys.stderr.write("no complete trio (%s,%s,%s) found in this vcf file.\n"%(father_index, mother_index, child_index))
			sys.exit(3)
		
		depthOutputWriter = outputDStruc.depthOutputWriter
		summaryOutputWriter = outputDStruc.summaryOutputWriter
		
		windowKey2data = {}	#(chr,window No.) as key, [no_of_inconsistent, no_of_total, ratio] as value
		frequencyKey2data = {}	#(chr, frequency) as key, [no_of_inconsistent, no_of_total, ratio] as value
		depth2data = {}	#depth is key, [no_of_inconsistent, no_of_total, ratio] as value
		#for summary statistic
		no_of_inconsistent = 0.
		no_of_total = 0.
		chr_set = set()
		minStart = None
		for vcfRecord in vcfFile.parseIter():
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
				fa_depth = fa_genotypeData.get('DP',None)
				mo_genotype = mo_genotypeData['GT']
				mo_depth = mo_genotypeData.get('DP', None)
				child_genotype = child_genotypeData['GT']
				child_depth = child_genotypeData.get('DP',None)
				
				if self.homoOnly and ( self.isGenotypeHet(fa_genotype) or \
						self.isGenotypeHet(mo_genotype) or \
						self.isGenotypeHet(child_genotype) ):
					#2011-9-28 ignore loci with het call in one of the trio if self.homoOnly
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

				frequency = vcfRecord.getAAF()
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
				depthOutputWriter.writerow([self.trio_set_str, chr, pos, fa_depth, mo_depth, child_depth, isInconsistent])
		
		del depthOutputWriter
		
		inconsistency = no_of_inconsistent/no_of_total
		sys.stderr.write("inconsistent rate: %s/%s=%s \n"%(no_of_inconsistent, no_of_total, inconsistency))
		chr_ls = list(chr_set)
		summaryOutputWriter.writerow([self.trio_set_str, ','.join(chr_ls), minStart, self.refSize, \
									no_of_inconsistent, no_of_total, inconsistency])
		del summaryOutputWriter
		
		self.outputWindowKey2Data(outputDStruc.windowOutputFname, trio_set_str=self.trio_set_str, windowKey2data=windowKey2data, \
								header=None, windowSize=self.windowSize)
		
		self.outputFrequencyKey2Data(outputDStruc.frequencyOutputFname, trio_set_str=self.trio_set_str, frequencyKey2data=frequencyKey2data)
		


	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		if NextGenSeq.isVCFFileEmpty(self.inputFname, checkContent=True):
			sys.stderr.write("Input %s doesn't exist or no variants in it.\n"%(self.inputFname))
			#make sure some output files will exist for downstream jobs.
			self.openOutputFiles(self.outputFnamePrefix, self.windowSize)
			sys.exit(0)
		
		vcfFile = VCFFile(inputFname=self.inputFname, minDepth=self.minDepth)
		trio_col_index_data = self.findTrioIndex(vcfFile.sample_id2index, self.trio_id)
		father_index = trio_col_index_data.father_index
		mother_index = trio_col_index_data.mother_index
		child_index = trio_col_index_data.child_index
		outputDStruc = self.openOutputFiles(self.outputFnamePrefix, self.windowSize)
		if (father_index==-1 and mother_index!=-1) or (father_index!=-1 and mother_index==-1):
			#one parent is missing. it's duo.
			self._calculateForDuo(vcfFile, outputDStruc=outputDStruc, trio_col_index_data=trio_col_index_data)
		else:
			self._calculateForTrio(vcfFile, outputDStruc=outputDStruc, trio_col_index_data=trio_col_index_data)
		
		"""
		try:
			pass
		except ValueError:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			#make sure some output files will exist for downstream jobs.
			self.openOutputFiles(self.outputFnamePrefix, self.windowSize)
		"""

if __name__ == '__main__':
	main_class = CalculateTrioInconsistency
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
