#!/usr/bin/env python
"""
Examples:
	# call one genome (nothing because no polymorphic sites)
	%s -i /tmp/vs_top150Contigs_by_aln_Contig0.bam -n 1 -o /tmp/vs_top150Contigs_by_aln_Contig0.call
	
	# call genotype from 8 genomes
	%s -n 8 -i script/vervet/data/8_genome_vs_Contig0.RG.bam -o script/vervet/data/8_genome_vs_Contig0.RG.call
	
Description:
	2011-7-12
		A multi-sample genotype caller based entirely on coverage of reads.
		It invokes VariantDiscovery.discoverHetsFromBAM() from vervet/src/misc.py
		sam/bam file has to be indexed beforehand.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions


class GenotypeCallByCoverage(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('refFastaFname', 1, ): [None, 'e', 1, 'the fasta file containing reference sequences.'],\
						('numberOfReadGroups', 1, int): [None, 'n', 1, 'number of read groups/genomes in the inputFname', ],\
						('minMinorAlleleCoverage', 1, int): [3, '', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, int): [7, '', 1, 'maximum read depth for the minor allele of a heterozygous call', ],\
						('maxNoOfReadsForGenotypingError', 1, int): [1, '', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, int): [20, '', 1, 'maximum read depth for one base to be considered'],\
						('maxMajorAlleleCoverage', 1, int): [10, '', 1, 'maximum read depth'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, int): [3, '', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	@classmethod
	def getIndividual2ColIndex(cls, header, col_name2index, sampleStartingColumn=9):
		"""
		2011-3-4
			called by discoverHetsFromVCF
		"""
		sys.stderr.write("Finding all individuals ...")
		no_of_cols = len(header)
		individual_name2col_index = {}	#individual's column name -> an opened file handler to store genetic data
		counter = 0
		for i in xrange(sampleStartingColumn, no_of_cols):
			individualName = header[i]
			col_index = col_name2index.get(individualName)
			if not individualName:	#ignore empty column
				continue
			if individualName[:-4]=='.bam':
				individualCode = individualName[:-4]	#get rid of .bam
			else:
				individualCode = individualName
			individual_name2col_index[individualCode] = col_index
			counter += 1
		sys.stderr.write("%s individuals added. Done.\n"%(counter))
		return individual_name2col_index
	
	@classmethod
	def discoverHetsFromVCF(cls, inputFname, outputFname, minMinorAlleleCoverage=4, VCFOutputType=1, maxMinorAlleleCoverage=8,\
						maxNoOfReads=30, minNoOfReads=2, \
						maxNoOfReadsForGenotypingError=1, maxMajorAlleleCoverage=30, maxNoOfReadsForAllSamples=1000,\
						nt_set = set(['a','c','g','t','A','C','G','T'])):
		"""
		2011-3-24
			add maxMinorAlleleCoverage
			Even a heterozygote's MAC is within [minMinorAlleleCoverage, maxMinorAlleleCoverage], it could still be
				a homozygous SNP.
		2011-3-4
			VCF output by GATK has a different format
			argument VCFOutputType
				1: output by samtools's vcfutils.pl
				2: output by GATK
		2011-1-6
			inputFname is VCF output by "vcfutils.pl varFilter" of samtools
		"""
		import csv
		from pymodule.utils import runLocalCommand, getColName2IndexFromHeader
		sys.stderr.write("Looking for heterozygous SNPs in %s (%s<=MAC<=%s).\n"%(os.path.basename(inputFname), \
																		minMinorAlleleCoverage, maxMinorAlleleCoverage))
		reader =csv.reader(open(inputFname), delimiter='\t')
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['sample', 'snp_id', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		moreHeader = ['GQ', 'GL', 'SB', 'QD', 'sndHighestGL', 'deltaGL']
		#['AF', 'AC','AN', 'Dels', 'HRun', 'HaplotypeScore','MQ0', 'QD']	#2011-3-4 useless
		if VCFOutputType==2:
			header += moreHeader
		writer.writerow(header)
		individual_name2col_index = None
		col_name2index = None
		counter = 0
		real_counter = 0
		for row in reader:
			if row[0] =='#CHROM':
				row[0] = 'CHROM'	#discard the #
				header = row
				col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
				individual_name2col_index = cls.getIndividual2ColIndex(header, col_name2index)
				continue
			elif row[0][0]=='#':	#2011-3-4
				continue
			chr = row[0][3:]
			pos = row[1]
			quality = row[5]
			
			outputHet= False
			
			info = row[7]
			info_ls = info.split(';')
			info_tag2value = {}
			for info in info_ls:
				try:
					tag, value = info.split('=')
				except:
					sys.stderr.write("Error in splitting %s by =.\n"%info)
					continue
				info_tag2value[tag] = value
			
			if VCFOutputType==2:	#2011-3-4
				format_column = row[col_name2index['FORMAT']]
				format_column_ls = format_column.split(':')
				format_column_name2index = getColName2IndexFromHeader(format_column_ls)
				for individual_name, individual_col_index in individual_name2col_index.iteritems():
					genotype_data = row[individual_col_index]
					genotype_data_ls = genotype_data.split(':')
					genotype_call = genotype_data_ls[format_column_name2index['GT']]
					genotype_quality_index = format_column_name2index.get('GQ')
					if genotype_quality_index is None:
						genotype_quality_index = format_column_name2index.get('DP')
					if len(genotype_data_ls)<len(format_column_name2index):
						continue
					genotype_quality = genotype_data_ls[genotype_quality_index]
					GL_index = format_column_name2index.get('GL')
					if genotype_call=='0/1':	#heterozygous
						GL_list = genotype_data_ls[GL_index]
						GL_list = GL_list.split(',')
						GL_list = map(float, GL_list)
						GL = GL_list[1]
						sndHighestGL = max([GL_list[0], GL_list[2]])
						deltaGL = GL-sndHighestGL
						AD = genotype_data_ls[format_column_name2index.get('AD')]
						AD = map(int, AD.split(','))
						minDP4 = min(AD)
						if minDP4<=maxMinorAlleleCoverage and minDP4>=minMinorAlleleCoverage:
							DP4_ratio = float(AD[0])/AD[1]
							data_row = [individual_name, 'chr%s:%s'%(chr, pos), chr, pos, quality, \
									genotype_data_ls[format_column_name2index.get('DP')], minDP4, DP4_ratio,\
									info_tag2value.get('MQ'), genotype_quality, GL,\
									info_tag2value.get('SB'), info_tag2value.get('QD'), sndHighestGL, deltaGL]
							#for i in range(3, len(moreHeader)):
							#	info_tag = moreHeader[i]
							#	data_row.append(info_tag2value.get(info_tag))
							writer.writerow(data_row)
							real_counter += 1
			elif VCFOutputType==1:
				sample_id = row[8]
				for tag in info_tag2value.keys():
					value = info_tag2value.get(tag)
					if tag=='DP4':
						tag = 'DP4_ratio'
						value = value.split(',')
						value = map(int, value)
						no_of_ref_allele = sum(value[0:2])
						no_of_non_ref_allele = sum(value[2:])
						MAC = min(no_of_ref_allele, no_of_non_ref_allele)
						if MAC<=maxMinorAlleleCoverage and MAC>=minMinorAlleleCoverage:
							outputHet = True
							value = float(no_of_ref_allele)/no_of_non_ref_allele
							info_tag2value['minDP4'] = min(no_of_ref_allele, no_of_non_ref_allele)
						else:
							value = None
						info_tag2value[tag] = value
				if outputHet:
					real_counter += 1
					output_row = [sample_id, 'chr%s:%s'%(chr, pos), chr, pos, quality, info_tag2value.get('DP'), \
								info_tag2value.get('minDP4'), info_tag2value.get('DP4_ratio'), info_tag2value.get('MQ')]
					writer.writerow(output_row)
			counter += 1
			if counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%("\x08"*80, counter, real_counter))
		del reader, writer
		sys.stderr.write("%s\t%s\t%s.\n"%("\x08"*80, counter, real_counter))
	"""
		#2011-1-6
		inputFname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.vcf'
		outputFname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=4)
		
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_p1'
		inputFname = '%s.raw.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		outputFname = '%s_min%s.raw.hets'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage)
		sys.exit(0)
		
		#2011-3-4
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		outputFname = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage, VCFOutputType=2)
		sys.exit(0)
		
		#2011-3-4
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		maxMinorAlleleCoverage=7
		outputFname = '%s_minMAC%s_maxMAC%s.hets'%(common_prefix, minMinorAlleleCoverage, maxMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage,\
			VCFOutputType=2, maxMinorAlleleCoverage=maxMinorAlleleCoverage)

		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet/data/eQTL summary.txt')
		outputType=3
		outputFname = '%s_min%s.outputType%s.tsv'%(common_prefix, minMinorAlleleCoverage, outputType)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname, outputType=outputType)
		sys.exit(0)
		
		#2011-3-24
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		#common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/454_illu_6_sub_vs_1MbBAC.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		maxMinorAlleleCoverage=7
		outputFname = '%s_minMAC%s_maxMAC%s.hets'%(common_prefix, minMinorAlleleCoverage, maxMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage,\
			VCFOutputType=1, maxMinorAlleleCoverage=maxMinorAlleleCoverage)
		sys.exit(2)
		
	"""
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		from vervet.src.misc import VariantDiscovery
		maxNoOfReadsForAllSamples = self.numberOfReadGroups*self.maxNoOfReads*self.maxNoOfReadsMultiSampleMultiplier
		VariantDiscovery.discoverHetsFromBAM(self.inputFname, self.outputFname, \
						refFastaFname=self.refFastaFname,\
						maxNoOfReads=self.maxNoOfReads, minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
						maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
						maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
						maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
						maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples)

if __name__ == '__main__':
	main_class = GenotypeCallByCoverage
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()