#!/usr/bin/env python
"""
Examples:
	# enter debug mode
	misc.py -b
	
	# change the db connection setting 
	misc.py -z localhost -u yh
	
Description:
	It's a file with all sorts of classes that haven't gone on to be independent.
	
	All those classes are manually chosen to be run in Main.run().
"""
import os, sys, numpy

#2007-03-05 common codes to initiate database connection
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	#sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	#sys.path.insert(0, os.path.expanduser('~/lib/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment


class VariantDiscovery(object):
	"""
	2011-1-6
	"""
	def __init__(self):
		pass
	
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
	def discoverHetsFromVCF(cls, input_fname, output_fname, minNumberOfReads=4, VCFOutputType=1):
		"""
		2011-3-4
			VCF output by GATK has a different format
			argument VCFOutputType
				1: output by samtools's vcfutils.pl
				2: output by GATK
		2011-1-6
			input_fname is VCF output by "vcfutils.pl varFilter" of samtools
		"""
		import csv
		from pymodule.utils import runLocalCommand, getColName2IndexFromHeader
		
		reader =csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['sample', 'snp_id', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		moreHeader = ['GQ', 'GL', 'SB']
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
					genotype_quality = genotype_data_ls[genotype_quality_index]
					GL_index = format_column_name2index.get('GL')
					if genotype_call=='0/1':	#heterozygous
						GL = genotype_data_ls[GL_index]
						GL = GL.split(',')
						GL = GL[1]
						
						AD = genotype_data_ls[format_column_name2index.get('AD')]
						AD = map(int, AD.split(','))
						if AD[0]>=minNumberOfReads and AD[1]>=minNumberOfReads:	#satisfy read support
							minDP4 = min(AD)
							DP4_ratio = float(AD[0])/AD[1]
							data_row = [individual_name, 'chr%s:%s'%(chr, pos), chr, pos, quality, \
									genotype_data_ls[format_column_name2index.get('DP')], minDP4, DP4_ratio,\
									info_tag2value.get('MQ'), genotype_quality, GL,\
									info_tag2value.get('SB')]
							for i in range(3, len(moreHeader)):
								info_tag = moreHeader[i]
								data_row.append(info_tag2value.get(info_tag))
							writer.writerow(data_row)
							real_counter += 1
			elif VCFOutputType==1:
				sample_id = row[8]
				for tag, value in info_tag2value.iteritems():
					if tag=='DP4':
						tag = 'DP4_ratio'
						value = value.split(',')
						value = map(int, value)
						no_of_ref_allele = sum(value[0:2])
						no_of_non_ref_allele = sum(value[2:])
						if no_of_ref_allele>=minNumberOfReads and no_of_non_ref_allele>=minNumberOfReads:
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
		input_fname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.vcf'
		output_fname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		VariantDiscovery.discoverHetsFromVCF(input_fname, output_fname, minNumberOfReads=4)
		
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_p1'
		input_fname = '%s.raw.vcf'%(common_prefix)
		minNumberOfReads=2
		output_fname = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHetsFromVCF(input_fname, output_fname, minNumberOfReads=minNumberOfReads)
		sys.exit(0)
		
		#2011-3-4
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		input_fname = '%s.vcf'%(common_prefix)
		minNumberOfReads=3
		output_fname = '%s_min%s.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHetsFromVCF(input_fname, output_fname, minNumberOfReads=minNumberOfReads, VCFOutputType=2)
		sys.exit(0)
	"""
	
	@classmethod
	def checkOverlapping(cls, myVariantFile, jessicaVariantFname, output_fname):
		"""
		2011-1-6
		"""
		import csv
		sys.stderr.write("Reading in Jessica's variants ... ")
		reader = csv.reader(open(jessicaVariantFname), delimiter='\t')
		reader.next()	#skip header
		snp_id2work_status = {}	
		chr_pos2data_ls = {}	#map chr_pos to [snp_id, work_status]
		for row in reader:
			if row[3]=='':	# the genotyping work/failure section. it comes first
				snp_id = row[0]
				genotyped = int(row[1])
				status = int(row[2])
				if genotyped==1:
					snp_id2work_status[snp_id] = status
				
			else:	# the snp_id, chromosome, position part. 2nd part. at this stage, snp_id2work_status should be ready.
				snp_id = row[1]
				chr_pos = row[2]
				if snp_id in snp_id2work_status:
					chr_pos2data_ls[chr_pos] = [snp_id, snp_id2work_status.get(snp_id)]
		del reader
		sys.stderr.write("%s genotyped jessica variants. Done.\n"%(len(chr_pos2data_ls)))
		
		sys.stderr.write("Reading in my variants ... ")
		reader = csv.reader(open(myVariantFile), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['snp_id', 'chr_pos|string', 'status', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		moreHeader = ['GQ', 'GL', 'SB']
		#['AF', 'AC','AN', 'Dels', 'HRun', 'HaplotypeScore','MQ0', 'QD']
		header += moreHeader
		writer.writerow(header)
		
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter = 0
		for row in reader:
			chr_pos = row[col_name2index['snp_id']]
			pos = row[col_name2index['pos']]
			if chr_pos in chr_pos2data_ls:
				snp_id, work_status = chr_pos2data_ls[chr_pos]
				data_row = [snp_id, chr_pos, work_status] + row[2:]
				writer.writerow(data_row)
				real_counter += 1
			
			counter += 1
		del writer
		sys.stderr.write("%s/%s of mine overlapped with jessica's. Done.\n"%(real_counter, counter))
	"""
	
		#2011-1-6
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_ANXA_p1'
		input_fname = '%s.raw.vcf'%(common_prefix)
		minNumberOfReads=3
		output_fname = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHetsFromVCF(input_fname, output_fname, minNumberOfReads=minNumberOfReads)
		sys.exit(0)
		
		#2011-1-6
		myVariantFile = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		myVariantFile = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		output_fname = '%s_min%s.overlap.tsv'%(common_prefix, minNumberOfReads)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, output_fname)
		sys.exit(0)
		
	"""
	
	class DrawHistogramOfChosenBWAOutputScore(object):
		"""
		2011-2-8
			
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include scoreType, plotType, mapq_ls, score_ls, C_ls
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-2-8
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			param_obj.real_counter += 1
			
			score = None
			if self.scoreType==1:
				score = read.mapq
			elif self.scoreType in [2,3]:
				if self.scoreType==2:
					no_of_bases = read.alen
				elif self.scoreType==3:
					no_of_bases = read.rlen
				for tag in read.tags:
					if tag[0]=='AS':
						score = tag[1]	#'AS'
						score = score/float(no_of_bases) 	#divide the alignment score by the read length
						break
			if score is not None:
				self.mapq_ls.append(read.mapq)
				self.score_ls.append(score)
				if self.plotType==2 or self.plotType==1:
					self.C_ls.append(1)
				elif self.plotType==3:
					self.C_ls.append(no_of_bases)	#read.rlen is different from read.alen
			
	
	@classmethod
	def drawHistogramOfChosenBWAOutputScore(cls, inputFname, outputFname, scoreType=1, plotType=2):
		"""
		2011-2-1
			scoreType
				1. mapping quality
				2. alignment score per aligned read length
				3. alignment score per read length
			plotType
				1: 1D histogram
				2: 2D histogram (scoreType has to be =2, color according to how many reads)
				3: 2D histogram (color according to median read length)
		"""
		import os,sys
		import pysam
		
		sys.stderr.write("Draw histogram of (scoreType=%s, plotType=%s) from %s ...\n"%\
						(scoreType, plotType, inputFname))
		mapq_ls = []
		score_ls = []
		C_ls = []
		processor = cls.DrawHistogramOfChosenBWAOutputScore(scoreType=scoreType, \
									plotType=plotType, mapq_ls=mapq_ls, C_ls=C_ls, score_ls=score_ls)
		samfile = pysam.Samfile(inputFname, "rb" )
		cls.traverseBamFile(samfile, processor=processor)
		
		sys.stderr.write("Done.\n")
		if scoreType==1:
			xlabel_1D = "read map quality"
		elif scoreType==2:
			xlabel_1D = 'alignment score per aligned read base'
		elif scoreType==3:
			xlabel_1D = 'alignment score per read base'
		title='%s'%(os.path.split(inputFname)[1])
		if plotType in [2, 3]:
			from variation.src.misc import CNV
			if plotType==2:
				reduce_C_function = CNV.logSum
				colorBarLabel='log(count)'
			elif plotType==3:
				import numpy
				reduce_C_function = numpy.mean
				if scoreType==2:
					colorBarLabel='median aligned read length'
				elif scoreType==3:
					colorBarLabel='median read length'
			CNV.drawHexbin(processor.mapq_ls, processor.score_ls, processor.C_ls, fig_fname=outputFname, gridsize=20, \
								title=title, \
								xlabel = 'read map quality', \
								ylabel = xlabel_1D,\
								colorBarLabel=colorBarLabel, reduce_C_function= reduce_C_function)
		elif plotType==1:
			import pylab
			pylab.hist(processor.score_ls, 20, log=True)
			pylab.title(title)
			pylab.xlabel(xlabel_1D)
			pylab.savefig(outputFname, dpi=200)
		
	"""
		# 2011-2-1
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_ref_1MbBAC.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType=1
		outputFname = os.path.expanduser("%s.score%s.hist.png"%(inputPrefix, scoreType))
		VariantDiscovery.drawHistogramOfChosenBWAOutputScore(inputFname, outputFname)
		sys.exit(0)
		
		# 2011-2-1
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC_c30.F4")
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19.chr18_26.7M_27.8M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19_20101230.chr9_124M_124.2M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 3
		plotType = 3
		outputFname = os.path.expanduser("%s.score%s.plot%s.hist.png"%(inputPrefix, scoreType, plotType))
		VariantDiscovery.drawHistogramOfChosenBWAOutputScore(inputFname, outputFname, scoreType=scoreType, plotType=plotType)
		sys.exit(0)
		
		
		#2011-2-24
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_hg19/454_vs_hg19.3eQTL")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 2
		plotType = 2
		outputFname = os.path.expanduser("%s.score%s.plot%s.hist.png"%(inputPrefix, scoreType, plotType))
		VariantDiscovery.drawHistogramOfChosenBWAOutputScore(inputFname, outputFname, scoreType=scoreType, plotType=plotType)
		sys.exit(0)
	"""
	
	class FilterReadsByASPerAlignedBaseAndMapQ(object):
		"""
		2011-2-8
			This class filters bwa output bam file based on two criteria.
				read.mapq>=minMapQ
				score >= minPerBaseAS
			only for bwasw, which has AS (alignment score) tag.
			
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include bamOutputF, scoreType, minMapQ, minPerBaseAS
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-2-8
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			param_obj.real_counter += 1
			score = None
			if self.scoreType==1:
				score = read.mapq
			elif self.scoreType in [2,3]:
				if self.scoreType==2:
					no_of_bases = read.alen
				else:
					# self.scoreType==3:
					no_of_bases = read.rlen
				for tag in read.tags:
					if tag[0]=='AS':
						score = tag[1]	#'AS'
						score = score/float(no_of_bases) 	#divide the alignment score by the read length
						break
			if score is not None:
				if read.mapq>=self.minMapQ and score>=self.minPerBaseAS:
					self.bamOutputF.write(read)
			
			
	@classmethod
	def filterReadsByASPerAlignedBaseAndMapQ(cls, inputFname, outputFname=None, minPerBaseAS=0.5, minMapQ=125, scoreType=2):
		"""
		2011-2-8
			scoreType
				2: minPerBaseAS is alignment score per aligned base.
				3: minPerBaseAS is alignment score per read base.
		"""
		import os, sys
		import pysam
		samfile = pysam.Samfile(inputFname, "rb" )
		bamOutputF = pysam.Samfile(outputFname, 'wb', template=samfile)
		sys.stderr.write("Filter reads from %s (minPerBaseAS>=%s, minMapQ>=%s) ...\n"%(inputFname, minPerBaseAS, minMapQ))
		processor = cls.FilterReadsByASPerAlignedBaseAndMapQ(bamOutputF=bamOutputF, \
									scoreType=scoreType, minMapQ=minMapQ, minPerBaseAS=minPerBaseAS)
		cls.traverseBamFile(samfile, processor=processor)
		
		bamOutputF.close()
	"""
		#2011-2-8
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC.F4")
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19.chr18_26.7M_27.8M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19_20101230.chr9_124M_124.2M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 3
		minPerBaseAS = 0.5
		minMapQ = 125
		outputFname = os.path.expanduser("%s.minPerBaseAS%s.minMapQ%s.score%s.bam"%(inputPrefix, minPerBaseAS, minMapQ, scoreType))
		
		VariantDiscovery.filterReadsByASPerAlignedBaseAndMapQ(inputFname, outputFname, minPerBaseAS=minPerBaseAS, minMapQ=minMapQ,\
				scoreType=scoreType)
		sys.exit(0)
		
		#2011-2-18
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_hg19/454_vs_hg19.3eQTL")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_hg19/454_vs_hg19.3eQTL")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 2
		minPerBaseAS = 0.4
		minMapQ = 125
		outputFname = os.path.expanduser("%s.minPerBaseAS%s.minMapQ%s.score%s.bam"%(inputPrefix, minPerBaseAS, minMapQ, scoreType))
		
		VariantDiscovery.filterReadsByASPerAlignedBaseAndMapQ(inputFname, outputFname, minPerBaseAS=minPerBaseAS, minMapQ=minMapQ,\
				scoreType=scoreType)
		sys.exit(0)
	"""
	@classmethod
	def traverseBamFile(cls, samfile, processor=None):
		"""
		2011-2-8
			a traverser used by other functions
		"""
		import os,sys
		import pysam
		from pymodule import PassingData
		it = samfile.fetch()
		counter = 0
		real_counter = 0
		qname2count = {}
		param_obj = PassingData(real_counter=real_counter, counter=counter, qname2count=qname2count)
		for read in it:
			counter += 1
			processor.run(read, param_obj=param_obj)
				
			if counter%10000==0:
				sys.stderr.write("%s\t%s"%('\x08'*80, counter))
		processor.qname2count = param_obj.qname2count	#2011-2-9 pass it to the processor
		max_redundant_read_count = max(param_obj.qname2count.values())
		sys.stderr.write(" %s unique reads among %s mapped reads, max redundant read count=%s. Done.\n"%\
						(len(param_obj.qname2count), param_obj.real_counter, max_redundant_read_count))
		del samfile
	"""
		
	"""
	
	
	class SelectReadsAlignedInMultiplePlaces(object):
		"""
		2011-2-8
			This class select reads that appear in more than minAlignmentOccurrence alignments.
				
			
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include minAlignmentOccurrence, bamOutputF
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-2-8
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			
			
			if read.qname not in self.qname2count:
				sys.stderr.write("Skip. Read %s not in qname2count.\n"%(read.qname))
				return
			count = self.qname2count.get(read.qname)
			if count>=self.minAlignmentOccurrence:
				param_obj.real_counter += 1
				self.bamOutputF.write(read)
	
	@classmethod
	def selectReadsAlignedInMultiplePlaces(cls, inputFname, outputFname=None, minAlignmentOccurrence=2):
		"""
		2011-2-8
		"""
		import os, sys
		import pysam
		samfile = pysam.Samfile(inputFname, "rb" )
		bamOutputF = pysam.Samfile(outputFname, 'wb', template=samfile)
		sys.stderr.write("Select multi-alined reads from %s...\n"%(inputFname, ))
		
		processorRecord = cls.DrawHistogramOfChosenBWAOutputScore(scoreType=2, \
									plotType=1, mapq_ls=[], C_ls=[], score_ls=[])
		# this processor is only used to get qname2count
		cls.traverseBamFile(samfile, processor=processorRecord)
		
		samfile.seek(0)
		
		processor = cls.SelectReadsAlignedInMultiplePlaces(bamOutputF=bamOutputF, \
											minAlignmentOccurrence=minAlignmentOccurrence, \
											qname2count=processorRecord.qname2count)
		cls.traverseBamFile(samfile, processor=processor)
		
		bamOutputF.close()
	"""
		#2011-2-9
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC_c30.F4")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19.chr18_26.7M_27.8M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19_20101230.chr9_124M_124.2M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		plotType = 1
		minAlignmentOccurrence = 2
		outputFname = os.path.expanduser("%s.minAlnOccu%s.bam"%(inputPrefix, minAlignmentOccurrence))
		VariantDiscovery.selectReadsAlignedInMultiplePlaces(inputFname, outputFname, minAlignmentOccurrence=minAlignmentOccurrence)
		sys.exit(0)

		#2011-2-9
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC.F4.minPerBaseAS0.75.minMapQ150.score2")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19.chr18_26.7M_27.8M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/fineMapping/454_vs_hg19_20101230.chr9_124M_124.2M")
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		plotType = 1
		minAlignmentOccurrence = 2
		outputFname = os.path.expanduser("%s.minAlnOccu%s.bam"%(inputPrefix, minAlignmentOccurrence))
		VariantDiscovery.selectReadsAlignedInMultiplePlaces(inputFname, outputFname, minAlignmentOccurrence=minAlignmentOccurrence)
		sys.exit(0)
	"""
	
	@classmethod
	def discoverHetsFromBAM(cls, inputFname, outputFname, minNumberOfReads=4, monomorphicDiameter=100):
		"""
		2011-2-18
			discover Hets from BAM files based on coverage of either allele and monomorphic span
			not tested and not finished.
		"""
		import pysam
		samfile = pysam.Samfile("ex1.bam", "rb" )
		current_locus = None	# record of polymorphic loci
		previous_locus = None
		candidate_locus = None
		good_polymorphic_loci = []
		for pileupcolumn in samfile.pileup():
			print
			print 'coverage at base %s %s = %s' % (pileupcolumn.tid, pileupcolumn.pos , pileupcolumn.n)
			current_locus = (pileupcolumn.tid, pileupcolumn.pos)
			base2count = {}
			for pileupread in pileupcolumn.pileups:
				base = pileupread.alignment.seq[pileupread.qpos]
				if base not in base2count:
					base2count[base] = 0
				base2count[base] += 1
				#print '\tbase in read %s = %s' % (pileupread.alignment.qname, base)
			if len(base2count)>=2:
				if previous_locus!=None and previous_locus[0]==current_locus[0]:
					gap = current_locus[1]-previous_locus[1]
					if gap>=monomorphicDiameter:
						if candidate_locus is not None and candidate_locus==previous_locus:
							#prior candidate locus is in proper distance. there's no polymorphic locus in between.
							good_polymorphic_loci.append(candidate_locus)
						candidate_locus = current_locus
					else:
						candidate_locus = None
				previous_locus = current_locus
		samfile.close()
		
	
class Main(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 0, ): ['', 'i', 1, 'common input file.', ],\
							('output_fname', 0, ): ['', 'o', 1, 'common output file', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		This class is the entry to all others.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self,):
		if self.debug:	# 2010-4-18 enter debug mode "~/.../variation/misc.py -b"
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		#import Stock_250kDB
		#db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
		#				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		#db_250k.setup(create_tables=False)
		#self.db_250k = db_250k
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		#curs = conn.cursor()
		
		#2011-3-4
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		input_fname = '%s.vcf'%(common_prefix)
		minNumberOfReads=2
		output_fname = '%s_min%s.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHetsFromVCF(input_fname, output_fname, minNumberOfReads=minNumberOfReads, VCFOutputType=2)
		#sys.exit(0)
		
		#2011-3-4
		minNumberOfReads =2
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		#common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		myVariantFile = '%s_min%s.hets'%(common_prefix, minNumberOfReads)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		output_fname = '%s_min%s.overlap.tsv'%(common_prefix, minNumberOfReads)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, output_fname)
		sys.exit(0)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Main
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
