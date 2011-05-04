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
else:	#32bit
	#sys.path.insert(0, os.path.expanduser('~/lib/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment



class PhastCons(object):
	"""
	2011-4-6
	"""
	def __init__(self):
		pass
	
	@classmethod
	def drawPhastConsScoreHist(cls, inputFname, outputFname, samplingFrequency=0.001):
		"""
		2011-4-6
		"""
		sys.stderr.write("Drawing histogram of phastCons score (sampling %s) from %s ... \n"%\
						(samplingFrequency, os.path.basename(inputFname)))
		import gzip, random, csv, re
		if inputFname[-2:]=='gz':
			inf = gzip.open(inputFname, 'rb')
		else:
			inf = open(inputFname)
		#reader = csv.reader(inf)
		
		header_pt = re.compile(r'^[a-zA-Z]')
		
		counter = 0
		real_counter = 0
		score_ls = []
		for line in inf:
			counter += 1
			if not header_pt.search(line):
				if random.random()<=samplingFrequency:
					score_ls.append(float(line.strip()))
					real_counter += 1
			if counter%50000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
		
		if len(score_ls)>100:
			title = 'phastCons hist of %s'%(os.path.basename(inputFname))
			xlabel_1D = 'phastCons score'
			import pylab
			pylab.clf()
			pylab.hist(score_ls, 30, log=True)
			pylab.title(title)
			pylab.xlabel(xlabel_1D)
			pylab.savefig(outputFname, dpi=200)
	"""
		
		#2011-4-6
		common_prefix = os.path.expanduser('/Network/Data/UCSC/phastCons46way/primates/chr1.phastCons46way.primates')
		common_prefix = os.path.expanduser('/Network/Data/UCSC/phastCons46way/primates/chrUn_gl000212.phastCons46way.primates')
		common_prefix = os.path.expanduser('/Network/Data/UCSC/phyloP46way/primates//chr22.phyloP46way.primate')
		inputFname = '%s.wigFix.gz'%(common_prefix)
		samplingFrequency=0.1
		outputFname = '%s.sample%s.hist.png'%(common_prefix, samplingFrequency)
		PhastCons.drawPhastConsScoreHist(inputFname, outputFname, samplingFrequency=samplingFrequency)
		sys.exit(0)
	"""
	
	@classmethod
	def addLocusToDB(cls, db_vervet, chromosome, start, stop, locus_method_id=None, ref_ind_seq_id=None):
		"""
		2011-4-6
			first check database to see if same locus exists in db
		"""
		import VervetDB
		locus_method = VervetDB.LocusMethod.get(locus_method_id)
		db_entry = VervetDB.Locus.query.filter_by(chromosome=chromosome).filter_by(start=start).filter_by(stop=stop).\
			filter_by(ref_ind_seq_id=ref_ind_seq_id).first()
		if db_entry:
			locus_method_id_ls = [locus_method.id for locus_method in db_entry.locus_method_ls]
			if locus_method_id not in locus_method_id_ls:
				db_entry.locus_method_ls.append(locus_method)
				db_vervet.session.add(db_entry)
		else:
			db_entry = VervetDB.Locus(chromosome=chromosome, start=start, stop=stop, \
					ref_ind_seq_id=ref_ind_seq_id)
			db_entry.locus_method_ls.append(locus_method)
			db_vervet.session.add(db_entry)
		return db_entry
	
	@classmethod
	def addLocusScoreToDB(cls, db_vervet, locus=None, score_method_id=None, score=None):
		"""
		2011-4-6
			
		"""
		import VervetDB
		if locus.id:
			db_entry = VervetDB.LocusScore.query.filter_by(score_method_id).filter_by(locus_id=locus.id).first()
		else:
			db_entry = None
		if not db_entry:
			db_entry = VervetDB.LocusScore(score=score, score_method_id=score_method_id)
			db_entry.locus = locus
			db_vervet.session.add(db_entry)
		return db_entry
		
	@classmethod
	def saveHighPhastConsLoci(cls, db_vervet, inputFname, minPhastConsScore=0.95, maxPhastConsScore=1.0, locus_method_id=1,\
							ref_ind_seq_id=9, score_method_id=1, minSegmentLength=1000, commit=False):
		"""
		2011-4-6
			choose segments in which every score is [minPhastConsScore, maxPhastConsScore] and length>minSegmentLength.
				It will save them into db (Locus, LocusScore).
			A better way to define segments is through segmentation algorithm (like GADA).
			
			phastCons scores are from UCSC genome browser.
			
		"""
		sys.stderr.write("Saving segments of high phastCons score (>=%s) from %s into db ... \n"%\
						(minPhastConsScore, os.path.basename(inputFname)))
		import VervetDB
		import gzip, random, re, numpy
		if inputFname[-2:]=='gz':
			inf = gzip.open(inputFname, 'rb')
		else:
			inf = open(inputFname)
		
		
		#lines beginned with letters, are headers 
		header_pt = re.compile(r'^[a-zA-Z]')
		#pattern to discover relevant coordinates
		chr_start_step_pt = re.compile(r'fixedStep chrom=(?P<chr>\w+) start=(?P<start>\d+) step=(?P<step>\d+)')
		
		
		counter = 0
		real_counter = 0
		score_ls = []
		
		chr = None
		start = None
		step = None
		conserved_segment_start = None
		conserved_segment_score_ls = []
		for line in inf:
			counter += 1
			if header_pt.search(line):
				pt_result = chr_start_step_pt.search(line)
				if pt_result:
					chr = pt_result.group('chr')
					if chr[:3]=='chr':	#get rid of the initial "chr" if it exists
						chr = chr[3:]
					start = int(pt_result.group('start'))
					step = int(pt_result.group('step'))
				else:	#can't parse the coordinates of this line. ignore it at this moment. (rather than input error into db)
					sys.stderr.write("Warning: no chr, start, step could be parsed out of %s.\n"%(line))
					chr = None
					start = None
					step = None
					conserved_segment_start = None
					conserved_segment_score_ls = []
			elif chr is not None and start is not None and step is not None:
				start += step
				score = float(line.strip())
				if score>=minPhastConsScore and score<=maxPhastConsScore:
					if conserved_segment_start is None:
						conserved_segment_start = start
						conserved_segment_score_ls = []
					conserved_segment_score_ls.append(score)
				else:
					if conserved_segment_start is not None:
						# found the end of the current conserved segment
						conserved_segment_stop = start - step	#stop of the segment is the prior position
						segment_length = conserved_segment_stop - conserved_segment_start + 1
						if segment_length>=minSegmentLength:
							locus = cls.addLocusToDB(db_vervet, chr, conserved_segment_start, conserved_segment_stop, \
											locus_method_id=locus_method_id, ref_ind_seq_id=ref_ind_seq_id)
							median_score = float(numpy.median(conserved_segment_score_ls))
								#without float(). type of median_score is numpy.float64, which would cause sqlalchemy "can't adapt programming error". 
							locus_score = cls.addLocusScoreToDB(db_vervet, locus=locus, score_method_id=score_method_id, \
												score=median_score)
							real_counter += 1
							if real_counter%1000==0 and real_counter>0:
								db_vervet.session.flush()
					#reset
					conserved_segment_start = None
					conserved_segment_score_ls = []
					
			if counter%50000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("%s%s\t%s\n"%('\x08'*80, counter, real_counter))
		
		db_vervet.session.flush()
		if commit:
			db_vervet.session.commit()
		
	"""
		#2011-4-6
		inputFname
		minPhastConsScore=0.95
		locus_method_id=1
		ref_ind_seq_id=9
		score_method_id=1
		minSegmentLength=1000
		commit=False
		PhastCons.saveHighPhastConsLoci(db_vervet, inputFname, minPhastConsScore=minPhastConsScore, \
			locus_method_id=locus_method_id, ref_ind_seq_id=ref_ind_seq_id, score_method_id=score_method_id, \
			minSegmentLength=minSegmentLength, commit=commit)
		sys.exit(0)
		
		#2011-4-6
		inputDir = '/Network/Data/UCSC/phastCons46way/primates/'
		for inputFname in os.listdir(inputDir):
			if inputFname[-9:]=='wigFix.gz' or inputFname[-6:]=='wigFix':
				#common_prefix = os.path.expanduser('/Network/Data/UCSC/phastCons46way/primates/chr22.phastCons46way.primates')
				#inputFname = '%s.wigFix.gz'%(common_prefix)
				inputFname = os.path.join(inputDir, inputFname)
				minPhastConsScore=0.4
				maxPhastConsScore=0.6
				locus_method_id=1
				ref_ind_seq_id=9
				score_method_id=1
				minSegmentLength=500
				commit = True
				try:
					PhastCons.saveHighPhastConsLoci(db_vervet, inputFname, minPhastConsScore=minPhastConsScore, \
												locus_method_id=locus_method_id, ref_ind_seq_id=ref_ind_seq_id, \
												score_method_id=score_method_id, \
												minSegmentLength=minSegmentLength, commit=commit)
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
		sys.exit(0)
	"""

def sortCMPBySecondTupleValue(a, b):
	"""
	2011-3-29
		a and b are list or tuple
	"""
	return cmp(a[1], b[1])


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
	def discoverHetsFromVCF(cls, inputFname, outputFname, minMinorAlleleCoverage=4, VCFOutputType=1, maxMinorAlleleCoverage=8,\
						maxNoOfReads=30):
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
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		maxMinorAlleleCoverage=7
		outputFname = '%s_minMAC%s_maxMAC%s.hets'%(common_prefix, minMinorAlleleCoverage, maxMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage,\
			VCFOutputType=2, maxMinorAlleleCoverage=maxMinorAlleleCoverage)

		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		outputType=3
		outputFname = '%s_min%s.outputType%s.tsv'%(common_prefix, minMinorAlleleCoverage, outputType)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname, outputType=outputType)
		sys.exit(0)
		
		#2011-3-24
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.D100'
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.D100'
		#common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/454_illu_6_sub_vs_1MbBAC.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		maxMinorAlleleCoverage=7
		outputFname = '%s_minMAC%s_maxMAC%s.hets'%(common_prefix, minMinorAlleleCoverage, maxMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage,\
			VCFOutputType=1, maxMinorAlleleCoverage=maxMinorAlleleCoverage)
		sys.exit(2)
		
	"""
	
	@classmethod
	def checkOverlapping(cls, myVariantFile, jessicaVariantFname, outputFname, outputType=1):
		"""
		2011-3-8
			add argument outputType
				1: output overlapping
				2: output variants only in myVariantFile
				3: output variants only in jessicaVariantFname
		2011-3-6
			deal with the varying number of columns in myVariantFile
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
		myVariantFileHeader = reader.next()
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		if outputType in (1,2):	#2011-3-8
			header = ['snp_id', 'chr_pos|string', 'status',] + myVariantFileHeader[2:]
		elif outputType==3:	#2011-3-8
			header = ['snp_id', 'chr_pos|string', 'status',]
		writer.writerow(header)
		
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(myVariantFileHeader)
		counter = 0
		real_counter = 0
		chr_pos_set_MyVariantFile = set()
		for row in reader:
			locusChosen = False
			chr_pos = row[col_name2index['snp_id']]
			chr_pos_set_MyVariantFile.add(chr_pos)
			pos = row[col_name2index['pos']]
			if outputType==1 and chr_pos in chr_pos2data_ls:	#2011-3-8
				locusChosen = True
			elif outputType==2 and chr_pos not in chr_pos2data_ls:	#2011-3-8
				locusChosen = True
			if locusChosen:
				snp_id, work_status = chr_pos2data_ls[chr_pos]
				
				data_row = [snp_id, chr_pos, work_status] + row[2:]
				writer.writerow(data_row)
				real_counter += 1
			
			counter += 1
		
		if outputType==3:	#2011-3-8
			for chr_pos, data_ls in chr_pos2data_ls.iteritems():
				if chr_pos not in chr_pos_set_MyVariantFile:
					snp_id, work_status = chr_pos2data_ls[chr_pos]
					data_row = [snp_id, chr_pos, work_status]
					writer.writerow(data_row)
					real_counter += 1
		del writer
		sys.stderr.write("%s/%s of mine overlapped with jessica's. Done.\n"%(real_counter, counter))
	"""
	
		#2011-1-6
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_ANXA_p1'
		inputFname = '%s.raw.vcf'%(common_prefix)
		minMinorAlleleCoverage=3
		outputFname = '%s_min%s.raw.hets'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage)
		sys.exit(0)
		
		#2011-1-6
		myVariantFile = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		myVariantFile = '%s_min%s.raw.hets'%(common_prefix, minMinorAlleleCoverage)
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		outputFname = '%s_min%s.overlap.tsv'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname)
		sys.exit(0)
		
		#2011-3-4
		minMinorAlleleCoverage =2
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		#common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		outputFname = '%s_min%s.overlap.tsv'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname)
		sys.exit(0)
		
		
		#2011-3-4
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/454_illu_6_sub_vs_1MbBAC.GATK')
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		outputFname = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage, VCFOutputType=2)
		sys.exit(0)
		
		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		outputType=3
		outputFname = '%s_min%s.outputType%s.tsv'%(common_prefix, minMinorAlleleCoverage, outputType)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname, outputType=outputType)
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
		cls.traverseBamByRead(samfile, processor=processor)
		
		sys.stderr.write("Done.\n")
		if scoreType==1:
			xlabel_1D = "read map quality"
		elif scoreType==2:
			xlabel_1D = 'alignment score per aligned read base'
		elif scoreType==3:
			xlabel_1D = 'alignment score per read base'
		title='%s'%(os.path.split(inputFname)[1])
		if plotType in [2, 3]:
			from pymodule import yh_matplotlib
			if plotType==2:
				reduce_C_function = yh_matplotlib.logSum
				colorBarLabel='log(count)'
			elif plotType==3:
				import numpy
				reduce_C_function = numpy.mean
				if scoreType==2:
					colorBarLabel='median aligned read length'
				elif scoreType==3:
					colorBarLabel='median read length'
			yh_matplotlib.drawHexbin(processor.mapq_ls, processor.score_ls, processor.C_ls, fig_fname=outputFname, gridsize=20, \
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
	
	class PairEndBWAOutputScoreTraverser(object):
		"""
		2011-3-23
			
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include scoreType, plotType, mapq_ls, score_ls, C_ls
			"""
			self.real_counter = 0
			self.readId2ScoreLs = {}	#2011-3-23
			
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-3-23
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			PE_id = read.qname[:-2]	#The last two characters (#0 or #1) are used to differentiate the PE pairs.
			if PE_id not in self.readId2ScoreLs:
				self.readId2ScoreLs[PE_id] = []
			
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
				self.readId2ScoreLs[PE_id].append(score)
				
				self.mapq_ls.append(read.mapq)
				self.score_ls.append(score)
				if self.plotType==2 or self.plotType==1:
					self.C_ls.append(1)
				elif self.plotType==3:
					self.C_ls.append(no_of_bases)	#read.rlen is different from read.alen
	
	
	@classmethod
	def drawHistogramOfPairEndBWAOutputScore(cls, inputFname, outputFnamePrefix, scoreType=1, plotType=2):
		"""
		2011-3-23
			2D histogram of mapq of mate1 vs mapq of mate2
			
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
		processor = cls.PairEndBWAOutputScoreTraverser(scoreType=scoreType, \
									plotType=plotType, mapq_ls=mapq_ls, C_ls=C_ls, score_ls=score_ls)
		samfile = pysam.Samfile(inputFname, "rb" )
		cls.traverseBamByRead(samfile, processor=processor)
		
		small_mapq_ls = []
		big_mapq_ls = []
		singleton_mapq_ls = []
		no_of_multi_mapped_reads = 0
		counter = 0
		for PE_id, score_ls in processor.readId2ScoreLs.iteritems():
			counter += 1
			if len(score_ls)==1:
				singleton_mapq_ls.append(score_ls[0])
			elif len(score_ls)==2:
				score_ls.sort()
				small_mapq_ls.append(score_ls[0])
				big_mapq_ls.append(score_ls[1])
			else:
				no_of_multi_mapped_reads += 1
		
		sys.stderr.write("%s singletons, %s multi-mapped, %s total PE pairs. Done.\n"%\
						(len(singleton_mapq_ls), no_of_multi_mapped_reads, counter))
		if scoreType==1:
			xlabel_1D = "read map quality of mate1"
		elif scoreType==2:
			xlabel_1D = 'alignment score per aligned read base'
		elif scoreType==3:
			xlabel_1D = 'alignment score per read base'
		title='%s'%(os.path.split(inputFname)[1])
		"""
		if plotType in [2, 3]:
			from pymodule import yh_matplotlib
			if plotType==2:
				reduce_C_function = yh_matplotlib.logSum
				colorBarLabel='log(count)'
			elif plotType==3:
				import numpy
				reduce_C_function = numpy.mean
				if scoreType==2:
					colorBarLabel='median aligned read length'
				elif scoreType==3:
					colorBarLabel='median read length'
		"""
		from pymodule import yh_matplotlib
		reduce_C_function = yh_matplotlib.logSum
		colorBarLabel = 'log(count)'
		yh_matplotlib.drawHexbin(small_mapq_ls, big_mapq_ls, [1]*len(small_mapq_ls), fig_fname='%s_2D_hist_PE_mapq.png'%(outputFnamePrefix), \
					gridsize=20, title=title, \
					xlabel = 'read map quality', ylabel = xlabel_1D,\
					colorBarLabel=colorBarLabel, reduce_C_function= reduce_C_function)
		
		import pylab
		pylab.clf()
		pylab.hist(singleton_mapq_ls, 20, log=True)
		pylab.title(title)
		pylab.xlabel(xlabel_1D)
		pylab.savefig('%s_singleton_mapq_hist.png'%outputFnamePrefix, dpi=200)
	
	"""
		#2011-3-23
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		plotType = 1
		outputFnamePrefix = os.path.expanduser("%s.score%s.plot%s"%(inputPrefix, scoreType, plotType))
		VariantDiscovery.drawHistogramOfPairEndBWAOutputScore(inputFname, outputFnamePrefix, scoreType=scoreType, plotType=plotType)
		sys.exit(0)
		
		#2011-3-24
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		plotType = 1
		outputFnamePrefix = os.path.expanduser("%s.score%s.plot%s"%(inputPrefix, scoreType, plotType))
		VariantDiscovery.drawHistogramOfPairEndBWAOutputScore(inputFname, outputFnamePrefix, scoreType=scoreType, plotType=plotType)
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
				if (self.minMapQ is None or read.mapq>=self.minMapQ) and (self.minPerBaseAS is None or score>=self.minPerBaseAS):
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
		cls.traverseBamByRead(samfile, processor=processor)
		
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
		
		#2011-3-7
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 2
		minPerBaseAS = 0.4
		minMapQ = 125
		outputFname = os.path.expanduser("%s.minPerBaseAS%s.minMapQ%s.score%s.bam"%(inputPrefix, minPerBaseAS, minMapQ, scoreType))
		
		VariantDiscovery.filterReadsByASPerAlignedBaseAndMapQ(inputFname, outputFname, minPerBaseAS=minPerBaseAS, minMapQ=minMapQ,\
				scoreType=scoreType)
		sys.exit(0)
		
		#2011-3-24 filter short-read (no AS score)
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		minPerBaseAS = None
		minMapQ = 40
		outputFname = os.path.expanduser("%s.minPerBaseAS%s.minMapQ%s.score%s.bam"%(inputPrefix, minPerBaseAS, minMapQ, scoreType))
		
		VariantDiscovery.filterReadsByASPerAlignedBaseAndMapQ(inputFname, outputFname, minPerBaseAS=minPerBaseAS, minMapQ=minMapQ,\
				scoreType=scoreType)
		sys.exit(0)
	"""
	@classmethod
	def traverseBamByRead(cls, samfile, processor=None):
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
		sys.stderr.write("\n %s unique reads among %s mapped reads, max redundant read count=%s. Done.\n"%\
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
		cls.traverseBamByRead(samfile, processor=processorRecord)
		
		samfile.seek(0)
		
		processor = cls.SelectReadsAlignedInMultiplePlaces(bamOutputF=bamOutputF, \
											minAlignmentOccurrence=minAlignmentOccurrence, \
											qname2count=processorRecord.qname2count)
		cls.traverseBamByRead(samfile, processor=processor)
		
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
	def addCountToDictionaryByKey(cls, dictionary, key):
		"""
		2011-3-24
		"""
		if key not in dictionary:
			dictionary[key] = 0
		dictionary[key] += 1
	
	@classmethod
	def reportValueOfDictionaryByKeyLs(cls, dictionary, key_ls, title=None):
		"""
		2011-3-24
		"""
		if title:
			sys.stderr.write('%s\n'%title)
		for key in key_ls:
			value = dictionary.get(key)
			sys.stderr.write("\t%s: %s\n"%(key, value))

	@classmethod
	def discoverHetsFromBAM(cls, inputFname, outputFname, monomorphicDiameter=100, \
						maxNoOfReads=300, minNoOfReads=2, minMinorAlleleCoverage=3, maxMinorAlleleCoverage=7,\
						maxNoOfReadsForGenotypingError=1, maxMajorAlleleCoverage=30, maxNoOfReadsForAllSamples=1000,\
						nt_set = set(['a','c','g','t','A','C','G','T'])):
		"""
		2011-3-24
			the BAM file needs to support RG tag for all its reads.
		2011-2-18
			discover Hets from BAM files based on coverage of either allele and monomorphic span
			not tested and not finished.
		"""
		import pysam, csv
		sys.stderr.write("Looking for heterozygous SNPs in %s (%s<=MinorAC<=%s), maxNoOfReads=%s, \
				maxNoOfReadsForGenotypingError=%s, maxMajorAlleleCoverage=%s, maxNoOfReadsForAllSamples=%s.\n"%\
					(os.path.basename(inputFname), minMinorAlleleCoverage, maxMinorAlleleCoverage ,\
					maxNoOfReads, maxNoOfReadsForGenotypingError, maxMajorAlleleCoverage, maxNoOfReadsForAllSamples))
		samfile = pysam.Samfile(inputFname, "rb" )
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		#header = ['RG', 'chr', 'pos', 'MinorAlleleCoverage', 'MajorAlleleCoverage']
		#writer.writerow(header)
		current_locus = None	# record of polymorphic loci
		previous_locus = None
		candidate_locus = None
		good_polymorphic_loci = []
		read_group2no_of_snps_with_trialleles = {}
		read_group2no_of_snps_with_quad_alleles = {}
		read_group2no_of_snps_with_penta_alleles = {}
		read_group2no_of_good_hets = {}
		read_group2no_of_good_tris = {}
		counter = 0
		real_counter = 0
		
		read_group2col_index = {}
		locus_id2row_index = {}
		data_matrix = []
		
		for pileupcolumn in samfile.pileup():
			#print
			#print 'coverage at base %s %s = %s'%(pileupcolumn.tid, pileupcolumn.pos , pileupcolumn.n)
			counter += 1
			
			current_locus = '%s_%s'%(pileupcolumn.tid, pileupcolumn.pos+1)
			read_group2base2count = {}
			read_group2depth = {}
			if pileupcolumn.n<=maxNoOfReadsForAllSamples:
				for pileupread in pileupcolumn.pileups:
					read_group = None
					# find the read group
					for tag in pileupread.alignment.tags:
						if tag[0]=='RG':
							tag_value = tag[1]
							if tag_value.find('sorted')==-1:	# sometimes one read has >1 RGs, take the one without 'sorted'
								read_group = tag_value
								break
					if read_group is None:
						sys.stderr.write("This read (tags:%s) has no non-sorted-embedded RG. Exit.\n"%(repr(pileupread.alignment.tags)))
						sys.exit(3)
					if read_group not in read_group2base2count:
						read_group2base2count[read_group] = {}
						read_group2depth[read_group] = 0
					if read_group not in read_group2col_index:
						read_group2col_index[read_group] = len(read_group2col_index)
					
					read_group2depth[read_group] += 1
					base = pileupread.alignment.seq[pileupread.qpos]
					base2count = read_group2base2count.get(read_group)
					if base in nt_set:	#make sure it's a nucleotide
						if base not in base2count:
							base2count[base] = 0
						base2count[base] += 1
					#print '\tbase in read %s = %s' % (pileupread.alignment.qname, base)
				data_row = ['NA']*len(read_group2col_index)
				
				found_one_het = False	#2011 flag to see if any het in all samples is called at this locus.
				allele2count = {}	#2011-3-29
				for read_group, base2count in read_group2base2count.iteritems():
					depth = read_group2depth.get(read_group)
					col_index = read_group2col_index.get(read_group)
					
					if depth>maxNoOfReads:	#2011-3-29 skip. coverage too high.
						continue
					allele = 'NA'	#default
					if len(base2count)>=2:
						item_ls = base2count.items()
						item_ls.sort(cmp=sortCMPBySecondTupleValue)
						
						if len(item_ls)==3:
							cls.addCountToDictionaryByKey(read_group2no_of_snps_with_trialleles, read_group)
							if item_ls[0][1]>maxNoOfReadsForGenotypingError:
								continue
						elif len(item_ls)==4:
							cls.addCountToDictionaryByKey(read_group2no_of_snps_with_quad_alleles, read_group)
							if item_ls[1][1]>maxNoOfReadsForGenotypingError:	# because sorted, count_ls[0] < count_ls[1]
								continue
						elif len(item_ls)>4:	#shouldn't happen. but maybe deletion/insertion + 4 nucleotides
							cls.addCountToDictionaryByKey(read_group2no_of_snps_with_penta_alleles, read_group)
							continue
						MinorAllele = item_ls[-2][0]
						MinorAC = item_ls[-2][1]
						
						MajorAllele = item_ls[-1][0]
						MajorAC = item_ls[-1][1]
						if MinorAC>=minMinorAlleleCoverage and MinorAC<=maxMinorAlleleCoverage and MajorAC<=maxMajorAlleleCoverage:
							real_counter += 1
							found_one_het = True
							#pysam position is 0-based.
							allele = min(MinorAllele, MajorAllele) + max(MinorAllele, MajorAllele)
							#data_row = [read_group, pileupcolumn.tid, pileupcolumn.pos+1, MinorAC, MajorAC]
							#writer.writerow(data_row)
							if len(item_ls)>2:
								cls.addCountToDictionaryByKey(read_group2no_of_good_tris, read_group)
							else:
								cls.addCountToDictionaryByKey(read_group2no_of_good_hets, read_group)
						elif MinorAC<=maxNoOfReadsForGenotypingError:	#2011-3-29 it's homozygous with the major allele
							allele = MajorAllele+MajorAllele
						else:
							continue
					elif len(base2count)==1:
						base = base2count.keys()[0]
						count = base2count.get(base)
						if count>=minMinorAlleleCoverage:
							allele = '%s%s'%(base, base)
					
					data_row[col_index] = allele
					if allele!='NA':
						if allele not in allele2count:
							allele2count[allele] = 0
						allele2count[allele] += 1
				if len(allele2count)>1:	#polymorphic across samples
					locus_id2row_index[current_locus] = len(locus_id2row_index)
					data_matrix.append(data_row)
			if counter%1000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*80, counter, real_counter))
				"""
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
				"""
		samfile.close()
		
		read_group_col_index_ls = read_group2col_index.items()
		read_group_col_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		header = ['locus_id', 'locus_id']+[row[0] for row in read_group_col_index_ls]
		writer.writerow(header)
		
		locus_id_and_row_index_ls = locus_id2row_index.items()
		locus_id_and_row_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		for i in xrange(len(locus_id_and_row_index_ls)):
			locus_id, row_index = locus_id_and_row_index_ls[i]
			data_row = data_matrix[i]
			for j in xrange(len(data_row), len(read_group_col_index_ls)):
				data_row.append('NA')
			writer.writerow([locus_id, locus_id] + data_row)
		del writer
		
		unique_read_group_ls = read_group2col_index.keys()
		unique_read_group_ls.sort()
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_good_hets, unique_read_group_ls, title="No of good hets")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_good_tris, unique_read_group_ls, title="No of good SNPs with tri-or-more alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_trialleles, unique_read_group_ls, title="No of SNPs with tri alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_quad_alleles, unique_read_group_ls, title="No of SNPs with 4 alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_penta_alleles, unique_read_group_ls, title="No of SNPs with 5-or-more alleles")

	
	"""
		#2011-3-24
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/454_illu_6_sub_vs_1MbBAC')
		inputFname = '%s.bam'%(common_prefix)
		minMinorAlleleCoverage=3
		maxMinorAlleleCoverage=7
		maxNoOfReadsForGenotypingError=1
		maxNoOfReads = 30
		maxMajorAlleleCoverage=10
		outputFname = '%s_minMAC%s_maxMAC%s_maxNoOfReadsForGenotypingError%s.maxCoverage%s.maxMajorAC%s.hets'%(common_prefix, \
								minMinorAlleleCoverage, maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError, maxNoOfReads,\
								maxMajorAlleleCoverage)
		VariantDiscovery.discoverHetsFromBAM(inputFname, outputFname, monomorphicDiameter=100, \
						maxNoOfReads=maxNoOfReads, minNoOfReads=None, minMinorAlleleCoverage=minMinorAlleleCoverage, \
						maxMinorAlleleCoverage=maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError=maxNoOfReadsForGenotypingError,\
						maxMajorAlleleCoverage=maxMajorAlleleCoverage)
		sys.exit(2)
	"""
	
	@classmethod
	def calculatePairwiseDistanceOutOfSNPXStrainMatrix(cls, inputFname, outputFname, convertHetero2NA=False,
													max_NA_rate=0.4, min_MAF=0.2):
		"""
		2011-4-7 output the pairwise distance as matrix
		2011-3-30
			add argument convertHetero2NA, max_NA_rate, min_MAF
		2011-3-29
			inputFname format:
				first row is sample id.
				first and 2nd column are the same locus id.
			genotypes are in alphabets.
		"""
		from pymodule.SNP import SNPData, transposeSNPData
		oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
							input_alphabet=1, turn_into_integer=1)
		snpData = transposeSNPData(oldSNPData)
		if convertHetero2NA:
			snpData = SNPData.convertHetero2NA(snpData)
		
		snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
		snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
		
		# add outputFname to function below to output the row pairwise distance
		row_id2pairwise_dist_ls = snpData.calRowPairwiseDist(assumeBiAllelic=True, outputFname=outputFname)
	
		#2011-4-7 output the pairwise distance as matrix
		import csv, numpy
		no_of_rows = len(row_id2pairwise_dist_ls)
		row_id_ls = row_id2pairwise_dist_ls.keys()
		row_id_ls.sort()
		row_id2index = {}
		for row_id in row_id_ls:
			row_id2index[row_id] = len(row_id2index)
		
		data_matrix = numpy.zeros([no_of_rows, no_of_rows], dtype=numpy.float)
		writer = csv.writer(open(outputFname, 'a'), delimiter='\t')
		for row_id in row_id_ls:
			pairwise_dist_ls = row_id2pairwise_dist_ls.get(row_id)
			for dist in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = dist[:4]
				i = row_id2index.get(row_id)
				j = row_id2index.get(row_id2)
				data_matrix[i][j] = mismatch_rate
				data_matrix[j][i] = mismatch_rate
		header = [''] + row_id_ls
		writer.writerow(header)
		for i in range(no_of_rows):
			row_id = row_id_ls[i]
			data_row = [row_id] + list(data_matrix[i])
			writer.writerow(data_row)
		del writer
		
	"""
	
		#2011-3-24
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet-web/data/454_illu_6_sub_vs_1MbBAC')
		inputFname = '%s.bam'%(common_prefix)
		minMinorAlleleCoverage=3
		maxMinorAlleleCoverage=7
		maxNoOfReadsForGenotypingError=1
		maxNoOfReads = 20
		maxMajorAlleleCoverage=10
		maxNoOfReadsForAllSamples = 1000
		outputFname = '%s_minMAC%s_maxMAC%s_maxNoOfReadsForGenotypingError%s.maxCoverage%s.maxMajorAC%s.maxNoOfReadsForAllSamples%s.snps'%(common_prefix, \
								minMinorAlleleCoverage, maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError, maxNoOfReads,\
								maxMajorAlleleCoverage, maxNoOfReadsForAllSamples)
		VariantDiscovery.discoverHetsFromBAM(inputFname, outputFname, monomorphicDiameter=100, \
						maxNoOfReads=maxNoOfReads, minNoOfReads=None, minMinorAlleleCoverage=minMinorAlleleCoverage, \
						maxMinorAlleleCoverage=maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError=maxNoOfReadsForGenotypingError,\
						maxMajorAlleleCoverage=maxMajorAlleleCoverage, maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples)
		#sys.exit(2)
		
		#2011-3-29
		inputFname = outputFname
		max_NA_rate = 0.4
		min_het_frequency = 0.2
		outputFnamePrefix = '%s.max_NA_rate%s.min_het_frequency%s'%(os.path.splitext(inputFname)[0], max_NA_rate, min_het_frequency)
		VariantDiscovery.studyHeterozygousCallOutOfSNPXStrainMatrix(inputFname, outputFnamePrefix,\
												max_NA_rate=max_NA_rate, min_het_frequency=min_het_frequency)
		sys.exit(0)
		
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0.2
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, \
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
	"""
	@classmethod
	def studyHeterozygousCallOutOfSNPXStrainMatrix(cls, inputFname, outputFnamePrefix,
										max_NA_rate=0.4, min_het_frequency=0.2, report=True):
		"""
		2011-3-31 add argument max_NA_rate=0.4, min_het_frequency=0.2
			filter out loci with low frequency heterozygous alleles
			turn all homozygous to missing
			calculate pairwise distance based on heterozygous calls only. 
		2011-3-29
			inputFname format:
				first row is sample id.
				first and 2nd column are the same locus id.
			genotypes are in alphabets.
		"""
		from pymodule.SNP import SNPData, transposeSNPData
		oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
							input_alphabet=1, turn_into_integer=1)
		snpData = transposeSNPData(oldSNPData)
		
		snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
		
		#snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
		
		# filter out loci with low frequency heterozygous alleles
		# turn all homozygous to missing
		# calculate pairwise distance based on heterozygous calls only. 
		no_of_rows = len(snpData.data_matrix)
		no_of_cols = len(snpData.data_matrix[0])
		allele2count_ls = []
		col_id_to_be_kept_ls = []
		no_of_SNPs_with_more_than_one_distinct_hets = 0
		for j in range(no_of_cols):
			col_id = snpData.col_id_ls[j]
			allele2count_ls.append({})
			for i in range(no_of_rows):
				allele = snpData.data_matrix[i][j]
				if allele>4:
					if allele not in allele2count_ls[j]:
						allele2count_ls[j][allele] = 0
					allele2count_ls[j][allele] += 1
					#if cls.report and allele_index>1:
					#sys.stderr.write("%s (more than 2) alleles at SNP %s (id=%s).\n"%((allele_index+1), j, snpData.col_id_ls[j]))
				else:
					snpData.data_matrix[i][j] = 0	#treat it as NA
			if len(allele2count_ls[j])>1:
				no_of_SNPs_with_more_than_one_distinct_hets += 1
				if report:
					sys.stderr.write("Warning: more than one hets at SNP %s (id=%s), %s.\n"%(j, snpData.col_id_ls[j], \
																					repr(allele2count_ls[j])))
			elif len(allele2count_ls[j])==1:
				MAF = min(allele2count_ls[j].values())/float(no_of_rows)
				if MAF>=min_het_frequency:
					col_id_to_be_kept_ls.append(col_id)
		
		snpData = SNPData.keepColsByColID(snpData, col_id_to_be_kept_ls)
		snpData.no_of_cols_removed = no_of_cols - len(col_id_to_be_kept_ls)
		sys.stderr.write("%s columns filtered by min_het_frequency and %s SNPs with >1 different hets. Done.\n"%\
						(snpData.no_of_cols_removed, no_of_SNPs_with_more_than_one_distinct_hets))		
		
		outputFname = '%s.Het.tsv'%(outputFnamePrefix)
		snpData.tofile(outputFname)
		
		outputFname = '%s.pairwiseDistBasedOnHet.tsv'%(outputFnamePrefix)
		row_id2pairwise_dist = snpData.calRowPairwiseDist(assumeBiAllelic=False, outputFname=outputFname)
		
		# calculate how many heterozygous calls each one has
		#row_id2fractionData[row_id1] = [mismatch_rate, no_of_mismatches, no_of_non_NA_pairs]
		row_id2fractionData = snpData.calFractionOfLociCarryingNonRefAllelePerRow(ref_allele=0,)
		for row_id , fractionData in row_id2fractionData.iteritems():
			print row_id, fractionData
		
		locus_id2het_count = {}
		for col_id,col_index in snpData.col_id2col_index.iteritems():
			data_row = snpData.data_matrix[:,col_index]
			het_count = 0
			for i in xrange(len(snpData.row_id_ls)):
				if data_row[i]>=5:
					het_count += 1
			locus_id2het_count[col_id] = het_count
		
		import pylab
		pylab.clf()
		het_count_ls = locus_id2het_count.values()
		pylab.hist(het_count_ls, 20, log=True)
		pylab.xlabel('Number of Hets in %s genomes, max_NA_rate=%s, min_het_frequency=%s.'%(len(snpData.row_id_ls), \
													max_NA_rate, min_het_frequency))
		pylab.ylabel('log10(count)')
		
		title='%s'%(os.path.split(inputFname)[1])
		pylab.title(title)
		outputFname = '%s.HetCount.png'%(outputFnamePrefix)
		pylab.savefig(outputFname, dpi=200)
		
		"""
		#2011-3-29
		inputFname = 
		outputFname = 
		VariantDiscovery.studyHeterozygousCallOutOfSNPXStrainMatrix(inputFname, outputFname)
		sys.exit(0)
		"""
	
	@classmethod
	def drawCoverageHistFromBAM(cls, inputFname, outputFname, minMinorAlleleCoverage=4, monomorphicDiameter=100, \
						maxNoOfReads=None, minNoOfReads=2):
		"""
		2011-3-24
			draw histogram of coverage
		"""
		import pysam, math
		samfile = pysam.Samfile(inputFname, "rb" )
		good_polymorphic_loci = []
		coverage_ls = []
		counter = 0
		real_counter = 0
		for pileupcolumn in samfile.pileup():
			counter += 1
			if pileupcolumn.n>0:
				coverage = math.log10(pileupcolumn.n)
			else:
				coverage = -1
			coverage_ls.append(coverage)
			if counter%10000==0:
				sys.stderr.write("%s\t%s"%('\x08'*80, counter))
		
		import numpy
		log10_median = numpy.median(coverage_ls)
		if log10_median>0:
			median_coverage = math.pow(10, log10_median)
		else:
			median_coverage = 0
		sys.stderr.write("\n median coverage of %s bases: %s.\n"%(counter, median_coverage))
		samfile.close()
		import pylab
		pylab.clf()
		title='%s'%(os.path.split(inputFname)[1])
		pylab.hist(coverage_ls, 40, log=True)
		pylab.title(title)
		pylab.xlabel('log10(base-level coverage)')
		pylab.savefig(outputFname, dpi=200)
	
	"""
		
		#2011-3-24
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4")
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4.realigned.sorted")
		inputPrefix = os.path.expanduser("/Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln.F4.minPerBaseASNone.minMapQ40.score1")
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script//vervet-web/data/ref/illumina/ref_illu_vs_1MbBAC_by_aln")
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script//vervet-web/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL")
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_BAC/454_vs_ref_1MbBAC.F4'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_BAC/454_vs_ref_1MbBAC.F4.minPerBaseAS0.5.minMapQ125.score2'
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		outputFname = os.path.expanduser("%s.coverage.hist.png"%(inputPrefix))
		VariantDiscovery.drawCoverageHistFromBAM(inputFname, outputFname)
		sys.exit(2)
	"""


class DBVervet(object):
	"""
	2011-4-27
		class to hold functions related to vervetdb
	"""
	def __init__(self):
		pass
	
	@classmethod
	def putCarribean2011IntoDB(cls, db_vervet, inputFname):
		"""
		2011-4-27
		"""
	
	@classmethod
	def filterValue(cls, value, data_type=None, NA_str_set=set(["", "NA", "N/A", 'n/a'])):
		"""
		2011-4-28
			adapted from variation.src.misc class DBGenome
		"""
		if value in NA_str_set:
			value = None
		if value is not None and data_type is not None:
			value = data_type(value)
		return value
	
	@classmethod
	def put2010SouthAfricanCollectionIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="VSA", \
										collector_name='Christopher A. Schmitt'):
		"""
		2011-4-28
		"""
		sys.stderr.write("Putting 2010 south african collection (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		reader = csv.reader(open(inputFname,), delimiter='\t')
		header = reader.next()
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Monkey ID")
		site_index = col_name2index.get("Location")
		gps_index = col_name2index.get("GPS")
		sex_index = col_name2index.get("Sex")
		age_index = col_name2index.get("Dental Age Category")
		age_cas_index = col_name2index.get("CAS Dental Age")
		phenotype_start_index = col_name2index.get("Age/Dentition present")
		collection_date_index = col_name2index.get("Date")
		
		gps_pattern = re.compile(r'(?P<lat_direction>[SN])(?P<lat_hr>\d+) (?P<lat_min>\d+) (?P<lat_sec>[\d.]+) (?P<lon_direction>[EW])(?P<lon_hr>\d+) (?P<lon_min>\d+) (?P<lon_sec>[\d.]+)')
		any_character_pattern = re.compile(r'[a-zA-Z]')
		collector = db_vervet.getUser(collector_name)
		
		for row in reader:
			monkey_id = monkeyIDPrefix + row[monkey_id_index]
			site_str = row[site_index]
			gps = row[gps_index]
			sex = row[sex_index]
			age = row[age_index]
			age = cls.filterValue(age, int)
			age_cas = row[age_cas_index]
			age_cas = cls.filterValue(age_cas, int)
			collection_date = row[collection_date_index]
			collection_date = datetime.strptime(collection_date, '%Y/%m/%d')	#2008/03/17
			
			gps_pattern_search = gps_pattern.search(gps)
			if gps_pattern_search:
				lat_direction = gps_pattern_search.group("lat_direction")
				lat_hr = float(gps_pattern_search.group('lat_hr'))
				lat_min = float(gps_pattern_search.group('lat_min'))
				lat_sec = float(gps_pattern_search.group('lat_sec'))
				latitude = lat_hr + lat_min/60. + lat_sec/3600.
				if lat_direction=='S':
					latitude = -latitude	#"-" because it's south
				
				lon_direction = gps_pattern_search.group("lon_direction")
				lon_hr = float(gps_pattern_search.group('lon_hr'))
				lon_min = float(gps_pattern_search.group('lon_min'))
				lon_sec = float(gps_pattern_search.group('lon_sec'))
				longitude = lon_hr + lon_min/60. + lon_sec/3600.
				if lon_direction=='W':
					longitude = -longitude
			else:
				sys.stderr.write("Error: gps %s is not parsable.\n"%(gps))
				sys.exit(0)
			
			site_name_ls = site_str.split(",")
			strip_func = lambda x: x.strip()
			site_name_ls = map(strip_func, site_name_ls)
			if len(site_name_ls)==3:
				city, province, country = site_name_ls
				description = None
			elif len(site_name_ls)==4:
				description, city, province, country = site_name_ls
			else:
				sys.stderr.write("Error: site %s is neither 3-entry nor 4-entry.\n"%site_str)
				sys.exit(0)
			site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex[0], age=age, age_cas=age_cas, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector)
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name not in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:
						if any_character_pattern.search(row[i]):	#ignore any column with any character in it.
							# should be number only
							continue
						value = cls.filterValue(row[i], data_type=float)
						comment =None
					else:
						value = None
						comment = row[i]
					db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
	
	"""
		#2011-4-28
		inputFname = os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/2010 SA.tsv")
		DBVervet.put2010SouthAfricanCollectionIntoDB(db_vervet, inputFname)
		sys.exit(0)
	"""
	
	@classmethod
	def put2010AfricanCollectionIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", \
										collector_name=''):
		"""
		2011-4-28
		"""
		sys.stderr.write("Putting 2010 African collection (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		# 2011-4-29 a handy function to strip blanks around strings
		strip_func = lambda x: x.strip()
		
		reader = csv.reader(open(inputFname,), delimiter='\t')
		
		sys.stderr.write("\t Getting data for all sites and its code ...")
		header = reader.next()
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		
		#2011-4-29 first get the location infomation into db.
		code_index = col_name2index.get("CODE")
		site_index = col_name2index.get("CITY")
		gps_index = col_name2index.get("COORDINATES")
		collector_index = col_name2index.get("Collection")
		# to match "S 140 56.564 E 0250 54.361  Altitude 1074 meters", "N 4 22 E 18 35"
		gps_pattern = re.compile(r"(?P<lat_direction>[SN]) (?P<lat_hr>\d+) (?P<lat_min>[\d.]+) +(?P<lon_direction>[EW]) (?P<lon_hr>\d+) (?P<lon_min>[\d.]+) *((Altitude (?P<alt>[\d.]+) meters)|)")
		# a dictionary which maps code to site and collector
		code2data = {}
		for row in reader:
			if row[0]=="Unique ID":
				#the actual monkey individual collection starts from here
				break
			code = row[code_index]
			if code not in code2data:
				code2data[code] = PassingData()
			site_str = row[site_index]
			gps = row[gps_index]
			if gps:
				gps_pattern_search = gps_pattern.search(gps)
				if gps_pattern_search:
					lat_direction = gps_pattern_search.group("lat_direction")
					lat_hr = float(gps_pattern_search.group('lat_hr'))
					lat_min = float(gps_pattern_search.group('lat_min'))
					latitude = lat_hr + lat_min/60.
					if lat_direction=='S':
						latitude = -latitude	#"-" because it's south
					
					lon_direction = gps_pattern_search.group("lon_direction")
					lon_hr = float(gps_pattern_search.group('lon_hr'))
					lon_min = float(gps_pattern_search.group('lon_min'))
					longitude = lon_hr + lon_min/60.
					if lon_direction=='W':
						longitude = -longitude
					
					altitude = gps_pattern_search.group("alt")
					if altitude:
						altitude = float(altitude)
					else:
						altitude = None
				else:
					lat_lon_ls = gps.split(',')
					lat_lon_ls = map(strip_func, lat_lon_ls)
					lat_lon_ls = map(float, lat_lon_ls)
					latitude, longitude = lat_lon_ls[:2]
					altitude = None
					#sys.stderr.write("Error: gps %s is not parsable.\n"%(gps))
					#sys.exit(0)
			else:
				longitude = None
				latitude = None
				altitude = None
			
			city = "Unknown"
			description = None
			province = None
			country = None
			if site_str!='?':
				site_name_ls = site_str.split(",")
				site_name_ls = map(strip_func, site_name_ls)
				if len(site_name_ls)==2:
					city, country = site_name_ls[:2]
				elif len(site_name_ls)==1:
					country = site_name_ls[0]
				else:
					sys.stderr.write("Error: site %s is neither 1-entry nor 2-entry.\n"%site_str)
					sys.exit(0)
				
				site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country,\
								latitude=latitude, longitude=longitude, altitude=altitude)
			else:
				site = None
			collector_name = row[collector_index]
			collector = db_vervet.getUser(collector_name)
			code2data[code].site = site
			code2data[code].collector = collector
		sys.stderr.write("%s sites found.\n"%(len(code2data)))
		
		sys.stderr.write("\t Putting individuals into db ...\n")
		header = row
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Animal ID")
		sex_index = col_name2index.get("Sex")
		age_index = col_name2index.get("Approximate Age")
		location_code_index = col_name2index.get("Location *")
		collection_date_index = col_name2index.get("Collection Date")
		phenotype_start_index = col_name2index.get("Weight (Kg)")
		species_index = col_name2index.get("Species (C=Chlorocebus)")
		
		any_character_pattern = re.compile(r'[a-zA-Z]')
		
		# to match either "June 23-25" or "2009/7/16"
		collection_date_pattern = re.compile(r'((?P<Month>\w+) (?P<Day>\d+)-\d+)|((?P<year>\d+)/(?P<month>\d+)/(?P<day>\d+))')
		
		# to match "Yearling", "Adult", "15 yrs", "5 years"
		age_pattern = re.compile(r"([a-z A-Z]+)|((?P<age>\d+) ((yrs)|(years)))")
		counter = 0
		real_counter = 0
		species2tax_id = {'sabaeus':60711, 'aethiops':101841, 'pygerythrus':460674, 'cynosurus':460675, 'tantalus':60712,\
						'cynosuros':460675, 'pygerytherus':460674 }	#typos
		for row in reader:
			counter += 1
			
			monkey_id = monkeyIDPrefix + row[monkey_id_index]
			location_code = row[location_code_index]
			sex = row[sex_index]
			species = row[species_index].strip()
			if species and species[-1]=='?':	#some name has a trailing ?
				species = species[:-1]
			if species in species2tax_id:
				tax_id = species2tax_id.get(species)
			else:
				if species!='??' and species!='?' and species!='':
					sys.stderr.write("Error, species %s didn't find its tax_id.\n"%species)
					sys.exit(2)
				tax_id = None
			age = row[age_index]
			age_pattern_search_result = age_pattern.search(age)
			if age_pattern_search_result:
				if age_pattern_search_result.group('age'):
					age = int(age_pattern_search_result.group('age'))
					approx_age_group_at_collection = None
				else:
					approx_age_group_at_collection = row[age_index]
					age = None
			else:
				approx_age_group_at_collection = None
				age = None
			collection_date = row[collection_date_index]
			collection_date_pattern_search_result = collection_date_pattern.search(collection_date)
			if collection_date_pattern_search_result:
				if collection_date_pattern_search_result.group('Month'):
					year = "2009"
					month = collection_date_pattern_search_result.group('Month')
					day = collection_date_pattern_search_result.group('Day')
					collection_date = datetime.strptime("%s %s %s"%(year, month, day), '%Y %B %d')	#2008 June 9
				else:
					collection_date = datetime.strptime(collection_date, '%Y/%m/%d')	#2008/03/17
			else:
				collection_date = None
			if location_code not in code2data:
				sys.stderr.write("Error: location code %s nowhere in the location table.\n"%(location_code))
				sys.exit(1)
			site = code2data.get(location_code).site
			collector = code2data.get(location_code).collector
			if site is not None:
				latitude = site.latitude
				longitude = site.longitude
				altitude = site.altitude
			else:
				latitude = None
				longitude = None
				altitude = None
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex[0], age=age, latitude=latitude,\
								longitude=longitude, altitude=altitude, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, \
								approx_age_group_at_collection = approx_age_group_at_collection, tax_id=tax_id)
			
			#phenotype
			value = row[phenotype_start_index]
			if value and not any_character_pattern.search(value):	#ignore any column with any character in it.
					# should be number only:
				phenotype_name = "Weight"
				value = cls.filterValue(value, data_type=float)
				comment =None
				db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
									individual_id=individual.id, comment=comment, collector_name=collector.realname)
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("%s individuals. Done.\n"%(counter))
	"""
		#2011-4-29
		inputFname = os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/Africa 20l09.tsv")
		DBVervet.put2010AfricanCollectionIntoDB(db_vervet, inputFname)
		sys.exit(0)
	"""
	
class Main(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
							('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('inputFname', 0, ): ['', 'i', 1, 'common input file.', ],\
							('outputFname', 0, ): ['', 'o', 1, 'common output file', ],\
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
		
		
		
		import VervetDB
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		#curs = conn.cursor()
		
		#2011-4-29
		inputFname = os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/Africa 20l09.tsv")
		DBVervet.put2010AfricanCollectionIntoDB(db_vervet, inputFname)
		sys.exit(0)
		
		
		#2011-4-7
		inputFname = os.path.expanduser('~/script/vervet-web/data/topConservedHG19_random_100loci_vervetSNP.tsv')
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0.2
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, \
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
		
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Main
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
