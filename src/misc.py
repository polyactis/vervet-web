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

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement



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
	def putFinishedFileIntoReplicaCatalog(cls, inputDir=None, site='uschpc',cutOffTime="2011 09 10 20:00", outputFname=None):
		"""
		2011-9-11
			for a half-finished pegasus workflow (inside inputDir) 
				regard any file before certain timestamp as finished and replica-catalog-eligible
			
			The outputFname contains the replica entries.
			
			inputDir has to be in absolute path. no relative path.
		"""
		from pymodule.utils import getAllFiles
		inputFiles = []
		inputDir = os.path.abspath(inputDir)
		sys.stderr.write("Finding all files in %s ..."%inputDir)
		getAllFiles(inputDir, inputFiles)
		sys.stderr.write("%s files found.\n"%(len(inputFiles)))
		from datetime import datetime
		
		cutOffTime = datetime.strptime(cutOffTime, "%Y %m %d %H:%M")
		outf = open(outputFname, 'w')
		no_of_files = 0
		counter = 0
		for fname in inputFiles:
			counter += 1
			statinfo = os.lstat(fname)	#lstat is same as stat except it doesn't follow sym link
			if datetime.fromtimestamp(statinfo.st_mtime)<=cutOffTime:	#st_mtime is unix time in seconds
				no_of_files += 1
				dirNameLength = len(inputDir)
				relativePath = fname[dirNameLength:]	#"/" should be excluded
				if relativePath[0]=="/":	#but in case
					relativePath = relativePath[1:]
				outf.write('%s file://%s pool="%s"\n'%(relativePath, fname, site))
		del outf
		sys.stderr.write("%s/%s files before cutOffTime %s.\n"%(no_of_files, counter, cutOffTime))
		
	
	"""
		#2011-9-11
		inputDir = os.path.expanduser("~/pg_work/crocea/pegasus/ShortRead2AlignmentPipeline/20110909T231932-0700/")
		site='uschpc'
		cutOffTime="2011 09 10 20:00"
		outputFname = ""
		VariantDiscovery.putFinishedFileIntoReplicaCatalog(inputDir=inputDir, site=site,cutOffTime=cutOffTime, \
			outputFname=outputFname)
		sys.exit(0)
		
		#inputDir = os.path.expanduser("~/pg_work/crocea/pegasus/ShortRead2AlignmentPipeline/20110909T231932-0700/")
		site='uschpc'
		cutOffTime="2011 09 10 20:00"
		outputFname = os.path.expanduser("~/replica.txt")
		VariantDiscovery.putFinishedFileIntoReplicaCatalog(inputDir=self.input, site=site,cutOffTime=cutOffTime, \
			outputFname=outputFname)
		sys.exit(0)
		
		
	"""
	
	
	@classmethod
	def moveFinishedBamIntoTargetFolder(cls, inputDirLs=None, outputDir=None, targetFolderSizeThresholdInMB=1000000,\
									timeGapInMinutes=20):
		"""
		2011-9-12
			targetFolderSizeThreshold is in MB
			As far as the size of outputDir is below targetFolderSizeThreshold, the move will continue. 
		"""
		from pymodule.utils import runLocalCommand
		targetFolderSizeThreshold = targetFolderSizeThresholdInMB*1000
		commandline = "du -s %s |awk '{print $1}'"%(outputDir)
		return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
		folderSizeInKB = int(return_data.stdout_content.strip())
		
		from datetime import datetime, timedelta
		import os, sys
		timeGap = timedelta(minutes=timeGapInMinutes)
		
		no_of_files = 0
		counter = 0
		sleepInSnds=5
		inputFnameLs = []
		while folderSizeInKB<targetFolderSizeThreshold:
			timeToBreak = False
			for inputDir in inputDirLs:
				fnameLs = os.listdir(inputDir)
				for fname in fnameLs:
					fname_ext = os.path.splitext(fname)[1]
					if fname_ext=='.bam':
						counter += 1
						bam_abs_path = os.path.join(inputDir, fname)
						bai_fname = '%s.bai'%(fname)
						bai_fname_abs_path = os.path.join(inputDir, bai_fname)
						
						inputFnameLs.append((bam_abs_path, bai_fname_abs_path))
			sys.stderr.write("%s pairs of potential bam/bai files.\n"%(len(inputFnameLs)))
			
			for bam_abs_path, bai_fname_abs_path in inputFnameLs:
				if folderSizeInKB<targetFolderSizeThreshold:
					if os.path.isfile(bai_fname_abs_path):
						current_time = datetime.now()
						statinfo = os.stat(bai_fname_abs_path)	#stat is same as lstat except it does follow sym link
						bai_mtime = datetime.fromtimestamp(statinfo.st_mtime)
						
						if (current_time-bai_mtime)>timeGap:
							no_of_files += 1
							#move both files
							commandline = 'mv %s %s'%(bam_abs_path, outputDir)
							runLocalCommand(commandline, report_stderr=True, report_stdout=True)
							
							commandline = 'mv %s %s'%(bai_fname_abs_path, outputDir)
							runLocalCommand(commandline, report_stderr=True, report_stdout=True)
							
							sys.stderr.write("Sleeping %s seconds..."%sleepInSnds)
							os.system("sleep %s"%sleepInSnds)
							sys.stderr.write("Wakeup.\n")
							#update folderSizeInKB
							commandline = "du -s %s |awk '{print $1}'"%(outputDir)
							return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
							folderSizeInKB = int(return_data.stdout_content.strip())
				else:	#folderSizeInKB exceeded threshold set. break
					timeToBreak = True
					break	#only break out of the forloop
			if timeToBreak:
				break	#break the while loop
		
		sys.stderr.write("%s/%s files moved. current folder size %s KB.\n"%(no_of_files, counter, folderSizeInKB))
		
	"""
		#2011-9-12
		inputDirLs = [os.path.expanduser("~/pg_work/crocea/pegasus/ShortRead2AlignmentPipeline/20110913T192749-0700/"),]
		
		#os.path.expanduser("~/pg_work/crocea/pegasus/ShortRead2AlignmentPipeline/20110910T002453-0700/"),
		#os.path.expanduser("~/pg_work/crocea/pegasus/ShortRead2AlignmentPipeline/20110910T105558-0700/")
		outputDir = os.path.expanduser("~/NetworkData/vervet/db/individual_alignment/")
		VariantDiscovery.moveFinishedBamIntoTargetFolder(inputDirLs=inputDirLs, outputDir=outputDir, \
														targetFolderSizeThresholdInMB=2000000)
		sys.exit(0)
	"""
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
		jessicaVariantFname = os.path.expanduser('~/script/vervet/data/eQTL summary.txt')
		outputFname = '%s_min%s.overlap.tsv'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname)
		sys.exit(0)
		
		#2011-3-4
		minMinorAlleleCoverage =2
		common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.GATK'
		#common_prefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK'
		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet/data/eQTL summary.txt')
		outputFname = '%s_min%s.overlap.tsv'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, outputFname)
		sys.exit(0)
		
		
		#2011-3-4
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/454_illu_6_sub_vs_1MbBAC.GATK')
		#common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL.RG.GATK')
		inputFname = '%s.vcf'%(common_prefix)
		minMinorAlleleCoverage=2
		outputFname = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		VariantDiscovery.discoverHetsFromVCF(inputFname, outputFname, minMinorAlleleCoverage=minMinorAlleleCoverage, VCFOutputType=2)
		sys.exit(0)
		
		myVariantFile = '%s_min%s.hets'%(common_prefix, minMinorAlleleCoverage)
		
		jessicaVariantFname = os.path.expanduser('~/script/vervet/data/eQTL summary.txt')
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
		
		#2011-7-1
		inputPrefix = os.path.expanduser("~/script/vervet/data/ref/BAC/vsTop150Contigs_by_bwasw/last_updated_vervet_GSS_seqs.sorted")
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
			self.exitAfterNumberOfReads = keywords.get("exitAfterNumberOfReads", None)
			
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
			if self.exitAfterNumberOfReads and param_obj.counter>=self.exitAfterNumberOfReads:
				return 1	#a non-zero or non-None return code will stop traverseBamByRead()
	
	
	@classmethod
	def drawHistogramOfPairEndBWAOutputScore(cls, inputFname, outputFnamePrefix, scoreType=1, plotType=2,
								exitAfterNumberOfReads=None):
		"""
		2011-11-3
			add argument exitAfterNumberOfReads
		2011-3-23
			2D histogram of mapq of mate1 vs mapq of mate2
			
			scoreType
				1. mapping quality
				2. alignment score per aligned read length
				3. alignment score per read length
			plotType is useless here. it's always a 2D histogram.
				2D histogram (each hexagon is colored according to how many reads)
		"""
		import os,sys
		import pysam
		
		sys.stderr.write("Draw histogram of (scoreType=%s, plotType=%s) from %s ...\n"%\
						(scoreType, plotType, inputFname))
		mapq_ls = []
		score_ls = []
		C_ls = []
		processor = cls.PairEndBWAOutputScoreTraverser(scoreType=scoreType, \
					plotType=plotType, mapq_ls=mapq_ls, C_ls=C_ls, score_ls=score_ls, \
					exitAfterNumberOfReads=exitAfterNumberOfReads)
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
				small_mapq_ls.append(score_ls[1])
				big_mapq_ls.append(score_ls[1])
				big_mapq_ls.append(score_ls[0])
			else:
				no_of_multi_mapped_reads += 1
		
		sys.stderr.write("%s singletons, %s multi-mapped, %s total PE pairs. Done.\n"%\
						(len(singleton_mapq_ls), no_of_multi_mapped_reads, counter))
		if scoreType==1:
			xlabel_1D = "read map quality of mate2"
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
		pylab.xlabel("mapq of singleton reads")
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
	
	class FilterAlignmentByReferenceIDs(object):
		"""
		2011-7-8
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
			2011-7-8
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			score = None
			target_name = param_obj.samfile.getrname(read.tid)
			if target_name in self.referenceIDSet:
				param_obj.real_counter += 1
				if self.readGroup:	#2011-7-11 add a read group
					read.tags = read.tags + [("RG", self.readGroup)]
				self.bamOutputF.write(read)
				return False	#like None, break the bam file reading
			else:
				return False	#if return True, means breaking once the target is no longer in referenceIDSet.
			
			
	@classmethod
	def filterAlignmentByReferenceIDs(cls, inputFname, outputFname=None, referenceIDSet = set(['Contig0', 'Contig1', 'Contig2']), \
									readGroup="", platform='LS454'):
		"""
		2011-7-8
		
		"""
		import os, sys
		import pysam, copy
		samfile = pysam.Samfile(inputFname, "rb" )
		header = copy.deepcopy(samfile.header)
		if readGroup:
			if "RG" not in header:
				header['RG'] = []
			header['RG'].append({'ID':readGroup, 'PL':platform, 'LB':platform, 'SM':readGroup})
		bamOutputF = pysam.Samfile(outputFname, 'wb', header=header)	# template=samfile)
		sys.stderr.write("Retain reads from %s only from these references %s ...\n"%(inputFname, referenceIDSet))
		processor = cls.FilterAlignmentByReferenceIDs(bamOutputF=bamOutputF, referenceIDSet=referenceIDSet, readGroup=readGroup)
		cls.traverseBamByRead(samfile, processor=processor)
		bamOutputF.close()
	
	"""
		#2011-7-8
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC.F4")
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/subspecies/aethiops/vs_top150Contigs_by_aln")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		referenceIDSet = set(['Contig0', 'Contig1', 'Contig2'])
		outputFname = os.path.expanduser("%s.Contig0_1_2.bam"%(inputPrefix))
		VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet, \
			readGroup='454_vs_Contig012')
		sys.exit(0)
		
		#2011-7-8
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/ref/454/vs_top150Contigs_by_bwasw")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		referenceIDSet = set(['Contig0',])
		outputFname = os.path.expanduser("%s.Contig0.bam"%(inputPrefix))
		VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet)
		sys.exit(0)
		
		for subspecies in ['cynosurus',  'pygerythrus']:
			inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/subspecies/%s/vs_top150Contigs_by_aln"%(subspecies))
			inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
			referenceIDSet = set(['Contig0', 'Contig1', 'Contig2'])
			outputFname = os.path.expanduser("%s.Contig0.bam"%(inputPrefix))
			VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet)
		sys.exit(0)
	"""
	
	@classmethod
	def filterAlignmentToOnlyTopContigs(cls, inputFname, outputFname=None, topNumber=156, readGroup=""):
		"""
		2011-7-8
		
		"""
		referenceIDSet = set([])
		
		#2011-7-7 output all the BACs in order
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_genome.setup(create_tables=False)
		query = GenomeDB.AnnotAssembly.query.filter_by(tax_id=60711).filter_by(sequence_type_id=9).order_by(GenomeDB.AnnotAssembly.stop)
		for row in query:
			referenceIDSet.add(row.chromosome)
			if len(referenceIDSet)>=topNumber:
				break
		
		#2011-7-8
		#inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC.F4")
		#inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/subspecies/aethiops/vs_top150Contigs_by_aln")
		#inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		#outputFname = os.path.expanduser("%s.top%sContigs.bam"%(inputPrefix, topNumber))
		VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet, readGroup=readGroup)
		sys.exit(0)
		
		"""
		#2011-7-10
		inputPrefix = os.path.expanduser("/Network/Data/vervet/ref/454_vs_BAC/454_vs_ref_1MbBAC.F4")
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/subspecies/aethiops/vs_top150Contigs_by_aln")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		outputFname = os.path.expanduser("%s.top%sContigs.bam"%(inputPrefix, topNumber))
		topNumber = 156
		VariantDiscovery.filterAlignmentToOnlyTopContigs(inputFname, outputFname, topNumber=topNumber)
		sys.exit(0)
		"""
	
	@classmethod
	def traverseBamByRead(cls, samfile, processor=None):
		"""
		2011-7-10
			add samfile to param_obj
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
		param_obj = PassingData(real_counter=real_counter, counter=counter, qname2count=qname2count, samfile=samfile)
		for read in it:
			counter += 1
			exitCode = processor.run(read, param_obj=param_obj)
				
			if counter%10000==0:
				sys.stderr.write("%s\t%s\t\t%s"%('\x08'*80, param_obj.counter, param_obj.real_counter))
			if exitCode:	#2011-7-8
				break
		processor.qname2count = param_obj.qname2count	#2011-2-9 pass it to the processor
		max_redundant_read_count = max(param_obj.qname2count.values())
		sys.stderr.write("\n %s unique reads among %s mapped reads, max redundant read count=%s. Done.\n"%\
						(len(param_obj.qname2count), param_obj.real_counter, max_redundant_read_count))
		del samfile
	"""
		
	"""
	
	class GetBACEndHitsOnContigs(object):
		"""
		2011-6-28
			
		"""
		def __init__(self, **keywords):
			"""
			2011-7-1
				keywords include min_alen_rlen_ratio 
			"""
			self.real_counter = 0
			self.no_of_aligned_bases = 0
			self.no_of_mismatches = 0
			self.no_of_insertions = 0
			self.no_of_deletions = 0
			self.BAC_id2contig_pos_ls = {}
			self.qname2match_percentage = {} 
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-6-28
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			read_name_split = read.qname.split('-')	#full ID looks like "DU835782 MUGQ_CH252-11A14-Sp6". its pair looks like "DU835788 MUGQ_CH252-11A14-T7"
			read_name_split[0] = '_'.join(read_name_split[0].split('_')[1:])
			read_id = '-'.join(read_name_split[:2])
			alen_rlen_ratio = float(read.alen)/float(read.rlen)
			if alen_rlen_ratio<self.min_alen_rlen_ratio:
				return
			
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			
			if read_id not in self.BAC_id2contig_pos_ls:
				self.BAC_id2contig_pos_ls[read_id] = []
			self.BAC_id2contig_pos_ls[read_id].append((read_name_split[-1], read.tid, read.pos))	#store the id of the two ends first
			
			if read.qname not in self.qname2match_percentage:
				cigar_code2no_of_bases = {}
				for cigar_tuple in read.cigar: #presented as a list of tuples (operation,length).
					#For example, the tuple [ (0,3), (1,5), (0,2) ] refers to an alignment with 3 matches, 5 insertions and another 2 matches.
					# CIGAR operation MIDNSHP=X =>012345678
					cigar_code = cigar_tuple[0]
					if cigar_code not in cigar_code2no_of_bases:
						cigar_code2no_of_bases[cigar_code] = 0
					cigar_code2no_of_bases[cigar_code] += cigar_tuple[1]
				no_of_matches = cigar_code2no_of_bases[0]	#0 represents match.
				self.qname2match_percentage[read.qname] = no_of_matches/float(read.alen)
				self.no_of_aligned_bases += read.alen
				self.no_of_insertions += cigar_code2no_of_bases.get(1, 0)
				self.no_of_deletions += cigar_code2no_of_bases.get(2, 0)
				no_of_mismatches = read.alen - no_of_matches
				if 8 in cigar_code2no_of_bases:
					if cigar_code2no_of_bases[8]!=no_of_matches:
						sys.stderr.write("No of mismatches in cigar (%s) doesn't match the number (%s) calculated.\n"%\
										(cigar_code2no_of_bases[8], no_of_mismatches))
						import pdb
						pdb.set_trace()
					self.no_of_mismatches += cigar_code2no_of_bases[8]
				else:
					self.no_of_mismatches += no_of_mismatches
			param_obj.real_counter += 1
	
	
	@classmethod
	def getBACEndHitsOnContigs(cls, inputFname, outputFnamePrefix=None, min_alen_rlen_ratio=0.9, minAlignmentOccurrence=2):
		"""
		2011-6-28
		"""
		import os, sys
		import pysam
		samfile = pysam.Samfile(inputFname, "rb" )
		#bamOutputF = pysam.Samfile(outputFname, 'wb', template=samfile)
		sys.stderr.write("Get where the BAC ends are aligned on the contigs from %s...\n"%(inputFname, ))
		
		
		processor = cls.GetBACEndHitsOnContigs(min_alen_rlen_ratio=min_alen_rlen_ratio)
		cls.traverseBamByRead(samfile, processor=processor)
		sys.stderr.write("One or both ends of %s BACs have been aligned.\n"%(len(processor.BAC_id2contig_pos_ls)))
		
		from pymodule import yh_matplotlib
		import math
		same_contig_BAC_end_pair_dist_ls = []
		contig_id_set = set()
		on_different_contig_BAC_end_contig_pos_ls = []
		for BAC_id, contig_pos_ls in processor.BAC_id2contig_pos_ls.iteritems():
			for contig_pos_item in contig_pos_ls:
				contig_id_set.add(contig_pos_item[0])
			if len(contig_pos_ls) == 2 and contig_pos_ls[0][0]!=contig_pos_ls[1][0]:	#two pairs and they have different end ids.
				if contig_pos_ls[0][1] == contig_pos_ls[1][1]:	#on the same target contig
					dist = abs( contig_pos_ls[0][2] - contig_pos_ls[1][2])
					if dist>0:
						dist = math.log10(dist)
					else:
						dist = -1
					same_contig_BAC_end_pair_dist_ls.append(dist)
				else:
					on_different_contig_BAC_end_contig_pos_ls.append(BAC_id)
		sys.stderr.write("%s BACs have their ends on the same contig. %s on different contigs.\n"%(len(same_contig_BAC_end_pair_dist_ls), \
																					len(on_different_contig_BAC_end_contig_pos_ls)))
		mismatch_rate = float(processor.no_of_mismatches)/processor.no_of_aligned_bases
		insertion_rate = float(processor.no_of_insertions)/processor.no_of_aligned_bases
		deletion_rate = float(processor.no_of_deletions)/processor.no_of_aligned_bases
		
		sys.stderr.write("Mismatch rate: %.4f.\n"%(mismatch_rate))
		sys.stderr.write("Insertion rate: %.4f.\n"%(insertion_rate))
		sys.stderr.write("Deletion rate: %.4f.\n"%(deletion_rate))
		
		if len(same_contig_BAC_end_pair_dist_ls)>10:
			yh_matplotlib.drawHist(same_contig_BAC_end_pair_dist_ls, title="distance of (%s) BAC-end pairs on the same contig"%\
							(len(same_contig_BAC_end_pair_dist_ls)), \
							xlabel_1D="distance of 2 BAC ends", outputFname='%s_same_contig_BAC_end_dist.png'%(outputFnamePrefix), \
							min_no_of_data_points=10, needLog=True)
		
		match_percentage_ls = processor.qname2match_percentage.values()
		yh_matplotlib.drawHist(match_percentage_ls, title="histogram of percentage of matches of %s BAC ends"%(len(match_percentage_ls)), \
							xlabel_1D="match percentage", outputFname='%s_match_percentage.png'%(outputFnamePrefix), \
							min_no_of_data_points=50, needLog=True)
		#bamOutputF.close()
		return contig_id_set
	"""
		#2011-7-1
		#inputPrefix = os.path.expanduser("~/script/vervet/data/ref/BAC/vsTop150Contigs_by_bwasw/last_updated_vervet_GSS_seqs.sorted")
		inputPrefix = os.path.expanduser("~/script/vervet/data/ref/BAC/vs_MinSize200Scaffolds_by_bwasw/last_updated_vervet_GSS_seqs.sorted")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		
		min_alen_rlen_ratio=0.95
		outputFnamePrefix= os.path.expanduser("%s.minAlenRatio%s"%(inputPrefix, min_alen_rlen_ratio))
		VariantDiscovery.getBACEndHitsOnContigs(inputFname=inputFname, outputFnamePrefix=outputFnamePrefix, min_alen_rlen_ratio=min_alen_rlen_ratio)
		sys.exit(0)
		
		#2011-7-1
		#inputPrefix = os.path.expanduser("~/script/vervet/data/ref/BAC/vsTop150Contigs_by_bwasw/last_updated_vervet_GSS_seqs.sorted")
		inputPrefix = os.path.expanduser("~/script/vervet/data/ref/BAC/vs_top150Scaffolds_by_bwasw_bwtsw/last_updated_vervet_GSS_seqs.sorted")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		
		min_alen_rlen_ratio=0.95
		outputFnamePrefix= os.path.expanduser("%s.minAlenRatio%s"%(inputPrefix, min_alen_rlen_ratio))
		contig_id_set = VariantDiscovery.getBACEndHitsOnContigs(inputFname=inputFname, outputFnamePrefix=outputFnamePrefix, min_alen_rlen_ratio=min_alen_rlen_ratio)
		sys.exit(0)
		
		
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
	def discoverHetsFromBAM(cls, inputFname, outputFname, refFastaFname=None, monomorphicDiameter=100, \
						maxNoOfReads=300, minNoOfReads=2, minMinorAlleleCoverage=3, maxMinorAlleleCoverage=7,\
						maxNoOfReadsForGenotypingError=1, maxMajorAlleleCoverage=30, maxNoOfReadsForAllSamples=1000,\
						nt_set = set(['a','c','g','t','A','C','G','T'])):
		"""
		2011-7-18
			add argument refFastaFname
			ref is at column 0. no read_group should bear the name "ref".
		2011-7-8
			it discovers homozygous SNPs as well.
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
		
		read_group2col_index = {'ref':0}	#ref is at column 0. "ref" must not be equal to any read_group.
		locus_id2row_index = {}
		data_matrix = []
		
		tid2refName = {}	#dictionary storing the target references which have SNP calls
		for pileupcolumn in samfile.pileup():
			#print
			#print 'coverage at base %s %s = %s'%(pileupcolumn.tid, pileupcolumn.pos , pileupcolumn.n)
			counter += 1
			
			current_locus = '%s_%s'%(pileupcolumn.tid, pileupcolumn.pos+1)
			if pileupcolumn.tid not in tid2refName:
				tid2refName[pileupcolumn.tid] = samfile.getrname(pileupcolumn.tid)
			
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
					if pileupread.qpos<0 or pileupread.qpos>=len(pileupread.alignment.seq):	#2011-7-13 need to investigate what happens here??
						continue	#
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
		
		#2011-7-18 read in the reference sequences in order to find out the ref base
		refNameSet = set(tid2refName.values())
		tid2Seq = {}
		from Bio import SeqIO
		handle = open(refFastaFname, "rU")
		for record in SeqIO.parse(handle, "fasta"):
			contig_id = record.id.split()[0]
			if contig_id in refNameSet:
				tid = samfile.gettid(contig_id)
				tid2Seq[tid] = record.seq
			if len(tid2Seq)>=len(refNameSet):	#enough data, exit.
				break
		handle.close()
		samfile.close()
		
		# output the matrix in the end
		read_group_col_index_ls = read_group2col_index.items()
		read_group_col_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		header = ['locus_id', 'locus_id']+[row[0] for row in read_group_col_index_ls]
		writer.writerow(header)
		
		locus_id_and_row_index_ls = locus_id2row_index.items()
		locus_id_and_row_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		for i in xrange(len(locus_id_and_row_index_ls)):
			locus_id, row_index = locus_id_and_row_index_ls[i]
			locus_id_split = locus_id.split('_')[:2]
			tid, pos = map(int, locus_id_split)[:2]
			refSeq = tid2Seq[tid]
			refBase = refSeq[pos-1]
			data_row = data_matrix[i]
			data_row[0] = refBase	#2011-7-18
			# if data_row is shorter than read_group_col_index_ls, add "NA" to fill it up
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
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/454_illu_6_sub_vs_1MbBAC')
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
	def calculatePairwiseDistanceOutOfSNPXStrainMatrix(cls, inputFname, outputFname, snpFnameToDoFiltering=None, \
											convertHetero2NA=False,
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
		if snpFnameToDoFiltering:
			oldSNPData = SNPData(input_fname=snpFnameToDoFiltering, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			
			snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
			if min_MAF>0 and min_MAF<0.5:
				snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
			col_id_ls = snpData.col_id_ls
			
			
			oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			snpData = SNPData.keepColsByColID(snpData, col_id_ls)
			
		else:
			oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
								input_alphabet=1, turn_into_integer=1)
			snpData = transposeSNPData(oldSNPData)
			if convertHetero2NA:
				snpData = SNPData.convertHetero2NA(snpData)
			
			snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
			if min_MAF>0 and min_MAF<0.5:
				snpData = snpData.removeColsByMAF(snpData, min_MAF=min_MAF)
		
		# add outputFname to function below to output the row pairwise distance
		row_id2pairwise_dist_ls = snpData.calRowPairwiseDist(assumeBiAllelic=True, outputFname=outputFname)
		cls.outputRow_id2pairwise_dist_lsInMatrix(row_id2pairwise_dist_ls, outputFname)
	
	"""
		#2011-4-7
		inputFname = os.path.expanduser('~/script/vervet/data/topConservedHG19_random_100loci_vervetSNP.tsv')
		
		inputFname = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC_minMAC3_maxMAC7_maxNoOfReadsForGenotypingError1.maxCoverage20.maxMajorAC10.maxNoOfReadsForAllSamples1000.snps')
		snpFnameToDoFiltering = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC_minMAC3_maxMAC7_maxNoOfReadsForGenotypingError1.maxCoverage20.maxMajorAC10.maxNoOfReadsForAllSamples1000_no_carribean.snps')
		snpFnameToDoFiltering = None
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, snpFnameToDoFiltering=snpFnameToDoFiltering,\
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
		
		#2011-7-8
		inputPrefix = os.path.expanduser("/usr/local/vervetData/ref/454/vs_MinSize200Scaffolds_by_bwasw")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		referenceIDSet = set(['Contig0',])
		outputFname = os.path.expanduser("%s.Contig0.bam"%(inputPrefix))
		VariantDiscovery.filterAlignmentByReferenceIDs(inputFname, outputFname, referenceIDSet=referenceIDSet, readGroup='454_vs_Contig0')
		sys.exit(0)
		
				
		#2011-4-7
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110712T234818-0700/8_genomes_vs_top156References_call.tsv'
		inputFname = '/usr/local/vervetData/vervetPipeline/outputs/call/vervet_path2.call'
		inputFname = '/usr/local/vervetData/vervetPipeline/7GenomeVsVervet1MbAsRef/7_genomes_vs_top1References_call.tsv'
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110714T013718-0700/call/Contig1.call'
		inputFname = '/usr/local/vervetData/vervetPipeline/top2Contigs8Genome/8_genomes_vs_top2References_call.tsv'
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110714T015458-0700/call/Contig0.call'
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110714T015458-0700/call/Contig1.call'
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110712T234818-0700/call/Contig110.call'
		inputFname = '/usr/local/vervetData/vervetPipeline/work/outputs/crocea/pegasus/AlignmentToCallPipeline/20110719T011659-0700/call/Contig0.call'
		inputFname = os.path.expanduser("~/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC.GATK.call")
		snpFnameToDoFiltering = None
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, snpFnameToDoFiltering=snpFnameToDoFiltering,\
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
		
		
		#2011-3-24
		common_prefix = os.path.expanduser('~/script/vervet/data/8_genome_vs_Contig0.RG')
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
		
		#2011-4-7
		inputFname = '%s_minMAC%s_maxMAC%s_maxNoOfReadsForGenotypingError%s.maxCoverage%s.maxMajorAC%s.hets'%(common_prefix, \
								minMinorAlleleCoverage, maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError, maxNoOfReads,\
								maxMajorAlleleCoverage)
		snpFnameToDoFiltering = None
		convertHetero2NA = True
		max_NA_rate = 0.4
		min_MAF = 0
		outputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(inputFname)[0], \
												convertHetero2NA, min_MAF, max_NA_rate)
		VariantDiscovery.calculatePairwiseDistanceOutOfSNPXStrainMatrix(inputFname, outputFname, snpFnameToDoFiltering=snpFnameToDoFiltering,\
							convertHetero2NA=convertHetero2NA, min_MAF=min_MAF, max_NA_rate=max_NA_rate)
		sys.exit(0)
	"""
	
	@classmethod
	def outputRow_id2pairwise_dist_lsInMatrix(cls, row_id2pairwise_dist_ls, outputFname):
		"""
		2011-5-14
			split out of calculatePairwiseDistanceOutOfSNPXStrainMatrix()
			it appends to outputFname, rather than overwriting.
		"""
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
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/454_illu_6_sub_vs_1MbBAC')
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
		row_id2pairwise_dist_ls = snpData.calRowPairwiseDist(assumeBiAllelic=False, outputFname=outputFname)
		
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
		inputFname = outputFname
		max_NA_rate = 0.4
		min_het_frequency = 0.2
		outputFnamePrefix = '%s.max_NA_rate%s.min_het_frequency%s'%(os.path.splitext(inputFname)[0], max_NA_rate, min_het_frequency)
		VariantDiscovery.studyHeterozygousCallOutOfSNPXStrainMatrix(inputFname, outputFnamePrefix,\
												max_NA_rate=max_NA_rate, min_het_frequency=min_het_frequency)
		sys.exit(0)
		"""
	
	@classmethod
	def filterSNPXStrainMatrixToRetainHetsCommonInSelectedIndividuals(cls, inputFname, outputFnamePrefix,
										max_NA_rate=0.3, selected_accession_id_ls=['VRC_ref_GA', 'Barbados_GA'], report=True):
		"""
		2011-5-14
			Steps:
				1. filter out SNPs with too many NAs
				2. retain loci that are have the same heterozygous call among all individuals in selected_accession_id_ls
				3. calculate pairwise distance ()
				4. output data
		"""
		from pymodule.SNP import SNPData, transposeSNPData
		oldSNPData = SNPData(input_fname=inputFname, turn_into_array=1, ignore_2nd_column=1,\
							input_alphabet=1, turn_into_integer=1)
		snpData = transposeSNPData(oldSNPData)
		
		snpData = snpData.removeColsByNARate(snpData, max_NA_rate=max_NA_rate)
		
		"""
		# filter out any non-Het calls
		sys.stderr.write("Marking all non-heterozygous calls as missing")
		import numpy as num
		newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=snpData.col_id_ls)
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		no_of_hets = 0
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				if snpData.data_matrix[i][j]>4:
					newSnpData.data_matrix[i][j] = snpData.data_matrix[i][j]
					no_of_hets += 1
		sys.stderr.write("%s heterozygous calls. Done.\n"%no_of_hets)
		snpData = newSnpData
		"""
		
		# retain loci at which selected_accession_id_ls have the same heterozygous call
		
		selected_row_index_ls = []
		for accession_id in selected_accession_id_ls:
			selected_row_index_ls.append(snpData.row_id2row_index[accession_id])
		
		no_of_rows = len(snpData.data_matrix)
		no_of_cols = len(snpData.data_matrix[0])
		allele2count_ls = []
		col_id_to_be_kept_ls = []
		no_of_SNPs_with_more_than_one_distinct_hets = 0
		for j in range(no_of_cols):
			col_id = snpData.col_id_ls[j]
			allele2count_ls.append({})
			for i in selected_row_index_ls:
				allele = snpData.data_matrix[i][j]
				if allele>4:
					if allele not in allele2count_ls[j]:
						allele2count_ls[j][allele] = 0
					allele2count_ls[j][allele] += 1
					#if cls.report and allele_index>1:
					#sys.stderr.write("%s (more than 2) alleles at SNP %s (id=%s).\n"%((allele_index+1), j, snpData.col_id_ls[j]))
			if len(allele2count_ls[j])>1:
				no_of_SNPs_with_more_than_one_distinct_hets += 1
				if report:
					sys.stderr.write("Warning: more than one hets at SNP %s (id=%s), %s.\n"%(j, snpData.col_id_ls[j], \
																					repr(allele2count_ls[j])))
			elif len(allele2count_ls[j])==1:	#only one type of heterozygous call
				numberOfAccessionWithThatHet = sum(allele2count_ls[j].values())
				if numberOfAccessionWithThatHet == len(selected_row_index_ls):	# all selected accessions bear the same call.
					col_id_to_be_kept_ls.append(col_id)
		
		snpData = SNPData.keepColsByColID(snpData, col_id_to_be_kept_ls)
		snpData.no_of_cols_removed = no_of_cols - len(col_id_to_be_kept_ls)
		sys.stderr.write("%s columns filtered out by Het-sharing and %s SNPs with >1 different hets. Done.\n"%\
						(snpData.no_of_cols_removed, no_of_SNPs_with_more_than_one_distinct_hets))		
		
		outputFnamePrefix = '%s.maxNARate%s.HetSharedIn_%s'%(outputFnamePrefix, max_NA_rate, '_'.join(selected_accession_id_ls))
		snpData.tofile('%s.tsv'%outputFnamePrefix)
		
		outputFname = '%s.pairwiseDistBasedOnHet.tsv'%(outputFnamePrefix)
		row_id2pairwise_dist_ls = snpData.calRowPairwiseDist(assumeBiAllelic=False, outputFname=outputFname)
		cls.outputRow_id2pairwise_dist_lsInMatrix(row_id2pairwise_dist_ls, outputFname)
		
		"""
		#2011-5-14
		common_prefix = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC')
		#inputFname = '%s.bam'%(common_prefix)
		minMinorAlleleCoverage=3
		maxMinorAlleleCoverage=7
		maxNoOfReadsForGenotypingError=1
		maxNoOfReads = 20
		maxMajorAlleleCoverage=10
		maxNoOfReadsForAllSamples = 1000
		inputFname = '%s_minMAC%s_maxMAC%s_maxNoOfReadsForGenotypingError%s.maxCoverage%s.maxMajorAC%s.maxNoOfReadsForAllSamples%s.snps'%(common_prefix, \
								minMinorAlleleCoverage, maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError, maxNoOfReads,\
								maxMajorAlleleCoverage, maxNoOfReadsForAllSamples)
		
		max_NA_rate = 0.4
		outputFnamePrefix = os.path.splitext(inputFname)[0]
		VariantDiscovery.filterSNPXStrainMatrixToRetainHetsCommonInSelectedIndividuals(inputFname, outputFnamePrefix,
										max_NA_rate=max_NA_rate, selected_accession_id_ls=['VRC_ref_GA', 'Barbados_GA'], report=True)
		sys.exit(0)
		
		
		#2011-5-14
		inputFolder = os.path.expanduser('~/script/vervet/data/1MbBAC_as_ref/')
		#inputFname = '%s.bam'%(common_prefix)
		minMinorAlleleCoverage=3
		maxMinorAlleleCoverage=7
		maxNoOfReadsForGenotypingError=1
		maxNoOfReads = 20
		maxMajorAlleleCoverage=10
		maxNoOfReadsForAllSamples = 1000
		inputFname = '454_illu_6_sub_vs_1MbBAC_minMAC%s_maxMAC%s_maxNoOfReadsForGenotypingError%s.maxCoverage%s.maxMajorAC%s.maxNoOfReadsForAllSamples%s.snps'%\
				(minMinorAlleleCoverage, maxMinorAlleleCoverage, maxNoOfReadsForGenotypingError, maxNoOfReads,\
				maxMajorAlleleCoverage, maxNoOfReadsForAllSamples)
		inputFname = os.path.join(inputFolder, inputFname)
		
		inputFolder = os.path.expanduser('~/script/vervet/data/hg19_as_ref/')
		inputFname = os.path.join(inputFolder, '8_alns_12345687_1000s_loci_phastCons_0.9_1_minMAC3_maxMAC7_maxNoOfReadsForGenotypingError1.maxCoverage20.maxMajorAC10.maxNoOfReadsForAllSamples480.snps')
		max_NA_rate = 0.4
		outputFnamePrefix = os.path.splitext(inputFname)[0]
		VariantDiscovery.filterSNPXStrainMatrixToRetainHetsCommonInSelectedIndividuals(inputFname, outputFnamePrefix,
										max_NA_rate=max_NA_rate, selected_accession_id_ls=['VRC_ref_GA', 'Barbados_GA'], report=True)
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
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script//vervet/data/ref/illumina/ref_illu_vs_1MbBAC_by_aln")
		inputPrefix = os.path.expanduser("~/mnt/hoffman2_home/script//vervet/data/ref/illumina/ref_illu_vs_hg9_by_aln.3eQTL")
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_hg19/454_vs_hg19.3eQTL'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_BAC/454_vs_ref_1MbBAC.F4'
		inputPrefix = '/Network/Data/vervet/ref/454/454_vs_BAC/454_vs_ref_1MbBAC.F4.minPerBaseAS0.5.minMapQ125.score2'
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		outputFname = os.path.expanduser("%s.coverage.hist.png"%(inputPrefix))
		VariantDiscovery.drawCoverageHistFromBAM(inputFname, outputFname)
		sys.exit(2)
	"""

	@classmethod
	def PCAContigByDistVectorFromOtherGenomes(cls, inputDir, outputFnamePrefix):
		"""
		2011-8-2
			Run PCA on the vervet ref contigs based on distance vector of other genomes to these contigs
			
			input file name looks like this:
			
				workflow_8GenomeVsTop156Contigs_GATK/call/Contig9.pairwiseDist.convertHetero2NATrue.minMAF0.maxNA0.4.tsv
			
			output the data in a matrix fashion that the web MotionChartAppMCPanel app would recognize 
		"""
		vectorData = cls.readContigDistVector(inputDir)

		# run PCA
		sys.stderr.write("Carrying out contig-wise PCA ...")
		phenotypePCA_fname = '%s_VRefContig.tsv'%outputFnamePrefix
		phenotypePCA_writer = csv.writer(open(phenotypePCA_fname, 'w'), delimiter='\t')
		
		import pca_module
		from pymodule.PCA import PCA
		#T, P, explained_var = pca_module.PCA_svd(phenData_trans.data_matrix, standardize=True)
		T, P, explained_var = PCA.eig(vectorData.data_matrix, normalize=False)	#normalize=True causes missing value in the covariance matrix
		# get the category information for each phenotype
		header = ['ContigID', 'DummyTime', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		phenotypePCA_writer.writerow(header)
		for i in range(len(vectorData.row_id_ls)):
			row_id = vectorData.row_id_ls[i]
			data_row = [row_id, '2011'] + list(T[i,0:6])
			phenotypePCA_writer.writerow(data_row)
		del phenotypePCA_writer
		sys.stderr.write("Done.\n")
		
	"""
		#2011-8-2
		inputDir = '/usr/local/vervetData/vervetPipeline/workflow_8GenomeVsTop156Contigs_GATK/call/'
		outputFnamePrefix = '/usr/local/vervetData/vervetPipeline/workflow_8GenomeVsTop156Contigs_GATK/contigPCAByDistVector'
		VariantDiscovery.PCAContigByDistVectorFromOtherGenomes(inputDir, outputFnamePrefix)
		sys.exit(3)
	"""
	
	@classmethod
	def readContigDistVector(cls, inputDir, refID='ref'):
		"""
		2011-10-18
			the xlabel will now distinguish the VRC_ref_454 and VRC_ref_GA
		2011-8-2
			used by PCAContigByDistVectorFromOtherGenomes() and drawContigByDistVectorFromOtherGenomes()
		"""
		from pymodule import SNPData
		import numpy, csv, os, sys, re
		contig_id_pattern = re.compile(r'Contig(\d+).*')
		row_id_ls = []
		col_id_ls = []
		data_matrix = []
		for fname in os.listdir(inputDir):
			#sys.stderr.write("%s ..."%fname)
			inputFname = os.path.join(inputDir, fname)
			contig_id_pattern_sr = contig_id_pattern.search(inputFname)
			if contig_id_pattern_sr:
				contig_id = contig_id_pattern_sr.group(1)
			else:
				contig_id = os.path.splitext(os.path.split(inputFname)[1])[0]
			row_id_ls.append(contig_id)
			reader = csv.reader(open(inputFname, ), delimiter='\t')
			matrixStart = False
			for row in reader:
				if row[0]=='':
					if not col_id_ls:
						col_id_ls = row[1:]
					matrixStart = True
				if matrixStart and row[0]==refID:
					data_row = row[1:]
					data_row = map(float, data_row)
					data_matrix.append(data_row)
			#sys.stderr.write("Done.\n")
		data_matrix = numpy.array(data_matrix)
		take1stSplit = lambda x: x.split("_GA_vs_")[0]
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			newXLabel = take1stSplit(col_id)
			if newXLabel[:7] == 'VRC_ref':
				newXLabel = col_id.split('_vs_')[0]	#a different split for VRC_ref to retain the sequencer
			col_id_ls[i] = newXLabel
		vectorData = SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
		sys.stderr.write("%s rows and %s columns.\n"%(data_matrix.shape[0], data_matrix.shape[1]))
		return vectorData
	
	@classmethod
	def reorderGenomeDistVector(cls, data_row, subspeciesName2index=None, subspeciesNameLs=[],\
							target_ID_order=['ref', 'VRC_ref_454', 'VRC_ref_GA', 'Barbados']):
		"""
		2011-10-13
			arrange data entries in data_row so that target_ID2index is met first and then everything else
		"""
		new_data_row = []
		indexUsedSet = set()
		new_subspeciesName2index = {}
		for target_ID in target_ID_order:
			indexOfTarget = subspeciesName2index.get(target_ID)
			if indexOfTarget is not None:
				indexUsedSet.add(indexOfTarget)
				new_subspeciesName2index[target_ID] = len(new_data_row)
				new_data_row.append(data_row[indexOfTarget])
		
		for i in xrange(len(data_row)):
			if i not in indexUsedSet:
				subspeciesName = subspeciesNameLs[i]
				new_subspeciesName2index[subspeciesName] = len(new_data_row)
				new_data_row.append(data_row[i])
		return new_data_row, new_subspeciesName2index
	
	@classmethod
	def drawContigByDistVectorFromOtherGenomes(cls, inputDir, outputFnamePrefix, refID="sabaeus_GA_vs_top156Contigs",\
											subspeciesName='sabaeus'):
		"""
		2011-8-2
			This program draws the distance vector of other genomes to each contig.
			
			The inputDir contains output by CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py
		"""
		vectorData = cls.readContigDistVector(inputDir, refID=refID)
		
		# run PCA
		sys.stderr.write("Drawing all distance vectors ...")
		outputFname = '%s.png'%outputFnamePrefix
		
		import pylab
		pylab.clf()
		target_ID_order=['ref', 'VRC_ref_454', 'VRC_ref_GA', 'Barbados']
		xlabel_ls, new_subspeciesName2index = cls.reorderGenomeDistVector(vectorData.col_id_ls, \
									subspeciesName2index=vectorData.col_id2col_index,\
									subspeciesNameLs=vectorData.col_id_ls, target_ID_order=target_ID_order)
		#skip the column for the subspeciesName itself.
		indexOfrefID = new_subspeciesName2index[subspeciesName]
		xlabel_ls = xlabel_ls[0:indexOfrefID] + xlabel_ls[indexOfrefID+1:]
		#xlabel_ls = vectorData.col_id_ls[1:3] + [vectorData.col_id_ls[0]] + vectorData.col_id_ls[3:6] + vectorData.col_id_ls[7:9]
		
		for i in range(len(xlabel_ls)):
			xlabel = xlabel_ls[i]
			if xlabel[:7] == 'VRC_ref':
				xlabel = xlabel[4:]
			xlabel_ls[i] = xlabel
		
		for data_row in vectorData.data_matrix:
			data_row = list(data_row)
			#new_data_row = data_row[1:3] + [data_row[0]] + data_row[3:6] + data_row[7:9]
			#skip the column for the subspeciesName itself.
			new_data_row, new_subspeciesName2index = cls.reorderGenomeDistVector(data_row, \
									subspeciesName2index=vectorData.col_id2col_index,\
									subspeciesNameLs=vectorData.col_id_ls, target_ID_order=target_ID_order)
			new_data_row = new_data_row[0:indexOfrefID] + new_data_row[indexOfrefID+1:]
			pylab.plot(range(len(new_data_row)), new_data_row)
		pylab.title("Distance vector from %s genomes to %s of %s contigs"%(len(xlabel_ls), subspeciesName, len(vectorData.row_id_ls)))
		pylab.xticks(range(len(xlabel_ls)), xlabel_ls)
		pylab.savefig(outputFname, dpi=200)
		sys.stderr.write("Done.\n")

	"""
		#2011-8-2
		inputDir = '/Network/Data/vervet/vervetPipeline/8GenomeVsTop156Contigs_GATK_all_bases/pairwiseDistMatrix/'
		outputFnamePrefix = '/Network/Data/vervet/vervetPipeline/8GenomeVsTop156Contigs_GATK_all_bases_ContigByDistVectorFrom8Genomes'
		VariantDiscovery.drawContigByDistVectorFromOtherGenomes(inputDir, outputFnamePrefix, refID='ref')
		sys.exit(0)
		
	"""

	@classmethod
	def compareContigByDistVectorFromTwoDifferentRuns(cls, inputDir1, inputDir2, outputFnamePrefix, partOfTitle=""):
		"""
		2011-8-26
			This program draws the distance vector of other genomes to each contig.
			
			The inputDir contains output by CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py
		"""
		vectorData1 = cls.readContigDistVector(inputDir1)
		vectorData2 = cls.readContigDistVector(inputDir2)
		import numpy	#numpy.linalg.norm(a-b)	#euclidean distance
		
		from scipy import spatial
		
		#>>> spatial.distance.correlation
		import rpy
		cor_ls = []
		index_ls = [0,] + range(3,6) + range(7,9)
		for row_id, row_index in vectorData1.row_id2row_index.iteritems():
			if row_id in vectorData2.row_id2row_index:
				row_index2 = vectorData2.row_id2row_index[row_id]
				data_row = vectorData1.data_matrix[row_index]
				vector1 = data_row[index_ls]
				data_row = vectorData2.data_matrix[row_index2]
				vector2 = data_row[index_ls]
				#cor = spatial.distance.correlation(vector1, vector2)
				cor = rpy.r.cor(vector1, vector2)
				cor_ls.append(cor)
		
		from pymodule import yh_matplotlib
		yh_matplotlib.drawHist(cor_ls, title="correlation of distance-2-contig vectors %s"%(partOfTitle), \
							xlabel_1D="correlation", xticks=None, outputFname="%s.png"%outputFnamePrefix,\
							min_no_of_data_points=50, needLog=False, dpi=200,)
		
	"""
		
	"""
	
	
	@classmethod
	def getContigTsTvList(cls, inputDir, suffix='.TiTv.100000.log'):
		"""
		2011-9-23
			used by drawContigTsTvHistogram
		"""
		from pymodule import SNPData
		import csv, os, sys, re
		contig_id_pattern = re.compile(r'Contig(\d+).*')
		row_id_ls = []
		TsTvLs = []
		for fname in os.listdir(inputDir):
			if fname.find(suffix)!=-1:
				
				sys.stderr.write("%s ..."%fname)
				inputFname = os.path.join(inputDir, fname)
				contig_id_pattern_sr = contig_id_pattern.search(inputFname)
				if contig_id_pattern_sr:
					contig_id = contig_id_pattern_sr.group(1)
				else:
					contig_id = os.path.splitext(os.path.split(inputFname)[1])[0]
				row_id_ls.append(contig_id)
				inf = open(inputFname)
				for line in inf:
					if line.find('Ts/Tv ratio: ')==0:
						row = line.strip().split(":")
						TsTv = row[-1].strip()
						TsTv = float(TsTv)
						TsTvLs.append(TsTv)
				sys.stderr.write("Done.\n")
		
		return TsTvLs
	
	@classmethod
	def drawContigTsTvHistogram(cls, inputDir, outputFname):
		"""
		2011-9-23
		"""
		TsTvLs = cls.getContigTsTvList(inputDir, suffix='.TiTv.100000.log')
		from pymodule.yh_matplotlib import drawHist
		drawHist(TsTvLs, title="hist of TsTv of %s contigs"%(len(TsTvLs)), \
				xlabel_1D="TsTv", xticks=None, outputFname=outputFname, min_no_of_data_points=20, needLog=False, \
				dpi=200)
	
	"""
		#2011-9-23
		inputDir = 
		outputFname
		VariantDiscovery.drawHist(inputDir, outputFname)
		sys.exit(0)
	"""
	
	@classmethod
	def getContig2Locus2Frequency(cls, inputDir, VCFOutputType=2):
		"""
		2011-9-21
		
		"""
		from pymodule import SNPData
		import numpy, csv, os, sys, re, gzip
		from pymodule.utils import runLocalCommand, getColName2IndexFromHeader
		from GenotypeCallByCoverage import GenotypeCallByCoverage
		contig_id_pattern = re.compile(r'Contig(\d+).*')
		contig2locus2frequency = {}
		for fname in os.listdir(inputDir):
			if fname[-6:]!='vcf.gz' and fname[-3:]!='vcf':
				continue
			sys.stderr.write("%s ..."%fname)
			inputFname = os.path.join(inputDir, fname)
			contig_id_pattern_sr = contig_id_pattern.search(inputFname)
			if contig_id_pattern_sr:
				contig_id = contig_id_pattern_sr.group(1)
			else:
				contig_id = os.path.splitext(os.path.split(inputFname)[1])[0]
			if contig_id not in contig2locus2frequency:
				contig2locus2frequency[contig_id] = {}
			
			if inputFname[-3:]=='.gz':
				import gzip
				inf = gzip.open(inputFname)
			else:
				inf = open(inputFname)
			reader = csv.reader(inf, delimiter='\t')
			
			counter = 0
			real_counter = 0
			
			
			for row in reader:
				if row[0] =='#CHROM':
					row[0] = 'CHROM'	#discard the #
					header = row
					col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
					#individual_name2col_index = GenotypeCallByCoverage.getIndividual2ColIndex(header, col_name2index)
					continue
				elif row[0][0]=='#':	#2011-3-4
					continue
				"""
				if chr_number_pattern.search(row[0]):
					chr = chr_number_pattern.search(row[0]).group(1)
				elif chr_pure_number_pattern.search(row[0]):
					chr = chr_pure_number_pattern.search(row[0]).group(1)
				else:
					sys.stderr.write("Couldn't parse the chromosome number/character from %s.\n Exit.\n"%(row[0]))
					sys.exit(4)
				"""
				chr = row[0]
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
						#sys.stderr.write("Error in splitting %s by =.\n"%info)	###Error in splitting DS by =.
						continue
					info_tag2value[tag] = value
				
				current_locus = '%s_%s'%(chr, pos)
				refBase = row[col_name2index['REF']]
				altBase = row[col_name2index['ALT']]
				if "AF1" in info_tag2value:
					AF1 = info_tag2value.get("AF1")
				else:
					AF1 = info_tag2value.get("AF")
				if AF1:
					AF1 = float(AF1)
					contig2locus2frequency[contig_id][current_locus] = AF1
			
			sys.stderr.write("%s loci. Done.\n"%(len(contig2locus2frequency[contig_id])))
		return contig2locus2frequency
	
	@classmethod
	def isLocusPolymorphicBasedOnAAF(cls, frequency, frequencyDelta=0.001):
		"""
		2011-9-28
		"""
		distanceTo1 = abs(frequency-1.0)
		distanceTo0 = abs(frequency-0.0)
		if distanceTo0<frequencyDelta or distanceTo1<frequencyDelta:
			return False
		else:
			return True
	
	@classmethod
	def compareAlternativeAlleleFrequencyOfTwoCallSets(cls, inputDir1, inputDir2, outputDir, min_no_of_data_points=10, \
													xlabel='', ylabel='', dpi=300, frequencyDelta=0.001):
		"""
		2011-9-21
			purpose is to compare alternative allele frequency and get to know what portion is polymorphic.
		
			argument frequencyDelta is used to judge whether one frequency is close to 0 or 1, which essentially means
				these loci are not polymorphic.
		"""
		contig2locus2frequency_1 = cls.getContig2Locus2Frequency(inputDir1)
		contig2locus2frequency_2 = cls.getContig2Locus2Frequency(inputDir2)
		
		import os, sys
		from pymodule import yh_matplotlib
		if not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		for contig in contig2locus2frequency_1:
			if contig in contig2locus2frequency_2:
				locus2frequency_1 = contig2locus2frequency_1.get(contig)
				locus2frequency_2 = contig2locus2frequency_2.get(contig)
				outputFname = os.path.join(outputDir, '%s_AF1.png'%contig)
				x_ls = []
				y_ls = []
				AAF_ls_1 = []
				AAF_ls_2 = []
				locus_set = set(locus2frequency_1.keys())|set(locus2frequency_2.keys())
				shared_locus_set = set(locus2frequency_1.keys())&set(locus2frequency_2.keys())
				no_of_loci_Dir1 = len(locus2frequency_1.keys())
				no_of_loci_Dir2 = len(locus2frequency_2.keys())
				no_of_total_loci = len(locus_set)
				no_of_shared_loci = no_of_loci_Dir1 + no_of_loci_Dir2 - no_of_total_loci
				
				shared_polymorphic_locus_set = set()
				for locus in locus_set:
					frequency_1 = locus2frequency_1.get(locus)
					frequency_2 = locus2frequency_2.get(locus)
					if frequency_1 is None:
						frequency_1 = 0
					if frequency_2 is None:
						frequency_2 = 0
					locusPolymorphicInPop1 = cls.isLocusPolymorphicBasedOnAAF(frequency_1, frequencyDelta=frequencyDelta)
					locusPolymorphicInPop2 = cls.isLocusPolymorphicBasedOnAAF(frequency_2, frequencyDelta=frequencyDelta)
					if locusPolymorphicInPop1 and locusPolymorphicInPop2:
						x_ls.append(frequency_1)
						y_ls.append(frequency_2)
						shared_polymorphic_locus_set.add(locus)
					elif locusPolymorphicInPop1:
						AAF_ls_1.append(frequency_1)
					elif locusPolymorphicInPop2:
						AAF_ls_2.append(frequency_2)
				
				no_of_shared_loci = len(shared_polymorphic_locus_set)
				no_of_loci_Dir1 = len(AAF_ls_1) + no_of_shared_loci
				no_of_loci_Dir2 = len(AAF_ls_2) + no_of_shared_loci
				no_of_total_loci = no_of_loci_Dir1 + no_of_loci_Dir2 - no_of_shared_loci
				
				
				import pylab
				pylab.clf()
				no_of_data_points = len(x_ls)
				if no_of_data_points>=min_no_of_data_points:
					pylab.plot(x_ls, y_ls, ".",)
					pylab.title('Contig%s %s shared & %s total loci'%(contig, no_of_shared_loci, no_of_total_loci))
					if xlabel:
						pylab.xlabel("%s loci in %s"%(no_of_loci_Dir1, xlabel))
					if ylabel:
						pylab.ylabel("%s loci in %s"%(no_of_loci_Dir2, ylabel))
					pylab.savefig(outputFname, dpi=dpi)
				
				outputFname = os.path.join(outputDir, '%s_AAF_hist_SNPs_unique_in_%s.png'%(contig, xlabel))
				yh_matplotlib.drawHist(AAF_ls_1, title="%s SNPs on contig%s unique in %s"%(no_of_loci_Dir1- no_of_shared_loci, contig, xlabel), \
									xlabel_1D="Alternative Allele Frequency", xticks=None, outputFname=outputFname, min_no_of_data_points=20, needLog=False, \
									dpi=200)
				outputFname = os.path.join(outputDir, '%s_AAF_hist_SNPs_unique_in_%s.png'%(contig, ylabel))
				yh_matplotlib.drawHist(AAF_ls_2, title="%s SNPs on contig%s unique in %s"%(no_of_loci_Dir2- no_of_shared_loci, contig, ylabel), \
									xlabel_1D="Alternative Allele Frequency", xticks=None, outputFname=outputFname, min_no_of_data_points=20, needLog=False, \
									dpi=200)
	
	"""
		#2011-9-21
		inputDir1 = "/Network/Data/vervet/vervetPipeline/AlignmentToCallPipeline_10VWP_627_629_650_656_vs_524_top_156Contigs_uschpc/call/"
		inputDir2 = "/Network/Data/vervet/vervetPipeline/AlignmentToCallPipeline_552_554_557_605_615_vs_524_top_156Contigs_uschpc_bugfix/call/"
		outputDir = "/Network/Data/vervet/vervetPipeline/10VWP_vs_552_554_557_605_615/"
		VariantDiscovery.compareAlternativeAlleleFrequencyOfTwoCallSets(inputDir1, inputDir2, \
				outputDir, min_no_of_data_points=10, \
				xlabel="VWP", ylabel='VRC', dpi=200)
		sys.exit(0)
		
		#2011-9-21
		inputDir1 = "/Network/Data/vervet/vervetPipeline/AlignmentToCallPipeline_4_8_vs_524_top_156Contigs_uschpc/call/"
		inputDir2 = "/Network/Data/vervet/vervetPipeline/AlignmentToCallPipeline_552_554_557_605_615_vs_524_top_156Contigs_uschpc_bugfix/call/"
		outputDir = "/Network/Data/vervet/vervetPipeline/5AfricanSubspecies_vs_552_554_557_605_615/"
		VariantDiscovery.compareAlternativeAlleleFrequencyOfTwoCallSets(inputDir1, inputDir2, \
				outputDir, min_no_of_data_points=10, \
				xlabel="5AfricanSubspecies", ylabel='VRC', dpi=200)
		sys.exit(0)
		
		#2011-9-21
		callFolder = 'call'
		callFolder = 'samtools'
		#callFolder = 'gatk'
		workDir1 = 'AlignmentToCallPipeline_10VWP_627_629_650_656_vs_524_top_156Contigs_condorpool_20110922T1837/call'
		xlabel = "VWP"
		workDir1 = "AlignmentToCallPipeline_5StKitts_isq_121-123_125_129_vs_524_top_156Contigs_condor_20110926T1323/%s"%(callFolder)
		xlabel = "5StKitts_%s"%callFolder
		workDir1 = "AlignmentToCallPipeline_5Nevis_isq_124_126-128_130_vs_524_top_156Contigs_condor_20110926T1325/%s"%(callFolder)
		xlabel="5Nevis_%s"%(callFolder)
		#workDir1 = "AlignmentToCallPipeline_555_559_589_vs_524_top_156Contigs_uschpc_bugfix/call"
		#xlabel = "VRC_555_559-589"
		inputDir1 = "/Network/Data/vervet/vervetPipeline/%s/"%(workDir1)
		
		workDir2 = "AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/%s"%(callFolder)
		ylabel = 'AllVRC_Barbados_%s'%(callFolder)
		#workDir2 = "AlignmentToCallPipeline_5Nevis_isq_124_126-128_130_vs_524_top_156Contigs_condor_20110926T1325/%s"%(callFolder)
		#ylabel="5Nevis_%s"%(callFolder)
		workDir2 = "AlignmentToCallPipeline_isq_id_22_31_33_36_43_48_52_56_57_64_68_70_71_81_83_vs_524_top_156Contigs_condor_20110923T1526/%s"%(callFolder)
		ylabel = "15DistantVRC_%s"%(callFolder)
		#workDir2 = "AlignmentToCallPipeline_558_616_656_vs_524_top_156Contigs_uschpc/call"
		#ylabel="VRC_558_616-656"
		inputDir2 = "/Network/Data/vervet/vervetPipeline/%s/"%(workDir2)
		
		outputDir = "/Network/Data/vervet/vervetPipeline/AAF_%s_vs_%s_polymorphic_only/"%(xlabel, ylabel)
		VariantDiscovery.compareAlternativeAlleleFrequencyOfTwoCallSets(inputDir1, inputDir2, \
				outputDir, min_no_of_data_points=10, \
				xlabel=xlabel, ylabel=ylabel, dpi=200)
		sys.exit(0)
	"""
	@classmethod
	def countHomoHetCallsForEachSampleFromVCF(cls, inputFname, outputFnamePrefix):
		"""
		2011-11-2
			given a VCF file, count the number of homo-ref, homo-alt, het calls
			
		"""
		sys.stderr.write("Count the number of homozygous-ref/alt & het from %s .\n"%(inputFname))
		from pymodule.VCFFile import VCFFile
		vcfFile = VCFFile(inputFname=inputFname)
		
		sampleID2data = {}	#key is sampleID, value is a list of 3 numbers. 'NoOfHomoRef', 'NoOfHomoAlt', 'NoOfHet'
		
		no_of_total = 0.
		minStart = None
		for vcfRecord in vcfFile.parseIter():
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			refBase = vcfRecord.data_row[0].get("GT")[0]
			
			for sample_id, sample_index in vcfFile.sample_id2index.iteritems():
				if sample_id=='ref':	#ignore the reference
					continue
				if sample_id not in sampleID2data:
					sampleID2data[sample_id] = [0, 0, 0]
				if not vcfRecord.data_row[sample_index]:	#None for this sample
					continue
				callForThisSample = vcfRecord.data_row[sample_index].get('GT')
				if not callForThisSample or callForThisSample=='NA':
					continue
				if callForThisSample[0]==refBase and callForThisSample[1]==refBase:
					#homozygous reference allele
					sampleID2data[sample_id][0]+=1
				elif callForThisSample[0]==callForThisSample[1] and callForThisSample[0]!=refBase:
					#homozygous alternative allele
					sampleID2data[sample_id][1]+=1
				elif callForThisSample[0]!=callForThisSample[1]:
					sampleID2data[sample_id][2]+=1
			
		outputFname = "%s.homoHetCountPerSample.tsv"%(outputFnamePrefix)
		import csv
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['sampleID','NoOfHomoRef', 'NoOfHomoAlt', 'NoOfHet'])
		sampleIDLs = sampleID2data.keys()
		sampleIDLs.sort()
		for sampleID in sampleIDLs:
			count_data = sampleID2data.get(sampleID)
			writer.writerow([sampleID] + count_data)
		del writer
		sys.stderr.write("Done.\n")
	"""
		#2011-11-2
		workflowName = 'AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top156Contigs_condor_20111101T2316'
		inputFname = '/Network/Data/vervet/vervetPipeline/%s/call/Contig0.vcf.gz'%(workflowName)
		outputFnamePrefix = '/tmp/%s_Contig0'%(workflowName)
		VariantDiscovery.countHomoHetCallsForEachSampleFromVCF(inputFname, outputFnamePrefix)
		sys.exit(0)
	"""

	class TallyFilteredSNPs(object):
		"""
		2011-11-21
			class to count how many genotypes are masked as missing, how many loci get filtered at particular steps.
		"""
		def __init__(self):
			pass
		
		@classmethod
		def findSitesFilteredOrMaskedByFilterVCFByDepth(cls, \
			FilterVCFByDepthStdoutFnamePattern="FilterVCF_4HighCovVRC_isq_15_18_vs_524_top7559Contigs_multi_sample_gatk_inter_minMAC2_condor.2011.11.17T1309/merge_workflow-FilterVCFByDepth-1.0_PID2_ID*out.000", \
			grepPattern="call\/Contig", sampleSize=4, reportPeriod=50):
			"""
			2011-11-21
				get the stats from the pegasus workflow stdout files
				
				the "call" before "\/Contig" is important to differentiate between two VCF input folders in FilterVCFPipeline.py
			"""
			#first get all files that match FilterVCFByDepthStdoutFnamePattern and contains the grepPattern
			from pymodule import runLocalCommand
			commandline = "grep -l %s %s"%(grepPattern, FilterVCFByDepthStdoutFnamePattern)
			#-l is to let grep to output filenames only
			return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=False)
			stdout_content = return_data.stdout_content
			import cStringIO
			import re
			noOfProcessedLociPattern = re.compile(r'FilterVCFByDepth - (\d+) loci processed.')
			noOfRetainedLociPattern = re.compile(r'FilterVCFByDepth - (\d+) loci retained.')
			noOfTotalGenotypesPattern = re.compile(r'FilterVCFByDepth - Number of total genotypes: (\d+)')
			noOfOutputGenotypesPattern = re.compile(r'FilterVCFByDepth - Number of outputted genotypes: (\d+)')
			noOfNonBiAllelicLociPattern = re.compile(r'Number of non-bi-allelic loci after GQ and Depth filter: (\d+)')
			noOfProcessedLoci = 0
			noOfRetainedLoci = 0
			noOfTotalGenotypes = 0
			noOfOutputGenotypes = 0
			noOfNonBiAllelicLoci = 0
			infWithFilenames = cStringIO.StringIO(stdout_content)
			counter = 0
			for line in infWithFilenames:
				counter += 1
				inputFname = line.strip()
				if counter%reportPeriod==0:
					sys.stderr.write("File %s: %s"%(counter, inputFname))
				inf = open(inputFname)
				for line in inf:
					noOfProcessedLoci = cls.handleNoOfSitesRemovedPattern(line, noOfProcessedLociPattern, \
															no_of_sites_removed_already=noOfProcessedLoci)
					noOfRetainedLoci = cls.handleNoOfSitesRemovedPattern(line, noOfRetainedLociPattern, \
															no_of_sites_removed_already=noOfRetainedLoci)
					noOfTotalGenotypes = cls.handleNoOfSitesRemovedPattern(line, noOfTotalGenotypesPattern, \
															no_of_sites_removed_already=noOfTotalGenotypes)
					noOfOutputGenotypes = cls.handleNoOfSitesRemovedPattern(line, noOfOutputGenotypesPattern, \
															no_of_sites_removed_already=noOfOutputGenotypes)
					noOfNonBiAllelicLoci = cls.handleNoOfSitesRemovedPattern(line, noOfNonBiAllelicLociPattern, \
															no_of_sites_removed_already=noOfNonBiAllelicLoci)
				del inf
				if counter%reportPeriod==0:
					sys.stderr.write(" .\n")
			import csv
			writer = csv.writer(sys.stdout, delimiter='\t')
			writer.writerow(["noOfProcessedLoci", noOfProcessedLoci])
			writer.writerow(["noOfRetainedLoci", noOfRetainedLoci])
			writer.writerow(["noOfNonBiAllelicLociAfterGQ&DepthFilter", noOfNonBiAllelicLoci])
			writer.writerow(["noOfTotalGenotypes", noOfTotalGenotypes])
			writer.writerow(["noOfGenotypesPassingFilter", noOfOutputGenotypes])
			no_of_masked_genotypes_in_output = noOfTotalGenotypes - noOfOutputGenotypes - (noOfProcessedLoci-noOfRetainedLoci)*sampleSize
			#no_of_masked_genotypes_in_output = noOfRetainedLoci*sampleSize - noOfOutputGenotypes	#same output
			writer.writerow(["no_of_masked_genotypes_in_output", no_of_masked_genotypes_in_output])
			writer.writerow(["no_of_non_missing_genotypes_in_output", noOfRetainedLoci*sampleSize-no_of_masked_genotypes_in_output])
			
		"""
		#2011-11-21
		
		inputDir='/Network/Data/vervet/vervetPipeline/work/'
		FilterVCFByDepthStdoutFnamePattern="FilterVCF_4HighCovVRC_isq_15_18_vs_524_top7559Contigs_multi_sample_gatk_inter_minMAC2_condor.2011.11.17T1309/merge_workflow-FilterVCFByDepth-1.0_PID2_ID*out.000"
		FilterVCFByDepthStdoutFnamePattern = os.path.join(inputDir, FilterVCFByDepthStdoutFnamePattern)
		VariantDiscovery.TallyFilteredSNPs.findSitesFilteredOrMaskedByFilterVCFByDepth(FilterVCFByDepthStdoutFnamePattern=FilterVCFByDepthStdoutFnamePattern, \
			sampleSize=4)
		sys.exit(0)
		
		#2011.12.9
		inputDir='/Network/Data/vervet/vervetPipeline/work/'
		workflowName = 'FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25_2011.12.1T1155'
		FilterVCFByDepthStdoutFnamePattern="%s/merge_workflow-FilterVCFByDepth*out.000"%(workflowName)
		FilterVCFByDepthStdoutFnamePattern = os.path.join(inputDir, FilterVCFByDepthStdoutFnamePattern)
		VariantDiscovery.TallyFilteredSNPs.findSitesFilteredOrMaskedByFilterVCFByDepth(FilterVCFByDepthStdoutFnamePattern=FilterVCFByDepthStdoutFnamePattern, \
			sampleSize=101, grepPattern='call\/Contig')
		sys.exit(0)
		"""
	
		@classmethod
		def countNoOfGATKSAMtoolsPerfectMatchSites(cls, inputDir, ):
			"""
			2011-11-21
				go into the folder which contains output left by CalculateSNPMismatchRateOfTwoVCF.py
					count the lines for each file
				and add them up
			"""
			from pymodule import runLocalCommand
			counter = 0
			no_of_PM_sites = 0
			for filename in os.listdir(inputDir):
				counter += 1
				sys.stderr.write("%s %s "%(counter, filename))
				inputFname = os.path.join(inputDir, filename)
				commandline = "wc -l %s|awk -F ' ' '{print $1-1}'"%inputFname
				return_data = runLocalCommand(commandline, report_stderr=False, report_stdout=False)
				stdout_content = return_data.stdout_content
				no_of_PM_sites += int(stdout_content.strip())-1
				sys.stderr.write(" .\n")
			print "no_of_PM_sites\t%s"%no_of_PM_sites
			
		"""
		#2011-11-21
		inputDir='/Network/Data/vervet/vervetPipeline/'
		subFolder = 'FilterVCF_4HighCovVRC_isq_15_18_vs_524_top7559Contigs_multi_sample_gatk_inter_minMAC2_condor.2011.11.17T1309/SNPMismatchStat/'
		inputDir = os.path.join(inputDir, subFolder)
		VariantDiscovery.TallyFilteredSNPs.countNoOfGATKSAMtoolsPerfectMatchSites(inputDir)
		sys.exit(0)
		"""
		
		@classmethod
		def handleNoOfSitesRemovedPattern(cls, line, pattern, no_of_sites_removed_already=0):
			"""
			2011.12.9
			"""
			searchResult = pattern.search(line)
			if searchResult:
				no_of_sites_removed_already += int(searchResult.group(1))
			return no_of_sites_removed_already
		
		@classmethod
		def countTotalNoOfSitesFilteredByVCFtools(cls, inputDir, reportPeriod=50, fileSuffix='filter_by_vcftools.log'):
			"""
			2011-11-21
				get the stats from the pegasus workflow stdout files
			"""
			sys.stderr.write("Counting sites filtered by vcftools from log files in %s ..."%(inputDir))
			import re
			searchPattern = re.compile(r'After filtering, kept (\d+) out of a possible (\d+) Sites')
			noOfSitesRemovedByGivenPositionsPattern = re.compile(r'Filtering sites by Positions file... (\d+) sites removed.')
			noOfSitesRemovedByAFAndCallRatePattern = re.compile(r'Filtering sites by allele frequency and call rate... (\d+) sites removed.')
			noOfSitesRemovedByACAndMissCountPattern = re.compile(r'Filtering sites by allele count and missing data... (\d+) sites removed.')
			noOfProcessedLoci = 0
			noOfRetainedLoci = 0
			noOfSitesRemovedByGivenPositions = 0
			noOfSitesRemovedByAFAndCallRate = 0
			noOfSitesRemovedByACAndMissCount = 0
			counter = 0
			for filename in os.listdir(inputDir):
				counter += 1
				if filename.find(fileSuffix)==-1:
					continue
				if counter%reportPeriod==0:
					sys.stderr.write("%s %s "%(counter, filename))
				inputFname = os.path.join(inputDir, filename)
				inf = open(inputFname)
				for line in inf:
					searchResult = searchPattern.search(line)
					
					if searchResult:
						noOfRetainedLoci += int(searchResult.group(1))
						noOfProcessedLoci += int(searchResult.group(2))
					noOfSitesRemovedByGivenPositions = cls.handleNoOfSitesRemovedPattern(line, noOfSitesRemovedByGivenPositionsPattern, \
															no_of_sites_removed_already=noOfSitesRemovedByGivenPositions)
					noOfSitesRemovedByAFAndCallRate = cls.handleNoOfSitesRemovedPattern(line, noOfSitesRemovedByAFAndCallRatePattern, \
															no_of_sites_removed_already=noOfSitesRemovedByAFAndCallRate)
					noOfSitesRemovedByACAndMissCount = cls.handleNoOfSitesRemovedPattern(line, noOfSitesRemovedByACAndMissCountPattern, \
															no_of_sites_removed_already=noOfSitesRemovedByACAndMissCount)
				del inf
				if counter%reportPeriod==0:
					sys.stderr.write(" .\n")
			import csv
			writer = csv.writer(sys.stdout, delimiter='\t')
			writer.writerow(["noOfProcessedLoci", noOfProcessedLoci])
			writer.writerow(["noOfRetainedLoci", noOfRetainedLoci])
			writer.writerow(["noOfSitesRemovedByGivenPositions", noOfSitesRemovedByGivenPositions])
			writer.writerow(["noOfSitesRemovedByAFAndCallRate", noOfSitesRemovedByAFAndCallRate])
			writer.writerow(["noOfSitesRemovedByACAndMissCount", noOfSitesRemovedByACAndMissCount])
			
		"""
		#2011-11-21
		
		inputDir='/Network/Data/vervet/vervetPipeline/work/outputs'
		subFolder="FilterVCF_4HighCovVRC_isq_15_18_vs_524_top7559Contigs_multi_sample_gatk_inter_minMAC2_condor.2011.11.17T1309/call_vcftoolsFilter"
		inputDir = os.path.join(inputDir, subFolder)
		VariantDiscovery.TallyFilteredSNPs.countTotalNoOfSitesFilteredByVCFtools(inputDir=inputDir)
		sys.exit(0)
		
		#2011.12.9
		inputDir='/Network/Data/vervet/vervetPipeline/scratch'
		workflowName = 'FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25_2011.12.1T1155'
		subFolder="%s/call_vcftoolsFilter"%(workflowName)
		inputDir = os.path.join(inputDir, subFolder)
		VariantDiscovery.TallyFilteredSNPs.countTotalNoOfSitesFilteredByVCFtools(inputDir=inputDir)
		sys.exit(0)
		"""
		
	@classmethod
	def drawAACHistogram(cls, AAC_ls = [], count_2D_ls =[], outputFnamePrefix=None, color_ls=['r', 'g', 'b', 'y','c', 'k'], **kwargs):
		"""
		2011-11-22
			count_2D_ls is a list of list. Each included list is the number of variants at that particular AAC.
			
		"""
		sys.stderr.write("Drawing barChart of %s bars, %s data series to %s..."%(len(AAC_ls), len(count_2D_ls), outputFnamePrefix))
		import pylab
		pylab.clf()
		rects_ls = []
		width = 0.75/len(count_2D_ls)
		for i in xrange(len(count_2D_ls)):
			c = color_ls[i%len(color_ls)]
			ind = numpy.array(AAC_ls) + width * i
			rects = pylab.bar(ind, count_2D_ls[i], width=width, bottom=0, color=c, log=False, **kwargs)
			rects_ls.append(rects)
		
		if len(rects_ls)==2:
			pylab.legend( (rects_ls[0][0], rects_ls[1][0]), ('version 0', 'version 1') )
		
		"""
		pylab.title(title)
		if xlabel_1D is not None:
			pylab.xlabel(xlabel_1D)
		if xticks:
			pylab.xticks(x_ls, xticks)
		"""
		pylab.savefig('%s.png'%outputFnamePrefix, dpi=200)
		pylab.savefig('%s.svg'%outputFnamePrefix, dpi=200)
		
		sys.stderr.write("Done.\n")
		
	"""
		#2011-11-22
		AAC_ls = range(1,9)
		# data from the 4 HC monkeys
		count_2D_ls = [[2499739, 1969155, 1567859, 1270833, 933173, 738506, 524289, 433315], \
					[0, 1245574, 949732, 718876, 546747, 430943, 0, 0]]
		outputFnamePrefix = '/tmp/AAC_distribution_of_4HCMonkeys'
		VariantDiscovery.drawAACHistogram(AAC_ls = AAC_ls, count_2D_ls =count_2D_ls, outputFnamePrefix=outputFnamePrefix, \
				color_ls=['r', 'g', 'b', 'y','c', 'k'])
		sys.exit(0)
	"""
	
	@classmethod
	def drawAACHistogramFromAACFile(cls, inputFname=None, outputFnamePrefix=None, color_ls=['r', 'g', 'b', 'y','c', 'k'], **kwargs):
		"""
		2011-12-7
			read data from inputFname (AAC_tally.tsv from CalculateVCFStatPipeline.py)
			and then call drawAACHistogram()
		"""
		import csv
		reader = csv.reader(open(inputFname), delimiter='\t')
		AAC2NoOfLoci = {}
		header = reader.next()
		for row in reader:
			AAC = int(row[0])
			noOfLoci = int(float(row[1]))
			if AAC not in AAC2NoOfLoci:
				AAC2NoOfLoci[AAC]= 0
			AAC2NoOfLoci[AAC] += noOfLoci
		sys.stderr.write("Found %s different AAC counts.\n"%(len(AAC2NoOfLoci)))
		
		AAC_NoOfLoci_ls = AAC2NoOfLoci.items()
		AAC_NoOfLoci_ls.sort()
		AAC_ls = [row[0] for row in AAC_NoOfLoci_ls]
		noOfLoci_ls = [row[1] for row in AAC_NoOfLoci_ls]
		
		VariantDiscovery.drawAACHistogram(AAC_ls = AAC_ls, count_2D_ls =[noOfLoci_ls], outputFnamePrefix=outputFnamePrefix, \
				color_ls=['r', 'g', 'b', 'y','c', 'k'])
		
	"""
		#2011-12.7
		workflowName = 'VCFStat_LowPass_top7559Contigs_no12eVarFilter_inter_100kbwindow_2011.12.6T2342'
		workflowName = 'VCFStat_FilterVCF_LowPass_top7559Contigs_no12eVarFilter_inter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25_100kbwindow.2011.12.6T2345'
		inputFname='/Network/Data/vervet/vervetPipeline/%s/AAC_tally.tsv'%(workflowName)
		outputFnamePrefix = '//Network/Data/vervet/vervetPipeline/AAC_distribution_of_101LCMonkeys-version1'
		VariantDiscovery.drawAACHistogramFromAACFile(inputFname=inputFname, outputFnamePrefix=outputFnamePrefix, \
				color_ls=['r', 'g', 'b', 'y','c', 'k'])
		sys.exit(0)
		
	"""

# 2011-4-29 a handy function to strip blanks around strings
strip_func = lambda x: x.strip()


class DBVervet(object):
	"""
	2011-4-27
		class to hold functions related to vervetdb
	"""
	def __init__(self):
		pass
	

	import re
	#any_character_pattern = re.compile(r'[^0-9.]')
	any_character_pattern = re.compile(r'[a-zA-Z\[\]]')	#[] is also a character that needs to be ignored.
	
	
	@classmethod
	def putAlignmentIntoDB(cls, db_vervet, inputFname, individual_code, \
						sequencer='GA', sequence_type='short-read', sequence_format='fastq', \
						ref_individual_sequence_id=10, \
						alignment_method_name='bwa-short-read', alignment_format='bam', createSymbolicLink=False):
		"""
		2011-7-11
			add argument createSymbolicLink
		2011-5-6
		"""
		sys.stderr.write("Putting alignment file %s into db ..."%(inputFname))
		db_vervet.session.begin()
		individual_alignment = db_vervet.getAlignment(individual_code=individual_code, path_to_original_alignment=inputFname, \
					sequencer=sequencer, sequence_type=sequence_type, sequence_format=sequence_format, \
					ref_individual_sequence_id=ref_individual_sequence_id, \
					alignment_method_name=alignment_method_name, alignment_format=alignment_format,\
					createSymbolicLink=createSymbolicLink)
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
	
	"""
		#2011-5-9
		inputFname = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/ref/illumina/ref_illu_vs_1MbBAC_by_aln.bam')
		individual_code = 'VRC_ref'
		DBVervet.putAlignmentIntoDB(db_vervet, inputFname, individual_code, ref_individual_sequence_id=9)
		sys.exit(4)
		
		#2011-5-9
		inputFname = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/ref/illumina/ref_illu_vs_1MbBAC_by_aln.bam')
		for individual_name in ['aethiops', 'Barbados', 'cynosurus', 'pygerythrus', 'sabaeus', 'tantalus']:
			inputFname = os.path.expanduser('~/mnt/hoffman2_home/script/vervet/data/subspecies/%s/%s_vs_1MbBAC_by_aln.bam'%\
									(individual_name, individual_name))
			individual_code = individual_name
			DBVervet.putAlignmentIntoDB(db_vervet, inputFname, individual_code, \
						sequencer='GA', sequence_type='short-read', sequence_format='fastq', \
						ref_individual_sequence_id=9, \
						alignment_method_name='bwa-short-read', alignment_format='bam', createSymbolicLink=False)
		sys.exit(0)
		
		#2011-7-11
		inputFname = os.path.expanduser('/usr/local/vervetData/ref/illumina/vs_MinSize200Scaffolds_by_aln.bam')
		individual_code = 'VRC_ref'
		DBVervet.putAlignmentIntoDB(db_vervet, inputFname, individual_code, \
						sequencer='GA', sequence_type='short-read', sequence_format='fastq', \
						ref_individual_sequence_id=120, \
						alignment_method_name='bwa-short-read', alignment_format='bam', createSymbolicLink=True)
		inputFname = os.path.expanduser('/usr/local/vervetData/ref/454/vs_MinSize200Scaffolds_by_bwasw.bam')
		individual_code = 'VRC_ref'
		DBVervet.putAlignmentIntoDB(db_vervet, inputFname, individual_code, \
						sequencer='454', sequence_type='short-read', sequence_format='fastq', \
						ref_individual_sequence_id=120, \
						alignment_method_name='bwa-long-read', alignment_format='bam', createSymbolicLink=True)
		
		
		for individual_name in ['aethiops', 'Barbados', 'cynosurus', 'pygerythrus', 'sabaeus', 'tantalus']:
			inputFname = os.path.expanduser('/usr/local/vervetData/subspecies/%s/vs_MinSize200Scaffolds_by_aln.bam'%\
									(individual_name))
			individual_code = individual_name
			DBVervet.putAlignmentIntoDB(db_vervet, inputFname, individual_code, \
						sequencer='GA', sequence_type='short-read', sequence_format='fastq', \
						ref_individual_sequence_id=120, \
						alignment_method_name='bwa-short-read', alignment_format='bam', createSymbolicLink=True)
		sys.exit(0)
	"""
	
	@classmethod
	def indexBamAlignmentFilesInDB(cls, db_vervet, samtools_path=os.path.expanduser("~/bin/samtools")):
		"""
		2011-5-9
		"""
		sys.stderr.write("Generating index files for all un-indexed bam alignment files in db ... \n")
		import VervetDB
		query = VervetDB.IndividualAlignment.query.filter_by(format='bam')
		from pymodule.utils import runLocalCommand
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			sys.stderr.write('%s: %s'%(counter, row.path))
			bam_pathname = os.path.join(db_vervet.data_dir, row.path)
			bam_index_pathname = os.path.join(db_vervet.data_dir, '%s.bai'%(row.path))
			if not os.path.isfile(bam_index_pathname):	#index doesn't exist
				commandline = '%s index %s'%(samtools_path, bam_pathname)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				real_counter += 1
			sys.stderr.write('.\n')
		sys.stderr.write("Done.\n")
	
	"""
		#2011-5-9
		DBVervet.indexBamAlignmentFilesInDB(db_vervet)
		sys.exit(0)
		
		
	"""
	
	@classmethod
	def pokeBamReadGroupPresence(cls, db_vervet, samtools_path=os.path.expanduser("~/bin/samtools"), dataDir=None, commit=False):
		"""
		2011-5-9
		"""
		sys.stderr.write("checking which bam file in individual_alignment table has read group or not ... \n")
		import VervetDB
		if dataDir is None:
			dataDir = db_vervet.data_dir
		session = db_vervet.session
		session.begin()
		import pysam
		query = VervetDB.IndividualAlignment.query.filter_by(format='bam')
		counter = 0
		real_counter = 0
		for row in query:
			counter += 1
			sys.stderr.write('%s: %s'%(counter, row.path))
			
			bam_pathname = os.path.join(dataDir, row.path)
			if os.path.isfile(bam_pathname):
				samfile = pysam.Samfile(bam_pathname, "rb" )
				if "RG" not in samfile.header:
					read_group_added=0
				else:
					read_group_added=1
				del samfile
			else:
				read_group_added=None	#file doesn't exist
			if row.read_group_added!=read_group_added:
				row.read_group_added = read_group_added
				real_counter += 1
				session.add(row)
			sys.stderr.write('.\n')
		if commit:
			session.commit()
		else:
			session.rollback()
		sys.stderr.write("%s/%s bams has read_group_added changed. Done.\n"%(real_counter, counter))
	
	"""
		#2011-9-15
		dataDir = os.path.expanduser("~/mnt/hpc-cmb_home/NetworkData/vervet/db/")
		DBVervet.pokeBamReadGroupPresence(db_vervet, samtools_path=os.path.expanduser("~/bin/samtools"), dataDir=dataDir, commit=True)
		sys.exit(0)
		
		
	"""
	
	@classmethod
	def moveLostIndividualSeqFileIntoCorrespondingFolder(cls, db_vervet, dataDir=None, inputDir=None, outputDir=None):
		"""
		2011-9-7
			dataDir is where the vanilla copy of db-affiliated storage is, used to generate the fname2folder dict
			inputDir is where "lost" individual_sequence files are.
			outputDir is where all "lost" files should be moved into. new folder will be created to house "lost files".
			
			If files in respective folder already exist, no move happens.
		"""
		sys.stderr.write("Moving lost individual_sequence files in %s into corresponding folders in %s ... \n"%(inputDir, outputDir))
		import VervetDB
		query = VervetDB.IndividualSequence.query.filter(VervetDB.IndividualSequence.id>=167)
		if dataDir is None:
			dataDir = db_vervet.data_dir
		fname2folder = {}
		no_of_folders = 0
		for row in query:
			if row.path:
				abs_path = os.path.join(dataDir, row.path)
				if os.path.isdir(abs_path):
					files = os.listdir(abs_path)
					folder = os.path.basename(abs_path)
					no_of_folders  += 1
					for fname in files:
						if fname  in fname2folder:
							sys.stderr.write("Error: %s already exists in this folder %s, here %s , 2nd folder.\n"%\
											(fname, fname2folder.get(fname), folder))
						else:
							fname2folder[fname] = folder
		sys.stderr.write("%s files in %s folders.\n"%(len(fname2folder), no_of_folders))
		
		from pymodule import runLocalCommand
		lostFiles = os.listdir(inputDir)
		no_of_files_moved = 0
		for lostFname in lostFiles:
			if lostFname in fname2folder:
				folder = fname2folder.get(lostFname)
				src_abs_path = os.path.join(inputDir, lostFname)
				dst_folder_abs_path = os.path.join(outputDir, folder)
				if not os.path.exists(dst_folder_abs_path):
					os.makedirs(dst_folder_abs_path)
				dst_abs_path = os.path.join(dst_folder_abs_path, lostFname)
				if not os.path.exists(dst_abs_path):
					commandline = "mv %s %s"%(src_abs_path, dst_folder_abs_path)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
					no_of_files_moved += 1
				else:
					sys.stderr.write("%s already exists in %s. ignored.\n"%(lostFname, os.path.basename(dst_folder_abs_path)))
		sys.stderr.write("Moved %s files.\n"%no_of_files_moved)
	
	"""
		#2011-9-7
		inputDir=os.path.expanduser('~/mnt/hpc-cmb_home/NetworkData/vervet/db/individual_sequence/')
		DBVervet.moveLostIndividualSeqFileIntoCorrespondingFolder(db_vervet, dataDir=None, \
			inputDir=inputDir, outputDir=inputDir)
		sys.exit(0)
		
	"""
	
	
	@classmethod
	def putPedigreeIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", \
						collector_name='Nelson B. Freimer', tax_id=60711, \
						city='VRC', country_name='United States of America', latitude=36.090, longitude=-80.270,\
						default_collection_date = '1970'):
		"""
		2011-4-27
		"""
		sys.stderr.write("Putting Pedigree monkeys (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr
		import VervetDB
		
		reader = csv.reader(open(inputFname,), delimiter='\t')
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("ID")
		sex_index = col_name2index.get("Sex")
		mother_id_index = col_name2index.get("Mo")
		father_id_index = col_name2index.get("Fa")
		counter = 0
		
		collector = db_vervet.getUser(collector_name)
		monkey_id2sex = {}
		monkey_id2mother_father = {}
		for row in reader:
			monkey_id = row[monkey_id_index].strip()
			if not monkey_id:	#skip if row doesn't have monkey_id
				continue
			sex = row[sex_index]
			if sex=='1':
				sex = 'M'
			else:
				sex= 'F'
			if monkey_id not in monkey_id2sex:
				monkey_id2sex[monkey_id] = sex
			
			mother_id = row[mother_id_index]
			father_id = row[father_id_index]
			
			if len(monkey_id)>4:	#monkeys born in the colony
				birthdate = monkey_id[:4]
				collection_date = birthdate
				birthdate = datetime.strptime(birthdate, '%Y')	#2010
			else:	#monkeys born outside
				collection_date = default_collection_date
				birthdate = None
			
			collection_date = datetime.strptime(collection_date, '%Y')	#2010
			
			site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country_name,\
								latitude=latitude, longitude=longitude, altitude=None)
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=None, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, \
								tax_id=tax_id, birthdate=birthdate)
			
			if father_id!='0':
				relationship_type_name = 'Father2Child'
				father = VervetDB.Individual.query.filter_by(ucla_id=father_id).filter_by(tax_id=tax_id).first()
				if not father:
					sys.stderr.write("father %s not in db yet.\n"%(father_id))
					sys.exit(3)
				db_vervet.getInd2Ind(individual1=father, individual2=individual, relationship_type_name=relationship_type_name)
			if mother_id!='0':
				relationship_type_name = 'Mother2Child'
				mother = VervetDB.Individual.query.filter_by(ucla_id=mother_id).filter_by(tax_id=tax_id).first()
				if not mother:
					sys.stderr.write("mother %s not in db yet.\n"%(mother_id))
					sys.exit(3)
				db_vervet.getInd2Ind(individual1=mother, individual2=individual, relationship_type_name=relationship_type_name)
		del reader
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
	
	"""
		#2011-5-6
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/VRC_pedigree1_v6 1.tsv")
		DBVervet.putPedigreeIntoDB(db_vervet, inputFname, )
		sys.exit(3)
		
		#2011-5-6
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/VRC_pedigree1_v6 1.tsv")
		DBVervet.putPedigreeIntoDB(db_vervet, inputFname )
		sys.exit(3)
		
	"""
	
	
	@classmethod
	def putSequencedPedigreeMonkeysIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", \
						collector_name='Nelson B. Freimer', tax_id=60711, \
						city='VRC', country_name='United States of America', latitude=36.090, longitude=-80.270,\
						default_collection_date = '1970'):
		"""
		2011-5-7
			record a IndividualSequence entry for each sequenced monkey, (add a coverage column)
		"""
		sys.stderr.write("Putting sequenced monkeys (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
		import VervetDB
		
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Individual Name")
		tissue_index = col_name2index.get("Tissue Name")
		coverage_index = col_name2index.get("Targeted sequencing coverage")
		sex_index = col_name2index.get("Gender Name")
		counter = 0
		
		collector = db_vervet.getUser(collector_name)
		for row in reader:
			monkey_id = row[monkey_id_index].strip()
			if not monkey_id:	#skip if row doesn't have monkey_id
				continue
			sex = row[sex_index]
			
			if len(monkey_id)>4:	#monkeys born in the colony
				birthdate = monkey_id[:4]
				collection_date = birthdate
				birthdate = datetime.strptime(birthdate, '%Y')	#2010
			else:	#monkeys born outside
				collection_date = default_collection_date
				birthdate = None
			
			collection_date = datetime.strptime(collection_date, '%Y')	#2010
			
			site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country_name,\
								latitude=latitude, longitude=longitude, altitude=None)
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=None, latitude=None,\
								longitude=None, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, \
								tax_id=tax_id, birthdate=birthdate)
			tissue_name = row[tissue_index]
			coverage = float(row[coverage_index][:-1])	#remove the final "x"
			individual_sequence = db_vervet.getIndividualSequence(individual_id=individual.id, sequencer='GA', \
						sequence_type='short-read',\
						sequence_format='fastq', path_to_original_sequence=None, tissue_name=tissue_name, coverage=coverage)
		
		del reader
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
		
	"""
		#2011-5-9
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/VRC/09 Avid chip informaiton.tsv")
		DBVervet.putSequencedPedigreeMonkeysIntoDB(db_vervet, inputFname)
		sys.exit(4)
		
	"""
	
	@classmethod
	def putVRCFoundersIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", \
						collector_name='Nelson B. Freimer', tax_id=60711, \
						city='VRC', country_name='United States of America', latitude=36.090, longitude=-80.270,\
						default_collection_date = '1970'):
		"""
		2011-5-7
			mark individuals who are founders of VRC in table Individual
		"""
		sys.stderr.write("Marking monkeys that are founders of VRC ...\n")
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
		import VervetDB
		
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("AnimalID")
		sex_index = col_name2index.get("Sex")
		counter = 0
		
		collector = db_vervet.getUser(collector_name)
		for row in reader:
			monkey_id = row[monkey_id_index].strip()
			if not monkey_id:	#skip if row doesn't have monkey_id
				continue
			sex = row[sex_index]
			
			if len(monkey_id)>4:	#monkeys born in the colony
				birthdate = monkey_id[:4]
				collection_date = birthdate
				birthdate = datetime.strptime(birthdate, '%Y')	#2010
			else:	#monkeys born outside
				collection_date = default_collection_date
				birthdate = None
			
			collection_date = datetime.strptime(collection_date, '%Y')	#2010
			
			site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country_name,\
								latitude=latitude, longitude=longitude, altitude=None)
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=None, latitude=None,\
								longitude=None, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, \
								tax_id=tax_id, birthdate=birthdate, vrc_founder=True)
			if individual.vrc_founder!=True:
				individual.vrc_founder = True
				db_vervet.session.add(individual)
				db_vervet.session.flush()
		del reader
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
		
	"""
		
		#2011-5-9
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/VRC/vrc founders.csv")
		DBVervet.putVRCFoundersIntoDB(db_vervet, inputFname)
		sys.exit(4)
		
		
	"""
	
	@classmethod
	def putSIVDataIntoDB(cls, db_vervet, inputFname, collector_name='University of Pittsburgh', country_name='South Africa',\
						tax_id=460674):
		"""
		2011-5-10
		
		"""
		sys.stderr.write("Putting SIV data from %s into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
		import VervetDB
		
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("UCLA ID")
		tissue_id_index = col_name2index.get("Unique ID")
		SIV_infection_index = col_name2index.get("SIV")
		STLV_infection_index = col_name2index.get("STLV")
		site_index = col_name2index.get("Chris List_Location")
		age_index = col_name2index.get("Age")
		sex_index = col_name2index.get("Sex")
		counter = 0
		
		collector = db_vervet.getUser(collector_name)
		for row in reader:
			monkey_id = row[monkey_id_index].strip()
			if not monkey_id:	#skip if row doesn't have monkey_id
				continue
			sex = row[sex_index]
			
			age = row[age_index].strip()
			age = cls.filterValue(age, data_type=int)
			site_str = row[site_index]
			
			individual = VervetDB.Individual.query.filter_by(ucla_id=monkey_id).first()
			if individual is None:
				sys.stderr.write("Error: Monkey %s not present in db.\n"%(monkey_id))
				
				city = None
				description = None
				province = None
			
				site_name_ls = site_str.split(",")
				site_name_ls = map(strip_func, site_name_ls)
				if len(site_name_ls)==2:
					city, province = site_name_ls
					description = None
				elif len(site_name_ls)==3:
					description, city, province = site_name_ls
				else:
					sys.stderr.write("Error: site %s is neither 3-entry nor 4-entry.\n"%site_str)
					sys.exit(0)
				site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country_name)
				
				individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=age, age_cas=None, latitude=site.latitude,\
									longitude=site.longitude, altitude=None, ucla_id=monkey_id, site=site, \
									collection_date=None, collector=collector, tax_id=tax_id)
			
			for i in [SIV_infection_index, STLV_infection_index]:
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					value = row[i].strip()
					if value == '-':
						value = 0
					elif value == '+':
						value = 1
					else:
						sys.stderr.write("phenotype value %s not parsable.\n"%(value))
						sys.exit(5)
					tissue_id = row[tissue_id_index]
					comment =tissue_id
					db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
		del reader
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
		
	"""
		
		#2011-5-9
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/SIV/South Africa AGMs data.tsv")
		DBVervet.putSIVDataIntoDB(db_vervet, inputFname)
		sys.exit(4)
		
		
	"""
	
	@classmethod
	def putSIVStrainClusterInfoIntoDB(cls, db_vervet, inputFname, collector_name='University of Pittsburgh'):
		"""
		2011-5-10
		"""
		sys.stderr.write("Putting SIV strain cluster info data from %s into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
		import VervetDB
		
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		tissue_id_index = col_name2index.get("Tissue ID")
		cluster_id_index = col_name2index.get("Cluster ID")
		counter = 0
		
		for row in reader:
			if not row:
				continue
			tissue_id = row[tissue_id_index].strip()
			cluster_id = row[cluster_id_index].strip()
			if not tissue_id:
				continue
			
			phenotype_query = VervetDB.Phenotype.query.filter_by(comment=tissue_id)
			individual_id_set = set()
			for phenotype in phenotype_query:
				individual_id_set.add(phenotype.individual.id)
			if len(individual_id_set)==0:
				sys.stderr.write("no individual could be found given tissue id %s.\n"%(tissue_id))
				sys.exit(4)
			elif len(individual_id_set)>1:
				sys.stderr.write("%s individuals were found to have this tissue id %s.\n"%(len(individual_id_set), tissue_id))
				sys.exit(5)
			else:
				individual = VervetDB.Individual.get(individual_id_set.pop())
				phenotype_name = 'SIVStrainCluster'
				value = int(cluster_id)
				
				comment =tissue_id
				db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
									individual_id=individual.id, comment=comment, collector_name=collector_name)
		del reader
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
		
	"""
		
		#2011-5-9
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/SIV/SIV_strain_onSouthAfricanVervets.csv")
		DBVervet.putSIVStrainClusterInfoIntoDB(db_vervet, inputFname)
		sys.exit(4)
		
		
		
		
	"""
	
	@classmethod
	def handleMonkeyID(cls, monkey_id, monkeyIDPrefix="", no_of_digits_in_id=5):
		"""
		2011-5-5
			convenient function to
			1. add extra "0"s in front of the pure-number form of monkey_id
			2. add monkeyIDPrefix
		"""
		monkey_id = str(monkey_id)
		monkey_id_len = len(monkey_id)
		monkey_id = '0'*(no_of_digits_in_id-monkey_id_len) + monkey_id	#add "0" in front of it
		monkey_id = monkeyIDPrefix + monkey_id	#add the prefix
		return monkey_id
	
	@classmethod
	def put2011StKittsIntoDB(cls, db_vervet, inputFname, gpsDataInputFname, monkeyIDPrefix="VWP", \
							collector_name='Nelson B. Freimer', tax_id=60711):
		"""
		2011-4-27
		"""
		sys.stderr.write("Putting 2011 St. Kitts collection (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule import PassingData
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr
		
		sys.stderr.write("\t Getting gps&site data from %s ..."%(gpsDataInputFname))
		reader = csv.reader(open(gpsDataInputFname,), delimiter='\t')
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		id_range_index = col_name2index.get("ID Range")
		site_index = col_name2index.get("Location Name")
		latitude_index = col_name2index.get("Latitude")
		longitude_index = col_name2index.get("Longitude")
		altitude_index = col_name2index.get("Altitude")
		
		collection_date_index = col_name2index.get("Collection Date")
		
		#2011-4-29 first get the location infomation into db.
		# a dictionary which maps code to site and collection_date
		code2data = {}
		counter = 0
		for row in reader:
			if not row[collection_date_index]:	#skip. no date means no data
				continue
			id_range = row[id_range_index]
			id_range = getListOutOfStr(id_range, data_type=int)
			
			site_str = row[site_index]
			latitude = float(row[latitude_index].strip())
			longitude = float(row[longitude_index].strip())
			
			altitude = row[altitude_index].strip()
			if altitude:
				altitude = float(altitude)
			else:
				altitude = None
			
			city = "Unknown"
			description = None
			province = None
			country = None
			if site_str and site_str!='?':
				site_name_ls = site_str.split(",")
				site_name_ls = map(strip_func, site_name_ls)
				if len(site_name_ls)==2:
					city, country = site_name_ls[:2]
					if country =='St. Kitts':
						country = "Saint Kitts"
				else:
					sys.stderr.write("Error: site %s is not in 2-entry form.\n"%site_str)
					sys.exit(0)
				
				site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country,\
								latitude=latitude, longitude=longitude, altitude=altitude)
			else:
				sys.stderr.write("No site info for these IDs: %s.\n"%(repr(row[id_range_index])))
				sys.exit(0)
			for monkey_id in id_range:
				monkey_id = cls.handleMonkeyID(monkey_id, monkeyIDPrefix=monkeyIDPrefix, no_of_digits_in_id=5)
				if monkey_id not in code2data:
					code2data[monkey_id] = PassingData()
				code2data[monkey_id].site = site
			counter += 1
		del reader
		sys.stderr.write("%s sites found.\n"%(counter))
		
		sys.stderr.write("Reading individual data from %s ..."%(inputFname))
		reader = csv.reader(open(inputFname,), delimiter='\t')
		header = reader.next()
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Monkey ID")
		sex_index = col_name2index.get("Sex")
		#age_index = col_name2index.get("Age Category")
		#age_cas_index = col_name2index.get("Age Category (CAS)")
		collection_date_index = col_name2index.get("Date")
		phenotype_start_index = col_name2index.get("Age/Dentition present")
		
		
		collector = db_vervet.getUser(collector_name)
		
		for row in reader:
			monkey_id = row[monkey_id_index].strip()
			if not monkey_id:	#skip if row doesn't have monkey_id
				continue
			monkey_id = cls.handleMonkeyID(monkey_id, monkeyIDPrefix=monkeyIDPrefix, no_of_digits_in_id=5)
			sex = row[sex_index]
			"""
			age = row[age_index]
			age = cls.filterValue(age, int)
			age_cas = row[age_cas_index]
			age_cas = cls.filterValue(age_cas, int)
			"""
			collection_date = row[collection_date_index].strip()
			if collection_date:
				collection_date = datetime.strptime(collection_date, '%Y/%m/%d')	#2010/1/13
			else:	#2011-5-5 skip if Date is not available (there are some phantom monkeys)
				continue
			
			site = code2data.get(monkey_id).site
			if site is not None:
				latitude = site.latitude
				longitude = site.longitude
				altitude = site.altitude
			else:
				latitude = None
				longitude = None
				altitude = None
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=None, latitude=latitude,\
								longitude=longitude, altitude=altitude, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, \
								tax_id=tax_id)
			
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:	#non-numeric phenotypes
						value = None
						comment = row[i]
					else:
						value_str, no_of_replacements = cls.any_character_pattern.subn('', row[i].strip())[:2]
						value = cls.filterValue(value_str, data_type=float)
						comment =None
					if value or comment:	#one of them has to be something
						db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
	
	"""
		# 2011-5-5
		inputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/Animal Data Sheets (Narek)ajj_refined.tsv")
		gpsDataInputFname =  os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/Date_GPS_SK2011_to Yu 4-15-11_refined.tsv")
		
		DBVervet.put2011StKittsIntoDB(db_vervet, inputFname, gpsDataInputFname, monkeyIDPrefix="VWP", \
								collector_name='Nelson B. Freimer', tax_id=60711)
		sys.exit(0)
	"""
	
	@classmethod
	def put2010StKittsIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", collector_name='Nelson B. Freimer', \
							country_name='Saint Kitts', tax_id=60711):
		"""
		2011-4-27
		"""
		
		sys.stderr.write("Putting 2010 St. Kitts collection (%s) into db ...\n"%(inputFname))
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
		age_index = col_name2index.get("Age Category")
		age_cas_index = col_name2index.get("Age Category (CAS)")
		collection_date_index = col_name2index.get("Date")
		phenotype_start_index = col_name2index.get("Weight")
		
		gps_pattern = re.compile(r'(?P<lat_direction>[SN])(?P<lat_hr>\d+) +(?P<lat_min>[\d.]+) +(?P<lon_direction>[EW])(?P<lon_hr>\d+) +(?P<lon_min>[\d.]+)')
		number_subtraction_pattern = re.compile(r'(?P<number1>[\d.]+)-(?P<number2>[\d.]+)')
		
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
			collection_date = row[collection_date_index].strip()
			if collection_date:
				collection_date = datetime.strptime(collection_date, '%Y/%m/%d')	#2010/1/13
			else:
				collection_date = None
			
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
			else:
				sys.stderr.write("Error: gps %s is not parsable.\n"%(gps))
				sys.exit(0)
			
			city = site_str.strip()
			site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country_name, latitude=latitude,\
									longitude=longitude)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=age, age_cas=age_cas, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector, tax_id=tax_id)
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:	#non-numeric phenotypes
						value = None
						comment = row[i]
					else:
						if cls.any_character_pattern.search(row[i]):	#ignore any column with any character in it.
							# should be number only
							continue
						number_subtraction_pattern_result = number_subtraction_pattern.search(row[i])
						if number_subtraction_pattern_result:
							number1 = float(number_subtraction_pattern_result.group('number1'))
							number2 = float(number_subtraction_pattern_result.group('number2'))
							value = number1 - number2
						else:
							value = cls.filterValue(row[i], data_type=float)
						comment =None
					if value or comment:
						db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("Done.\n")
	
	"""
		
		#2011-5-5
		inputFname = os.path.expanduser("~/mnt/banyan/mnt/win/vervet-analysis/Yu/SK Jan 2010_no_SKBRF.tsv")
		DBVervet.put2010StKittsIntoDB(db_vervet, inputFname)
		sys.exit(0)
		
	"""
	
	
	@classmethod
	def put2010StKittsColonyIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", collector_name='St. Kitts Colony', \
							country_name='Saint Kitts', tax_id=60711):
		"""
		2011-9-8
			input is a spreasheet prepared by Qiao based on paper slides from Ania
		"""
		
		sys.stderr.write("Putting 2010 St. Kitts Colony Samples (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule.utils import getColName2IndexFromHeader, figureOutDelimiter
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("AnimalID")
		site_index = col_name2index.get("ORIGIN")
		gps_index = col_name2index.get("GPS position")
		sex_index = col_name2index.get("Sex")
		age_index = col_name2index.get("Approx. Age (years)")
		collection_date_index = col_name2index.get("CollectionDate")
		whatIsIt_index = col_name2index.get("WhatIsIt")
		comment_index = col_name2index.get("Comment")
		originalORIGIN_index = col_name2index.get("OriginalORIGIN")
		
		gps_pattern = re.compile(r'(?P<lat>[-\d.]+), +(?P<lon>[-\d.]+)')
		
		collector = db_vervet.getUser(collector_name)
		no_of_monkeys = 0
		for row in reader:
			monkey_id = monkeyIDPrefix + row[monkey_id_index]
			site_str = row[site_index]
			gps = row[gps_index]
			sex = row[sex_index]
			if len(row)>age_index:
				age = row[age_index]
				age = cls.filterValue(age, float)
			else:
				age = None
			#collection_date = row[collection_date_index].strip()
			collection_date = "04 2010"	#original text is "April 2010 SK"
			if collection_date:
				collection_date = datetime.strptime(collection_date, '%m %Y')
			else:
				collection_date = None
			if gps:
				gps_pattern_search = gps_pattern.search(gps)
				if gps_pattern_search:
					latitude = float(gps_pattern_search.group('lat'))
					longitude = float(gps_pattern_search.group('lon'))
				else:
					sys.stderr.write("Error: gps %s is not parsable.\n"%(gps))
					sys.exit(0)
			else:
				latitude = None
				longitude = None
			
			if len(row)>comment_index:
				comment = row[comment_index]
			else:
				comment = ""
			if len(row)>whatIsIt_index:
				whatIsIt = row[whatIsIt_index]
			else:
				whatIsIt = ""
			if len(row)>originalORIGIN_index:
				originalORIGIN = row[originalORIGIN_index]
			else:
				originalORIGIN = ""
			if originalORIGIN:
				originalORIGIN = "originalORIGIN: %s"%(originalORIGIN)
			
			comment_ls = [whatIsIt, comment, originalORIGIN]
			comment = "; ".join(comment_ls)
			
			city = site_str.strip()
			site = db_vervet.getSite(description=None, city=city, stateprovince=None, country_name=country_name, latitude=latitude,\
									longitude=longitude)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex, age=age, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=None, site=site, \
								collection_date=collection_date, collector=collector, tax_id=tax_id, comment=comment)
			no_of_monkeys += 1
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("%s monkeys. Done.\n"%(no_of_monkeys))
	
	"""
		
		#2011-9-8
		inputFname = os.path.expanduser("/tmp/sk blood rna summary gps.csv")
		DBVervet.put2010StKittsColonyIntoDB(db_vervet, inputFname)
		sys.exit(0)
		
	"""
	
	@classmethod
	def drawPedigree(cls, db_vervet, outputFnamePrefix=None, baseNodeSize=40):
		"""
		2011-5-6
			
		"""
		sys.stderr.write("Getting pedigree from db ...")
		import networkx as nx
		DG=nx.DiGraph()
		
		import VervetDB
		
		for row in VervetDB.Ind2Ind.query:
			DG.add_edge(row.individual1_id, row.individual2_id)
		
		
		sys.stderr.write("%s nodes. %s edges. %s connected components. Done.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
																	nx.number_connected_components(DG.to_undirected())))
		
		from pymodule import yh_matplotlib
		sys.stderr.write("Plotting out degree histogram ....")
		out_degree_ls = []
		for n in DG.nodes_iter():
			out_degree_ls.append(DG.out_degree(n))
		outputFname = '%s_outDegreeHist.png'%(outputFnamePrefix)
		yh_matplotlib.drawHist(out_degree_ls, title="Histogram of no. of children per monkey", xlabel_1D="no. of children", \
							outputFname=outputFname, min_no_of_data_points=50, needLog=True)
		
		sys.stderr.write("Assigning each node with different size/color ...")
		
		sex2node_property_list = {}		# in node_property_list, each entry is (node, size, color)
		#size depends on whether it's deep-sequenced (30X), 4 for 30X, 8 for the REF, 2 for all others. 
		#color depends on it's sequenced or not
		for v in DG:
			individual = VervetDB.Individual.get(v)
			sex = individual.sex
			if individual.sex not in sex2node_property_list:
				sex2node_property_list[sex] = []
			individual_sequence = VervetDB.IndividualSequence.query.filter_by(sequencer='GA').\
				filter_by(individual_id=individual.id).first()
			node_color = 0.8	#color for the sequenced
			
			
			if individual_sequence is None:
				node_size =baseNodeSize
				node_color = 0	#not sequenced in a different color
			elif individual_sequence.coverage<5:
				node_size = baseNodeSize*12
			elif individual_sequence.coverage==30:
				node_size = baseNodeSize*36
			elif individual_sequence.individual_id==1:	#the reference
				node_size = baseNodeSize*108
			else:
				node_size = baseNodeSize
			
			if individual.vrc_founder:
				node_color = 0.25
				if individual_sequence is None:
					node_size = baseNodeSize*12
			
			sex2node_property_list[sex].append((v, node_size, node_color))
		sys.stderr.write("Done.\n")
		
		import pylab, numpy
		#nx.draw_circular(DG,with_labels=False, alpha=0.5)
		pylab.clf()
		pylab.axis("off")
		pylab.figure(figsize=(100, 60))
		layout = 'dot'
		pos = nx.graphviz_layout(DG, prog=layout)
		nx.draw_networkx_edges(DG, pos, alpha=0.9, width=0.8)
		
		import matplotlib as mpl
		"""
		phenotype_cmap = mpl.cm.jet
		max_phenotype = 1
		min_phenotype = 0
		phenotype_gap = max_phenotype - min_phenotype
		phenotype_jitter = phenotype_gap/10.
		phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
		axe_x_offset4 = 0
		axe_y_offset1 = 0.9
		axe_width4 = 0.1
		axe_height4 = 0.2
		#axe_map_phenotype_legend = pylab.axes([axe_x_offset4+0.02, axe_y_offset1, axe_width4-0.02, axe_height4], frameon=False)
		#cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
		#							norm=phenotype_norm,
		#							orientation='horizontal')
		#cb.set_label('Legend for the color')
		"""
		def drawGraphNodes(G, pos, sex2node_property_list):
			import matplotlib as mpl
			for sex, node_property_ls in sex2node_property_list.iteritems():
				if sex=='M':
					node_shape = 's'
				else:
					node_shape = 'o'
				node_list = []
				node_color_list = []
				node_size_list = []
				for node_property in node_property_ls:
					node, node_size, node_color = node_property[:3]
					node_list.append(node)
					node_size_list.append(node_size)
					node_color_list.append(node_color)
				node_size_ar = numpy.array(node_size_list)
				if sex=='M':
					node_size_ar = node_size_ar*1.5	#by default, the square is smaller than a circle icon.
				
				nx.draw_networkx_nodes(DG, pos, nodelist=node_list, node_color=node_color_list, node_size=node_size_ar, \
									node_shape=node_shape, alpha=1, width=0, linewidths=0, cmap =mpl.cm.jet, vmin=0, vmax=1.0)
		
		drawGraphNodes(DG, pos, sex2node_property_list)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=100)
		
		
		
		pylab.clf()
		pylab.axis("off")
		pylab.figure(figsize=(60, 60))
		G = DG.to_undirected()
		pos = nx.spectral_layout(G)
		nx.draw_networkx_edges(G, pos, alpha=0.9, width=0.5)
		drawGraphNodes(G, pos, sex2node_property_list)
		#nx.draw_networkx_nodes(G, pos, node_color='r', node_size=baseNodeSize, alpha=0.8, width=0, linewidths=0)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_spectral_layout_graph.png'%(outputFnamePrefix), dpi=100)
		"""
		layout = 'twopi'
		pylab.clf()
		nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=300)
		
		layout = 'circo'
		pylab.clf()
		nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=300)
		"""
	
	"""
		
		#2011-5-6
		outputFnamePrefix = os.path.expanduser('~/script/vervet/data/pedigree')
		DBVervet.drawPedigree(db_vervet, outputFnamePrefix=outputFnamePrefix)
		sys.exit(0)
		
		
		
	"""
	
	@classmethod
	def calculateAndDrawPairwiseDistanceInPedigree(cls, db_vervet, inputFname, outputFnamePrefix=None, baseNodeSize=40):
		"""
		2011-5-6
			
		"""
		sys.stderr.write("Getting pedigree from db ...")
		import networkx as nx
		DG=nx.Graph()
		
		import VervetDB
		
		for row in VervetDB.Ind2Ind.query:
			DG.add_edge(row.individual1_id, row.individual2_id)
		sys.stderr.write("%s edges. %s nodes.\n"%(len(DG.edges()), len(DG.nodes()) ))
		
		import csv
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		monkey_id_ls = []
		monkey_seq_id_ls = []
		header = reader.next()
		for row in reader:
			if row and row[0]:
				monkey = VervetDB.Individual.query.filter_by(code=row[0]).first()
				if monkey:
					ind_seq_ls = VervetDB.IndividualSequence.query.filter_by(individual_id=monkey.id).filter_by(sequencer='GA').\
						filter_by(filtered=0).all()
					if len(ind_seq_ls)>1:
						sys.stderr.write("monkey %s has %s sequences.\n"%(monkey.code, len(ind_seq_ls)))
					elif len(ind_seq_ls)==1:
						ind_seq = ind_seq_ls[0]
						monkey_seq_id_ls.append(ind_seq.id)
					monkey_id_ls.append(monkey.id)
		sys.stderr.write("%s monkeys from %s.\n"%(len(monkey_id_ls), inputFname))
		
		monkey_seq_id_ls.sort()
		monkey_seq_id_ls = map(repr, monkey_seq_id_ls)
		print "individual sequence id: " + ",".join(monkey_seq_id_ls)
		
		monkey_id_ls.sort()
		print "individual id: ", monkey_id_ls
		shortest_distance_ls = []
		for i in range(len(monkey_id_ls)):
			for j in range(i+1, len(monkey_id_ls)):
				dist = nx.shortest_path_length(DG,source=monkey_id_ls[i],target=monkey_id_ls[j])
				shortest_distance_ls.append(dist)
		
		from pymodule.yh_matplotlib import drawHist
		outputFname = "%s_%smonkeys.png"%(outputFnamePrefix, len(monkey_id_ls))
		drawHist(shortest_distance_ls, title="hist of shortest dist of %s monkeys"%(len(monkey_id_ls)), \
				xlabel_1D="shortest distance", xticks=None, outputFname=outputFname, min_no_of_data_points=20, needLog=False, \
				dpi=200)
	
	"""
		#2011-9-22
		inputFname = '/tmp/LeastRelatedVervets.csv'
		outputFnamePrefix = '/tmp/LeastRelatedVervets_pedigree_distance'
		DBVervet.calculateAndDrawPairwiseDistanceInPedigree(db_vervet, inputFname, outputFnamePrefix=outputFnamePrefix)
		sys.exit(0)
	"""
	@classmethod
	def filterValue(cls, value, data_type=None, NA_str_set=set(["", "NA", "N/A", 'n/a'])):
		"""
		2011-5-5
			strip the value first.
		2011-4-28
			adapted from variation.src.misc class DBGenome
		"""
		value = value.strip()
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
			site_name_ls = map(strip_func, site_name_ls)
			if len(site_name_ls)==3:
				city, province, country = site_name_ls
				description = None
			elif len(site_name_ls)==4:
				description, city, province, country = site_name_ls
			else:
				sys.stderr.write("Error: site %s is neither 3-entry nor 4-entry.\n"%site_str)
				sys.exit(0)
			site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country, \
									latitude=latitude, longitude=longitude)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex[0], age=age, age_cas=age_cas, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector)
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name not in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:
						if cls.any_character_pattern.search(row[i]):	#ignore any column with any character in it.
							# should be number only
							continue
						value = cls.filterValue(row[i], data_type=float)
						comment =None
					else:
						value = None
						comment = row[i]
					if value or comment:
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
	def put2011SouthAfricanCollectionIntoDB(cls, db_vervet, inputFname, monkeyIDPrefix="", \
										collector_name='Christopher A. Schmitt', province="East Cape", country="South Africa"):
		"""
		2011-10-7
			collection from east cape
		"""
		sys.stderr.write("Putting 2011 south african collection (%s) into db ...\n"%(inputFname))
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule.utils import getColName2IndexFromHeader, figureOutDelimiter
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Monkey ID")
		site_index = col_name2index.get("Location")
		gps_index = col_name2index.get("GPS")
		sex_index = col_name2index.get("Sex")
		age_index = col_name2index.get("Dental Age Category")
		microchip_index = col_name2index.get("Microchip ID")
		collection_date_index = col_name2index.get("Date")
		phenotype_start_index = col_name2index.get("diagnostic physical features")
		
		collector = db_vervet.getUser(collector_name)
		counter = 0
		no_of_phenotype_entries = 0
		for row in reader:
			monkey_id = monkeyIDPrefix + row[monkey_id_index]
			site_str = row[site_index]
			gps = row[gps_index]
			sex = row[sex_index]
			age = row[age_index]
			age = cls.filterValue(age, int)
			collection_date = row[collection_date_index]
			collection_date = datetime.strptime(collection_date, '%m/%d/%y')	#6/20/11
			
			latitude, longitude = gps.split(',')[:2]
			latitude = float(latitude)
			longitude = float(longitude)
			
			microchip_id = row[microchip_index].strip()
			
			description = site_str.strip()
			city = None
			site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country, \
									latitude=latitude, longitude=longitude,\
									altitude=None)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex[0], age=age, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector,tax_id=None, \
								birthdate=None, vrc_founder=None, comment=None, microchip_id=microchip_id)
			counter += 1
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name not in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:
						if cls.any_character_pattern.search(row[i]):	#ignore any column with any character in it.
							# should be number only
							continue
						value = cls.filterValue(row[i], data_type=float)
						comment =None
					else:
						value = None
						comment = row[i]
					if value or comment:
						db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
						no_of_phenotype_entries += 1
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("%s individuals, %s phenotype entries. Done.\n"%(counter, no_of_phenotype_entries))
	
	"""
		#2011-10-7
		inputFname = os.path.expanduser("/tmp/2011SouthAfrica_cas30sep2011.csv")
		DBVervet.put2011SouthAfricanCollectionIntoDB(db_vervet, inputFname)
		sys.exit(0)
	"""
	
	@classmethod
	def put2011GambiaCollectionIntoDB(cls, db_vervet, inputFname, gpsDataInputFname, monkeyIDPrefix="", \
						collector_name='TeamGambia', province=None, country="Gambia"):
		"""
		2011-10-7
			collection from Gambia
		"""
		db_vervet.session.begin()
		import csv, re
		from datetime import datetime
		from pymodule.utils import getColName2IndexFromHeader, figureOutDelimiter
		
		sys.stderr.write("\t Getting gps&site data from %s ..."%(gpsDataInputFname))
		reader = csv.reader(open(gpsDataInputFname,), delimiter=figureOutDelimiter(gpsDataInputFname))
		header = reader.next()
		
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		site_index = col_name2index.get("LocationID")
		gps_index = col_name2index.get("GPS")
		gps_pattern = re.compile(r"(?P<lat_direction>[SN]) (?P<lat_hr>\d+) (?P<lat_min>[\.\d]+)' (?P<lon_direction>[EW]) (?P<lon_hr>\d+) (?P<lon_min>[\.\d]+)'")

		#2011-4-29 first get the location infomation into db.
		# a dictionary which maps code to site and collection_date
		site_str2gps = {}
		counter = 0
		for row in reader:
			gps = row[gps_index]
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
			elif gps:
				gps_ls = gps.split(',')
				gps_ls = map(float, gps_ls)
				latitude, longitude = gps_ls[:2]
			else:
				sys.stderr.write("Error: gps %s is not parsable.\n"%(gps))
				sys.exit(0)
			site_str = row[site_index].strip()
			city = None
			description = site_str
			#site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country,\
			#					latitude=latitude, longitude=longitude)
			counter += 1
			site_str2gps[site_str] = (latitude, longitude)
		del reader
		sys.stderr.write("%s sites found.\n"%(counter))
		
		sys.stderr.write("Putting 2011 Gambia collection (%s) into db ...\n"%(inputFname))
		
		reader = csv.reader(open(inputFname,), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("Monkey ID")
		site_index = col_name2index.get("Location")
		gps_index = col_name2index.get("GPS")
		sex_index = col_name2index.get("Sex")
		age_index = col_name2index.get("Dental Age Category")
		age_cas_index = col_name2index.get("CAS Dental Age")
		
		
		microchip_index = col_name2index.get("Microchip ID")
		collection_date_index = col_name2index.get("Date")
		phenotype_start_index = col_name2index.get("diagnostic physical features")
		
		collector = db_vervet.getUser(collector_name)
		counter = 0
		no_of_phenotype_entries = 0
		for row in reader:
			monkey_id = monkeyIDPrefix + row[monkey_id_index]
			site_str = row[site_index].strip()
			sex = row[sex_index]
			age = row[age_index]
			age = cls.filterValue(age, int)
			age_cas = row[age_cas_index]
			age_cas = cls.filterValue(age_cas, int)
			collection_date = row[collection_date_index]
			collection_date = datetime.strptime(collection_date, '%m/%d/%y')	#6/20/11
			
			microchip_id = row[microchip_index].strip()
			
			gps = site_str2gps.get(site_str)
			if gps:
				latitude, longitude = gps[:2]
			else:
				latitude, longitude = None, None
			description = site_str
			city = None
			site = db_vervet.getSite(description=description, city=city, stateprovince=province, country_name=country, \
									latitude=latitude, longitude=longitude,\
									altitude=None)
			
			individual = db_vervet.getIndividual(code=monkey_id, sex=sex[0], age=age, age_cas=age_cas, latitude=latitude,\
								longitude=longitude, altitude=None, ucla_id=monkey_id, site=site, \
								collection_date=collection_date, collector=collector,tax_id=None, \
								birthdate=None, vrc_founder=None, comment=None, microchip_id=microchip_id)
			counter += 1
			for i in range(phenotype_start_index, len(row)):
				if header[i] and row[i]:
					phenotype_name = header[i].strip()
					if phenotype_name not in ['Age/Dentition present', 'diagnostic physical features', 'Additional comments']:
						if cls.any_character_pattern.search(row[i]):	#ignore any column with any character in it.
							# should be number only
							continue
						value = cls.filterValue(row[i], data_type=float)
						comment =None
					else:
						value = None
						comment = row[i]
					if value or comment:
						db_vervet.getPhenotype(phenotype_name=phenotype_name, value=value, replicate=None, \
										individual_id=individual.id, comment=comment, collector_name=collector_name)
						no_of_phenotype_entries += 1
		db_vervet.session.flush()
		db_vervet.session.commit()
		sys.stderr.write("%s individuals, %s phenotype entries. Done.\n"%(counter, no_of_phenotype_entries))
	
	"""
		#2011-10-7
		inputFname = os.path.expanduser("/tmp/gambia 2011 animal datasheet08012011cas.csv")
		gpsDataInputFname = os.path.expanduser("/tmp/Gambia2011Location_GPS_AJJ.csv")
		DBVervet.put2011GambiaCollectionIntoDB(db_vervet, inputFname, gpsDataInputFname)
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
		
		
		# to match either "June 23-25" or "2009/7/16"
		collection_date_pattern = re.compile(r'((?P<Month>\w+) (?P<Day>\d+)-\d+)|((?P<year>\d+)/(?P<month>\d+)/(?P<day>\d+))')
		
		# to match "Yearling", "Adult", "15 yrs", "5 years"
		age_pattern = re.compile(r"([a-z A-Z]+)|((?P<age>\d+) ((yrs)|(years)))")
		counter = 0
		real_counter = 0
		species2tax_id = {'sabaeus':60711, 'aethiops':101841, 'pygerythrus':460674, 'cynosurus':460675, 'tantalus':60712,\
						'cynosuros':460675, 'pygerythrus':460674 }	#typos
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
			if value and not cls.any_character_pattern.search(value):	#ignore any column with any character in it.
					# should be number only:
				phenotype_name = "Weight"
				value = cls.filterValue(value, data_type=float)
				comment =None
				if value or comment:
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
	@classmethod
	def putDOCWalkerResultsIntoDB(cls, db_vervet, inputFname, commit=True):
		"""
		2011-11-28
			inputFname is output of InspectAlignmentPipeline.py.
			
			this function updates IndividualAlignment.pass_qc_read_base_count and median_depth, mean_depth
		"""
		import VervetDB
		import csv
		from pymodule import PassingData
		db_vervet.session.begin()
		reader = csv.reader(open(inputFname,), delimiter='\t')
		
		sys.stderr.write("\t putting coverage data from %s into db ..."%(inputFname))
		header = reader.next()
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		sample_id_index = col_name2index.get("sample_id")
		total_base_count_index = col_name2index.get('total')
		mean_depth_index = col_name2index.get("mean")
		median_depth_index = col_name2index.get("granular_median")
		counter = 0
		for row in reader:
			sample_id = row[sample_id_index]
			if sample_id=='Total':	#ignore rows with this as sample id
				continue
			alignment_id = int(sample_id.split("_")[0])
			total_base_count = int(row[total_base_count_index])
			mean_depth = float(row[mean_depth_index])
			median_depth = float(row[median_depth_index])
			individual_alignment = VervetDB.IndividualAlignment.get(alignment_id)
			individual_alignment.pass_qc_read_base_count = total_base_count
			individual_alignment.mean_depth = mean_depth
			individual_alignment.median_depth = median_depth
			db_vervet.session.add(individual_alignment)
			db_vervet.session.flush()
			counter += 1
		if commit:
			db_vervet.session.commit()
		sys.stderr.write("%s alignments updated. Done.\n"%(counter))
	
	"""
		#2011-11-28
		inputFname = '/Network/Data/vervet/vervetPipeline/InspectRefSeq524WholeAlignment_2011.11.25T2311/DepthOfCoverage.tsv'
		DBVervet.putDOCWalkerResultsIntoDB(db_vervet, inputFname)
		sys.exit(0)
	"""
	
class VervetGenome(object):
	"""
	2011-6-27
		
	"""
	def __init__(self):
		pass
	
	@classmethod
	def getContigID2SizeFromAGPFile(cls, contigAGPFname):
		"""
		2011-6-27
			a helper function used by other functions
		"""
		sys.stderr.write("Getting information about all contigs from %s ..."%(contigAGPFname))
		contig_id2size = {}
		import csv
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(contigAGPFname, 'r'), delimiter=figureOutDelimiter(contigAGPFname))
		for row in reader:
			contig_id = row[0]
			if contig_id not in contig_id2size:
				contig_id2size[contig_id] = 0
			contig_end_pos = int(row[2])
			if contig_end_pos>contig_id2size[contig_id]:
				contig_id2size[contig_id] = contig_end_pos
			
		del reader
		sys.stderr.write("%s contigs. Done.\n"%(len(contig_id2size)))
		return contig_id2size
		
	
	@classmethod
	def drawHistogramOfNFractionInContig(cls, contigAGPFname=None, outputFnamePrefix=None, topNumber=150, contig_id_set=None):
		"""
		2011-6-27
			In contigAGPFname, each contig is a supercontig composed of a couple of small contigs.
			format:
				Contig0 1       1865    1       W       Contig0.1       1       1865    +
				Contig0 1866    1941    2       N       76      fragment        yes
				Contig0 1942    3526    3       W       Contig0.2       1       1585    +
				Contig0 3527    4532    4       N       1006    fragment        yes
			
		"""
		contig_id2size = cls.getContigID2SizeFromAGPFile(contigAGPFname)
		
		size_contig_id_ls = []
		for contig_id, contig_size in contig_id2size.iteritems():
			size_contig_id_ls.append((contig_size, contig_id))
		size_contig_id_ls.sort()
		
		selected_contig_id_set = set()
		for i in xrange(topNumber):
			contig_id = size_contig_id_ls[-i-1][1]
			if len(size_contig_id_ls)<=i:
				break
			if contig_id_set and contig_id not in contig_id_set:
				break
			selected_contig_id_set.add(contig_id)
		sys.stderr.write("%s contigs selected.\n"%(len(selected_contig_id_set)))
		
		sys.stderr.write("Getting base-data for N-percentage in all contigs from %s ..."%(contigAGPFname))
		contig_id2size_tuple = {}
		import csv,  math
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(contigAGPFname, 'r'), delimiter=figureOutDelimiter(contigAGPFname))
		for row in reader:
			contig_id = row[0]
			base_type = row[4]
			if contig_id not in selected_contig_id_set:
				continue
			if contig_id not in contig_id2size_tuple:
				contig_id2size_tuple[contig_id] = [0,0]	#(#bases that are not N, #bases that are N)
			contig_span = int(row[2])-int(row[1]) + 1
			if base_type=='N':
				contig_id2size_tuple[contig_id][1] += contig_span
			else:
				contig_id2size_tuple[contig_id][0] += contig_span
			
		del reader
		sys.stderr.write("%s contigs. Done.\n"%(len(contig_id2size_tuple)))
	
		sys.stderr.write("Calculating N percentage for %s contigs ..."%(len(selected_contig_id_set)))
		N_percentage_ls = []
		scaffoldSizeLs = []
		for contig_id in selected_contig_id_set:
			if contig_id in contig_id2size_tuple:
				W_size, N_size = contig_id2size_tuple.get(contig_id)
				total_size = W_size + N_size
				if total_size>0:
					N_percentage = N_size/float(total_size)
					N_percentage_ls.append(N_percentage)
					scaffoldSizeLs.append(math.log10(total_size))
		
		sys.stderr.write("Drawing histogram of contig N-percentage ...")
		"""
		contig_size_ls.sort()
		for i in xrange(len(contig_size_ls)):
			contig_size_ls[i] = math.log10(contig_size_ls[i])
		"""
		from pymodule import yh_matplotlib
		yh_matplotlib.drawHist(N_percentage_ls, title="histogram of N-percentage of top %s contigs"%(len(N_percentage_ls)), \
							xlabel_1D="N-percentage", outputFname='%s.png'%(outputFnamePrefix), min_no_of_data_points=50, needLog=True)
		sys.stderr.write("Done.\n")
		
		fig_fname = '%s_2D_size_hist.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(N_percentage_ls, scaffoldSizeLs, [1]*len(scaffoldSizeLs), fig_fname=fig_fname, \
				gridsize=20, title="%s contigs in N-fraction X size space"%(len(scaffoldSizeLs)), \
				xlabel="N-fraction", ylabel="log10(size)",\
				colorBarLabel="log10(count)", reduce_C_function=yh_matplotlib.logSum, dpi=200)
	
	
	"""
		#2011-7-5
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		topNumber=200000
		outputFnamePrefix = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/top%ssupercontigs_N_percentage"%topNumber)
		VervetGenome.drawHistogramOfNFractionInContig(contigAGPFname=contigAGPFname, outputFnamePrefix=outputFnamePrefix, \
													topNumber=topNumber)
		sys.exit(3)
		
		contig_id_set = VariantDiscovery.getBACEndHitsOnContigs(inputFname=inputFname, outputFnamePrefix=outputFnamePrefix, min_alen_rlen_ratio=min_alen_rlen_ratio)
		sys.exit(0)
		
		
		#2011-7-5
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		topNumber=485000
		outputFnamePrefix = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/MinSize200ScaffoldsWithBACEndHits_N_percentage"%topNumber)
		VervetGenome.drawHistogramOfNFractionInContig(contigAGPFname=contigAGPFname, outputFnamePrefix=outputFnamePrefix, \
													topNumber=topNumber, contig_id_set=contig_id_set)
		sys.exit(0)
		
		
		#2011-7-5
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		topNumber=485000
		outputFnamePrefix = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/MinSize200ScaffoldsWithBACEndHits_N_percentage"%topNumber)
		VervetGenome.drawHistogramOfNFractionInContig(contigAGPFname=contigAGPFname, outputFnamePrefix=outputFnamePrefix, \
													topNumber=topNumber, contig_id_set=contig_id_set)
		sys.exit(0)
		
	"""
	
	@classmethod
	def drawContigSizeHistogram(cls, contigAGPFname=None, outputFnamePrefix=None):
		"""
		2011-6-27
			In contigAGPFname, each contig is a supercontig composed of a couple of small contigs.
			format:
				Contig0 1       1865    1       W       Contig0.1       1       1865    +
				Contig0 1866    1941    2       N       76      fragment        yes
				Contig0 1942    3526    3       W       Contig0.2       1       1585    +
				Contig0 3527    4532    4       N       1006    fragment        yes
			
		"""
		contig_id2size = cls.getContigID2SizeFromAGPFile(contigAGPFname)
		sys.stderr.write("Drawing histogram of contig sizes ...")
		contig_size_ls = contig_id2size.values()
		contig_size_ls.sort()
		for i in xrange(len(contig_size_ls)):
			contig_size_ls[i] = math.log10(contig_size_ls[i])
		
		from pymodule import yh_matplotlib
		yh_matplotlib.drawHist(contig_size_ls, title="histogram of size of %s contigs"%(len(contig_size_ls)), \
							xlabel_1D="log10(size)", outputFname='%s.png'%(outputFnamePrefix), min_no_of_data_points=50, needLog=True)
		sys.stderr.write("Done.\n")
		
	"""
		#2011-6-27
		contigAGPFname = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/Draft_June_2011/supercontigs.agp")
		outputFnamePrefix = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigSizeHistogram")
		VervetGenome.drawContigSizeHistogram(contigAGPFname, outputFnamePrefix)
		sys.exit(2)
		
		#2011-6-27
		contigAGPFname = os.path.expanduser("~/mnt/hoffman2_home/script/vervet/data/Draft_June_2011/supercontigs.agp")
		outputFnamePrefix = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigSizeHistogram")
		VervetGenome.drawContigSizeHistogram(contigAGPFname, outputFnamePrefix)
		sys.exit(2)
		
		
		
	"""
	
	@classmethod
	def outputTopBigContigs(cls, contigAGPFname, contigFastaFname, outputFname=None, topNumber=150):
		"""
		2011-6-27
			
		"""
		contig_id2size = cls.getContigID2SizeFromAGPFile(contigAGPFname)
		size_contig_id_ls = []
		for contig_id, contig_size in contig_id2size.iteritems():
			size_contig_id_ls.append((contig_size, contig_id))
		size_contig_id_ls.sort()
		
		selected_contig_id_set = set()
		for i in xrange(topNumber):
			if len(size_contig_id_ls)<=i:
				break
			contig_id = size_contig_id_ls[-i-1][1]
			selected_contig_id_set.add(contig_id)
		sys.stderr.write("%s contigs selected.\n"%(len(selected_contig_id_set)))
		
		sys.stderr.write("Picking selected contigs and their sequences ...")
		sequences = []
		from Bio import SeqIO
		handle = open(contigFastaFname, "rU")
		for record in SeqIO.parse(handle, "fasta"):
			contig_id = record.id.split()[0]
			if contig_id in selected_contig_id_set:
				record.id = contig_id	#superfluous. record.id is 1st word of description.
				record.description = contig_id
				sequences.append(record)
			if len(sequences)>=len(selected_contig_id_set):
				break
		handle.close()
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Writing output to %s ..."%outputFname)
		output_handle = open(outputFname, "w")
		SeqIO.write(sequences, output_handle, "fasta")
		output_handle.close()
		sys.stderr.write("Done.\n")
		
	"""
		#2011-6-27
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs.agp")
		contigFastaFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs.fasta")
		topNumber=150
		outputFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/top%ssupercontigs.fasta"%topNumber)
		VervetGenome.outputTopBigContigs(contigAGPFname=contigAGPFname, contigFastaFname=contigFastaFname, \
				outputFname=outputFname, topNumber=topNumber)
		sys.exit(3)
		
		#2011-6-27
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs.agp")
		contigFastaFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs.fasta")
		topNumber=150
		outputFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/top%ssupercontigs.fasta"%topNumber)
		VervetGenome.outputTopBigContigs(contigAGPFname=contigAGPFname, contigFastaFname=contigFastaFname, \
				outputFname=outputFname, topNumber=topNumber)
		sys.exit(3)
		
	"""
	
	@classmethod
	def outputContigIDSizeInBED(cls, contigAGPFname, outputFname=None, minSize=2000):
		"""
		2011-11-3
		"""
		contig_id2size = cls.getContigID2SizeFromAGPFile(contigAGPFname)
		size_contig_id_ls = []
		for contig_id, contig_size in contig_id2size.iteritems():
			size_contig_id_ls.append((contig_size, contig_id))
		size_contig_id_ls.sort()
		size_contig_id_ls.reverse()	#big ones come up first
		import csv
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		for contig_size, contig_id in size_contig_id_ls:
			if contig_size>=minSize:
				writer.writerow([contig_id, 0, contig_size])
		del writer
		
	"""
		#2011-11-3
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		outputFname = os.path.expanduser("/tmp/ContigIDSize.bed")
		VervetGenome.outputContigIDSizeInBED(contigAGPFname, outputFname)
		sys.exit(0)
		
	"""
	
	@classmethod
	def outputContigsAboveCertainSize(cls, contigFastaFname, outputFname=None, minSize=200):
		"""
		2011-6-27
			
		"""
		
		sys.stderr.write("Picking selected contigs and their sequences ...")
		from Bio import SeqIO
		from pymodule import PassingData
		handle = open(contigFastaFname, "rU")
		input_seq_iterator = SeqIO.parse(handle, 'fasta')
		passingData = PassingData(input_seq_iterator =input_seq_iterator, counter=0, no_of_good_records=0)
		#setattr(input_seq_iterator, 'no_of_good_records', 0)
		#setattr(input_seq_iterator, 'counter', 0)
		def filterRecords(passingData, minSize):
			for record in passingData.input_seq_iterator:
				passingData.counter += 1
				if len(record)>=minSize:
					contig_id = record.id.split()[0]
					record.description = contig_id
					#record.id = contig_id
					passingData.no_of_good_records += 1
					yield record
		#seq_iterator = (record for record in input_seq_iterator if len(record.seq) >=minSize)
		seq_iterator = filterRecords(passingData, minSize)
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Writing output to %s ..."%outputFname)
		output_handle = open(outputFname, "w")
		SeqIO.write(seq_iterator, output_handle, "fasta")
		output_handle.close()
		sys.stderr.write("%s/%s records. Done.\n"%(passingData.no_of_good_records, passingData.counter))
		
	"""
		#2011-6-27
		contigFastaFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.fasta")
		minSize = 200
		outputFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/superContigsMinSize%s.fasta"%minSize)
		VervetGenome.outputContigsAboveCertainSize(contigFastaFname=contigFastaFname, \
				outputFname=outputFname, minSize=200)
		sys.exit(3)
		
		
			#2011-6-27
		contigFastaFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.fasta")
		minSize = 2000
		outputFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/superContigsMinSize%s.fasta"%minSize)
		VervetGenome.outputContigsAboveCertainSize(contigFastaFname=contigFastaFname, \
				outputFname=outputFname, minSize=minSize)
		sys.exit(3)
		
	"""
	
	class ContigFastaRecordFilter(object):
		"""
		2011-7-7
			a class used by FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles()
				as fastaRecordFilterHandler.
				
			The whole class is instantiated in outputTopNumberContigsIntoSeparateFiles().
		"""
		def __init__(self, **keywords):
			"""
			2011-7-7
				keywords include argument maxContigNumber
			"""
			# 2011-7-7
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, fastaTitle=""):
			"""
			2011-7-7
			"""
			import re
			p_contig_id = re.compile(r'Contig(\d+)')
			if p_contig_id.search(fastaTitle):
				contig_number = int(p_contig_id.search(fastaTitle).groups()[0])
				if contig_number<=self.maxContigNumber:
					return True
				else:
					return False
	
	@classmethod
	def outputTopNumberContigsIntoSeparateFiles(cls, contigFastaFname, outputDir=None, maxContigNumber=193):
		"""
		2011-7-7
			for symap
		"""
		contigFastaRecordFilter = VervetGenome.ContigFastaRecordFilter(maxContigNumber=maxContigNumber)
		from variation.src.misc import FileFormatExchange
		FileFormatExchange.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(contigFastaFname, outputDir, \
														fastaRecordFilterHandler=contigFastaRecordFilter)
		
	"""
		#2011-7-7
		contigFastaFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.fasta")
		outputDir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/vervetScaffolds/sequence/pseudo/")
		maxContigNumber = 193
		VervetGenome.outputTopNumberContigsIntoSeparateFiles(contigFastaFname, outputDir=outputDir, maxContigNumber=maxContigNumber)
		sys.exit(0)
		
	"""
	
	@classmethod
	def outputVervetMSSequenceHumanChr(cls, db_genome, vervetMSTemplateFname, vervetMSHg19CoordinatesFname, outputFname):
		"""
		2011-10-24
			vervetMSHg19CoordinatesFname was a custom track (microsatellites) downloaded from ucsc vervet browser (link from mcgill site).
				Its start is 0-based.
				Its end is 1-based.
		"""
		import csv, sys, os, gzip
		from pymodule.utils import figureOutDelimiter, getColName2IndexFromHeader
		sys.stderr.write("Getting human genome coordinates of vervet MS sequences from %s ...\n"%(vervetMSHg19CoordinatesFname))
		inf = gzip.open(vervetMSHg19CoordinatesFname, 'r')
		reader = csv.reader(inf, delimiter='\t')
		ms_id2hg19_coords = {}
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		ms_id_index = col_name2index.get('name')
		hg19_chr_index = col_name2index.get('#chrom')
		hg19_start_index = col_name2index.get('chromStart')
		hg19_end_index = col_name2index.get('chromEnd')
		for row in reader:
			ms_id = row[ms_id_index]
			if ms_id[-2]==".":	#get rid of version in D7S1820.0
				ms_id = ms_id[:-2]
			hg19_chr = row[hg19_chr_index][3:]	#3: is to get rid of "chr"
			hg19_start = int(row[hg19_start_index])+1	# 0-based start
			hg19_end = int(row[hg19_end_index])	#1-based end
			if ms_id in ms_id2hg19_coords:
				sys.stderr.write("error: %s already in ms_id2hg19_coords.\n"%(ms_id))
				sys.exit(3)
			ms_id2hg19_coords[ms_id] = [hg19_chr, hg19_start, hg19_end]
		del reader
		sys.stderr.write("%s microsatellites. Done.\n"%(len(ms_id2hg19_coords)))
		
		
		sys.stderr.write("Outputting vervet microsatellite sequences and human chromosomes to %s ...\n"%outputFname)
		
		delimiter = figureOutDelimiter(vervetMSTemplateFname)
		reader = csv.reader(open(vervetMSTemplateFname, 'r'), delimiter=delimiter)
		header =reader.next()
		header.append('hg19_start')
		header.append('hg19_stop')
		writer = csv.writer(open(outputFname, 'w'), delimiter=delimiter)
		writer.writerow(header)
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		ms_id_index = col_name2index.get('Marker')
		Vchrom_index = col_name2index.get("Vchrom")
		Hchromosome_index = col_name2index.get("Hchromosome")
		mapKoscM_index = col_name2index.get("mapKoscM")
		Mapped1000to1_index = col_name2index.get("Mapped1000to1")
		human_sequence_index = col_name2index.get("human sequence")
		hg19_start_index = col_name2index.get('hg19_start')
		hg19_end_index = col_name2index.get('hg19_stop')
		for row in reader:
			ms_id = row[ms_id_index]
			Vchrom = row[Vchrom_index]
			mapKoscM = row[mapKoscM_index]
			
			new_row = ['']*len(header)
			new_row[ms_id_index] = ms_id
			new_row[Vchrom_index] = Vchrom
			new_row[mapKoscM_index] = mapKoscM
			new_row[Mapped1000to1_index] = row[Mapped1000to1_index] 
			chr_start_stop = ms_id2hg19_coords.get(ms_id)
			if chr_start_stop is None:
				sys.stderr.write("Warning: %s not in %s. Output as it is.\n"%(ms_id, vervetMSHg19CoordinatesFname))
			else:
				chr, start, stop = chr_start_stop[:3]
				seq = db_genome.getSequenceSegment(tax_id=9606, chromosome=chr, start=start, stop=stop)
				new_row[Hchromosome_index] = chr
				new_row[human_sequence_index] = seq
				new_row[hg19_start_index] = start
				new_row[hg19_end_index] = stop
			writer.writerow(new_row)
		del reader, writer
		sys.stderr.write("Done.\n")
	"""
		#2011-10-25
		vervetMSTemplateFname = '/Network/Data/vervet/microsatellites/vervet microsatellite info.csv'
		vervetMSHg19CoordinatesFname = '/Network/Data/vervet/microsatellites/vervetMicrosatellitesOnHg19.gz'
		outputFname = '/Network/Data/vervet/microsatellites/vervetMicrosatelliteInfo.csv'
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		vervetMSTemplateFname = '/Network/Data/vervet/microsatellites/vervet microsatellite info.csv'
		vervetMSTemplateFname = '/Network/Data/vervet/microsatellites/vervetstrmap.csv'
		vervetMSHg19CoordinatesFname = '/Network/Data/vervet/microsatellites/vervetMicrosatellitesOnHg19.gz'
		outputFname = '/Network/Data/vervet/microsatellites/vervetStrMap.hsSeq.csv'
		VervetGenome.outputVervetMSSequenceHumanChr(db_genome, vervetMSTemplateFname, vervetMSHg19CoordinatesFname, outputFname)
		sys.exit(0)
	"""
	
	@classmethod
	def addMoreContigsAsChromosomeRecordsInGenomeDB(cls, db_genome, contigAGPFname, tax_id=60711,
												sequence_type_name='Scaffold', minSize=500, **keywords):
		"""
		2011-11-6
			only the name, start, stop are saved. no raw sequence inserted into db.
			This is just for programs that rely on this AnnotAssembly table to fetch information about these contigs.
		"""
		
		contig_id2size = cls.getContigID2SizeFromAGPFile(contigAGPFname)
		size_contig_id_ls = []
		for contig_id, contig_size in contig_id2size.iteritems():
			size_contig_id_ls.append((contig_size, contig_id))
		size_contig_id_ls.sort()
		size_contig_id_ls.reverse()	#big ones come up first
		
		sys.stderr.write("Adding contigs as chromosome records (no sequence) into genome db (minSize=%s) ...\n"%(minSize))
		counter = 0
		real_counter = 0
		from pymodule import GenomeDB
		
		db_genome.session.begin()
		sequence_type = GenomeDB.SequenceType.query.filter_by(type=sequence_type_name).first()
		for contig_size, contig_id in size_contig_id_ls:
			if contig_size>=minSize:
				counter += 1
				aa_attr_instance = GenomeDB.AnnotAssembly.query.filter_by(chromosome=contig_id).filter_by(tax_id=tax_id).\
					filter_by(start=1).filter_by(stop=contig_size).\
					filter_by(sequence_type_id=sequence_type.id).first()
				if not aa_attr_instance:
					real_counter += 1
					aa_attr_instance = GenomeDB.AnnotAssembly()
					aa_attr_instance.gi = None
					aa_attr_instance.acc_ver = None
					aa_attr_instance.tax_id = tax_id
					aa_attr_instance.chromosome = contig_id
					aa_attr_instance.start = 1
					aa_attr_instance.sequence_type_id = sequence_type.id
					aa_attr_instance.stop = contig_size
					db_genome.session.add(aa_attr_instance)
					db_genome.session.flush()
		db_genome.session.commit()
		sys.stderr.write("%s/%s contigs saved.\n"%(real_counter, counter))
		"""
		#2011-11-6 add more contig records into db
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		VervetGenome.addMoreContigsAsChromosomeRecordsInGenomeDB(db_genome, contigAGPFname, maxRank=1001, tax_id=60711,
												sequence_type_name='Scaffold', minSize=500)
		sys.exit(0)
		
				#2011-11-6 add more contig records into db
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		contigAGPFname = os.path.expanduser("~/script/vervet/data/Draft_June_2011/supercontigs/supercontigs.agp")
		VervetGenome.addMoreContigsAsChromosomeRecordsInGenomeDB(db_genome, contigAGPFname, maxRank=1001, tax_id=60711,
												sequence_type_name='Scaffold', minSize=1000)
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
							('input', 0, ): ['', 'i', 1, 'common input.', ],\
							('output', 0, ): ['', 'o', 1, 'common output', ],\
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
		
		#2011.12.9
		inputDir='/Network/Data/vervet/vervetPipeline/work/'
		workflowName = 'FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25_2011.12.1T1155'
		workflowName = ''
		workflowName = 'Filter_Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25.2011.12.9T0643'
		FilterVCFByDepthStdoutFnamePattern="%s/merge_workflow-FilterVCFByDepth*out.000"%(workflowName)
		FilterVCFByDepthStdoutFnamePattern = os.path.join(inputDir, FilterVCFByDepthStdoutFnamePattern)
		grepPattern='call_vcftoolsFilter\/Contig'
		#grepPattern='call_vcftoolsFilter_vcftoolsFilter\/Contig'
		VariantDiscovery.TallyFilteredSNPs.findSitesFilteredOrMaskedByFilterVCFByDepth(FilterVCFByDepthStdoutFnamePattern=FilterVCFByDepthStdoutFnamePattern, \
			sampleSize=101, grepPattern=grepPattern)
		sys.exit(0)
		
		#2011.12.9
		inputDir='/Network/Data/vervet/vervetPipeline/scratch/'
		workflowName = 'FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25_2011.12.1T1155'
		workflowName = 'Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4.2011.12.9T0505'
		workflowName = 'Filter_Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25.2011.12.9T0643'
		subFolder="%s/call_vcftoolsFilter"%(workflowName)
		subFolder="%s/call_vcftoolsFilter_vcftoolsFilter"%(workflowName)
		inputDir = os.path.join(inputDir, subFolder)
		fileSuffix='filter_by_vcftools.log'
		#fileSuffix='keepGivenSNP.log'
		VariantDiscovery.TallyFilteredSNPs.countTotalNoOfSitesFilteredByVCFtools(inputDir=inputDir, fileSuffix=fileSuffix)
		sys.exit(0)
		
		
		
		#2011-11-6
		inputDirLs = [os.path.expanduser("~/NetworkData/vervet/db/individual_alignment/"),]
		outputDir = os.path.expanduser("~/panfs/NetworkData/vervet/db/individual_alignment/")
		VariantDiscovery.moveFinishedBamIntoTargetFolder(inputDirLs=inputDirLs, outputDir=outputDir, \
														targetFolderSizeThresholdInMB=1300000, timeGapInMinutes=2400)
		sys.exit(0)
		
		#2011-3-24
		inputPrefix = os.path.expanduser("/Network/Data/vervet/db/individual_alignment/577_37_vs_524_by_2")
		inputFname = os.path.expanduser("%s.bam"%(inputPrefix))
		scoreType = 1
		plotType = 1
		outputFnamePrefix = os.path.expanduser("/tmp/%s.score%s.plot%s"%(os.path.basename(inputPrefix), scoreType, plotType))
		VariantDiscovery.drawHistogramOfPairEndBWAOutputScore(inputFname, outputFnamePrefix, scoreType=scoreType, \
						plotType=plotType, exitAfterNumberOfReads=500000)
		sys.exit(0)
		
		#2011-11-2
		workflowName = 'AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top156Contigs_condor_20111101T2316'
		inputFname = '/Network/Data/vervet/vervetPipeline/%s/call/Contig0.vcf.gz'%(workflowName)
		outputFnamePrefix = '/tmp/%s_Contig0'%(workflowName)
		VariantDiscovery.countHomoHetCallsForEachSampleFromVCF(inputFname, outputFnamePrefix)
		sys.exit(0)
		
		workflowName= '8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA_20111014T0043'
		#workflowName = '8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA_20111014T0059'
		inputDir= "/Network/Data/vervet/vervetPipeline/%s/pairwiseDistMatrix/"%(workflowName)
		for subspeciesName in ['ref', 'Barbados', 'VRC_ref_454', "VRC_ref_GA", 'aethiops', 'cynosurus', 'sabaeus', 'tantalus','pygerythrus',]:
			outputFnamePrefix = '/Network/Data/vervet/vervetPipeline/%s/GATK_all_bases_maxNA0_minMAF0_het2NA_DistVectorFrom8Genomes2%s'%\
				(workflowName, subspeciesName)
			if subspeciesName=='ref':
				refID = subspeciesName
			elif subspeciesName =='VRC_ref_GA' or subspeciesName=='VRC_ref_454':
				refID = '%s_vs_top156Contigs'%(subspeciesName)
			else:
				refID = '%s_GA_vs_top156Contigs'%(subspeciesName)
			VariantDiscovery.drawContigByDistVectorFromOtherGenomes(inputDir, outputFnamePrefix, refID=refID, subspeciesName=subspeciesName)
		sys.exit(0)
		
		
		#2011-9-22
		inputFname = '/tmp/LeastRelatedVervets.csv'
		outputFnamePrefix = '/tmp/LeastRelatedVervets_pedigree_distance'
		DBVervet.calculateAndDrawPairwiseDistanceInPedigree(db_vervet, inputFname, outputFnamePrefix=outputFnamePrefix)
		sys.exit(0)
		
		#2011-9-23
		workflowDir = 'AlignmentToCallPipeline_4_8_vs_524_top_156Contigs_uschpc'
		workflowDir = 'AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110920T2319'
		workflowDir = 'AlignmentToCallPipeline_10VWP_627_629_650_656_vs_524_top_156Contigs_condorpool_20110922T1837'
		inputDir = '/Network/Data/vervet/vervetPipeline/%s/samtools/'%workflowDir
		outputFname = '/Network/Data/vervet/vervetPipeline/%s/samtools_TsTv_hist.png'%(workflowDir)
		VariantDiscovery.drawContigTsTvHistogram(inputDir, outputFname)
		sys.exit(0)
		
		
		#2011-9-15
		dataDir = os.path.expanduser("~/mnt/hpc-cmb_home/NetworkData/vervet/db/")
		DBVervet.pokeBamReadGroupPresence(db_vervet, samtools_path=os.path.expanduser("~/bin/samtools"), dataDir=dataDir, commit=True)
		sys.exit(0)
		
		
		#2011-8-26
		
		inputDir1 = '/Network/Data/vervet/vervetPipeline/8GenomeVsTop156Contigs_GATK_all_bases/pairwiseDistMatrix/'
		inputDir2 = '/Network/Data/vervet/vervetPipeline/workflow_8GenomeVsTop156Contigs_GATK/call/'
		outputFnamePrefix = '/Network/Data/vervet/vervetPipeline/8GenomeVsTop156Contigs_GATK_ContigByDistVectorFrom8Genomes_all_bases_vs_variants_only'
		VariantDiscovery.compareContigByDistVectorFromTwoDifferentRuns(inputDir1, inputDir2, outputFnamePrefix, partOfTitle='all_sites vs variants only')
		sys.exit(3)
		
		#2011-8-2
		inputDir = '/Network/Data/vervet/vervetPipeline/workflow_8GenomeVsTop156Contigs_GATK/call/'
		outputFnamePrefix = '/Network/Data/vervet/vervetPipeline/workflow_8GenomeVsTop156Contigs_GATK/contigPCAByDistVector'
		VariantDiscovery.PCAContigByDistVectorFromOtherGenomes(inputDir, outputFnamePrefix)
		sys.exit(3)
		
		
		
		#2011-7-7 output all the BACs in order
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_genome.setup(create_tables=False)
		outputDir = os.path.expanduser("~/script/variation/bin/symap_3_3/data/pseudo/vervet176BAC/sequence/pseudo/")
		db_genome.outputGenomeSequence(tax_id=60711, sequence_type_id=10, fastaTitlePrefix='BAC', outputDir=outputDir, chunkSize=70)
		sys.exit(0)
		
		
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Main
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
