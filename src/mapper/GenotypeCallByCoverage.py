#!/usr/bin/env python
"""
Examples:
	# call one genome (nothing because no polymorphic sites)
	%s -i /tmp/vs_top150Contigs_by_aln_Contig0.bam -n 1 -o /tmp/vs_top150Contigs_by_aln_Contig0.call
		-y2 -e ...
	
	# call genotype from 8 genomes
	%s -n 8 -i script/vervet/data/8_genome_vs_Contig0.RG.bam -o script/vervet/data/8_genome_vs_Contig0.RG.call
		-y2 -e ...
	
	# 2011-7-21 call from GATK vcf file
	%s -n 8 -y1 -i ~/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC.GATK.vcf
		-o ~/script/vervet/data/1MbBAC_as_ref/454_illu_6_sub_vs_1MbBAC.GATK.call
		-e /Network/Data/vervet/db/individual_sequence/9_1Mb-BAC.fa
	
Description:
	2011-9-1
		Two functions:
		1. A multi-sample genotype caller based entirely on coverage of reads.
			sam/bam file has to be indexed beforehand.
		2. A coverage-based GATK-generated VCF file filter.
		
		For multi-read-group input, seqCoverageFname will provide coverage for each individual sequence.
			Each read group has the individual_sequence.id embedded as 2nd value if it's split by "_".
			Read group formula: alnID_isqID_individual.code_sequencer_vs_alignment.ref_ind_seq_id
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO, re, csv
from pymodule import ProcessOptions, figureOutDelimiter
from pymodule.utils import sortCMPBySecondTupleValue
from pymodule.VCFFile import VCFFile

class GenotypeCallByCoverage(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('refFastaFname', 0, ): [None, 'e', 1, 'the fasta file containing reference sequences.'],\
						('numberOfReadGroups', 1, int): [None, 'n', 1, 'number of read groups/genomes in the inputFname', ],\
						('minMinorAlleleCoverage', 1, float): [1/4., 'M', 1, 'minimum read depth multiplier for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, float): [3/4., 'A', 1, 'maximum read depth multiplier for the minor allele of a heterozygous call', ],\
						('maxMajorAlleleCoverage', 1, float): [7/8., 'a', 1, 'maximum read depth multiplier for the major allele of het call'],\
						('maxNoOfReadsForGenotypingError', 1, float): [1, 'x', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, float): [2, 'm', 1, 'maximum read depth multiplier for one base to be considered'],\
						('minNoOfReads', 1, float): [1/4., 'O', 1, 'minimum read depth multiplier for one base to be considered'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, float): [3, 'N', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						('defaultCoverage', 1, float): [5, 'f', 1, 'default coverage when coverage is not available for a read group'],\
						('minDepth', 0, float): [0, '', 1, 'minimum depth for a VCF call to regarded as non-missing, used in 3rd run_type', ],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('outputDelimiter', 1, ): ['\t', 'u', 1, 'delimiter in the csv output file.'],\
						("run_type", 1, int): [1, 'y', 1, '1: discoverFromVCF (output of GATK), 2: discoverFromBAM, 3: discoverFromVCFWithoutFilter'],\
						("site_type", 1, int): [1, 's', 1, '1: all sites, 2: variants only'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.discoverFuncDict = {1: self.discoverFromVCF, 2: self.discoverFromBAM, 3:self.discoverFromVCFWithoutFilter}
	
	def get_isqID2coverage(self, seqCoverageFname, defaultCoverage=None):
		"""
		2011-9-2
		"""
		sys.stderr.write("Fetching sequence coverage info from %s ..."%(seqCoverageFname))
		
		reader = csv.reader(open(seqCoverageFname, 'r'), delimiter=figureOutDelimiter(seqCoverageFname))
		isqID2coverage = {}
		header = reader.next()
		for row in reader:
			isqID = int(row[0])
			coverage = float(row[1])
			isqID2coverage[isqID] = coverage
		sys.stderr.write("%s entries.\n"%len(isqID2coverage))
		return isqID2coverage
	
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
	def discoverFromBAM(cls, inputFname, outputFname, refFastaFname=None, monomorphicDiameter=100, \
						maxNoOfReads=300, minNoOfReads=2, minMinorAlleleCoverage=3, maxMinorAlleleCoverage=7,\
						maxNoOfReadsForGenotypingError=1, maxMajorAlleleCoverage=30, maxNoOfReadsForAllSamples=1000,\
						nt_set = set(['a','c','g','t','A','C','G','T']), VCFOutputType=None, \
						outputDelimiter='\t',\
						isqID2coverage=None, defaultCoverage=10, report=0, site_type=1):
		"""
		2011-9-2
			additional arguments, isqID2coverage and defaultCoverage, are not used right now.
		2011-8-31
			add argument site_type
		2011-7-20
			copied from discoverHetsFromBAM() of vervet.src.misc
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
		refNameSet = set()	#reverse map of tid2refName
		for pileupcolumn in samfile.pileup():
			#print
			#print 'coverage at base %s %s = %s'%(pileupcolumn.tid, pileupcolumn.pos , pileupcolumn.n)
			counter += 1
			refName = samfile.getrname(pileupcolumn.tid)
			current_locus = '%s_%s'%(refName, pileupcolumn.pos+1)
			if pileupcolumn.tid not in tid2refName:
				tid_str = str(pileupcolumn.tid)
				tid2refName[tid_str] = refName
				refNameSet.add(refName)
			
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
				if len(allele2count)>site_type-1:	#polymorphic across samples
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
		cls.outputCallMatrix(data_matrix, refFastaFname, outputFname=outputFname, \
					refNameSet=refNameSet, read_group2col_index=read_group2col_index, \
					locus_id2row_index=locus_id2row_index, outputDelimiter=outputDelimiter)
		
		unique_read_group_ls = read_group2col_index.keys()
		unique_read_group_ls.sort()
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_good_hets, unique_read_group_ls, title="No of good hets")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_good_tris, unique_read_group_ls, title="No of good SNPs with tri-or-more alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_trialleles, unique_read_group_ls, title="No of SNPs with tri alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_quad_alleles, unique_read_group_ls, title="No of SNPs with 4 alleles")
		cls.reportValueOfDictionaryByKeyLs(read_group2no_of_snps_with_penta_alleles, unique_read_group_ls, title="No of SNPs with 5-or-more alleles")
	
	@classmethod
	def outputCallMatrix(cls, data_matrix=None, refFastaFname=None, outputFname=None, refNameSet=None, read_group2col_index=None, \
						locus_id2row_index=None, outputDelimiter='\t'):
		"""
		2012.8.20
			add argument outputDelimiter
		2012.5.8
			if refNameSet is empty or None, stop reading the refFastaFname and outputting the ref base
		2011-7-26
			replace arguments refName2tid & tid2refName with refNameSet
		2011-7-20
			Meaning of tid from refName2tid and tid2refName depends on whether it's discoverFromBAM() or discoverFromVCF().
				In the former, tid is a consecutive number ID used in bam file to track reference sequences.
				In the latter, tid is chromosome number (could be 1,2,3 or X, Y).
			So in either case, the type of tid has been casted to str.
		"""
		#2011-7-18 read in the reference sequences in order to find out the ref base
		refName2Seq = {}
		if refNameSet:	#2012.5.8 not empty
			from Bio import SeqIO
			handle = open(refFastaFname, "rU")
			for record in SeqIO.parse(handle, "fasta"):
				contig_id = record.id.split()[0]
				if contig_id in refNameSet:
					refName2Seq[contig_id] = record.seq
				if len(refName2Seq)>=len(refNameSet):	#enough data, exit.
					break
			handle.close()
		
		# output the matrix in the end
		#read_group2col_index.pop('ref', None)	#remove the "ref" item if "ref" is in read_group2col_index. None is for failsafe when "ref" is not present.
		read_group_col_index_ls = read_group2col_index.items()
		read_group_col_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		header = ['locus_id', 'locus_id']+[row[0] for row in read_group_col_index_ls]
		writer = csv.writer(open(outputFname, 'w'), delimiter=outputDelimiter)
		#header = ['RG', 'chr', 'pos', 'MinorAlleleCoverage', 'MajorAlleleCoverage']
		#writer.writerow(header)
		writer.writerow(header)
		
		locus_id_and_row_index_ls = locus_id2row_index.items()
		locus_id_and_row_index_ls.sort(cmp=sortCMPBySecondTupleValue)
		for i in xrange(len(locus_id_and_row_index_ls)):
			locus_id, row_index = locus_id_and_row_index_ls[i]
			data_row = data_matrix[i]
			if refNameSet:
				refName, pos = locus_id.split('_')[:2]
				refSeq = refName2Seq[refName]
				pos = int(pos)
				refBase = refSeq[pos-1]
				data_row[0] = refBase	#2011-7-18
			# if data_row is shorter than read_group_col_index_ls, add "NA" to fill it up
			for j in xrange(len(data_row), len(read_group_col_index_ls)):
				data_row.append('NA')
			writer.writerow([locus_id, locus_id] + data_row)
		del writer
		
	
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
	def discoverFromVCF(cls, inputFname, outputFname, refFastaFname=None, VCFOutputType=2, \
					minMinorAlleleCoverage=1/4., maxMinorAlleleCoverage=3/4.,\
					maxNoOfReads=2., minNoOfReads=1/4., \
					maxNoOfReadsForGenotypingError=1, maxMajorAlleleCoverage=7/8., maxNoOfReadsForAllSamples=1000,\
					nt_set = set(['a','c','g','t','A','C','G','T']), isqID2coverage=None, defaultCoverage=10, \
					outputDelimiter='\t',\
					report=0, site_type=1):
		"""
		2011-9-2
			add argument isqID2coverage, defaultCoverage
		2011-8-26
			add argument site_type
			function is also more robust against missing fields etc.
		2011-7-20
			copied from discoverHetsFromVCF() of vervet.src.misc
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
		
		
		read_group2col_index = {'ref':0}	#ref is at column 0. "ref" must not be equal to any read_group.
		read_group2coverage = {}	#2011-9-2
		locus_id2row_index = {}
		data_matrix = []
		
		tid2refName = {}	#dictionary storing the target references which have SNP calls
		refNameSet = set()
		"""
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['sample', 'snp_id', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		moreHeader = ['GQ', 'GL', 'SB', 'QD', 'sndHighestGL', 'deltaGL']
		#['AF', 'AC','AN', 'Dels', 'HRun', 'HaplotypeScore','MQ0', 'QD']	#2011-3-4 useless
		if VCFOutputType==2:
			header += moreHeader
		chr_pure_number_pattern = re.compile(r'[a-z_A-Z]+(\d+)')
		chr_number_pattern = re.compile(r'chr(\d+)')
		"""
		
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
			refNameSet.add(chr)
			
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
			if VCFOutputType==2:	#2011-3-4 GATK
				format_column = row[col_name2index['FORMAT']]
				format_column_ls = format_column.split(':')
				format_column_name2index = getColName2IndexFromHeader(format_column_ls)
				data_row = ['NA']*(len(individual_name2col_index)+1)	# extra 1 for the ref
				allele2count = {}
				for individual_name, individual_col_index in individual_name2col_index.iteritems():
					read_group = individual_name
					if read_group not in read_group2col_index:
						read_group2col_index[read_group] = len(read_group2col_index)
						#2011-9-2
						if isqID2coverage:
							try:
								isqID = read_group.split('_')[1]
								isqID = int(isqID)
								coverage = isqID2coverage.get(isqID, defaultCoverage)
							except:
								sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
								import traceback
								traceback.print_exc()
								sys.stderr.write("Coverage for %s not available. use default=%s.\n"%(read_group, defaultCoverage))
								coverage = defaultCoverage
						else:
							coverage = defaultCoverage
						read_group2coverage[read_group] = coverage
					
					coverage = read_group2coverage[read_group]
					genotype_data = row[individual_col_index]
					genotype_data_ls = genotype_data.split(':')
					genotype_call_index = format_column_name2index.get('GT')
					genotype_quality_index = format_column_name2index.get('GQ')
					if genotype_quality_index is None:
						genotype_quality_index = format_column_name2index.get('DP')
					depth_index = format_column_name2index.get("DP")
					#GL_index = format_column_name2index.get('GL')
					if len(genotype_data_ls)<len(format_column_name2index):
						continue
					if depth_index is None or genotype_call_index is None:
						continue
					#genotype_quality = genotype_data_ls[genotype_quality_index]
					genotype_call = genotype_data_ls[genotype_call_index]
					depth = int(genotype_data_ls[depth_index])
					if depth>maxNoOfReads*coverage or depth<minNoOfReads*coverage:	#2011-3-29 skip. coverage too high or too low
						continue
					allele = 'NA'
					if genotype_call=='0/1' or genotype_call =='1/0':	#heterozygous, the latter notation is never used though.
						"""
						GL_list = genotype_data_ls[GL_index]
						GL_list = GL_list.split(',')
						GL_list = map(float, GL_list)
						GL = GL_list[1]
						sndHighestGL = max([GL_list[0], GL_list[2]])
						deltaGL = GL-sndHighestGL
						"""
						AD = genotype_data_ls[format_column_name2index.get('AD')]
						AD = map(int, AD.split(','))
						minorAlleleCoverage = min(AD)
						majorAlleleCoverage = max(AD)
						
						if minorAlleleCoverage<=maxMinorAlleleCoverage*coverage and minorAlleleCoverage>=minMinorAlleleCoverage*coverage \
								and majorAlleleCoverage<=maxMajorAlleleCoverage*coverage:
							DP4_ratio = float(AD[0])/AD[1]
							allele = '%s%s'%(refBase, altBase)
							"""
							data_row = [individual_name, 'chr%s:%s'%(chr, pos), chr, pos, quality, \
									depth, minorAlleleCoverage, DP4_ratio,\
									info_tag2value.get('MQ'), genotype_quality, GL,\
									info_tag2value.get('SB'), info_tag2value.get('QD'), sndHighestGL, deltaGL]
							#for i in range(3, len(moreHeader)):
							#	info_tag = moreHeader[i]
							#	data_row.append(info_tag2value.get(info_tag))
							writer.writerow(data_row)
							"""
					elif genotype_call=='./.':	#missing
						continue
					elif genotype_call =='1/1':
						allele = '%s%s'%(altBase, altBase)
					elif genotype_call =='0/0':
						allele = '%s%s'%(refBase, refBase)
					col_index = read_group2col_index.get(read_group)
					data_row[col_index] = allele
					if allele!='NA':
						if allele not in allele2count:
							allele2count[allele] = 0
						allele2count[allele] += 1
				
				if len(allele2count)>site_type-1:	#whether polymorphic across samples or all sites in vcf
					real_counter += 1
					locus_id2row_index[current_locus] = len(locus_id2row_index)
					data_matrix.append(data_row)
			"""
			elif VCFOutputType==1:	#samtools. 2011-7-20 outdated.
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
			"""
			counter += 1
			if counter%2000==0 and report:
				sys.stderr.write("%s\t%s\t%s"%("\x08"*80, counter, real_counter))
		del reader
		
		cls.outputCallMatrix(data_matrix, refFastaFname, outputFname=outputFname, refNameSet=refNameSet, \
					read_group2col_index=read_group2col_index, \
					locus_id2row_index=locus_id2row_index, outputDelimiter=outputDelimiter)
		
		sys.stderr.write("%s\t%s\t%s.\n"%("\x08"*80, counter, real_counter))
	
	
	def discoverFromVCFWithoutFilter(self, inputFname=None, outputFname=None, **keywords):
		"""
		2012.9.11
			read minDepth from self.minDepth
		2012.9.5
			add minDepth=0 to VCFFile
		#2012.8.20 locus_id2row_index from VCFFile is using (chr, pos) as key, not chr_pos
			need a conversion in between
		2012.5.8
		"""
		vcfFile = VCFFile(inputFname=inputFname, minDepth=self.minDepth)
		vcfFile.parseFile()
		
		read_group2col_index = vcfFile.sample_id2index
		locus_id2row_index = vcfFile.locus_id2row_index
		#2012.8.20 locus_id2row_index from VCFFile is using (chr, pos) as key, not chr_pos
		new_locus_id2row_index = {}
		for locus_id, row_index  in locus_id2row_index.iteritems():
			new_locus_id = '%s_%s'%(locus_id[0], locus_id[1])
			new_locus_id2row_index[new_locus_id] = row_index
		locus_id2row_index = new_locus_id2row_index
		
		data_matrix = vcfFile.genotype_call_matrix
		
		self.outputCallMatrix(data_matrix, refFastaFname=None, outputFname=outputFname, refNameSet=None, \
					read_group2col_index=read_group2col_index, \
					locus_id2row_index=locus_id2row_index, outputDelimiter=self.outputDelimiter)
	
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
		
		if self.seqCoverageFname:
			isqID2coverage = self.get_isqID2coverage(self.seqCoverageFname)
		else:
			isqID2coverage = {}
		from vervet.src.misc import VariantDiscovery
		maxNoOfReadsForAllSamples = self.numberOfReadGroups*self.maxNoOfReads*self.defaultCoverage*self.maxNoOfReadsMultiSampleMultiplier
		if self.run_type in self.discoverFuncDict:
			self.discoverFuncDict[self.run_type](self.inputFname, self.outputFname, \
					refFastaFname=self.refFastaFname,\
					maxNoOfReads=self.maxNoOfReads, minNoOfReads=self.minNoOfReads,\
					minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
					maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
					maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
					maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
					maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples, VCFOutputType=2, \
					outputDelimiter=self.outputDelimiter,\
					isqID2coverage=isqID2coverage, defaultCoverage=self.defaultCoverage, \
					report=self.report, site_type=self.site_type)
		else:
			sys.stderr.write("Unsupported run_type %s.Exit.\n"%(self.run_type))
			sys.exit(5)
		"""
		if self.run_type==2:
			self.discoverFromBAM(self.inputFname, self.outputFname, \
						refFastaFname=self.refFastaFname,\
						maxNoOfReads=self.maxNoOfReads, minNoOfReads=self.minNoOfReads,\
						minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
						maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
						maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
						maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
						maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples, \
						report=self.report, site_type=self.site_type)
		elif self.run_type==1:
			self.discoverFromVCF(self.inputFname, self.outputFname, \
					refFastaFname=self.refFastaFname,\
					maxNoOfReads=self.maxNoOfReads, minNoOfReads=self.minNoOfReads,\
					minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
					maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
					maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
					maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
					maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples, VCFOutputType=2, \
					isqID2coverage=isqID2coverage, defaultCoverage=self.defaultCoverage, \
					report=self.report, site_type=self.site_type)
		"""
		

if __name__ == '__main__':
	main_class = GenotypeCallByCoverage
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()