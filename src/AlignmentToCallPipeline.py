#!/usr/bin/env python
"""
Examples:
	# scaffold (120) is used as reference in alignment. get genome sequence id from 1 to 8.
	%s -o workflow.xml -a 120 -i 1-8 -y2 -s2
	
	# 1Mb-BAC (9) is used as reference.
	%s -o workflow.xml -a 9 -i 1-4 -y2 -s2
	
	# 8 genomes versus top 156 contigs
	%s -o workflow_8GenomeVsTop156Contigs.xml -u yh  -a 128 -i 1-8 -N 156 -y2 -s2
	
	# 2011-7-21 use GATK + coverage filter, top 2 contigs
	~/script/vervet/src/AlignmentToCallPipeline.py -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8 -N 2 -y1  -s2
	
	# 2011-7-21 use GATK + coverage filter on hoffman2 and site_handler, top 2 contigs
	~/script/vervet/src/AlignmentToCallPipeline.py -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8
		-N 2 -y1 -l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db  -s2
	
	#2011-8-31 work 10 VWP and one VRC ref monkeys, variants only
	~/script/vervet/src/AlignmentToCallPipeline.py -a 9 -I 495,498-507 -u yh  
		-l condorpool -y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC.xml -s2
	
	~/script/vervet/src/AlignmentToCallPipeline.py -a 120 -I 34,38 -u yh -l hoffman2
	-y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC_hoffman2.xml  -s2 -e /u/home/eeskin/polyacti
	-t /u/home/eeskin/polyacti/NetworkData/vervet/db -N 4
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData
from Pegasus.DAX3 import *


class AlignmentToCallPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 0, ): ['1-3,4,5', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['1-3,4,5', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						('minMinorAlleleCoverage', 1, int): [3, '', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, int): [7, '', 1, 'maximum read depth for the minor allele of a heterozygous call', ],\
						('maxNoOfReadsForGenotypingError', 1, int): [1, '', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, int): [20, '', 1, 'maximum read depth for one base to be considered'],\
						('maxMajorAlleleCoverage', 1, int): [10, '', 1, 'maximum read depth'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, int): [3, '', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/vervet/bin/GenomeAnalysisTK", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["condorpool", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\
						("topNumberOfContigs", 1, int): [156, 'N', 1, 'number of contigs'],\
						("needFastaIndexJob", 0, int): [0, '', 1, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, '', 1, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						("site_type", 1, int): [1, 's', 1, '1: all genome sites, 2: variants only'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		if self.ind_aln_id_ls:
			self.ind_aln_id_ls = getListOutOfStr(self.ind_aln_id_ls, data_type=int)
		
		self.samtools_path = self.samtools_path%self.home_path
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	def getTopNumberOfContigs(self, topNumberOfContigs, tax_id=60711, sequence_type_id=9):
		"""
		2011-7-12
			get all the top contigs
		"""
		sys.stderr.write("Getting %s top big contigs ..."%(self.topNumberOfContigs))
		refNameSet = set([])
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		from sqlalchemy import desc
		query = GenomeDB.AnnotAssembly.query.filter_by(tax_id=tax_id).filter_by(sequence_type_id=sequence_type_id).order_by(desc(GenomeDB.AnnotAssembly.stop))
		for row in query:
			refNameSet.add(row.chromosome)
			if len(refNameSet)>=topNumberOfContigs:
				break
		sys.stderr.write("Done.\n")
		return refNameSet
	
	def getAlignments(self, aln_ref_ind_seq_id=None, ind_seq_id_ls=[], ind_aln_id_ls=[], aln_method_id=2):
		"""
		2011-8-31
			add argument aln_method_id
		2011-7-12
		
		"""
		alignmentLs = []
		TableClass = VervetDB.IndividualAlignment
		if ind_aln_id_ls:
			sys.stderr.write("Getting %s alignments ..."%(len(ind_aln_id_ls)) )
			query = TableClass.query.filter(TableClass.id.in_(ind_aln_id_ls))
		elif ind_seq_id_ls:
			sys.stderr.write("Getting all alignments for %s sequences with %s as reference ..."%(len(ind_seq_id_ls), aln_ref_ind_seq_id))
			query = TableClass.query.filter_by(ref_ind_seq_id=aln_ref_ind_seq_id).filter(TableClass.ind_seq_id.in_(ind_seq_id_ls))\
				.filter_by(aln_method_id=aln_method_id)
		else:
			sys.stderr.write("Both ind_seq_id_ls and ind_aln_id_ls are empty. no alignment to be fetched.\n")
			sys.exit(3)
		for row in query:
			if row.path:	#it's not None
				alignmentLs.append(row)
		sys.stderr.write("%s alignments Done.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def addRefFastaFileSplitJobs(self, workflow, refFastaF, selectAndSplitFasta, refNameLs, mkdir=None, samtools=None,
								java=None, createSequenceDictionaryJar=None,\
								site_handler=None, namespace='workflow', version='1.0',):
		"""
		2011-7-25
			split the whole fasta file into files, each containing one fasta record (from refNameLs)
			return the data
		"""
		
		# Add a mkdir job for the call directory.
		#letting numerou genotype call jobs detect&create this directory runs into race condition.
		mkDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
		fastaOutputDir = "fasta"
		mkDirJob.addArguments(fastaOutputDir)
		#mkDirJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
		workflow.addJob(mkDirJob)
		
		selectAndSplitFastaJob = Job(namespace=namespace, name=selectAndSplitFasta.name, version=version)
		selectAndSplitFastaJob.addArguments('-i', refFastaF, "-o", fastaOutputDir)
		#selectAndSplitFastaJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		workflow.addJob(selectAndSplitFastaJob)
		workflow.addDependency(parent=mkDirJob, child=selectAndSplitFastaJob)
		
		refName2jobDataLs = {}
		for refName in refNameLs:
			if refName not in refName2jobDataLs:
				refName2jobDataLs[refName] = []
			selectAndSplitFastaJob.addArguments(refName)
			
			fastaFname = os.path.join(fastaOutputDir, '%s.fasta'%(refName))
			fastaFile = File(fastaFname)
			
			# add the index job
			fai_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			fai_index_job.addArguments("faidx", fastaFname)
			fastaFAIIndexFname = os.path.join(fastaOutputDir, '%s.fasta.fai'%(refName))
			fastaFAIIndexFile = File(fastaFAIIndexFname)
			fai_index_job.uses(fastaFAIIndexFile, transfer=True, register=False, link=Link.OUTPUT)	#this file is input & output
			#fai_index_job.uses(output, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
			#fai_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(fai_index_job)
			workflow.addDependency(parent=selectAndSplitFastaJob, child=fai_index_job)
			
			# the add-read-group job
			createSeqDictJob = Job(namespace=namespace, name=java.name, version=version)
			dictFname = os.path.join(fastaOutputDir, '%s.dict'%(refName))
			#outputRGFname = '%s_%s.RG.bam'%(inputFileBaseNamePrefix, refName)
			dictFile = File(dictFname)
			createSeqDictJob.addArguments('-jar', createSequenceDictionaryJar, \
								"R=", fastaFile, 'O=', dictFile)
			createSeqDictJob.uses(dictFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			createSeqDictJob.uses(fastaFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			workflow.addJob(createSeqDictJob)
			workflow.addDependency(parent=selectAndSplitFastaJob, child=createSeqDictJob)
			
			
			refName2jobDataLs[refName] = [fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile]
		return PassingData(refName2jobDataLs=refName2jobDataLs, workflow=workflow)
	
	def addSelectAndSplitBamJobs(self, db_vervet, workflow, alignmentLs, site_handler, topNumberOfContigs, refNameLs, samtools=None, \
							java=None, addOrReplaceReadGroupsAndCleanSQHeaderJar=None, mkdir=None, namespace="workflow",\
							version="1.0", mkCallDirJob=None,\
							addRGAndCleanSQHeaderAlignment=None, dataDir=None):
		"""
		2011-7-14
			1. select reference out of whole-alignment
			2. index
			3. add read groups
			4. index
		"""
		sys.stderr.write("Adding Bam select and split jobs for %s alignments ..."%(len(alignmentLs)))
		refName2jobDataLs = {}
		for alignment in alignmentLs:
			# Add input file to the DAX-level replica catalog
			inputFname = os.path.join(dataDir, alignment.path)
			input = File(inputFname)
			input.addPFN(PFN("file://" + inputFname, site_handler))
			workflow.addFile(input)
			
			inputFileBaseNamePrefix = os.path.splitext(os.path.basename(alignment.path))[0]
			outputDir = inputFileBaseNamePrefix
			
			
			# add RG to this bam
			sequencer = alignment.ind_sequence.sequencer
			read_group = '%s_%s_%s_vs_%s'%(alignment.ind_seq_id, alignment.ind_sequence.individual.code, sequencer, \
										alignment.ref_ind_seq_id)
			if sequencer=='454':
				platform_id = 'LS454'
			elif sequencer=='GA':
				platform_id = 'ILLUMINA'
			else:
				platform_id = 'ILLUMINA'
			# Add a mkdir job
			mkdirJob = Job(namespace=namespace, name=mkdir.name, version=version)
			
			mkdirJob.addArguments(outputDir)
			#mkdir.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
			workflow.addJob(mkdirJob)
			workflow.addDependency(parent=mkCallDirJob, child=mkdirJob)
			
			for refName in refNameLs:
				if refName not in refName2jobDataLs:
					refName2jobDataLs[refName] = []
				
				selectRefJob = Job(namespace=namespace, name=samtools.name, version=version)
				
				outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix, refName))
				#outputFname = '%s_%s.bam'%(inputFileBaseNamePrefix, refName)
				output = File(outputFname)
				selectRefJob.addArguments('view', '-h', input, refName, "-o", output, "-b", "-u")	# -b -u forces uncompressed bam output
				#selectRefJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
				workflow.addJob(selectRefJob)
				workflow.addDependency(parent=mkdirJob, child=selectRefJob)
				
				# add the index job
				index_1_job = Job(namespace=namespace, name=samtools.name, version=version)
				index_1_job.addArguments("index", output)
				#index_1_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				#index_1_job.uses(output, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputFname)
				#index_1_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(index_1_job)
				workflow.addDependency(parent=selectRefJob, child=index_1_job)
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGAndCleanSQHeaderAlignment.name, version=version)
				addRGJob = Job(namespace=namespace, name=java.name, version=version)
				#2011-7-27 somehow AddOrReplaceReadGroupsAndCleanSQHeader.jar couldn't output a un-corrupted bam file. so sam first.
				outputRGSAMFname = os.path.join(outputDir, '%s_%s.RG.sam'%(inputFileBaseNamePrefix, refName))
				outputRGSAM = File(outputRGSAMFname)
				"""
				addRGJob.addArguments("-i", outputFname, '-o', outputRG, '-p', platform_id, '-a', read_group, \
									'-e', refName)
				"""
				addRGJob.addArguments('-jar', addOrReplaceReadGroupsAndCleanSQHeaderJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SQName=%s'%(refName))	#'SORT_ORDER=coordinate', (adding this is useless)
				addRGJob.uses(output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(outputRGSAM, transfer=False, register=True, link=Link.OUTPUT)	#time to discard them
				workflow.addJob(addRGJob)
				workflow.addDependency(parent=index_1_job, child=addRGJob)
				#output.addPFN(PFN("file://" + outputFname, site_handler))
				#selectAndSplitJob.uses(output, transfer=True, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
				
				samToBamJob = Job(namespace=namespace, name=samtools.name, version=version)
				outputRGFname = os.path.join(outputDir, '%s_%s.RG.bam'%(inputFileBaseNamePrefix, refName))
				outputRG = File(outputRGFname)
				samToBamJob.addArguments('view', '-F4', '-Sbh', "-o", outputRG, "-u", outputRGSAM)	# -b -u forces uncompressed bam output
				samToBamJob.uses(outputRGSAM, transfer=False, register=True, link=Link.INPUT)
				#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
				workflow.addJob(samToBamJob)
				workflow.addDependency(parent=addRGJob, child=samToBamJob)
				
				# add the index job
				samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
				samtools_index_job.addArguments("index", outputRG)
				#samtools_index_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				samtools_index_job.uses(outputRG, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputRGFname)
				samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(samtools_index_job)
				workflow.addDependency(parent=samToBamJob, child=samtools_index_job)
				refName2jobDataLs[refName].append((outputRG, samtools_index_job, bai_output))
		sys.stderr.write(".Done\n")
		return PassingData(refName2jobDataLs=refName2jobDataLs, workflow=workflow)
	
	def addMergeAlignmentAndGenotypeCallJobs(self, workflow, refName2jobDataLs, refNameLs, samtools, \
				java, createSequenceDictionaryJar=None, genotypeCallByCoverage=None, refFastaF=None, \
				namespace='workflow', version="1.0", callOutputDir = "call", genotypeCallerType=1,\
				mergeSamFilesJar=None, genomeAnalysisTKJar=None, calcula=None, refName2splitFastaJobDataLs=None, \
				needFastaIndexJob=False, \
				needFastaDictJob=False, site_type=1):
		"""
		2011-8-26
			add argument needFastaDictJob, needFastaIndexJob, site_type
		2011-7-26
			add refName2splitFastaJobDataLs, for GATK
		2011-7-14
			1. merge alignments based on same reference from different genomes into one
			2. run genotypeCallByCoverage on each merged alignment
		"""
		sys.stderr.write("Adding alignment merge and genotype call jobs for %s references ..."%(len(refName2jobDataLs)))
		if needFastaDictJob:	# the .dict file is required for GATK
			refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
			refFastaDictF = File(refFastaDictFname)
			#not os.path.isfile(refFastaDictFname) or 
			fastaDictJob = Job(namespace=namespace, name=java.name, version=version)
			fastaDictJob.addArguments('-jar', createSequenceDictionaryJar, \
					'REFERENCE=', refFastaF, 'OUTPUT=', refFastaDictF)
			fastaDictJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaDictJob.uses(refFastaDictF, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaDictJob)
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
			refFastaIndexF = File(refFastaIndexFname)
			#not os.path.isfile(refFastaIndexFname)
			fastaIndexJob = Job(namespace=namespace, name=samtools.name, version=version)
			fastaIndexJob.addArguments("faidx", refFastaF)
			fastaIndexJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaIndexJob.uses(refFastaIndexFname, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaIndexJob)
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		#workflow.addDependency(parent=fastaDictJob, child=fastaIndexJob)	#no dependency between these two jobs
		
		# add merge jobs for every reference
		refName2mergedBamCallJob = {}
		for refName, jobDataLs in refName2jobDataLs.iteritems():
			# add the index job
			picard_job = Job(namespace=namespace, name=java.name, version=version)
			outputFname = '%s.bam'%(refName)
			picard_output = File(outputFname)
			picard_job.addArguments('-jar', mergeSamFilesJar, \
				'USE_THREADING=true', 'SORT_ORDER=coordinate', 'ASSUME_SORTED=false', 'OUTPUT=', picard_output)
			#picard_job.uses(picard_output, transfer=False, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
			workflow.addJob(picard_job)
			
			for jobData in jobDataLs:
				bamFile, samtools_index_job, bamFileBai = jobData[:3]
				picard_job.addArguments('INPUT=', bamFile)
				picard_job.uses(bamFileBai, transfer=True, register=True, link=Link.INPUT)	#register them here to be deleted 
				picard_job.uses(bamFile, transfer=True, register=True, link=Link.INPUT)
				#this picard merge job depends on a bunch of prior samtools index jobs
				workflow.addDependency(parent=samtools_index_job, child=picard_job)
			
			# add the index job on the merged bam file
			samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			samtools_index_job.addArguments("index", picard_output)
			samtools_index_job.uses(picard_output, transfer=True, register=True, link=Link.OUTPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
			bai_output = File('%s.bai'%outputFname)
			samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(samtools_index_job)
			workflow.addDependency(parent=picard_job, child=samtools_index_job)
			
			
			if genotypeCallerType==2:
				coverageFilterParentJob = samtools_index_job
				coverageFilterParentOutput = picard_output
				coverageFilterParentOutput_bai = bai_output
			else:
				splitFastaJobDataLs = refName2splitFastaJobDataLs.get(refName)
				fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile = splitFastaJobDataLs[:5]
				
				gatk_job = Job(namespace=namespace, name=java.name, version=version)
				gatkOutputFname = os.path.join(callOutputDir, '%s.vcf'%(refName))
				gatk_output = File(gatkOutputFname)
				gatkIDXOutputFname = os.path.join(callOutputDir, '%s.vcf.idx'%(refName))
				gatkIDXOutput = File(gatkIDXOutputFname)
				
				gatk_job.addArguments('-jar', genomeAnalysisTKJar, \
					"-I", picard_output, "-R", fastaFile, "-T", "UnifiedGenotyper","--out", gatk_output,\
					'-U', '-S SILENT')
				if site_type==1:
					gatk_job.addArguments('--output_mode EMIT_ALL_SITES')	#2011-8-24 new GATK no longers ues "-all_bases"
				gatk_job.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
				gatk_job.uses(picard_output, transfer=False, register=True, link=Link.INPUT)
				
				gatk_job.uses(fastaFAIIndexFile, transfer=True, register=True, link=Link.INPUT)
				gatk_job.uses(fastaFile, transfer=True, register=True, link=Link.INPUT)
				gatk_job.uses(dictFile, transfer=False, register=True, link=Link.INPUT)
				
				gatk_job.uses(gatk_output, transfer=True, register=True, link=Link.OUTPUT)
				gatk_job.uses(gatkIDXOutput, transfer=True, register=True, link=Link.OUTPUT)
				
				
				gatk_job_max_memory = 1500	#in MB
				gatk_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%gatk_job_max_memory))
				gatk_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%gatk_job_max_memory))
				
				workflow.addJob(gatk_job)
				workflow.addDependency(parent=samtools_index_job, child=gatk_job)
				workflow.addDependency(parent=fai_index_job, child=gatk_job)
				workflow.addDependency(parent=createSeqDictJob, child=gatk_job)
				if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
					workflow.addDependency(parent=fastaIndexJob, child=gatk_job)
				if fastaDictJob:
					workflow.addDependency(parent=fastaDictJob, child=gatk_job)
				
				coverageFilterParentJob = gatk_job
				coverageFilterParentOutput = gatk_output
				coverageFilterParentOutput_bai = None
				
			
			#add the cover filter (or filter+call) job after index is done
			genotypeCallByCoverage_job = Job(namespace=namespace, name=genotypeCallByCoverage.name, version=version)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(refName))	#genotypeCallByCoverage_job would create directory "call".
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job.addArguments("-i", coverageFilterParentOutput, "-n", str(len(jobDataLs)), \
					"-o", genotypeCallOutput, '-e', refFastaF, '-y', str(genotypeCallerType), \
					'-a', repr(self.maxMajorAlleleCoverage), '-M', repr(self.minMinorAlleleCoverage), '-s', repr(self.site_type))
			if coverageFilterParentOutput_bai:
				genotypeCallByCoverage_job.uses(coverageFilterParentOutput_bai, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
			genotypeCallByCoverage_job.uses(coverageFilterParentOutput, transfer=True, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(genotypeCallOutput, transfer=True, register=True, link=Link.OUTPUT)
			job_max_memory = 5000	#in MB
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			workflow.addJob(genotypeCallByCoverage_job)
			workflow.addDependency(parent=coverageFilterParentJob, child=genotypeCallByCoverage_job)
			
			
			#add the pairwise distance matrix job after filter is done
			calcula_job = Job(namespace=namespace, name=calcula.name, version=version)
			
			convertHetero2NA = True
			max_NA_rate = 0.4
			min_MAF = 0
			calculaOutputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(genotypeCallOutputFname)[0], \
								convertHetero2NA, min_MAF, max_NA_rate)
			calculaOutput = File(calculaOutputFname)
			calcula_job.addArguments("-i", genotypeCallOutput, "-n", str(min_MAF), \
								"-o", calculaOutput, '-m', repr(max_NA_rate), '-c', str(1))
			calcula_job.uses(genotypeCallOutput, transfer=True, register=False, link=Link.INPUT)
			calcula_job.uses(calculaOutput, transfer=True, register=False, link=Link.OUTPUT)
			
			workflow.addJob(calcula_job)
			workflow.addDependency(parent=genotypeCallByCoverage_job, child=calcula_job)
			
			refName2mergedBamCallJob[refName] = [genotypeCallOutput, genotypeCallByCoverage_job]
		sys.stderr.write(".Done\n")
		return PassingData(refName2mergedBamCallJob=refName2mergedBamCallJob)
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		refNameSet = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refNameSet = set(['Contig149'])	#temporary when testing Contig149
		#refNameSet = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = list(refNameSet)
		refNameLsStr = ','.join(refNameLs)
		
		alignmentLs = self.getAlignments(self.aln_ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2)
		
		refSequence = VervetDB.IndividualSequence.get(self.aln_ref_ind_seq_id)
		
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		if self.needFastaDictJob and self.needFastaIndexJob:
			refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path
			# Add it into replica only when needed.
			refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
			workflow.addFile(refFastaF)
			# If it's not needed, assume the index is done and all relevant files are in absolute path.
			# and no replica transfer
		else:
			refFastaF = File(refFastaFname)
		# Create a abstract dag
		workflow = ADAG("AlignmentToCallPipeline")
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		#add the MergeSamFiles.jar file into workflow
		abs_path = os.path.join(self.picard_path, 'MergeSamFiles.jar')
		mergeSamFilesJar = File(abs_path)	#using abs_path avoids add this jar to every job as Link.INPUT
		mergeSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(mergeSamFilesJar)
		
		abs_path = os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar')
		genomeAnalysisTKJar = File(abs_path)
		genomeAnalysisTKJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(genomeAnalysisTKJar)
		
		abs_path = os.path.join(self.picard_path, 'CreateSequenceDictionary.jar')
		createSequenceDictionaryJar = File(abs_path)
		createSequenceDictionaryJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(createSequenceDictionaryJar)
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroupsAndCleanSQHeader.jar')
		addOrReplaceReadGroupsAndCleanSQHeaderJar = File(abs_path)
		addOrReplaceReadGroupsAndCleanSQHeaderJar.addPFN(PFN("file://" + abs_path, \
												site_handler))
		workflow.addFile(addOrReplaceReadGroupsAndCleanSQHeaderJar)
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		
		mkdir = Executable(namespace=namespace, name="mkdir", version=version, os=operatingSystem, arch=architecture, installed=True)
		mkdir.addPFN(PFN("file://" + '/bin/mkdir', site_handler))
		workflow.addExecutable(mkdir)
		
		selectAndSplit = Executable(namespace=namespace, name="SelectAndSplitAlignment", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		workflow.addExecutable(selectAndSplit)
		
		# 2011-9-1 not used anymore. too slow. replace with addOrReplaceReadGroupsAndCleanSQHeaderJar.
		addRGAndCleanSQHeaderAlignment = Executable(namespace=namespace, name="AddRGAndCleanSQHeaderAlignment", version=version, \
												os=operatingSystem, arch=architecture, installed=True)
		addRGAndCleanSQHeaderAlignment.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AddRGAndCleanSQHeaderAlignment.py"), site_handler))
		workflow.addExecutable(addRGAndCleanSQHeaderAlignment)
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "SelectAndSplitFastaRecords.py"), site_handler))
		workflow.addExecutable(selectAndSplitFasta)
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(java)
		
		genotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		workflow.addExecutable(genotypeCallByCoverage)
		
		mergeGenotypeMatrix = Executable(namespace=namespace, name="MergeGenotypeMatrix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		mergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "MergeGenotypeMatrix.py"), site_handler))
		workflow.addExecutable(mergeGenotypeMatrix)
		
		calcula = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), site_handler))
		workflow.addExecutable(calcula)
		
		
		# Add a mkdir job for the call directory.
		#letting numerou genotype call jobs detect&create this directory runs into race condition.
		mkCallDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
		callOutputDir = "call"
		mkCallDirJob.addArguments(callOutputDir)
		#mkdir.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
		workflow.addJob(mkCallDirJob)
		
		addSelectAndSplitBamReturnData = self.addSelectAndSplitBamJobs(db_vervet, workflow, alignmentLs, site_handler, \
							self.topNumberOfContigs, refNameLs, samtools=samtools, \
							java=java, addOrReplaceReadGroupsAndCleanSQHeaderJar=addOrReplaceReadGroupsAndCleanSQHeaderJar, \
							mkdir=mkdir, namespace=namespace, \
							version=version, mkCallDirJob=mkCallDirJob,\
							addRGAndCleanSQHeaderAlignment=addRGAndCleanSQHeaderAlignment, dataDir=self.dataDir)
		refName2jobDataLs = addSelectAndSplitBamReturnData.refName2jobDataLs
		
		returnData3 = self.addRefFastaFileSplitJobs(workflow, refFastaF, selectAndSplitFasta, refNameLs, mkdir=mkdir, samtools=samtools,
								java=java, createSequenceDictionaryJar=createSequenceDictionaryJar,\
								site_handler=site_handler, namespace=namespace, version=version)
		refName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		returnData2 = self.addMergeAlignmentAndGenotypeCallJobs(workflow, refName2jobDataLs, refNameLs, samtools, \
							java, createSequenceDictionaryJar=createSequenceDictionaryJar, genotypeCallByCoverage=genotypeCallByCoverage, \
							refFastaF=refFastaF, \
							namespace=namespace, version=version, callOutputDir = callOutputDir, \
							genotypeCallerType=self.genotypeCallerType, \
							mergeSamFilesJar=mergeSamFilesJar, genomeAnalysisTKJar=genomeAnalysisTKJar, calcula=calcula, \
							refName2splitFastaJobDataLs=refName2splitFastaJobDataLs, needFastaIndexJob=self.needFastaIndexJob, \
							needFastaDictJob=self.needFastaDictJob, site_type=self.site_type)
		
		refName2mergedBamCallJob =returnData2.refName2mergedBamCallJob
		
		#merge all genotype call files
		mergeGenotypeMatrix_job = Job(namespace=namespace, name=mergeGenotypeMatrix.name, version=version)
		finalCallOutputFname = '%s_genomes_vs_top%sReferences_call.tsv'%(len(alignmentLs), len(refNameLs))
		finalCallOutput = File(finalCallOutputFname)
		mergeGenotypeMatrix_job.addArguments("-o", finalCallOutput)
		mergeGenotypeMatrix_job.uses(finalCallOutput, transfer=True, link=Link.OUTPUT, register=True)
		workflow.addJob(mergeGenotypeMatrix_job)
		inputFnameLs = []
		for refName, callJobData in refName2mergedBamCallJob.iteritems():
			callJobOutput, callJob = callJobData[:2]
			#callJobOutput = callJob.used_files[1]
			
			mergeGenotypeMatrix_job.uses(callJobOutput, transfer=False, register=True, link=Link.INPUT)
			inputFnameLs.append(callJobOutput.name)
			workflow.addDependency(parent=callJob, child=mergeGenotypeMatrix_job)
		
			mergeGenotypeMatrix_job.addArguments(callJobOutput)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = AlignmentToCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
