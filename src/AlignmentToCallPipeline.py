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
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8 -N 2 -y1  -s2
	
	# 2011-7-21 use GATK + coverage filter on hoffman2 and site_handler, top 5 contigs
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8
		-N 5 -y1 -l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db  -s2
	
	#2011-8-31 work 10 VWP and one VRC ref monkeys, variants only
	%s -a 9 -I 495,498-507 -u yh  
		-l condorpool -y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC.xml -s2 -q /tmp/all_isq_coverage.tsv
	
	%s -a 120 -I 34,38 -u yh -l hoffman2
	-y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC_hoffman2.xml  -s2 -e /u/home/eeskin/polyacti
	-t /u/home/eeskin/polyacti/NetworkData/vervet/db -N 4 -q /tmp/all_isq_coverage.tsv
	
	#2011-9-14 top 25 contigs, variants only, run on uschpc cluster
	%s -I 559-656 -j uschpc -l uschpc -u yh -a 524 -s 2 -e /home/cmb-03/mn/yuhuang -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/
		-z 10.8.0.10 -o ./AlignmentToCallPipeline_559_656_vs_524_top_25Contigs_uschpc.xml
		-D ~/mnt/hpc-cmb_home/NetworkData/vervet/db/ -N25
	
	# 2011-11-4 run GATK/samtools on single-sample at a time, for 4 high-coverage VRC monkeys, top 804 contigs
	%s -a 524 -i 15-18 -u yh -l condorpool -j condorpool -s 2 -N 804
		-o AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_single_sample_condor.xml -z uclaOffice -n 2
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0],\
				sys.argv[0], sys.argv[0], sys.argv[0])

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *


class AlignmentToCallPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\
						("topNumberOfContigs", 1, int): [156, 'N', 1, 'number of contigs'],\
						("needFastaIndexJob", 0, int): [0, '', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, '', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						("site_type", 1, int): [2, 's', 1, '1: all genome sites, 2: variants only'],\
						("run_type", 1, int): [1, 'n', 1, '1: cross-sample calling, 2: single-sample one by one'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
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
		2011-9-13
			return refName2size instead of a set of ref names
		2011-7-12
			get all the top contigs
		"""
		sys.stderr.write("Getting %s top big contigs ..."%(self.topNumberOfContigs))
		refName2size = {}
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		from sqlalchemy import desc
		query = GenomeDB.AnnotAssembly.query.filter_by(tax_id=tax_id).filter_by(sequence_type_id=sequence_type_id).order_by(desc(GenomeDB.AnnotAssembly.stop))
		for row in query:
			refName2size[row.chromosome] = row.stop
			if len(refName2size)>=topNumberOfContigs:
				break
		sys.stderr.write("Done.\n")
		return refName2size
	
	def getAlignments(self, aln_ref_ind_seq_id=None, ind_seq_id_ls=[], ind_aln_id_ls=[], aln_method_id=2, dataDir=None):
		"""
		2011-9-16
			order each alignment by id. It is important because this is the order that gatk&samtools take input bams.
			#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
			# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		2011-9-13
			add argument dataDir, to filter out alignments that don't exist on file storage
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
		#order by TableClass.id is important because this is the order that gatk&samtools take input bams.
		#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
		# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		query = query.order_by(TableClass.id)
		for row in query:
			if row.path:	#it's not None
				abs_path = os.path.join(dataDir, row.path)
				if os.path.isfile(abs_path):
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
		workflow.depends(parent=mkDirJob, child=selectAndSplitFastaJob)
		
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
			workflow.depends(parent=selectAndSplitFastaJob, child=fai_index_job)
			
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
			workflow.depends(parent=selectAndSplitFastaJob, child=createSeqDictJob)
			
			
			refName2jobDataLs[refName] = [fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile]
		return PassingData(refName2jobDataLs=refName2jobDataLs, workflow=workflow)
	
	def addSelectAndSplitBamJobs(self, db_vervet, workflow, alignmentLs, site_handler, topNumberOfContigs, refNameLs, samtools=None, \
							java=None, addOrReplaceReadGroupsAndCleanSQHeaderJar=None, BuildBamIndexFilesJar=None, mkdir=None, namespace="workflow",\
							version="1.0", mkCallDirJob=None,\
							addRGExecutable=None, dataDir=None):
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
			read_group = '%s_%s_%s_%s_vs_%s'%(alignment.id, alignment.ind_seq_id, alignment.ind_sequence.individual.code, \
									sequencer, alignment.ref_ind_seq_id)
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
			workflow.depends(parent=mkCallDirJob, child=mkdirJob)
			
			for refName in refNameLs:
				if refName not in refName2jobDataLs:
					refName2jobDataLs[refName] = []
				#select reads that are aligned to one reference name
				selectRefJob = Job(namespace=namespace, name=samtools.name, version=version)
				
				outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix, refName))
				#outputFname = '%s_%s.bam'%(inputFileBaseNamePrefix, refName)
				output = File(outputFname)
				selectRefJob.addArguments('view', '-h', input, refName, "-o", output, "-b", "-u")	# -b -u forces uncompressed bam output
				#selectRefJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
				workflow.addJob(selectRefJob)
				workflow.depends(parent=mkdirJob, child=selectRefJob)
				
				# add the index job
				index_1_job = Job(namespace=namespace, name=samtools.name, version=version)
				index_1_job.addArguments("index", output)
				#index_1_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				index_1_job.uses(output, transfer=False, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputFname)
				index_1_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(index_1_job)
				workflow.depends(parent=selectRefJob, child=index_1_job)
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=java.name, version=version)
				#2011-7-27 somehow AddOrReplaceReadGroupsAndCleanSQHeader.jar couldn't output a un-corrupted bam file. so sam first.
				outputRGSAMFname = os.path.join(outputDir, '%s_%s.RG.sam'%(inputFileBaseNamePrefix, refName))
				outputRGSAM = File(outputRGSAMFname)
				"""
				tmpRGFname = os.path.join(outputDir, '%s_%s.RG.txt'%(inputFileBaseNamePrefix, refName))
				addRGJob.addArguments(read_group, platform_id, output, tmpRGFname, outputRGSAM)
				"""
				addRGJob.addArguments('-jar', addOrReplaceReadGroupsAndCleanSQHeaderJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SQName=%s'%(refName))	#'SORT_ORDER=coordinate', (adding this is useless)
				"""
				
				#addRGJob.addArguments('-jar', addOrReplaceReadGroupsJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM)	#'SORT_ORDER=coordinate', (adding this is useless)
				"""
				addRGJob.uses(output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(outputRGSAM, transfer=False, register=True, link=Link.OUTPUT)	#time to discard them
				workflow.addJob(addRGJob)
				workflow.depends(parent=index_1_job, child=addRGJob)
				#output.addPFN(PFN("file://" + outputFname, site_handler))
				#selectAndSplitJob.uses(output, transfer=True, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
				
				samToBamJob = Job(namespace=namespace, name=samtools.name, version=version)
				outputRGFname = os.path.join(outputDir, '%s_%s.RG.bam'%(inputFileBaseNamePrefix, refName))
				outputRG = File(outputRGFname)
				samToBamJob.addArguments('view', '-F4', '-Sbh', "-o", outputRG, "-u", outputRGSAM)	# -b -u forces uncompressed bam output
				samToBamJob.uses(outputRGSAM, transfer=False, register=True, link=Link.INPUT)
				#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
				workflow.addJob(samToBamJob)
				workflow.depends(parent=addRGJob, child=samToBamJob)
				
				# add the index job
				index_sam_job = Job(namespace=namespace, name=java.name, version=version)
				bai_output = File('%s.bai'%outputRGFname)
				index_sam_job.addArguments("-Xms128m", "-Xmx2500m", "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
								"INPUT=", outputRG, \
								"OUTPUT=", bai_output)
				index_sam_job.uses(outputRG, transfer=False, register=False, link=Link.INPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
				index_sam_job.uses(bai_output, transfer=False, register=False, link=Link.OUTPUT)
				
				"""
				samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
				samtools_index_job.addArguments("index", outputRG)
				#samtools_index_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				samtools_index_job.uses(outputRG, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputRGFname)
				samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				"""
				workflow.addJob(index_sam_job)
				workflow.depends(parent=samToBamJob, child=index_sam_job)
				refName2jobDataLs[refName].append((outputRG, index_sam_job, bai_output))
				
				"""
				#2011-9-1 temporary addition to make sure it's sorted for Vasily
				# input bam doesn't need to be indexed.
				outputRGSortSAMFnamePrefix = os.path.join(outputDir, '%s_%s.RG.sorted'%(inputFileBaseNamePrefix, refName))
				outputRGSortSAMFname = File('%s.bam'%(outputRGSortSAMFnamePrefix))
				sam_sort_job = Job(namespace=namespace, name=samtools.name, version=version)
				sam_sort_job.addArguments('sort', '-m', '2000000000', outputRG, outputRGSortSAMFnamePrefix)	#maxMemory is down to 2G
				sam_sort_job.uses(outputRG, transfer=False, register=False, link=Link.INPUT)
				sam_sort_job.uses(outputRGSortSAMFname, transfer=True, register=False, link=Link.OUTPUT)
				job_max_memory = 2000	#in MB
				sam_sort_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				sam_sort_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				workflow.addJob(sam_sort_job)
				workflow.depends(parent=samToBamJob, child=sam_sort_job)
				"""
				
		sys.stderr.write(".Done\n")
		return PassingData(refName2jobDataLs=refName2jobDataLs, workflow=workflow)
	
	def addMergeAlignmentAndGenotypeCallJobs(self, workflow, refName2jobDataLs, refNameLs, samtools, \
				java, createSequenceDictionaryJar=None, genotypeCallByCoverage=None, refFastaF=None, \
				namespace='workflow', version="1.0", callOutputDir = "call", genotypeCallerType=1,\
				mergeSamFilesJar=None, genomeAnalysisTKJar=None, calcula=None, refName2splitFastaJobDataLs=None, \
				seqCoverageF=None, needFastaIndexJob=False, \
				needFastaDictJob=False, site_type=1):
		"""
		2011-9-2
			add argument seqCoverageF
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
		#workflow.depends(parent=fastaDictJob, child=fastaIndexJob)	#no dependency between these two jobs
		
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
				workflow.depends(parent=samtools_index_job, child=picard_job)
			
			# add the index job on the merged bam file
			samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			samtools_index_job.addArguments("index", picard_output)
			samtools_index_job.uses(picard_output, transfer=True, register=True, link=Link.OUTPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
			bai_output = File('%s.bai'%outputFname)
			samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(samtools_index_job)
			workflow.depends(parent=picard_job, child=samtools_index_job)
			
			
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
					'-U', '-S SILENT', "-nt 4")
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
				workflow.depends(parent=samtools_index_job, child=gatk_job)
				workflow.depends(parent=fai_index_job, child=gatk_job)
				workflow.depends(parent=createSeqDictJob, child=gatk_job)
				if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
					workflow.depends(parent=fastaIndexJob, child=gatk_job)
				if fastaDictJob:
					workflow.depends(parent=fastaDictJob, child=gatk_job)
				
				coverageFilterParentJob = gatk_job
				coverageFilterParentOutput = gatk_output
				coverageFilterParentOutput_bai = None
				
			
			#add the cover filter (or filter+call) job after index is done
			genotypeCallByCoverage_job = Job(namespace=namespace, name=genotypeCallByCoverage.name, version=version)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(refName))	#genotypeCallByCoverage_job would create directory "call".
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job.addArguments("-i", coverageFilterParentOutput, "-n", str(len(jobDataLs)), \
					"-o", genotypeCallOutput, '-e', refFastaF, '-y', str(genotypeCallerType), \
					'-s', repr(self.site_type))
			if seqCoverageF:
				genotypeCallByCoverage_job.addArguments("-q", seqCoverageF)
				genotypeCallByCoverage_job.uses(seqCoverageF, transfer=True, register=True, link=Link.INPUT)
			if coverageFilterParentOutput_bai:
				genotypeCallByCoverage_job.uses(coverageFilterParentOutput_bai, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
			genotypeCallByCoverage_job.uses(coverageFilterParentOutput, transfer=True, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(genotypeCallOutput, transfer=True, register=True, link=Link.OUTPUT)
			if self.site_type==1:	#all sites require additional memory
				job_max_memory = 5000	#in MB
			else:	#variants only, less memory
				job_max_memory=1000	#in MB. 
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			workflow.addJob(genotypeCallByCoverage_job)
			workflow.depends(parent=coverageFilterParentJob, child=genotypeCallByCoverage_job)
			
			
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
			workflow.depends(parent=genotypeCallByCoverage_job, child=calcula_job)
			
			refName2mergedBamCallJob[refName] = [genotypeCallOutput, genotypeCallByCoverage_job]
		sys.stderr.write(".Done\n")
		return PassingData(refName2mergedBamCallJob=refName2mergedBamCallJob)
	
	def addAddRG2BamJobsAsNeeded(self, workflow, alignmentLs, site_handler, input_site_handler=None, \
							addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
							BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
							mv=None, namespace="workflow", version="1.0", \
							dataDir=None):
		"""
		2011-9-15
			add a read group only when the alignment doesn't have it according to db record
			DBVervet.pokeBamReadGroupPresence() from misc.py helps to fill in db records if it's unclear.
		2011-9-14
			The read-group adding jobs will have a "move" part that overwrites the original bam&bai if site_handler and input_site_handler is same.
			For those alignment files that don't need to. It doesn't matter. pegasus will transfer/symlink them.
		"""
		sys.stderr.write("Adding add RG 2 BAM jobs for %s alignments ..."%(len(alignmentLs)))
		job_max_memory = 3500	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		no_of_rg_jobs = 0
		alignmentId2RGJobDataLs = {}
		for alignment in alignmentLs:
			if alignment.read_group_added==1:
				inputFname = os.path.join(dataDir, alignment.path)
				input = File(alignment.path)	#relative path, induces symlinking or stage-in
				input.addPFN(PFN("file://" + inputFname, input_site_handler))
				workflow.addFile(input)
				bai_output = File('%s.bai'%alignment.path)
				bai_output.addPFN(PFN("file://" + "%s.bai"%inputFname, input_site_handler))
				workflow.addFile(bai_output)
				alignmentId2RGJobDataLs[alignment.id] = [None, input, bai_output]
			else:
				# Add input file to the DAX-level replica catalog
				inputFname = os.path.join(dataDir, alignment.path)
				input = File(alignment.path)
				input.addPFN(PFN("file://" + inputFname, input_site_handler))
				workflow.addFile(input)
				# add RG to this bam
				sequencer = alignment.ind_sequence.sequencer
				read_group = '%s_%s_%s_%s_vs_%s'%(alignment.id, alignment.ind_seq_id, alignment.ind_sequence.individual.code, \
										sequencer, alignment.ref_ind_seq_id)
				read_group = alignment.getReadGroup()	##2011-11-02
				if sequencer=='454':
					platform_id = 'LS454'
				elif sequencer=='GA':
					platform_id = 'ILLUMINA'
				else:
					platform_id = 'ILLUMINA'
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=addOrReplaceReadGroupsJava.name, version=version)
				inputFileBaseNamePrefix = os.path.splitext(os.path.basename(alignment.path))[0]
				outputRGSAMFname = '%s.RG.bam'%(inputFileBaseNamePrefix)
				outputRGSAM = File(outputRGSAMFname)
	
				addRGJob.addArguments(javaMemRequirement, '-jar', addOrReplaceReadGroupsJar, \
									"INPUT=", input,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SORT_ORDER=coordinate', "VALIDATION_STRINGENCY=LENIENT")
									#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
				addRGJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#input is in absolute path. so don't register here
				workflow.addJob(addRGJob)
				addRGJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				addRGJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				
				# add the index job to the input (needs to be re-indexed??)
				index_sam_job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
				
				if input_site_handler==site_handler:	#on the same site. overwrite the original file without RG
					mvJob = Job(namespace=namespace, name=mv.name, version=version)
					mvJob.addArguments(outputRGSAM, inputFname)	#watch, it's inputFname, not input. input is in relative path.
					#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
					workflow.addJob(mvJob)
					workflow.depends(parent=addRGJob, child=mvJob)
					bai_output = File('%s.bai'%inputFname)	#in absolute path, don't register it to the job
				else:
					##on different site, input for index should be outputRGSAM and register it as well
					mvJob = addRGJob
					input = outputRGSAM	
					addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
					bai_output = File('%s.bai'%outputRGSAMFname)
					index_sam_job.uses(bai_output, transfer=True, register=False, link=Link.OUTPUT)
				
				index_sam_job.addArguments(javaMemRequirement, "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
								"INPUT=", input, \
								"OUTPUT=", bai_output)
				# input is in relative path, either symlink to inputFname (same site) or outputRGSAM (different site)
				index_sam_job.uses(input, transfer=True, register=False, link=Link.INPUT)
				
				if input_site_handler==site_handler:	#on the same site. need to change bai_output's name to relative path and register it as well
					#for symlinking. this has to be done after addArguments() for index_sam_job is done.
					bai_output = File('%s.bai'%alignment.path)
					bai_output.addPFN(PFN("file://" + '%s.bai'%inputFname, input_site_handler))
					workflow.addFile(bai_output)
				
				index_sam_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				index_sam_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				
				workflow.addJob(index_sam_job)
				workflow.depends(parent=mvJob, child=index_sam_job)
				alignmentId2RGJobDataLs[alignment.id]= [index_sam_job, input, bai_output]
				no_of_rg_jobs += 1
		sys.stderr.write(" %s alignments need read-group addition. Done\n"%(no_of_rg_jobs))
		return PassingData(alignmentId2RGJobDataLs=alignmentId2RGJobDataLs)
	
	def addRefFastaJobDependency(self, workflow, job, refFastaF=None, fastaDictJob=None, refFastaDictF=None, fastaIndexJob = None, refFastaIndexF = None):
		"""
		2011-9-14
		"""
		if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
			workflow.depends(parent=fastaIndexJob, child=job)
			job.uses(refFastaIndexF, transfer=True, register=True, link=Link.INPUT)
		if fastaDictJob:
			workflow.depends(parent=fastaDictJob, child=job)
			job.uses(refFastaDictF, transfer=True, register=True, link=Link.INPUT)
		if fastaIndexJob or fastaDictJob:
			job.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
	
	def addVCFFormatConvertJob(self, workflow, vcf_convert=None, parentJob=None, inputF=None, outputF=None, \
							namespace=None, version=None, transferOutput=False):
		"""
		2011-11-4
		"""
		vcf_convert_job = Job(namespace=namespace, name=vcf_convert.name, version=version)
		vcf_convert_job.addArguments(inputF, outputF)
		vcf_convert_job.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		vcf_convert_job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(vcf_convert_job)
		workflow.depends(parent=parentJob, child=vcf_convert_job)
		return vcf_convert_job
	
	def addBGZIP_tabix_Job(self, workflow, bgzip_tabix=None, parentJob=None, inputF=None, outputF=None, \
							namespace=None, version=None, transferOutput=False):
		"""
		2011-11-4
		"""
		bgzip_tabix_job = Job(namespace=namespace, name=bgzip_tabix.name, version=version)
		tbi_F = File("%s.tbi"%outputF.name)
		bgzip_tabix_job.addArguments(inputF, outputF)
		bgzip_tabix_job.uses(inputF, transfer=False, register=True, link=Link.INPUT)
		bgzip_tabix_job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		bgzip_tabix_job.uses(tbi_F, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(bgzip_tabix_job)
		workflow.depends(parent=parentJob, child=bgzip_tabix_job)
		return bgzip_tabix_job
	
	def addGenotypeCallJobs(self, workflow, alignmentLs, refName2size, samtools=None, \
				genotyperJava=None,  genomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				mergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, CallVariantBySamtools=None,\
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None,\
				genotypeCallByCoverage=None, refFastaF=None, bamListF=None, \
				callOutputDirJob =None, gatkDirJob=None, samtoolsDirJob=None, unionDirJob=None, intersectionDirJob=None,\
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				seqCoverageF=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				chunkSize=2000000, site_type=1, dataDir=None):
		"""
		2011-9-22
			add argument concatGATK, concatSamtools.
		2011-9-15
			bamListF is now useless. samtools_job could accept variable-length list of bam input files
		2011-9-14
			argument chunkSize determines how many sites gatk/samtools works on at a time
		"""
		sys.stderr.write("Adding genotype call jobs for %s references ..."%(len(refName2size)))
		job_max_memory = 2000	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		no_of_gatk_threads = 3
		vcf_job_max_memory = 1000
		if needFastaDictJob:	# the .dict file is required for GATK
			refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
			refFastaDictF = File(refFastaDictFname)
			#not os.path.isfile(refFastaDictFname) or 
			fastaDictJob = Job(namespace=namespace, name=createSequenceDictionaryJava.name, version=version)
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
		
		returnData = self.addAddRG2BamJobsAsNeeded(workflow, alignmentLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		
		# add merge jobs for every reference
		refName2mergedBamCallJob = {}
		for refName, refSize in refName2size.iteritems():
			no_of_chunks = int(math.ceil(refSize/float(chunkSize)))
			wholeRefUnionOfIntersectionJob = Job(namespace=namespace, name=vcf_concat.name, version=version)
			callOutputFname = os.path.join(callOutputDirJob.folder, '%s.vcf.gz'%refName)
			callOutputF = File(callOutputFname)
			wholeRefUnionOfIntersectionJob.addArguments(callOutputF)
			
			wholeRefUnionOfIntersectionJob.uses(callOutputF, transfer=True, register=True, link=Link.OUTPUT)
			callOutput_tbi_F = File("%s.tbi"%callOutputFname)
			wholeRefUnionOfIntersectionJob.uses(callOutput_tbi_F, transfer=True, register=True, link=Link.OUTPUT)
			
			wholeRefUnionOfIntersectionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%vcf_job_max_memory))
			wholeRefUnionOfIntersectionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%vcf_job_max_memory))
			
			workflow.addJob(wholeRefUnionOfIntersectionJob)
			workflow.depends(parent=callOutputDirJob, child=wholeRefUnionOfIntersectionJob)
			
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			
			#2011-9-22 union of all GATK intervals for one contig
			gatkUnionJob = Job(namespace=namespace, name=concatGATK.name, version=version)
			gatkUnionOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf.gz'%refName)
			gatkUnionOutputF = File(gatkUnionOutputFname)
			gatkUnionJob.addArguments(gatkUnionOutputF)
			gatkUnionJob.uses(gatkUnionOutputF, transfer=True, register=True, link=Link.OUTPUT)
			gatkUnionOutput_tbi_F = File("%s.tbi"%gatkUnionOutputFname)
			gatkUnionJob.uses(gatkUnionOutput_tbi_F, transfer=True, register=True, link=Link.OUTPUT)
			
			gatkUnionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%vcf_job_max_memory))
			gatkUnionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%vcf_job_max_memory))
			workflow.addJob(gatkUnionJob)
			workflow.depends(parent=gatkDirJob, child=gatkUnionJob)
			
			#2011-9-22 union of all samtools intervals for one contig
			samtoolsUnionJob = Job(namespace=namespace, name=concatSamtools.name, version=version)
			samtoolsUnionOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf.gz'%refName)
			samtoolsUnionOutputF = File(samtoolsUnionOutputFname)
			samtoolsUnionJob.addArguments(samtoolsUnionOutputF)
			samtoolsUnionJob.uses(samtoolsUnionOutputF, transfer=True, register=True, link=Link.OUTPUT)
			samtoolsUnionOutput_tbi_F = File("%s.tbi"%samtoolsUnionOutputFname)
			samtoolsUnionJob.uses(samtoolsUnionOutput_tbi_F, transfer=True, register=True, link=Link.OUTPUT)
			
			samtoolsUnionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%vcf_job_max_memory))
			samtoolsUnionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%vcf_job_max_memory))
			workflow.addJob(samtoolsUnionJob)
			workflow.depends(parent=samtoolsDirJob, child=samtoolsUnionJob)
			
			for i in range(no_of_chunks):
				startPos = i*chunkSize + 1
				stopPos = min((i+1)*chunkSize, refSize)
				interval = "%s:%s-%s"%(refName, startPos, stopPos)
				vcfBaseFname = '%s_%s_%s'%(refName, startPos, stopPos)
				
				#GATK job
				gatk_job = Job(namespace=namespace, name=genotyperJava.name, version=version)
				gatkOutputFname = os.path.join(gatkDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
				gatkOutputF = File(gatkOutputFname)
				gatkIDXOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf.idx'%(vcfBaseFname))
				gatkIDXOutput = File(gatkIDXOutputFname)
				gatk_job.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "UnifiedGenotyper",\
					"-L", interval, "-mbq 20", "-mmq 30", "-R", refFastaF, "--out", gatkOutputF,\
					'-U', '-S SILENT', "-nt %s"%no_of_gatk_threads)
				if site_type==1:
					gatk_job.addArguments('--output_mode EMIT_ALL_SITES')	#2011-8-24 new GATK no longers ues "-all_bases"
				
				gatk_job.uses(gatkOutputF, transfer=False, register=True, link=Link.OUTPUT)
				gatk_job.uses(gatkIDXOutput, transfer=False, register=True, link=Link.OUTPUT)
				
				gatk_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				gatk_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				workflow.addJob(gatk_job)
				workflow.depends(parent=gatkDirJob, child=gatk_job)
				
				vcf_convert_gatkOutputF_job = Job(namespace=namespace, name=vcf_convert.name, version=version)
				vcf4_gatkOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf'%vcfBaseFname)
				vcf4_gatkOutputF = File(vcf4_gatkOutputFname)
				vcf_convert_gatkOutputF_job.addArguments(gatkOutputF, vcf4_gatkOutputF)
				vcf_convert_gatkOutputF_job.uses(gatkOutputF, transfer=False, register=True, link=Link.INPUT)
				vcf_convert_gatkOutputF_job.uses(vcf4_gatkOutputF, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(vcf_convert_gatkOutputF_job)
				workflow.depends(parent=gatk_job, child=vcf_convert_gatkOutputF_job)
				
				bgzip_tabix_gatkOutputF_job = Job(namespace=namespace, name=bgzip_tabix.name, version=version)
				gatkGzipOutputF = File("%s.gz"%vcf4_gatkOutputFname)
				gatkGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_gatkOutputFname)
				bgzip_tabix_gatkOutputF_job.addArguments(vcf4_gatkOutputF, gatkGzipOutputF)
				bgzip_tabix_gatkOutputF_job.uses(vcf4_gatkOutputF, transfer=False, register=True, link=Link.INPUT)
				bgzip_tabix_gatkOutputF_job.uses(gatkGzipOutputF, transfer=False, register=True, link=Link.OUTPUT)
				bgzip_tabix_gatkOutputF_job.uses(gatkGzipOutput_tbi_F, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(bgzip_tabix_gatkOutputF_job)
				workflow.depends(parent=vcf_convert_gatkOutputF_job, child=bgzip_tabix_gatkOutputF_job)
				
				#add this output to a GATK union job
				gatkUnionJob.addArguments(gatkGzipOutputF)
				gatkUnionJob.uses(gatkGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				gatkUnionJob.uses(gatkGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=bgzip_tabix_gatkOutputF_job, child=gatkUnionJob)
				
				#samtools part
				samtools_job = Job(namespace=namespace, name=CallVariantBySamtools.name, version=version)
				samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
				samtoolsOutputF = File(samtoolsOutputFname)
				indelVCFOutputFname = "%s.indel.vcf"%(samtoolsOutputFname)
				indelVCFOutputF = File(indelVCFOutputFname)
				
				samtools_job.addArguments(refFastaF, interval, samtoolsOutputF, repr(site_type))
				#samtools_job.uses(bamListF, transfer=True, register=True, link=Link.INPUT)
				samtools_job.uses(samtoolsOutputF, transfer=False, register=True, link=Link.OUTPUT)
				samtools_job.uses(indelVCFOutputF, transfer=False, register=True, link=Link.OUTPUT)
				
				samtools_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%vcf_job_max_memory))
				samtools_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%vcf_job_max_memory))
				workflow.addJob(samtools_job)
				workflow.depends(parent=samtoolsDirJob, child=samtools_job)
				
				#deal with samtools's indel VCF
				samtoolsIndelVCF4Fname = os.path.join(samtoolsDirJob.folder, '%s.indel.vcf'%vcfBaseFname)
				samtoolsIndelVCF4F = File(samtoolsIndelVCF4Fname)
				samtoolsIndelVCF4F_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=samtools_job, inputF=indelVCFOutputF, outputF=samtoolsIndelVCF4F, \
							namespace=namespace, version=version, transferOutput=False)
				
				samtoolsIndelGzipOutputF = File("%s.gz"%samtoolsIndelVCF4Fname)
				samtoolsIndelGzipOutput_tbi_F = File("%s.gz.tbi"%samtoolsIndelVCF4Fname)
				samtoolsIndelGzipOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=samtoolsIndelVCF4F_job, inputF=samtoolsIndelVCF4F, outputF=samtoolsIndelGzipOutputF, \
						namespace=namespace, version=version, transferOutput=True)
				
				#deal with samtools's snp VCF
				vcf4_samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf'%vcfBaseFname)
				vcf4_samtoolsOutputF = File(vcf4_samtoolsOutputFname)
				vcf_convert_samtoolsOutputF_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=samtools_job, inputF=samtoolsOutputF, outputF=vcf4_samtoolsOutputF, \
							namespace=namespace, version=version, transferOutput=False)
				
				samtoolsGzipOutputF = File("%s.gz"%vcf4_samtoolsOutputFname)
				samtoolsGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_samtoolsOutputFname)
				bgzip_tabix_samtoolsOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=vcf_convert_samtoolsOutputF_job, inputF=vcf4_samtoolsOutputF, outputF=samtoolsGzipOutputF, \
						namespace=namespace, version=version, transferOutput=False)
				
				#add this output to a samtools union job
				samtoolsUnionJob.addArguments(samtoolsGzipOutputF)
				samtoolsUnionJob.uses(samtoolsGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				samtoolsUnionJob.uses(samtoolsGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=bgzip_tabix_samtoolsOutputF_job, child=samtoolsUnionJob)
				
				
				# intersection of GATK and samtools
				isectJob = Job(namespace=namespace, name=vcf_isec.name, version=version)
				intersectionOutputFname = os.path.join(intersectionDirJob.folder, '%s.isec.vcf.gz'%vcfBaseFname)
				intersectionOutputF = File(intersectionOutputFname)
				intersectionOutput_tbi_F = File("%s.tbi"%intersectionOutputFname)
				
				#put samtools in front of GATK cuz the former outputs better-quality data based on inconsistent rates
				#	estimated from trios
				isectJob.addArguments(intersectionOutputF, samtoolsGzipOutputF, gatkGzipOutputF)
				isectJob.uses(samtoolsGzipOutputF, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(samtoolsGzipOutput_tbi_F, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(gatkGzipOutputF, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(gatkGzipOutput_tbi_F, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(intersectionOutputF, transfer=False, register=True, link=Link.OUTPUT)
				isectJob.uses(intersectionOutput_tbi_F, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(isectJob)
				workflow.depends(parent=intersectionDirJob, child=isectJob)
				workflow.depends(parent=bgzip_tabix_gatkOutputF_job, child=isectJob)
				workflow.depends(parent=bgzip_tabix_samtoolsOutputF_job, child=isectJob)
				
				isectJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%vcf_job_max_memory))
				isectJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%vcf_job_max_memory))
				#union of intersected-intervals from one contig
				wholeRefUnionOfIntersectionJob.addArguments(intersectionOutputF)
				wholeRefUnionOfIntersectionJob.uses(intersectionOutputF, transfer=False, register=True, link=Link.INPUT)
				wholeRefUnionOfIntersectionJob.uses(intersectionOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=isectJob, child=wholeRefUnionOfIntersectionJob)
				
				"""
				java -jar /path/to/dist/GenomeAnalysisTK.jar \
				  -T CombineVariants \
				  -R /path/to/reference.fasta \
				  -o /path/to/output.vcf \
				  -B:foo,VCF /path/to/variants1.vcf \
				  -B:bar,VCF /path/to/variants2.vcf \
				  -genotypeMergeOptions PRIORITIZE \
				
				java -Xmx2g -jar dist/GenomeAnalysisTK.jar -T SelectVariants 
				-R ~/Desktop/broadLocal/localData/human_g1k_v37.fasta 
				-L 1:1-1,000,000 -B:variant,VCF union.vcf 
				-select 'set == "Intersection"' -o intersect.vcf
				
				CombineVariants has problem working with samtools vcf output in which GT is not placed in the first position.
				
				genotypeUnionJob = Job(namespace=namespace, name=java.name, version=version)
				unionOutputFname = os.path.join(unionDirJob.folder, '%s.vcf'%vcfBaseFname)
				unionOutputF = File(unionOutputFname)
				genotypeUnionJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "CombineVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:gatk,VCF", gatkOutputF, "-B:samtools,VCF", samtoolsOutputF, \
					"--out", unionOutputF, \
					'-U', '-S LENIENT', "-genotypeMergeOptions PRIORITIZE")	#no support for multi-threading
				genotypeUnionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeUnionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				genotypeUnionJob.uses(gatkOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(samtoolsOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(unionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeUnionJob)
				workflow.depends(parent=unionDirJob, child=genotypeUnionJob)
				workflow.depends(parent=samtools_job, child=genotypeUnionJob)
				workflow.depends(parent=gatk_job, child=genotypeUnionJob)
				
				genotypeIntersectionJob = Job(namespace=namespace, name=java.name, version=version)
				intersectionOutputFname = os.path.join(intersectionDirJob.folder, '%s.vcf'%vcfBaseFname)
				intersectionOutputF = File(intersectionOutputFname)
				genotypeIntersectionJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "SelectVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:%s,VCF"%vcfBaseFname, unionOutputF, \
					"-select 'set == "Intersection"'", "--out", intersectionOutputF, \
					'-U', '-S LENIENT')	#no support for multi-threading
				genotypeIntersectionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeIntersectionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				genotypeIntersectionJob.uses(unionOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeIntersectionJob.uses(intersectionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeIntersectionJob)
				workflow.depends(parent=intersectionDirJob, child=genotypeIntersectionJob)
				workflow.depends(parent=genotypeUnionJob, child=genotypeIntersectionJob)
				
				#add this intersected output to the union job for the whole refName
				wholeRefUnionOfIntersectionJob.addArguments("-B:%s,VCF"%vcfBaseFname, intersectionOutputF)
				workflow.depends(parent=genotypeIntersectionJob, child=wholeRefUnionOfIntersectionJob)
				"""

				lisOfJobs = [gatk_job, samtools_job]
				for job in lisOfJobs:
					self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
				
				for alignment in alignmentLs:
					if alignment.id in alignmentId2RGJobDataLs:
						index_sam_job, bamF, bai_output = alignmentId2RGJobDataLs[alignment.id]
						gatk_job.addArguments("-I", bamF)
						samtools_job.addArguments(bamF)
						#it's either symlink or stage-in
						gatk_job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
						gatk_job.uses(bai_output, transfer=True, register=True, link=Link.INPUT)
						samtools_job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
						samtools_job.uses(bai_output, transfer=True, register=True, link=Link.INPUT)
						if index_sam_job is not None:
							workflow.depends(parent=index_sam_job, child=gatk_job)
							workflow.depends(parent=index_sam_job, child=samtools_job)
		
		sys.stderr.write(".Done\n")
		return None
	
	def outputSeqCoverage(self, outputFname, ind_seq_id_ls=[]):
		"""
		2011-9-2
		"""
		sys.stderr.write("Outputting sequence coverage to %s ..."%outputFname)
		import csv
		counter = 0
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['isq.id', 'coverage'])
		TableClass = VervetDB.IndividualSequence
		query = TableClass.query.filter(TableClass.coverage!=None)
		if ind_seq_id_ls:
			query = query.filter(TableClass.id.in_(ind_seq_id_ls))
		query = query.order_by(TableClass.id)
		for row in query:
			writer.writerow([row.id, row.coverage])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
	def outputListOfBam(self, alignmentLs, dataDir=None, outputFname=None, workflow=None, site_handler=None):
		"""
		2011-9-15
			deprecated. CallVariantBySamtools.sh doesn't need this file anymore.
			using relative path, alignments are symlinked or transferred
		2011-9-14
			samtools mpileup uses this file.
		"""
		sys.stderr.write("Outputting path to %s alignments to %s ..."%(len(alignmentLs), outputFname))
		outf = open(outputFname, 'w')
		for alignment in alignmentLs:
			#abs_path = os.path.join(dataDir, alignment.path)
			outf.write("%s\n"%alignment.path)	#using relative path, alignments are symlinked or transferred
		del outf
		bamListF = File(os.path.basename(outputFname))
		bamListF.addPFN(PFN("file://" + os.path.abspath(outputFname), site_handler))
		workflow.addFile(bamListF)
		sys.stderr.write("Done.\n")
		return bamListF
	
	def addMkDirJob(self, workflow, mkdir=None, outputDir=None, namespace=None, version=None):
		"""
		2011-9-14
		"""
		return yh_pegasus.addMkDirJob(workflow, mkdir=mkdir, outputDir=outputDir, namespace=namespace, version=version)
		"""
		# Add a mkdir job for any directory.
		mkDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
		mkDirJob.addArguments(outputDir)
		mkDirJob.folder = outputDir	#custom attribute
		workflow.addJob(mkDirJob)
		return mkDirJob
		"""
		
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
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		refName2size = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refName2size = set(['Contig149'])	#temporary when testing Contig149
		#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = refName2size.keys()
		
		alignmentLs = self.getAlignments(self.aln_ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2, dataDir=self.localDataDir)
		
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
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = ADAG(workflowName)
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		#add the MergeSamFiles.jar file into workflow
		abs_path = os.path.join(self.picard_path, 'MergeSamFiles.jar')
		mergeSamFilesJar = File(abs_path)	#using abs_path avoids add this jar to every job as Link.INPUT
		mergeSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(mergeSamFilesJar)
		
		abs_path = os.path.join(self.picard_path, 'BuildBamIndex.jar')
		BuildBamIndexFilesJar = File(abs_path)
		BuildBamIndexFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(BuildBamIndexFilesJar)
		
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
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroups.jar')
		addOrReplaceReadGroupsJar = File(abs_path)
		addOrReplaceReadGroupsJar.addPFN(PFN("file://" + abs_path, \
											site_handler))
		workflow.addFile(addOrReplaceReadGroupsJar)
		
		
		"""
		#2011-9-2
		self.outputSeqCoverage(self.seqCoverageFname)
		seqCoverageF = File(os.path.basename(self.seqCoverageFname))
		seqCoverageF.addPFN(PFN("file://" + os.path.abspath(self.seqCoverageFname), \
											self.input_site_handler))
		workflow.addFile(seqCoverageF)
		"""
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		clusters_size = 20
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mkdirWrap)
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		mv.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mv)
		
		selectAndSplit = Executable(namespace=namespace, name="SelectAndSplitAlignment", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		workflow.addExecutable(selectAndSplit)
		
		"""
		# 2011-9-1 not used anymore. too slow. replace with addOrReplaceReadGroupsAndCleanSQHeaderJar.
		addRGAndCleanSQHeaderAlignment = Executable(namespace=namespace, name="AddRGAndCleanSQHeaderAlignment", version=version, \
												os=operatingSystem, arch=architecture, installed=True)
		addRGAndCleanSQHeaderAlignment.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AddRGAndCleanSQHeaderAlignment.py"), site_handler))
		workflow.addExecutable(addRGAndCleanSQHeaderAlignment)
		"""
		# 2011-9-1 used
		addRGToBAM = Executable(namespace=namespace, name="addRGToBAM", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		addRGToBAM.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "addRGToBAM.sh"), site_handler))
		workflow.addExecutable(addRGToBAM)
		
		
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
		
		addOrReplaceReadGroupsJava = Executable(namespace=namespace, name="addOrReplaceReadGroupsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		addOrReplaceReadGroupsJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(addOrReplaceReadGroupsJava)
		
		genotyperJava = Executable(namespace=namespace, name="genotyperJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		genotyperJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(genotyperJava)
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(BuildBamIndexFilesJava)
		
		createSequenceDictionaryJava = Executable(namespace=namespace, name="createSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		createSequenceDictionaryJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(createSequenceDictionaryJava)
		
		
		CallVariantBySamtools = Executable(namespace=namespace, name="CallVariantBySamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CallVariantBySamtools.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "CallVariantBySamtools.sh"), site_handler))
		workflow.addExecutable(CallVariantBySamtools)
		
		genotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		workflow.addExecutable(genotypeCallByCoverage)
		
		mergeGenotypeMatrix = Executable(namespace=namespace, name="MergeGenotypeMatrix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		mergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "MergeGenotypeMatrix.py"), site_handler))
		workflow.addExecutable(mergeGenotypeMatrix)
		
		bgzip_tabix = Executable(namespace=namespace, name="bgzip_tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		bgzip_tabix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "bgzip_tabix.sh"), site_handler))
		bgzip_tabix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(bgzip_tabix)
		
		vcf_convert = Executable(namespace=namespace, name="vcf_convert", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		vcf_convert.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "vcf_convert.sh"), site_handler))
		vcf_convert.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_convert)
		
		vcf_isec = Executable(namespace=namespace, name="vcf_isec", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_isec.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "vcf_isec.sh"), site_handler))
		vcf_isec.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_isec)
		
		vcf_concat = Executable(namespace=namespace, name="vcf_concat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcf_concat.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "vcf_concat.sh"), site_handler))
		vcf_concat.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcf_concat)
		
		concatGATK = Executable(namespace=namespace, name="concatGATK", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		concatGATK.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "vcf_concat.sh"), site_handler))
		concatGATK.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(concatGATK)
		
		concatSamtools = Executable(namespace=namespace, name="concatSamtools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		concatSamtools.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "vcf_concat.sh"), site_handler))
		concatSamtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(concatSamtools)
		
		calcula = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), site_handler))
		workflow.addExecutable(calcula)
		
		#2011-11-4 for single-alignment-calling pipeline, adjust the folder name so that they are unique from each other
		if self.run_type==1:
			dirPrefix = ""
			# Add a mkdir job for the call directory.
			callOutputDir = "%scall"%(dirPrefix)
			callOutputDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=callOutputDir, namespace=namespace, version=version)
			gatkDir = "%sgatk"%(dirPrefix)
			gatkDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=gatkDir, namespace=namespace, version=version)
			samtoolsDir = "%ssamtools"%(dirPrefix)
			samtoolsDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=samtoolsDir, namespace=namespace, version=version)
			unionDir = "%sgatk_samtools_union"%(dirPrefix)
			unionDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=unionDir, namespace=namespace, version=version)
			intersectionDir = "%sgatk_samtools_intersection"%(dirPrefix)
			intersectionDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=intersectionDir, namespace=namespace, version=version)
			
			if self.site_handler=='condorpool':
				bamListF_site_handler = 'condorpool'
			else:	#every other site requires transfer from local
				bamListF_site_handler = 'local'
			#bamListF = self.outputListOfBam(alignmentLs, self.dataDir, self.bamListFname, workflow=workflow, site_handler=bamListF_site_handler)
			
			
			self.addGenotypeCallJobs(workflow, alignmentLs, refName2size, samtools=samtools, \
					genotyperJava=genotyperJava, genomeAnalysisTKJar=genomeAnalysisTKJar, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					createSequenceDictionaryJava=createSequenceDictionaryJava, createSequenceDictionaryJar=createSequenceDictionaryJar, \
					mergeSamFilesJar=mergeSamFilesJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, CallVariantBySamtools=CallVariantBySamtools,\
					bgzip_tabix=bgzip_tabix, vcf_convert=vcf_convert, vcf_isec=vcf_isec, vcf_concat=vcf_concat, \
					concatGATK=concatGATK, concatSamtools=concatSamtools,\
					genotypeCallByCoverage=genotypeCallByCoverage, refFastaF=refFastaF, bamListF=None, \
					callOutputDirJob =callOutputDirJob, gatkDirJob=gatkDirJob, samtoolsDirJob=samtoolsDirJob, unionDirJob=unionDirJob, intersectionDirJob=intersectionDirJob,\
					namespace=namespace, version=version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
					seqCoverageF=None, \
					needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
					chunkSize=2000000, site_type=self.site_type, dataDir=self.dataDir)
		elif self.run_type==2:
			for alignment in alignmentLs:
				dirPrefix = "%s"%(alignment.id)
				# Add a mkdir job for the call directory.
				callOutputDir = "%s_call"%(dirPrefix)
				callOutputDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=callOutputDir, namespace=namespace, version=version)
				gatkDir = "%s_gatk"%(dirPrefix)
				gatkDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=gatkDir, namespace=namespace, version=version)
				samtoolsDir = "%s_samtools"%(dirPrefix)
				samtoolsDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=samtoolsDir, namespace=namespace, version=version)
				unionDir = "%s_gatk_samtools_union"%(dirPrefix)
				unionDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=unionDir, namespace=namespace, version=version)
				intersectionDir = "%s_gatk_samtools_intersection"%(dirPrefix)
				intersectionDirJob = self.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=intersectionDir, namespace=namespace, version=version)
				
				if self.site_handler=='condorpool':
					bamListF_site_handler = 'condorpool'
				else:	#every other site requires transfer from local
					bamListF_site_handler = 'local'
				#bamListF = self.outputListOfBam(alignmentLs, self.dataDir, self.bamListFname, workflow=workflow, site_handler=bamListF_site_handler)
				
				
				self.addGenotypeCallJobs(workflow, [alignment], refName2size, samtools=samtools, \
						genotyperJava=genotyperJava, genomeAnalysisTKJar=genomeAnalysisTKJar, \
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
						createSequenceDictionaryJava=createSequenceDictionaryJava, createSequenceDictionaryJar=createSequenceDictionaryJar, \
						mergeSamFilesJar=mergeSamFilesJar, \
						BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
						mv=mv, CallVariantBySamtools=CallVariantBySamtools,\
						bgzip_tabix=bgzip_tabix, vcf_convert=vcf_convert, vcf_isec=vcf_isec, vcf_concat=vcf_concat, \
						concatGATK=concatGATK, concatSamtools=concatSamtools,\
						genotypeCallByCoverage=genotypeCallByCoverage, refFastaF=refFastaF, bamListF=None, \
						callOutputDirJob =callOutputDirJob, gatkDirJob=gatkDirJob, samtoolsDirJob=samtoolsDirJob, unionDirJob=unionDirJob, intersectionDirJob=intersectionDirJob,\
						namespace=namespace, version=version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						seqCoverageF=None, \
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						chunkSize=2000000, site_type=self.site_type, dataDir=self.dataDir)
				
		
		"""
		addSelectAndSplitBamReturnData = self.addSelectAndSplitBamJobs(db_vervet, workflow, alignmentLs, site_handler, \
							self.topNumberOfContigs, refNameLs, samtools=samtools, \
							java=java, addOrReplaceReadGroupsAndCleanSQHeaderJar=addOrReplaceReadGroupsAndCleanSQHeaderJar, \
							BuildBamIndexFilesJar=BuildBamIndexFilesJar,\
							mkdir=mkdirWrap, namespace=namespace, \
							version=version, mkCallDirJob=callOutputDirJob,\
							addRGExecutable=addRGToBAM, dataDir=self.dataDir)
		refName2jobDataLs = addSelectAndSplitBamReturnData.refName2jobDataLs
		
		returnData3 = self.addRefFastaFileSplitJobs(workflow, refFastaF, selectAndSplitFasta, refNameLs, mkdir=mkdirWrap, samtools=samtools,
								java=java, createSequenceDictionaryJar=createSequenceDictionaryJar,\
								site_handler=site_handler, namespace=namespace, version=version)
		refName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		returnData2 = self.addMergeAlignmentAndGenotypeCallJobs(workflow, refName2jobDataLs, refNameLs, samtools, \
							java, createSequenceDictionaryJar=createSequenceDictionaryJar, genotypeCallByCoverage=genotypeCallByCoverage, \
							refFastaF=refFastaF, \
							namespace=namespace, version=version, callOutputDir = callOutputDir, \
							genotypeCallerType=self.genotypeCallerType, \
							mergeSamFilesJar=mergeSamFilesJar, genomeAnalysisTKJar=genomeAnalysisTKJar, calcula=calcula, \
							refName2splitFastaJobDataLs=refName2splitFastaJobDataLs, seqCoverageF=seqCoverageF, \
							needFastaIndexJob=self.needFastaIndexJob, \
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
			workflow.depends(parent=callJob, child=mergeGenotypeMatrix_job)
		
			mergeGenotypeMatrix_job.addArguments(callJobOutput)
		"""
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = AlignmentToCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
