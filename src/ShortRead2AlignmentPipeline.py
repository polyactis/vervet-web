#!/usr/bin/env python
"""
Examples:
	# 2011-8-30 workflow on condor
	%s -i 165-167 -o ShortRead2AlignmentPipeline_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool -n1 -z dl324b-1.cmb.usc.edu -c
	
	# 2011-8-30 a workflow with 454 long-read and short-read PE 
	%s -i 165-167 -o ShortRead2AlignmentPipeline_isq_id_165_167_vs_9.xml -u yh -a 9
	-e /u/home/eeskin/polyacti -l hoffman2 -t /u/home/eeskin/polyacti/NetworkData/vervet/db -n1 -z dl324b-1.cmb.usc.edu -c
	
	# 2011-8-30 output a workflow to run alignments on hoffman2 (and add -D if your local db_vervet.data_dir is outdated.)
	%s -i 165-495 -o ShortRead2AlignmentPipeline_isq_id_165_495_vs_120.xml -u yh -a 120 -e /u/home/eeskin/polyacti 
		-l hoffman2 -t /u/home/eeskin/polyacti/NetworkData/vervet/db -n1 -z dl324b-1.cmb.usc.edu -c
		-D ~/mnt/hoffman2/u/home/eeskintmp2/polyacti/NetworkData/vervet/db/
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files.
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2AlignmentPipeline_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool -n0 -z dl324b-1.cmb.usc.edu -p secret  -c
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files.
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2AlignmentPipeline_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool -n0 -z dl324b-1.cmb.usc.edu -p secret  -c
		
	# 2011-8-30 a workflow to run on uschpc, with ref index job. Note the site_handler and input_site_handler.
	# to enable replica-transfer.
	%s -i 391-397,456,473,493
		-o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j local -l uschpc -n1 -e /home/rcf-47/yuhuang -z 10.8.0.10 -p secret  -c

	# 2011-8-30 a workflow to run on uschpc, with ref index job, and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler and "-t ..." to enable symlink
	%s -i 391-397,456,473,493
		-o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j uschpc -l uschpc -n1 -p secret -c -m 4
		-e /home/rcf-47/yuhuang -z 10.8.0.1 
		-t /home/rcf-47/yuhuang/NetworkData/vervet/db/
	
	# 2011-8-31 output the same workflow above but for condorpool
	%s -i 391-397,456,473,493, -o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
		-u yh -a 9 -j condorpool -l condorpool -n0 -z 10.8.0.10  -p secret  -c
	
	# 2011-8-31
	%s -i 167,176,178,182,183,207-211,391-397,456,473,493
		-o ShortRead2AlignmentPipeline_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
		-u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  -c
	
Description:
	2011-8-30
		a program which generates a pegasus workflow dag (xml file) which does the alignment for all available sequences.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
		It will also stage out every output file.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0])

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from Pegasus.DAX3 import *


class ShortRead2AlignmentPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number. must be non-empty if need ssh tunnel'],\
						('sshTunnelCredential', 0, ): ['', 's', 1, 'a ssh credential to allow machine to access db server. \
										polyacti@login3, yuhuang@hpc-login2. if empty or port is empty, no tunnel', ],\
						
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						('additionalArguments', 0, ): ["-q 20", '', 1, 'a string of additional arguments passed to aln, not bwasw, add double quote if space'],\
						("bwa_path", 1, ): ["%s/bin/bwa", '', 1, 'bwa binary'],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/vervet/bin/GenomeAnalysisTK", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						("alignment_method_name", 1, ): ["bwa-short-read", '', 1, 'alignment_method.short_name from db.'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local will induce file transfer and input_site_handler=hoffman2 induces symlink.'],\
						("needRefIndexJob", 0, int): [0, 'n', 1, 'need to add a reference index job by bwa?'],\
						('no_of_aln_threads', 1, int): [3, 'm', 1, 'number of threads during alignment'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('stageOutFinalOutput', 0, int):[0, 'O', 0, 'toggle to stage out final output (bam + bam.bai)'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
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
		
		self.bwa_path = self.bwa_path%self.home_path
		self.samtools_path = self.samtools_path%self.home_path
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	
	def addAllAlignmentJobs(self, db_vervet, individualSequenceID2FilePairLs=None, dataDir=None, \
					refSequence=None, refFastaF=None, refIndexJob=None,
					workflow=None, bwa=None, additionalArguments=None, samtools=None, mkdirWrap=None, mv=None,\
					alignment_method_name='bwa-short-read', alignment_format='bam',\
					namespace='workflow', version='1.0', stageOutFinalOutput=False,\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, no_of_aln_threads=3):
		"""
		2011-8-28
		"""
		sys.stderr.write("Adding alignment jobs for %s individual sequences ..."%(len(individualSequenceID2FilePairLs)))
		
		no_of_alignment_jobs = 0
		
		for ind_seq_id, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if individual_sequence is not None and individual_sequence.format=='fastq':
				
				tmpOutputDir = os.path.basename(individual_sequence.path)
				# add a mkdir job
				mkdirJob = self.addMKDIRJob(workflow, mkdirWrap=mkdirWrap, dirName=tmpOutputDir, \
										namespace=namespace, version=version)
				
				# get or add an alignment
				alignment_method = db_vervet.getAlignmentMethod(alignment_method_name)
				individual_alignment = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
										individual_sequence_id=ind_seq_id,\
							path_to_original_alignment=None, sequencer=individual_sequence.sequencer, \
							sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
							ref_individual_sequence_id=refSequence.id, \
							alignment_method_name=alignment_method_name, alignment_format=alignment_format,\
							individual_sequence_filtered=individual_sequence.filtered)
				if not individual_alignment.path:
					individual_alignment.path = db_vervet.constructRelativePathForIndividualAlignment(individual_alignment_id=individual_alignment.id, \
							individual_sequence_id=ind_seq_id, \
							ref_individual_sequence_id=refSequence.id, alignment_method=alignment_method, \
							alignment_format=alignment_format,)
					session.add(individual_alignment)
					session.flush()
				
				AlignmentJobAndOutputLs = []
				# for each pair, add an alignment job
				for filePair in FilePairLs:
					newFilePair = self.registerFileToWorkflow(filePair, workflow, dataDir=dataDir)
					alignmentJob, alignmentOutput = self.addAlignmentJob(workflow, newFilePair, dataDir=dataDir, refFastaF=refFastaF, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, refIndexJob=refIndexJob, mkdirJob=mkdirJob, \
						alignment_method=alignment_method, \
						outputDir=tmpOutputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						no_of_aln_threads=no_of_aln_threads)
					AlignmentJobAndOutputLs.append([alignmentJob, alignmentOutput])
					no_of_alignment_jobs += 1
				#finalBamFileName = os.path.join(dataDir, individual_alignment.path)
				finalBamFileName = os.path.basename(individual_alignment.path)	#relative path in the scratch
				finalBamFile = File(finalBamFileName)
				self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, finalBamFile=finalBamFile, \
									samtools=samtools, mv=mv, namespace=namespace, version=version, \
									stageOutFinalOutput=self.stageOutFinalOutput)
					
		sys.stderr.write("%s alignment jobs.\n"%(no_of_alignment_jobs))
	
	def registerFileToWorkflow(self, filePair, workflow, dataDir=None):
		'''
		2011-8-30
		'''
		
		newFilePair = []
		for fileRecord in filePair:
			relativePath = fileRecord[0]
			fastqF = File(relativePath)
			fastqF.addPFN(PFN("file://" + os.path.join(dataDir, relativePath), self.input_site_handler))
			workflow.addFile(fastqF)
			newFileRecord = [fastqF] + fileRecord[1:]
			newFilePair.append(newFileRecord)
		return newFilePair
	
	def addRefIndexJob(self, workflow, refFastaF=None, refSequence=None, bwa=None,\
					namespace='workflow', version='1.0'):
		"""
		2011-8-28
		"""
		bwa_index_job = Job(namespace=namespace, name=bwa.name, version=version)
		if refSequence.base_count is None:	#default	#maybe add a base count job here
			index_algorithm = 'bwtsw'
		elif refSequence.base_count<500000000:	#500 million
			index_algorithm = "is"
		else:
			index_algorithm = "bwtsw"
		bwa_index_job.addArguments("index", "-a", index_algorithm, refFastaF)
		bwa_index_job.uses(refFastaF, transfer=True, register=False, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
		workflow.addJob(bwa_index_job)
		return bwa_index_job
	
	def addMKDIRJob(self, workflow, mkdirWrap=None, dirName=None, namespace='workflow', version='1.0'):
		"""
		"""
		mkdirJob = Job(namespace=namespace, name=mkdirWrap.name, version=version)
		mkdirJob.addArguments(dirName)
		workflow.addJob(mkdirJob)
		return mkdirJob
	
	def addAlignmentJobPipe(self, workflow, filePair, dataDir=None, refFastaF=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3):
		"""
		2011-9-7
			version that uses shell(bash) pipes to skip files.
		"""
		#in MB, 2.5GB for one  aln, 5G for sampe (input is 4.5G gzipped fastq, versus 480K contigs (total size~3G)), 2G for sort
		# paired end
		pe_job_max_memory = 11000
		#pe_job_max_memory = 7000	#lie about the usage. on usc cluster, this is translated into min memory requirement.
		## single end
		se_job_max_memory = 9000
		#se_job_max_memory = 6500	#lie about the usage. on usc cluster, this is translated into min memory requirement.
		if len(filePair)==1:	#single end
			fileRecord = filePair[0]
			fastqF, format, sequence_type, sequencer = fileRecord[:4]
			relativePath = fastqF.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
			sortBamF = File('%s.bam'%(bam_output_fname_prefix))
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				alignmentJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				alignmentJob = Job(namespace=namespace, name=LongSEAlignmentByBWA.name, version=version)
			
			alignmentJob.addArguments(refFastaF, fastqF, bam_output_fname_prefix)
			alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
			alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%se_job_max_memory))
			alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%se_job_max_memory))
			#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
		
		elif len(filePair)==2:	#paired end
			fastqF1, format, sequence_type = filePair[0][:3]
			fastqF2, format, sequence_type = filePair[1][:3]
			relativePath = fastqF1.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			fname_prefix = fname_prefix[:-2]
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
			sortBamF = File('%s.bam'%(bam_output_fname_prefix))
			alignmentJob = Job(namespace=namespace, name=PEAlignmentByBWA.name, version=version)
			alignmentJob.addArguments(refFastaF, fastqF1, fastqF2, bam_output_fname_prefix)
			alignmentJob.uses(fastqF1, transfer=True, register=False, link=Link.INPUT)
			alignmentJob.uses(fastqF2, transfer=True, register=False, link=Link.INPUT)
			
			alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%pe_job_max_memory))
			alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%pe_job_max_memory))
			#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
			
		alignmentJob.uses(sortBamF, transfer=False, register=False, link=Link.OUTPUT)
		workflow.addJob(alignmentJob)
		if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
			alignmentJob.uses(refFastaF, transfer=True, register=False, link=Link.INPUT)
		if mkdirJob:
			workflow.addDependency(parent=mkdirJob, child=alignmentJob)
		return alignmentJob, sortBamF
	
	def addAlignmentJob(self, workflow, filePair, dataDir=None, refFastaF=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3):
		"""
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 3800	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		if len(filePair)==1:	#single end
			fileRecord = filePair[0]
			fastqF, format, sequence_type, sequencer = fileRecord[:4]
			relativePath = fastqF.name
			
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				alignmentJob.addArguments(alignment_method.command, additionalArguments,"-t %s"%no_of_aln_threads, \
										"-f", saiOutput, refFastaF, \
										fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, register=False, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=False, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				
				bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
				sortBamF = File('%s.bam'%(bam_output_fname_prefix))
				sai2samJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
				
				
				sai2samJob.addArguments(refFastaF, saiOutput, fastqF, bam_output_fname_prefix)
				sai2samJob.uses(saiOutput, transfer=False, register=False, link=Link.INPUT)
				sai2samJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%samse_job_max_memory))
				sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%samse_job_max_memory))
				
				workflow.addJob(sai2samJob)
				workflow.addDependency(parent=alignmentJob, child=sai2samJob)
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
				sortBamF = File('%s.bam'%(bam_output_fname_prefix))
				alignmentJob = Job(namespace=namespace, name=LongSEAlignmentByBWA.name, version=version)
				
				alignmentJob.addArguments(refFastaF, fastqF, bam_output_fname_prefix, repr(no_of_aln_threads))
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
			
			if refIndexJob:
				workflow.addDependency(parent=refIndexJob, child=alignmentJob)
			if mkdirJob:
				workflow.addDependency(parent=mkdirJob, child=alignmentJob)
			
		elif len(filePair)==2:	#paired end
			fastqF1, format, sequence_type = filePair[0][:3]
			#fastqF2, format, sequence_type = filePair[1][:3]
			relativePath = fastqF1.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			fname_prefix = fname_prefix[:-2]
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
			sortBamF = File('%s.bam'%(bam_output_fname_prefix))
			sai2samJob = Job(namespace=namespace, name=PEAlignmentByBWA.name, version=version)
			sai2samJob.addArguments(refFastaF)
			sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%sampe_job_max_memory))
			sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%sampe_job_max_memory))
			#sai2samJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
			if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
				sai2samJob.uses(refFastaF, register=False, link=Link.INPUT)
			workflow.addJob(sai2samJob)
			for fileRecord in filePair:
				fastqF, format, sequence_type = fileRecord[:3]
				relativePath = fastqF.name
				#relativePath, format, sequence_type = fileRecord[:3]
				tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments, "-t %s"%(no_of_aln_threads), \
										"-f", saiOutput, refFastaF, fastqF)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, register=False, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=False, link=Link.OUTPUT)
				workflow.addJob(alignmentJob)
				
				if refIndexJob:
					workflow.addDependency(parent=refIndexJob, child=alignmentJob)
				if mkdirJob:
					workflow.addDependency(parent=mkdirJob, child=alignmentJob)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(saiOutput, transfer=False, register=False, link=Link.INPUT)
				workflow.addDependency(parent=alignmentJob, child=sai2samJob)
				
			
			#add a pair of fastq files to sampe in the end
			for fileRecord in filePair:
				fastqF = fileRecord[0]
				sai2samJob.addArguments(fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
			sai2samJob.addArguments(bam_output_fname_prefix)
		
		sai2samJob.uses(sortBamF, transfer=False, register=False, link=Link.OUTPUT)
		return sai2samJob, sortBamF
	
	def addAlignmentJobNoPipe(self, workflow, filePair, dataDir=None, refFastaF=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3):
		"""
		2011-8-28
		"""
		aln_job_max_memory = 3000	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		if len(filePair)==1:	#single end
			fileRecord = filePair[0]
			fastqF, format, sequence_type, sequencer = fileRecord[:4]
			relativePath = fastqF.name
			
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments,"-t %s"%no_of_aln_threads, \
										"-f", saiOutput, refFastaF, \
										fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, register=False, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=False, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				
				sai2samJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				sai2samJob.addArguments("samse", "-f", outputFname, refFastaF, saiOutput, fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					sai2samJob.uses(refFastaF, register=False, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=False, register=False, link=Link.INPUT)
				sai2samJob.uses(alignmentSamF, transfer=False, register=False, link=Link.OUTPUT)
				job_max_memory = 2500	#in MB
				sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				
				
				workflow.addJob(sai2samJob)
				workflow.addDependency(parent=alignmentJob, child=sai2samJob)
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				alignmentJob.addArguments('bwasw', additionalArguments, "-t %s"%(no_of_aln_threads), "-f", alignmentSamF, refFastaF, \
									fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, register=False, link=Link.INPUT)
				alignmentJob.uses(alignmentSamF, transfer=False, register=False, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
				
			if refIndexJob:
				workflow.addDependency(parent=refIndexJob, child=alignmentJob)
			if mkdirJob:
				workflow.addDependency(parent=mkdirJob, child=alignmentJob)
			
		elif len(filePair)==2:	#paired end
			### run sampe to combine two paired-end results into one sam file
			"bwa sampe -a 1000 -P hsref.fa ga1.sai ga2.sai ga1.fq ga2.fq | gzip > ga.sam.gz"
			sai2samJob = Job(namespace=namespace, name=bwa.name, version=version)
			# -P of sampe speeds things up but requires 4-5G memory for a human-size genome
			# "-a 1000" means maximum insert size is 1000.
			sai2samJob.addArguments('sampe', "-P")
			job_max_memory = 6000	#in MB
			sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
			sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			fileRecord, fileRecord2 = filePair[:2]
			fastqF, format, sequence_type = fileRecord[:3]
			relativePath = fastqF.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			# remove the _1 or _2 in the end of the two paired-end filenames.
			fname_prefix = fname_prefix[:-2]
			samOutputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
			alignmentSamF = File(samOutputFname)
			sai2samJob.addArguments("-f", alignmentSamF)
			sai2samJob.addArguments(refFastaF)
			if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
				sai2samJob.uses(refFastaF, register=False, link=Link.INPUT)
			sai2samJob.uses(alignmentSamF, transfer=False, register=False, link=Link.OUTPUT)
			workflow.addJob(sai2samJob)
			for fileRecord in filePair:
				fastqF, format, sequence_type = fileRecord[:3]
				relativePath = fastqF.name
				#relativePath, format, sequence_type = fileRecord[:3]
				tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments, "-t %s"%(no_of_aln_threads), \
										"-f", saiOutput, refFastaF, fastqF)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				alignmentJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, register=False, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=False, link=Link.OUTPUT)
				workflow.addJob(alignmentJob)
				
				if refIndexJob:
					workflow.addDependency(parent=refIndexJob, child=alignmentJob)
				if mkdirJob:
					workflow.addDependency(parent=mkdirJob, child=alignmentJob)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(fastqF, transfer=True, register=False, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=False, register=False, link=Link.INPUT)
				workflow.addDependency(parent=alignmentJob, child=sai2samJob)
				
			
			#add a pair of fastq files to sampe in the end
			for fileRecord in filePair:
				fastqF = fileRecord[0]
				sai2samJob.addArguments(fastqF)
		
		
		## convert sam into bam and remove unmapped reads (or queries)
		sam_convert_job = Job(namespace=namespace, name=samtools.name, version=version)
		bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fname_prefix)))
		sam_convert_job.addArguments('view',  '-F', '4', '-bSh', '-o', bamOutputF, alignmentSamF)
		sam_convert_job.uses(alignmentSamF, transfer=False, register=False, link=Link.INPUT)
		sam_convert_job.uses(bamOutputF, transfer=False, register=False, link=Link.OUTPUT)
		workflow.addJob(sam_convert_job)
		workflow.addDependency(parent=sai2samJob, child=sam_convert_job)
		
		"""
		# 2010-2-4
			sort it so that it could be used for merge
		"""
		bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
		sam_sort_job = Job(namespace=namespace, name=samtools.name, version=version)
		sam_sort_job.addArguments('sort', '-m', '2000000000', bamOutputF, bam_output_fname_prefix)	#maxMemory is down to 2G
		sortBamF = File('%s.bam'%(bam_output_fname_prefix))
		sam_sort_job.uses(bamOutputF, transfer=False, register=False, link=Link.INPUT)
		sam_sort_job.uses(sortBamF, transfer=False, register=False, link=Link.OUTPUT)
		job_max_memory = 2000	#in MB
		sam_sort_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
		sam_sort_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
		workflow.addJob(sam_sort_job)
		workflow.addDependency(parent=sam_convert_job, child=sam_sort_job)
		
		return sam_sort_job, sortBamF
	
	def addAlignmentMergeJob(self, workflow, AlignmentJobAndOutputLs=[], finalBamFile=None,samtools=None,\
					mv=None, namespace='workflow', version='1.0', stageOutFinalOutput=False):
		"""
		2011-9-4
			add argument stageOutFinalOutput, default=False, which means leave the output files where they are
		2011-8-28
			merge alignment
			index it
		"""
		if len(AlignmentJobAndOutputLs)>1:
			samtools_merge_job = Job(namespace=namespace, name=samtools.name, version=version)
			samtools_merge_job.addArguments("merge", finalBamFile)
			if stageOutFinalOutput:
				#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
				samtools_merge_job.uses(finalBamFile, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(samtools_merge_job)
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				alignmentJob, alignmentOutput = AlignmentJobAndOutput[:2]
				samtools_merge_job.addArguments(alignmentOutput)
				samtools_merge_job.uses(alignmentOutput, transfer=False, register=False, link=Link.INPUT)
				workflow.addDependency(parent=alignmentJob, child=samtools_merge_job)
		else:	#one input file, no samtools merge. use "mv" to rename it instead
			alignmentJob, alignmentOutput = AlignmentJobAndOutputLs[0][:2]
			samtools_merge_job = Job(namespace=namespace, name=mv.name, version=version)
			samtools_merge_job.addArguments(alignmentOutput, finalBamFile)
			workflow.addDependency(parent=alignmentJob, child=samtools_merge_job)
			samtools_merge_job.uses(alignmentOutput, transfer=False, register=False, link=Link.INPUT)
			if stageOutFinalOutput:
				samtools_merge_job.uses(finalBamFile, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(samtools_merge_job)
		
		# add the index job on the merged bam file
		samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
		samtools_index_job.addArguments("index", finalBamFile)
		if stageOutFinalOutput:
			samtools_index_job.uses(finalBamFile, transfer=True, register=False, link=Link.INPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
			bai_output = File('%s.bai'%finalBamFile)
			samtools_index_job.uses(bai_output, transfer=True, register=False, link=Link.OUTPUT)
		workflow.addJob(samtools_index_job)
		workflow.addDependency(parent=samtools_merge_job, child=samtools_index_job)
	
	
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
		session = db_vervet.session
		session.begin()
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		# Create a abstract dag
		workflow = ADAG("ShortRead2AlignmentPipeline")
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		
		#add the MergeSamFiles.jar file into workflow
		abs_path = os.path.join(self.picard_path, 'MergeSamFiles.jar')
		mergeSamFilesJar = File(abs_path)
		mergeSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(mergeSamFilesJar)
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mkdirWrap.sh"), site_handler))
		workflow.addExecutable(mkdirWrap)
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		workflow.addExecutable(mv)
		
		bwa = Executable(namespace=namespace, name="bwa", version=version, os=operatingSystem, arch=architecture, installed=True)
		bwa.addPFN(PFN("file://" + self.bwa_path, site_handler))
		workflow.addExecutable(bwa)
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
		
		
		PEAlignmentByBWA = Executable(namespace=namespace, name="PEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		PEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "PEAlignmentByBWA.sh"), site_handler))
		workflow.addExecutable(PEAlignmentByBWA)
		
		ShortSEAlignmentByBWA = Executable(namespace=namespace, name="ShortSEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		ShortSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "ShortSEAlignmentByBWA.sh"), site_handler))
		workflow.addExecutable(ShortSEAlignmentByBWA)
		
		LongSEAlignmentByBWA = Executable(namespace=namespace, name="LongSEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		LongSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "LongSEAlignmentByBWA.sh"), site_handler))
		workflow.addExecutable(LongSEAlignmentByBWA)
		
		#must use db_vervet.data_dir.
		# If self.dataDir differs from db_vervet.data_dir, this program (must be run on submission host) won't find files.
		individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, dataDir=self.localDataDir)
		refSequence = VervetDB.IndividualSequence.get(self.aln_ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)

		if self.needRefIndexJob:
			refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path
			# Add it into replica only when needed.
			refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
			workflow.addFile(refFastaF)
			# If it's not needed, assume the index is done and all relevant files are in absolute path.
			# and no replica transfer
			refIndexJob = self.addRefIndexJob(workflow, refFastaF=refFastaF, refSequence=refSequence, bwa=bwa, \
											namespace=namespace, version=version)
		else:
			refFastaF = File(refFastaFname)	#absolute path if no index job and no replica transfer
			refIndexJob = None
		
		self.addAllAlignmentJobs(db_vervet, individualSequenceID2FilePairLs=individualSequenceID2FilePairLs, dataDir=self.dataDir,\
					refSequence=refSequence, refFastaF=refFastaF, refIndexJob=refIndexJob,
					workflow=workflow, bwa=bwa, additionalArguments=self.additionalArguments, samtools=samtools, \
					mkdirWrap=mkdirWrap, mv=mv, alignment_method_name=self.alignment_method_name, alignment_format='bam',\
					namespace=namespace, version=version, stageOutFinalOutput=self.stageOutFinalOutput,\
					PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
					LongSEAlignmentByBWA=LongSEAlignmentByBWA, no_of_aln_threads=self.no_of_aln_threads)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = ShortRead2AlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
