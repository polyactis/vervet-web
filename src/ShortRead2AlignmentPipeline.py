#!/usr/bin/env python
"""
Examples:
	# 2011-8-30 workflow on condor
	%s -i 165-167 -o ShortRead2AlignmentPipeline_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool -n1 -z dl324b-1.cmb.usc.edu -c
	
	# 2011-8-30 a workflow with 454 long-read and short-read PE. need a ref index job (-n1). 
	%s -i 165-167 -o ShortRead2AlignmentPipeline_isq_id_165_167_vs_9.xml -u yh -a 9
		-e /u/home/eeskin/polyacti -l hoffman2 -t /u/home/eeskin/polyacti/NetworkData/vervet/db -n1 -z dl324b-1.cmb.usc.edu -c
		-f /work/
	
	# 2011-8-30 output a workflow to run alignments on hoffman2's condor pool (-D changes localDataDir. -t changes dataDir.)
	# 2012.3.20 use /work/ or /u/scratch/p/polyacti/tmp as TMP_DIR for MarkDuplicates.jar (/tmp is too small for 30X genome)
	# 2012.5.4 cluster 10 alignment jobs (before merging) as a unit (-w 10), skip done alignment (-K)
	%s  -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-l hcondor -j hcondor 
		-z localhost -u yh -c
		-i 631-700 -o workflow/ShortRead2Alignment_Isq_631-700_vs_524_hcondor.xml  -a 524 
		-f /work/ -e /u/home/eeskin/polyacti  -w 10 -K
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files. no ref index job (-n0).
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
		-j local -l uschpc -n1 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10 -p secret  -c

	# 2011-8-30 a workflow to run on uschpc, NO ref index job (-n0), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler and "-t ..." to enable symlink
	%s -i 391-397,456,473,493
		-o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j uschpc -l uschpc -n0 -p secret -c -m 4
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		-t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -J /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
	
	# 2011-11-16 a workflow to run on uschpc, NO ref index job (-n0), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler. this will stage in all input and output (-O).
	%s -i 391-397,456,473,493
		-o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_local2usc.xml -u yh -a 9
		-j local -l uschpc -n0 -p secret -c -m 4
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		-J /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
	
	
	#2011-9-13 no ref index job, staging input files from localhost to uschpc, stage output files back to localhost
	# modify the refFastaF's path in xml manually
	%s -i 1-3 -o ShortRead2AlignmentPipeline_1_3_vs_524_local2uschpc.xml -u yh -a 524
		-j local -l uschpc -n0 -p secret -c -m 4 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
		-t /Network/Data/vervet/db/
	
	# 2011-8-31 output the same workflow above but for condorpool
	%s -i 391-397,456,473,493, -o ShortRead2AlignmentPipeline_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
		-u yh -a 9 -j condorpool -l condorpool -n0 -z 10.8.0.10  -p secret  -c
	
	# 2012-4-5 new alignment method, stampy (-g)
	%s -i 167,176,178,182,183,207-211,391-397,456,473,493
		-o ShortRead2AlignmentPipeline_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
		-u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  -c -g stampy
	
Description:
	2012.5.3
		a program which generates a pegasus workflow dag (xml file) which does the alignment for all available sequences.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
		It will also stage out every output file.
		Be careful about -R, only toggle it if you know every input individual_sequence_file is not empty.
			Empty read files would fail alignment jobs and thus no final alignment for a few indivdiuals.
		Use "-w 10" to cluster the alignment jobs if the input read file is small enough (~1Million reads for bwa, ~300K for stampy). 
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow


class ShortRead2AlignmentPipeline(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('noCheckEmptyReadFile', 0, int):[0, 'R', 0, "toggle to not check whether each read file is empty (if empty, exclude it). \
							If IndividualSequenceFile.read_count is null, it'll try to count them on the fly and take a lot of time.\
							however, only toggle it if you know every input individual_sequence_file is not empty. empty read file fails alignment jobs."],\
						('ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						('additionalArguments', 0, ): ["-q 20", '', 1, 'a string of additional arguments passed to aln, not bwasw, add double quote if space'],\
						("bwa_path", 1, ): ["%s/bin/bwa", '', 1, 'bwa binary'],\
						("stampy_path", 1, ): ["%s/bin/stampy.py", '', 1, 'path to stampy.py'],\
						("alignment_method_name", 1, ): ["bwa-short-read", 'g', 1, 'alignment_method.short_name from db.\
								used only when unable to guess based on individual_sequence.sequencer and individual_sequence.sequence_type'],\
						("needRefIndexJob", 0, int): [0, 'n', 1, 'need to add a reference index job by bwa?'],\
						('no_of_aln_threads', 1, int): [1, 'm', 1, 'number of threads during alignment'],\
						('cluster_size_for_aln_jobs', 1, int): [1, 'w', 1, 'job cluster size for bwa/PEAlignmentByBWA/LongSEAlignmentByBWA/addOrReplaceReadGroupsJava/SortSamFilesJava/samtools jobs'],\
						('notStageOutFinalOutput', 0, int):[0, 'O', 0, 'toggle to not stage out final output (bam + bam.bai)'],\
						('skipDoneAlignment', 0, int):[0, 'K', 0, 'skip alignment whose path is a valid file'],\
						("tmpDir", 1, ): ["/tmp/", 'f', 1, 'for MarkDuplicates.jar, default is /tmp/ but sometimes it is too small'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						})

	def __init__(self,  **keywords):
		"""
		2012.3.29
			default to stage out final output.
			Argument stageOutFinalOutput morphs into notStageOutFinalOutput.
		2011-7-11
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		
		self.bwa_path =  self.insertHomePath(self.bwa_path, self.home_path)
		self.stampy_path =  self.insertHomePath(self.stampy_path, self.home_path)
		
		if self.notStageOutFinalOutput:
			self.stageOutFinalOutput = False
		else:
			self.stageOutFinalOutput = True
		
		if self.noCheckEmptyReadFile:
			self.ignoreEmptyReadFile = False
		else:
			self.ignoreEmptyReadFile = True
	
	def addAllAlignmentJobs(self, db_vervet, individualSequenceID2FilePairLs=None, \
						dataDir=None, \
						isq_id2LibrarySplitOrder2FileLs=None,\
					refSequence=None, refFastaFList=None, refIndexJob=None,
					workflow=None, bwa=None, additionalArguments=None, samtools=None, mkdirWrap=None, mv=None,\
					java=None, mergeSamFilesJava=None, mergeSamFilesJar=None, \
					MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir='/tmp',\
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					SortSamFilesJava=None, SortSamFilesJar=None, \
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					alignment_method_name='bwa-short-read', alignment_format='bam',\
					namespace='workflow', version='1.0', stageOutFinalOutput=False,\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					no_of_aln_threads=3, stampyExecutable=None, skipDoneAlignment=False):
		"""
		2012.4.20
			bugfix, pass alignment_method.short_name instead of alignment_method_name to db_vervet.getAlignment()
				because alignment_method might be changed according to sequencer regardless of alignment_method_name.
		2012.4.12
			add argument skipDoneAlignment
		2012.4.5
			fetch the alignment_method directly based on the alignment_method_name, except for 454 sequences.
				officially merges "bwa-short-read-SR" (single-read) into "bwa-short-read"
		2012.2.24
			pass dataDir to db_vervet.getAlignment()
			add stampy part
		2011-9-15
			adjust alignment_method_name according to individual_sequence.sequencer and individual_sequence.sequence_type
			only when this is not possible, value of argument alignment_method_name is used.
		2011-9-14
			give different names to different java jobs according to jars
		2011-8-28
		"""
		sys.stderr.write("Adding alignment jobs for %s individual sequences ..."%(len(isq_id2LibrarySplitOrder2FileLs)))
		
		no_of_alignment_jobs = 0
		no_of_merging_jobs = 0
		for ind_seq_id, LibrarySplitOrder2FileLs in isq_id2LibrarySplitOrder2FileLs.iteritems():
		#for ind_seq_id, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if individual_sequence is not None and individual_sequence.format=='fastq':
				
				AlignmentJobAndOutputLs = []
				#2012.3.29 this folder will store the alignment output by the alignment jbos
				tmpOutputDir = os.path.basename(individual_sequence.path)
				# add a mkdir job
				mkdirJob = None
				
				for key, fileObjLs in LibrarySplitOrder2FileLs.iteritems():
					library, split_order = key[:2]
					# get or add an alignment
					if individual_sequence.sequencer=='454':
						alignment_method = db_vervet.getAlignmentMethod("bwa-long-read")
					else:
						alignment_method = db_vervet.getAlignmentMethod(alignment_method_name)
					"""
					#2012.4.5 have all this commented out
					elif individual_sequence.sequencer=='GA':
						if individual_sequence.sequence_type=='SR':	#single-end
							alignment_method = db_vervet.getAlignmentMethod("bwa-short-read-SR")
						else:	#default is PE
							alignment_method = db_vervet.getAlignmentMethod("bwa-short-read")
					"""
					individual_alignment = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
											individual_sequence_id=ind_seq_id,\
								path_to_original_alignment=None, sequencer=individual_sequence.sequencer, \
								sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
								ref_individual_sequence_id=refSequence.id, \
								alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
								individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
								dataDir=dataDir)	#read-group addition is part of pipeline
					if not individual_alignment.path:
						individual_alignment.path = individual_alignment.constructRelativePath()
						session.add(individual_alignment)
						session.flush()
					
					alignmentAbsPath= os.path.join(dataDir, individual_alignment.path)
					#2012.3.29	check if the alignment exists or not. if it already exists, no alignment jobs.
					if skipDoneAlignment and os.path.isfile(alignmentAbsPath):
						break
					if mkdirJob is None:	#now it's time to add the mkdirJob
						# add a mkdir job
						mkdirJob = self.addMKDIRJob(workflow, mkdirWrap=mkdirWrap, dirName=tmpOutputDir, \
											namespace=namespace, version=version)
					#newFilePair = self.registerFileToWorkflow(filePair, workflow, dataDir=dataDir)
					newFileObjLs = self.registerISQFileObjLsToWorkflow(fileObjLs, workflow)
					alignmentJob, alignmentOutput = self.addAlignmentJob(workflow, newFileObjLs, individual_alignment=individual_alignment, \
						dataDir=dataDir, refFastaFList=refFastaFList, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, refIndexJob=refIndexJob, mkdirJob=mkdirJob, \
						alignment_method=alignment_method, \
						outputDir=tmpOutputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
						LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
						no_of_aln_threads=no_of_aln_threads, stampyExecutable=stampyExecutable)
					AlignmentJobAndOutputLs.append([alignmentJob, alignmentOutput])
					no_of_alignment_jobs += 1
				if len(AlignmentJobAndOutputLs)>0:	#2012.3.29	merge alignment output only when there is something to merge!
					fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(individual_alignment.path))[0]
					
					mergedBamFile = File('%s_merged.bam'%(fname_prefix))
					alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
										outputBamFile=mergedBamFile, \
										samtools=samtools, java=java, \
										mergeSamFilesJava=mergeSamFilesJava, mergeSamFilesJar=mergeSamFilesJar, \
										BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
										mv=mv, namespace=namespace, version=version, \
										stageOutFinalOutput=False)
					#finalBamFileName = os.path.join(dataDir, individual_alignment.path)
					finalBamFileName = os.path.basename(individual_alignment.path)	#relative path in the scratch
					finalBamFile = File(finalBamFileName)
					
					markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[alignmentMergeJob, bamIndexJob], inputBamF=mergedBamFile, \
							inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
							MarkDuplicatesJava=MarkDuplicatesJava, MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
							BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
							namespace=namespace, version=version, stageOutFinalOutput=self.stageOutFinalOutput)
					no_of_merging_jobs += 1
		sys.stderr.write("%s alignment jobs; %s merge alignment jobs.\n"%(no_of_alignment_jobs, no_of_merging_jobs))
	
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
	
	def registerISQFileObjLsToWorkflow(self, fileObjLs, workflow, dataDir=None):
		'''
		2012-2.24
			similar to registerFileToWorkflow but for a different input
		'''
		
		newFilePair = []
		for fileObject in fileObjLs:
			relativePath = fileObject.db_entry.path
			fastqF = File(relativePath)
			fastqF.addPFN(PFN("file://" + fileObject.path, self.input_site_handler))
			workflow.addFile(fastqF)
			fileObject.fastqF = fastqF
			newFilePair.append(fileObject)
		return newFilePair
	
	def addStampyGenomeIndexHashJob(self, workflow, executable=None, refFastaFList=None, \
						parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
						extraDependentInputLs=[], \
						transferOutput=False, **keywords):
		"""
		2012.2.23
		"""
		stampyGenomeIndexJob= Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		refFastaFile = refFastaFList[0]
		
		genomeIndexFile = File('%s.stidx'%(refFastaFile.name))
		stampyGenomeIndexJob.addArguments("-G ", refFastaFile, refFastaFile)
		stampyGenomeIndexJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		#genomeIndexFile is output of this stampyGenomeIndexJob but ignore it as output otherwise it'll get auto-cleaned by pegasus
		#stampyGenomeIndexJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(stampyGenomeIndexJob)
		for input in extraDependentInputLs:
			stampyGenomeIndexJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=stampyGenomeIndexJob)
		yh_pegasus.setJobProperRequirement(stampyGenomeIndexJob, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		
		stampyGenomeHashJob = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		genomeHashFile = File('%s.sthash'%(refFastaFile.name))
		stampyGenomeHashJob.addArguments("-g", refFastaFile, "-H", refFastaFile)
		
		stampyGenomeHashJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		#genomeIndexFile is input to this stampyGenomeHashJob but mark it as output otherwise it'll get auto-cleaned by pegasus
		stampyGenomeHashJob.uses(genomeIndexFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		stampyGenomeHashJob.uses(genomeHashFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		stampyGenomeHashJob.outputLs = [genomeIndexFile, genomeHashFile]
		
		workflow.addJob(stampyGenomeHashJob)
		workflow.depends(parent=stampyGenomeIndexJob, child=stampyGenomeHashJob)
		yh_pegasus.setJobProperRequirement(stampyGenomeHashJob, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		return stampyGenomeHashJob
	
	def addRefIndexJob(self, workflow, refFastaFList=None, refSequence=None, bwa=None,\
					namespace='workflow', version='1.0', \
					bwaIndexFileSuffixLs = ['amb', 'ann', 'bwt', 'pac', 'sa', 'rbwt', 'rpac', 'rsa',]):
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
		for refFastaFile in refFastaFList:
			bwa_index_job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		bwa_index_job.outputLs = []
		for suffix in bwaIndexFileSuffixLs:
			file = File("%s.%s"%(refFastaFile.name, suffix))
			bwa_index_job.outputLs.append(file)
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
					no_of_aln_threads=3, **keywords):
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
			alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
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
			alignmentJob.uses(fastqF1, transfer=True, register=True, link=Link.INPUT)
			alignmentJob.uses(fastqF2, transfer=True, register=True, link=Link.INPUT)
			
			alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%pe_job_max_memory))
			alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%pe_job_max_memory))
			#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
			
		alignmentJob.uses(sortBamF, transfer=False, register=True, link=Link.OUTPUT)
		workflow.addJob(alignmentJob)
		if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
			alignmentJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
		if mkdirJob:
			workflow.depends(parent=mkdirJob, child=alignmentJob)
		return alignmentJob, sortBamF
	
	def addRefIndexJobAndItsOutputAsParent(self, workflow, refIndexJob, childJob=None):
		"""
		2012.2.24
		"""
		workflow.depends(parent=refIndexJob, child=childJob)
		for output in refIndexJob.outputLs:
			childJob.uses(output, transfer=True, register=True, link=Link.INPUT)
	
	def addStampyAlignmentJob(self, workflow, fileObjLs, individual_alignment=None, \
					dataDir=None, refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, \
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, stampyExecutable=None, stageOutFinalOutput=False, **keywords):
		"""
		2014.4.4
			update for stampy 1.0.17
				which added option --gatkcigarworkaround removing adjacent I/D events from CIGAR strings, which trips up GATK
		2012.3.5
			added "--overwrite" to stampy.py to overwrite partial output alignment file
			added two stampy options
				--baq                              (SAM format) Compute base-alignment quality (BAQ; BQ tag)
				--alignquals                       (SAM format) Compute posterior alignment probabilities (YQ tag)
		2012.2.26
			handle read files with quality_score_format="Illumina" (post 1.3 solexa).
			
		2012.2.24
			alignment job for stampy
			no_of_aln_threads is only for bwa.
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 6000	#in MB, bwa needs 3G. stampy needs 3G.
		#bwa: memory 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		aln_job_max_walltime= 4800	#80 hours, in minutes
		aln_job_max_walltime = 1320	#22 hours, because all reads are stored in chunks of 5-million-read files
		
		javaMemRequirement = "-Xms128m -Xmx%sm"%addRGJob_max_memory
		refFastaFile = refFastaFList[0]
		firstFileObject = fileObjLs[0]
		fastqF = firstFileObject.fastqF
		relativePath = fastqF.name
		fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
		outputSamFile = File('%s.sam'%(os.path.join(outputDir, fname_prefix)))
		
		alignmentJob = Job(namespace=workflow.namespace, name=stampyExecutable.name, version=workflow.version)
		# make sure to use ', rather than ", to wrap the bwaoptions. double-quote(") would disappear during xml translation.
		alignmentJob.addArguments(" --bwa=%s "%(yh_pegasus.getAbsPathOutOfExecutable(bwa)), \
					"--bwaoptions='%s -t%s %s' "%(additionalArguments, no_of_aln_threads, refFastaFile.name),  \
					"-g", refFastaFile, "-h", refFastaFile, "-o", outputSamFile, '--overwrite',\
					'--baq', '--alignquals', '--gatkcigarworkaround')
		#Added option --gatkcigarworkaround removing adjacent I/D events from CIGAR strings, which trips up GATK
		if firstFileObject.db_entry.quality_score_format=='Illumina':
			alignmentJob.addArguments("--solexa")
		alignmentJob.uses(outputSamFile, transfer=stageOutFinalOutput, register=True, link=Link.OUTPUT)
		alignmentJob.output = outputSamFile
		alignmentJob.fname_prefix = fname_prefix
		
		for refFastaFile in refFastaFList:
			alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
										no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
		workflow.addJob(alignmentJob)
		#add fastq files after "-M"
		alignmentJob.addArguments(" -M ")	# -M FILE[,FILE] Map fastq file(s).  Use FILE.recaldata for recalibration if available
		for fileObject in fileObjLs:
			fastqF = fileObject.fastqF
			relativePath = fastqF.name
			alignmentJob.addArguments(fastqF)
			alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
		if refIndexJob:
			self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=alignmentJob)
		if mkdirJob:
			workflow.depends(parent=mkdirJob, child=alignmentJob)
		return alignmentJob
	
	def addBWAAlignmentJob(self, workflow, fileObjectLs, individual_alignment=None, \
					dataDir=None, refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, **keywords):
		"""
		2012.4.5
			handle read files with quality_score_format="Illumina" (post 1.3 solexa).
			input fileObjectLs is different now.
		2012.2.23
			split out of addAlignmentJob()
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 2600	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		
		aln_job_max_walltime= 4800	#80 hours, in minutes
		aln_job_max_walltime = 1020	#22 hours, because all reads are stored in chunks of 5-million-read files
		
		javaMemRequirement = "-Xms128m -Xmx%sm"%addRGJob_max_memory
		if len(fileObjectLs)==1:	#single end
			fileObject = fileObjectLs[0]
			fastqF = fileObject.fastqF
			relativePath = fastqF.name
			sequence_type = fileObject.db_entry.individual_sequence.sequence_type
			sequencer = fileObject.db_entry.individual_sequence.sequencer
			
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				alignmentJob.addArguments(alignment_method.command, additionalArguments,"-t %s"%no_of_aln_threads, \
										"-f", saiOutput)
				if fileObject.db_entry.quality_score_format=='Illumina':	#2012.4.5
					alignmentJob.addArguments("-I")
				alignmentJob.addArguments(refFastaFList[0], fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				for refFastaFile in refFastaFList:
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
												no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
				
				workflow.addJob(alignmentJob)
				
				alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
				sai2samJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
				
				sai2samJob.addArguments(refFastaFList[0], saiOutput, fastqF, alignmentSamF)
				for refFastaFile in refFastaFList:
					sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				yh_pegasus.setJobProperRequirement(sai2samJob, job_max_memory=samse_job_max_memory)
				workflow.addJob(sai2samJob)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
				alignmentJob = Job(namespace=namespace, name=LongSEAlignmentByBWA.name, version=version)
				
				alignmentJob.addArguments(refFastaFList[0], fastqF, alignmentSamF, repr(no_of_aln_threads))
				for refFastaFile in refFastaFList:
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=bwasw_job_max_memory, \
												no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
			
			if refIndexJob:
				self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=alignmentJob)
				#workflow.depends(parent=refIndexJob, child=alignmentJob)
			if mkdirJob:
				workflow.depends(parent=mkdirJob, child=alignmentJob)
			
		elif len(fileObjectLs)==2:	#paired end
			fileObject = fileObjectLs[0]
			fastqF1 = fileObject.fastqF
			relativePath = fastqF1.name
			sequence_type = fileObject.db_entry.individual_sequence.sequence_type
			sequencer = fileObject.db_entry.individual_sequence.sequencer
			
			#fastqF1, format, sequence_type = fileObjectLs[0][:3]
			#fastqF2, format, sequence_type = fileObjectLs[1][:3]
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			fname_prefix = fname_prefix[:-2]
			alignmentSamF = File('%s.sam.gz'%(os.path.join(outputDir, fname_prefix)))
			sai2samJob = Job(namespace=namespace, name=PEAlignmentByBWA.name, version=version)
			sai2samJob.addArguments(refFastaFList[0])
			yh_pegasus.setJobProperRequirement(sai2samJob, job_max_memory=sampe_job_max_memory)
			for refFastaFile in refFastaFList:
				sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
			workflow.addJob(sai2samJob)
			for fileObject in fileObjectLs:
				fastqF = fileObject.fastqF
				relativePath = fastqF.name
				sequence_type = fileObject.db_entry.individual_sequence.sequence_type
				sequencer = fileObject.db_entry.individual_sequence.sequencer
				
				#fastqF, format, sequence_type = fileObject[:3]
				
				#relativePath, format, sequence_type = fileObject[:3]
				tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments, "-t %s"%(no_of_aln_threads), \
										"-f", saiOutput)
				if fileObject.db_entry.quality_score_format=='Illumina':	#2012.4.5
					alignmentJob.addArguments("-I")
				alignmentJob.addArguments(refFastaFList[0], fastqF)
				yh_pegasus.setJobProperRequirement(alignmentJob, job_max_memory=aln_job_max_memory, \
												no_of_cpus=no_of_aln_threads, max_walltime=aln_job_max_walltime)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				for refFastaFile in refFastaFList:
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(alignmentJob)
				
				if refIndexJob:
					self.addRefIndexJobAndItsOutputAsParent(workflow, refIndexJob, childJob=alignmentJob)
					#workflow.depends(parent=refIndexJob, child=alignmentJob)
				if mkdirJob:
					workflow.depends(parent=mkdirJob, child=alignmentJob)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
			
			#add a pair of fastq files to sampe in the end
			for fileObject in fileObjectLs:
				fastqF = fileObject.fastqF
				sai2samJob.addArguments(fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
			sai2samJob.addArguments(alignmentSamF)
		sai2samJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
		sai2samJob.output = alignmentSamF
		sai2samJob.fname_prefix = fname_prefix
		return sai2samJob
	
	def addAlignmentJob(self, workflow, fileObjLs, individual_alignment=None, \
					dataDir=None, refFastaFList=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					java=None, SortSamFilesJava=None, SortSamFilesJar=None,\
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					no_of_aln_threads=3, stampyExecutable=None, **keywords):
		"""
		2012.4.5
			choose which alignment program to use based on alignment_method.short_name
		2012.2.23
			split the BWA alignment part to addBWAAlignmentJob()
		2011-9-13
			add argument java & SortSamFilesJar
		2011-9-9
			two steps:
				1. aln doesn't use pipe, outputs to sai files.
				2. sampe/samse, convert, sort => connected through pipe
		"""
		aln_job_max_memory = 2600	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		bwasw_job_max_memory = 3800	#in MB, same ref, bwasw needs more memory
		samse_job_max_memory = 4500	#in MB
		sampe_job_max_memory = 6000	#in MB
		addRGJob_max_memory = 2500	#in MB
		
		aln_job_max_walltime= 4800	#80 hours, in minutes
		
		javaMemRequirement = "-Xms128m -Xmx%sm"%addRGJob_max_memory
		
		
		if alignment_method.short_name=='stampy':
			alignmentJob = self.addStampyAlignmentJob(workflow, fileObjLs, individual_alignment=individual_alignment, \
						dataDir=dataDir, refFastaFList=refFastaFList, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, \
						refIndexJob=refIndexJob, mkdirJob=mkdirJob, \
						alignment_method=alignment_method, \
						outputDir=outputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
						LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
						no_of_aln_threads=no_of_aln_threads, stampyExecutable=stampyExecutable,\
						stageOutFinalOutput=False)
		elif alignment_method.short_name.find('bwa')==0:
			alignmentJob = self.addBWAAlignmentJob(workflow, fileObjLs, individual_alignment=individual_alignment, \
						dataDir=dataDir, refFastaFList=refFastaFList, bwa=bwa, \
						additionalArguments=additionalArguments, samtools=samtools, refIndexJob=refIndexJob, mkdirJob=mkdirJob, \
						alignment_method=alignment_method, \
						outputDir=outputDir, namespace=namespace, version=version,\
						PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
						LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
						java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
						addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
						no_of_aln_threads=no_of_aln_threads)
		else:
			sys.stderr.write("Alignment method %s is not supported.\n"%(alignment_method_name))
			sys.exit(3)
		fname_prefix = alignmentJob.fname_prefix
		
		## convert sam into bam and remove unmapped reads (or queries)
		sam_convert_job = Job(namespace=namespace, name=samtools.name, version=version)
		bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fname_prefix)))
		sam_convert_job.addArguments('view', '-bSh', '-o', bamOutputF, alignmentJob.output)
		#'-F', '4', (remove un-mapped reads) is removed
		sam_convert_job.uses(alignmentJob.output, transfer=True, register=True, link=Link.INPUT)
		sam_convert_job.uses(bamOutputF, transfer=False, register=True, link=Link.OUTPUT)
		workflow.addJob(sam_convert_job)
		workflow.depends(parent=alignmentJob, child=sam_convert_job)
		
		#2011-9-14 add/replace read group in the bam file
		# add RG to this bam
		alignment=individual_alignment
		sequencer = alignment.ind_sequence.sequencer
		read_group = '%s_%s_%s_%s_vs_%s'%(alignment.id, alignment.ind_seq_id, alignment.ind_sequence.individual.code, \
								sequencer, alignment.ref_ind_seq_id)
		if sequencer=='454':
			platform_id = 'LS454'
		elif sequencer=='GA':
			platform_id = 'ILLUMINA'
		else:
			platform_id = 'ILLUMINA'
		# the add-read-group job
		#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
		addRGJob = Job(namespace=namespace, name=addOrReplaceReadGroupsJava.name, version=version)
		outputRGBAM = File(os.path.join(outputDir, "%s.RG.bam"%(fname_prefix)))
		addRGJob.addArguments(javaMemRequirement, '-jar', addOrReplaceReadGroupsJar, \
							"INPUT=", bamOutputF,\
							'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
							'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
							'OUTPUT=', outputRGBAM, "VALIDATION_STRINGENCY=LENIENT")	#not including 'SORT_ORDER=coordinate'
							#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
		addRGJob.uses(bamOutputF, transfer=True, register=True, link=Link.INPUT)
		addRGJob.uses(outputRGBAM, transfer=False, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(addRGJob, job_max_memory=addRGJob_max_memory, no_of_cpus=None)
		workflow.addJob(addRGJob)
		workflow.depends(parent=sam_convert_job, child=addRGJob)
		
		
		"""
		# 2010-2-4
			sort it so that it could be used for merge
		"""
		bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
		sortBamF = File('%s.bam'%(bam_output_fname_prefix))
		sort_sam_job = Job(namespace=namespace, name=SortSamFilesJava.name, version=version)
		sort_sam_job.addArguments(javaMemRequirement, '-jar', SortSamFilesJar, "SORT_ORDER=coordinate", "I=", outputRGBAM, \
								"O=", sortBamF, "VALIDATION_STRINGENCY=LENIENT")
		sort_sam_job.uses(outputRGBAM, transfer=True, register=True, link=Link.INPUT)
		sort_sam_job.uses(sortBamF, transfer=False, register=True, link=Link.OUTPUT)
		sort_sam_job.output = sortBamF
		#job_max_memory = 2500	#in MB, 2.5g resident memory, 6458Mb Virtual in one example 
		yh_pegasus.setJobProperRequirement(sort_sam_job, job_max_memory=addRGJob_max_memory, no_of_cpus=None)
		workflow.addJob(sort_sam_job)
		workflow.depends(parent=addRGJob, child=sort_sam_job)
		
		return sort_sam_job, sortBamF
	
	def addAlignmentJobNoPipe(self, workflow, filePair, dataDir=None, refFastaF=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3, **keywords):
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
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				
				sai2samJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				sai2samJob.addArguments("samse", "-f", outputFname, refFastaF, saiOutput, fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					sai2samJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
				job_max_memory = 2500	#in MB
				sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				
				
				workflow.addJob(sai2samJob)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				alignmentJob.addArguments('bwasw', additionalArguments, "-t %s"%(no_of_aln_threads), "-f", alignmentSamF, refFastaF, \
									fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
				
			if refIndexJob:
				workflow.depends(parent=refIndexJob, child=alignmentJob)
			if mkdirJob:
				workflow.depends(parent=mkdirJob, child=alignmentJob)
			
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
			sai2samJob.addProfile(Profile(Namespace.CONDOR, key="request_memory", value="%s"%job_max_memory))
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
				sai2samJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
			sai2samJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
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
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"%no_of_aln_threads))
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(alignmentJob)
				
				if refIndexJob:
					workflow.depends(parent=refIndexJob, child=alignmentJob)
				if mkdirJob:
					workflow.depends(parent=mkdirJob, child=alignmentJob)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
				
			
			#add a pair of fastq files to sampe in the end
			for fileRecord in filePair:
				fastqF = fileRecord[0]
				sai2samJob.addArguments(fastqF)
		
		
		## convert sam into bam and remove unmapped reads (or queries)
		sam_convert_job = Job(namespace=namespace, name=samtools.name, version=version)
		bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fname_prefix)))
		sam_convert_job.addArguments('view',  '-F', '4', '-bSh', '-o', bamOutputF, alignmentSamF)
		sam_convert_job.uses(alignmentSamF, transfer=True, register=True, link=Link.INPUT)
		sam_convert_job.uses(bamOutputF, transfer=False, register=True, link=Link.OUTPUT)
		workflow.addJob(sam_convert_job)
		workflow.depends(parent=sai2samJob, child=sam_convert_job)
		
		"""
		# 2010-2-4
			sort it so that it could be used for merge
		"""
		bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
		sam_sort_job = Job(namespace=namespace, name=samtools.name, version=version)
		sam_sort_job.addArguments('sort', '-m', '2000000000', bamOutputF, bam_output_fname_prefix)	#maxMemory is down to 2G
		sortBamF = File('%s.bam'%(bam_output_fname_prefix))
		sam_sort_job.uses(bamOutputF, transfer=True, register=True, link=Link.INPUT)
		sam_sort_job.uses(sortBamF, transfer=False, register=True, link=Link.OUTPUT)
		job_max_memory = 2000	#in MB
		sam_sort_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
		sam_sort_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
		workflow.addJob(sam_sort_job)
		workflow.depends(parent=sam_convert_job, child=sam_sort_job)
		
		return sam_sort_job, sortBamF
	
	def addAlignmentMergeBySAMtoolsJob(self, workflow, AlignmentJobAndOutputLs=[], outputBamFile=None,samtools=None,\
					java=None, mergeSamFilesJava=None, mergeSamFilesJar=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					mv=None, stageOutFinalOutput=False, job_max_memory=2500, **keywords):
		"""
		2012-3.22
			samtools merge version of addAlignmentMergeJob(), which uses picard's MergeSamFiles.jar
			**untested**
		"""
		javaMaxMemory=2500
		if len(AlignmentJobAndOutputLs)>1:
			merge_sam_job = Job(namespace=workflow.namespace, name=samtools.name, version=workflow.version)
			merge_sam_job.addArguments('-f', outputBamFile)		#'overwrite the output BAM if exist'
			merge_sam_job.uses(outputBamFile, transfer=stageOutFinalOutput, register=True, link=Link.OUTPUT)
			yh_pegasus.setJobProperRequirement(merge_sam_job, job_max_memory=job_max_memory)
			workflow.addJob(merge_sam_job)
			for AlignmentJobAndOutput in AlignmentJobAndOutputLs:
				alignmentJob, alignmentOutput = AlignmentJobAndOutput[:2]
				merge_sam_job.addArguments(alignmentOutput)
				merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=merge_sam_job)
		else:	#one input file, no samtools merge. use "mv" to rename it instead
			alignmentJob, alignmentOutput = AlignmentJobAndOutputLs[0][:2]
			merge_sam_job = Job(namespace=workflow.namespace, name=mv.name, version=workflow.version)
			merge_sam_job.addArguments(alignmentOutput, outputBamFile)
			workflow.depends(parent=alignmentJob, child=merge_sam_job)
			merge_sam_job.uses(alignmentOutput, transfer=True, register=True, link=Link.INPUT)
			merge_sam_job.uses(outputBamFile, transfer=stageOutFinalOutput, register=True, link=Link.OUTPUT)
			workflow.addJob(merge_sam_job)
		
		# add the index job on the merged bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, \
						BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
						inputBamF=outputBamFile, parentJobLs=[merge_sam_job], \
						stageOutFinalOutput=stageOutFinalOutput, javaMaxMemory=javaMaxMemory)
		return merge_sam_job, bamIndexJob
	
	def addMarkDupJob(self, workflow, parentJobLs=[], inputBamF=None, inputBaiF=None, outputBamFile=None,\
					MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir="/Network/Data/vervet/vervetPipeline/tmp/",\
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', stageOutFinalOutput=True, no_of_cpus=1):
		"""
		2012.3.21
			improve it
			set no_of_cpus=1 (was 2) to avoid thread problem in some linux kernels.
		#2011-11-10 duplicate-marking job
		"""
		MarkDupJobMaxMemory=4000
		MarkDupJob = Job(namespace=getattr(workflow, 'namespace', namespace), name=MarkDuplicatesJava.name, \
						version=getattr(workflow, 'version', version))
		bamFnamePrefix = os.path.splitext(outputBamFile.name)[0]
		MarkDupOutputF = outputBamFile
		MarkDupOutputMetricF = '%s.metric'%(bamFnamePrefix)
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=MarkDupJobMaxMemory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		MarkDupJob.addArguments(javaMemRequirement, '-jar', MarkDuplicatesJar, "MAX_FILE_HANDLES=500",\
			"VALIDATION_STRINGENCY=LENIENT", "ASSUME_SORTED=true", "INPUT=", inputBamF, \
			'OUTPUT=', MarkDupOutputF, "M=", MarkDupOutputMetricF, "MAX_RECORDS_IN_RAM=500000",\
			"TMP_DIR=%s"%tmpDir)
		MarkDupJob.uses(inputBaiF, transfer=True, register=True, link=Link.INPUT)
		MarkDupJob.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
		MarkDupJob.output = MarkDupOutputF
		MarkDupJob.MarkDupOutputMetricF = MarkDupOutputMetricF
		if stageOutFinalOutput:
			MarkDupJob.uses(MarkDupOutputF, transfer=True, register=True, link=Link.OUTPUT)
			MarkDupJob.uses(MarkDupOutputMetricF, transfer=True, register=True, link=Link.OUTPUT)
		else:
			MarkDupJob.uses(MarkDupOutputF, transfer=False, register=True, link=Link.OUTPUT)
			MarkDupJob.uses(MarkDupOutputMetricF, transfer=False, register=True, link=Link.OUTPUT)
			#pass	#don't register the files so leave them there
		workflow.addJob(MarkDupJob)
		yh_pegasus.setJobProperRequirement(MarkDupJob, job_max_memory=job_max_memory, no_of_cpus=no_of_cpus)
		for parentJob in parentJobLs:
			if parentJob:
				workflow.depends(parent=parentJob, child=MarkDupJob)
		
		
		# add the index job on the bam file
		bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, \
								BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
								inputBamF=MarkDupOutputF,\
								parentJobLs=[MarkDupJob], namespace=namespace, version=version,\
								stageOutFinalOutput=stageOutFinalOutput)
		return MarkDupJob, bamIndexJob
		
	@classmethod
	def addSAMtoolsCalmdJob(cls, workflow, samtoolsCalmd=None, inputBamF=None, \
					refFastaFList=None, outputBamF=None, \
					parentJob=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', stageOutFinalOutput=True,\
					**keywords):
		"""
		2011-11-20
			run "samtools calmd" on the input bam and index the output bam
		"""
		job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
		job.addArguments(inputBamF, refFastaFList[0], outputBamF)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=1000)
		job.uses(inputBamF, transfer=True, register=True, link=Link.INPUT)
		for refFastaFile in refFastaFList:
			job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		if stageOutFinalOutput:
			job.uses(outputBamF, transfer=True, register=True, link=Link.OUTPUT)
		else:
			pass	#don't register the files so leave them there
		workflow.addJob(job)
		if parentJob:
			workflow.depends(parent=parentJob, child=job)
			
		# add the index job on the bam
		return self.addBAMIndexJob(BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					inputBamF=outputBamF,\
					parentJobLs=[job], namespace=namespace, version=version,\
					stageOutFinalOutput=stageOutFinalOutput)
	
	@classmethod
	def addBaseQualRecalibrationJob(cls, workflow, samtoolsCalmd=None, inputBamF=None, \
					refFastaFList=None, outputBamF=None, \
					parentJob=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', stageOutFinalOutput=True,\
					**keywords):
		"""
		2011-11-20
			create a sub-workflow
			
			for each 2million bp interval
				split a small bam out of input bam (add 1kb on both ends of interval)
				index the small bam
				run RealignerTargetCreator
				run IndelRealigner
				(extract the exact interval out of bam if 1kb is added to both ends of input interval)
			merge all intervals back into one bam
			index the merged bam
		"""
	
	@classmethod
	def addLocalRealignmentJob(cls, workflow, samtoolsCalmd=None, inputBamF=None, \
					refFastaFList=None, outputBamF=None, \
					parentJob=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', stageOutFinalOutput=True,\
					**keywords):
		"""
		2011-11-20
			create a sub-workflow to do local re-alignment for one inputBAM
			for each 2million bp interval
				split a small bam out of input bam (add 1kb on both ends of interval)
				index the small bam
				run RealignerTargetCreator
				run IndelRealigner
				(extract the exact interval out of bam if 1kb is added to both ends of input interval)
			merge all intervals back into one bam
			index the merged bam
			
		"""
	
	def registerCustomJars(self, workflow=None):
		"""
		2012.1.9
		"""
		site_handler = self.site_handler
		
		abs_path = os.path.join(self.picard_path, 'SortSam.jar')
		SortSamFilesJar = File(abs_path)
		SortSamFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(SortSamFilesJar)
		workflow.SortSamFilesJar = SortSamFilesJar
		
		abs_path = os.path.join(self.picard_path, 'SamFormatConverter.jar')
		SamFormatConverterJar = File(abs_path)
		SamFormatConverterJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(SamFormatConverterJar)
		workflow.SamFormatConverterJar = SamFormatConverterJar
		
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.5.3
			add clusters.size profile for alignment specific jobs (self.cluster_size_for_aln_jobs)
		2012.1.3
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		stampy = Executable(namespace=namespace, name="stampy", version=version, os=operatingSystem, \
						arch=architecture, installed=True)
		stampy.addPFN(PFN("file://" + self.stampy_path, site_handler))
		stampy.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(stampy)
		workflow.stampy = stampy
		
		workflow.IndexMergedBamIndexJava = workflow.BuildBamIndexFilesJava #no clustering for this as this runs on merged/dup-marked big bams
		#workflow.BuildBamIndexFilesJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addOrReplaceReadGroupsJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.samtools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		
		
		bwa = Executable(namespace=namespace, name="bwa", version=version, os=operatingSystem, arch=architecture, installed=True)
		bwa.addPFN(PFN("file://" + self.bwa_path, site_handler))
		bwa.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(bwa)
		workflow.bwa = bwa
		
		PEAlignmentByBWA = Executable(namespace=namespace, name="PEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		PEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/PEAlignmentByBWA.sh"), site_handler))
		PEAlignmentByBWA.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(PEAlignmentByBWA)
		workflow.PEAlignmentByBWA = PEAlignmentByBWA
		
		ShortSEAlignmentByBWA = Executable(namespace=namespace, name="ShortSEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		ShortSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/ShortSEAlignmentByBWA.sh"), site_handler))
		ShortSEAlignmentByBWA.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(ShortSEAlignmentByBWA)
		workflow.ShortSEAlignmentByBWA = ShortSEAlignmentByBWA
		
		LongSEAlignmentByBWA = Executable(namespace=namespace, name="LongSEAlignmentByBWA.sh", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		LongSEAlignmentByBWA.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/LongSEAlignmentByBWA.sh"), site_handler))
		LongSEAlignmentByBWA.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(LongSEAlignmentByBWA)
		workflow.LongSEAlignmentByBWA = LongSEAlignmentByBWA
		
		mergeSamFilesJava = Executable(namespace=namespace, name="mergeSamFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		mergeSamFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		workflow.addExecutable(mergeSamFilesJava)
		workflow.mergeSamFilesJava = mergeSamFilesJava
		
		SortSamFilesJava = Executable(namespace=namespace, name="SortSamFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		SortSamFilesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		SortSamFilesJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.cluster_size_for_aln_jobs))
		workflow.addExecutable(SortSamFilesJava)
		workflow.SortSamFilesJava = SortSamFilesJava

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
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerCustomJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		#individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, dataDir=self.localDataDir)
		isq_id2LibrarySplitOrder2FileLs = db_vervet.getISQ_ID2LibrarySplitOrder2FileLs(self.ind_seq_id_ls, dataDir=self.dataDir, \
												filtered=None, ignoreEmptyReadFile=self.ignoreEmptyReadFile)
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		
		#2011-11-16 new way of registering reference fasta file. but still dont' want to trasnfer 7Gb of data
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]

		if self.needRefIndexJob:
			#refIndexJob = self.addRefIndexJob(workflow, refFastaFList=refFastaFList, refSequence=refSequence, bwa=workflow.bwa, \
			#								namespace=namespace, version=version)
			refIndexJob = self.addStampyGenomeIndexHashJob(workflow, executable=workflow.stampy, refFastaFList=refFastaFList, \
						parentJobLs=[], job_max_memory=3000, job_max_walltime = 200, \
						extraDependentInputLs=[], \
						transferOutput=True)
		else:
			refIndexJob = None
		
		self.addAllAlignmentJobs(db_vervet, individualSequenceID2FilePairLs=None, \
					isq_id2LibrarySplitOrder2FileLs = isq_id2LibrarySplitOrder2FileLs,\
					dataDir=self.dataDir,\
					refSequence=refSequence, refFastaFList=refFastaFList, refIndexJob=refIndexJob,
					workflow=workflow, bwa=workflow.bwa, additionalArguments=self.additionalArguments, \
					samtools=workflow.samtools, \
					mkdirWrap=workflow.mkdirWrap, mv=workflow.cp, \
					java=workflow.java, \
					mergeSamFilesJava=workflow.mergeSamFilesJava, mergeSamFilesJar=workflow.mergeSamFilesJar, \
					MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=self.tmpDir,\
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					SortSamFilesJava=workflow.SortSamFilesJava, SortSamFilesJar=workflow.SortSamFilesJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar,\
					alignment_method_name=self.alignment_method_name, alignment_format='bam',\
					namespace=workflow.namespace, version=workflow.version, stageOutFinalOutput=self.stageOutFinalOutput,\
					PEAlignmentByBWA=workflow.PEAlignmentByBWA, ShortSEAlignmentByBWA=workflow.ShortSEAlignmentByBWA, \
					LongSEAlignmentByBWA=workflow.LongSEAlignmentByBWA, no_of_aln_threads=self.no_of_aln_threads,\
					stampyExecutable=workflow.stampy, skipDoneAlignment=self.skipDoneAlignment)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = ShortRead2AlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
