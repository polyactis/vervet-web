#!/usr/bin/env python
"""
Examples:
	
	#2011-11-4 run on condorpool
	%s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice -o InspectTop804ContigsAlnRefSeq524Alignments.xml -I 552-661 -N 804
	
	#2011-11-5 run it on hoffman2, need ssh tunnel for db (-H)
	%s -a 524 -j hoffman2 -l hoffman2 -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_hoffman2.xml 
		-I 552-661 -e /u/home/eeskin/polyacti/ -m /u/home/eeskin/polyacti/NetworkData/ 
		-J /u/local/apps/java/jre1.6.0_23/bin/java -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /Network/Data/vervet/db/
		-H
	
	#2011-11-5 run on uschpc (input data is on uschpc), for each top contig as well
	%s -a 524 -j uschpc -l uschpc -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_uschpc.xml
		-I 552-661 -P -e /home/cmb-03/mn/yuhuang/ -m /home/cmb-03/mn/yuhuang/tmp/
		-J /usr/usc/jdk/default/bin/java -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -D /Network/Data/vervet/db/
	
	#2011-11-25 on hoffman2's condor pool, need ssh tunnel for db (-H)
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost -N 7559 -o InspectRefSeq524WholeAlignment.xml -C 30
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -J ~/bin/jdk/bin/java -H
	
	#2012.4.3 change tmpDir (-m) for AddOrReplaceReadGroups, no job clustering (-C 1)
	%s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice -o InspectAln1_To_661_RefSeq524Alignments.xml -I 1-661
		-m /Network/Data/vervet/vervetPipeline/tmp/ -C 1
	
	#2012.5.8 do perContig depth estimation (-P) and skip alignments with stats in db already (-s), need ssh tunnel for db (-H)
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost -N 7559 -o workflow/InspectAln1_To_1251_RefSeq524Alignments.xml
		-I 1-1251
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java 
		-P -s -H
	
Description:
	2012.3.21
		use samtools flagstat
	2011-11-4
		a pegasus workflow that inspects no-of-reads-aligned, inferred insert size and etc.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])	#, sys.argv[0], sys.argv[0]

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from AlignmentToCallPipeline import AlignmentToCallPipeline
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class InspectAlignmentPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("tmpDir", 1, ): ["/tmp/", 'm', 1, 'for MarkDuplicates.jar or AddOrReplaceReadGroups.jar, default is /tmp/ but sometimes too small'],\
						("needPerContigJob", 0, int): [0, 'P', 0, 'toggle to add DepthOfCoverage and VariousReadCount jobs for each contig.'],\
						("skipAlignmentWithStats", 0, int): [0, 's', 0, 'If an alignment has depth stats filled, not DOC job will be run. similar for flagstat job.'],\
						("fractionToSample", 0, float): [0.001, '', 1, 'fraction of loci to walk through for DepthOfCoverage walker.'],\
						("needFastaIndexJob", 0, int): [0, '', 0, 'toggle to add a reference index job by samtools'],\
						("needFastaDictJob", 0, int): [0, '', 0, 'toggle to add a reference dict job by picard CreateSequenceDictionary.jar'],\
						("sequence_filtered", 0, int): [None, '', 1, 'To filter alignments. None: whatever; 0: unfiltered sequences, 1: filtered sequences'],\
						("alignment_method_id", 0, int): [None, '', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						})

	def __init__(self, **keywords):
		"""
		2011-11-4
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		#AlignmentToCallPipeline.__init__(self, **keywords)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		if self.ind_aln_id_ls:
			self.ind_aln_id_ls = getListOutOfStr(self.ind_aln_id_ls, data_type=int)
		
	
	def addDepthOfCoverageJob(self, workflow, DOCWalkerJava=None, genomeAnalysisTKJar=None,\
							refFastaFList=None, bamF=None, baiF=None, DOCOutputFnamePrefix=None,\
							fractionToSample=None,\
							parentJobLs=[], job_max_memory = 1000, additionalArguments="", \
							transferOutput=False, minMappingQuality=30, minBaseQuality=20):
		"""
		2012.5.7
			no longer used, superceded by addSAMtoolsDepthJob()
		2012.4.17
			add --omitIntervalStatistics and --omitLocusTable to the walker
		2012.4.12
			add "--read_filter BadCigar" to GATK to avoid stopping because of malformed Cigar
				malformed: Read starting with deletion. Cigar: 1D65M299S 
		2012.4.3
			add argument fractionToSample
		2011-11-25
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		DOCJob = Job(namespace=workflow.namespace, name=DOCWalkerJava.name, version=workflow.version)
		DOCJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "DepthOfCoverage",\
			'-R', refFastaF, '-o', DOCOutputFnamePrefix, "-pt sample", "--omitDepthOutputAtEachBase",\
			"-mmq %s"%(minMappingQuality), "-mbq %s"%(minBaseQuality), "--read_filter BadCigar", \
			'--omitLocusTable', '--omitIntervalStatistics')
		DOCJob.addArguments("-I", bamF)
		if fractionToSample and fractionToSample>0 and fractionToSample<=1:
			DOCJob.addArguments("--fractionToSample %s"%(fractionToSample))
		if additionalArguments:
			DOCJob.addArguments(additionalArguments)
		#it's either symlink or stage-in
		DOCJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		DOCJob.uses(baiF, transfer=True, register=True, link=Link.INPUT)
		self.registerFilesAsInputToJob(DOCJob, refFastaFList)
		DOCJob.sample_summary_file = File('%s.sample_summary'%(DOCOutputFnamePrefix))
		DOCJob.sample_interval_summary_file =File('%s.sample_interval_summary'%(DOCOutputFnamePrefix)) 
		DOCJob.uses(DOCJob.sample_summary_file, transfer=transferOutput, register=True, link=Link.OUTPUT)
		DOCJob.uses(DOCJob.sample_interval_summary_file, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(DOCJob)
		yh_pegasus.setJobProperRequirement(DOCJob, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=DOCJob)
		return DOCJob
	
	def addSAMtoolsDepthJob(self, workflow, samtoolsDepth=None, samtools_path=None,\
							bamF=None, outputFile=None, baiF=None, \
							parentJobLs=[], job_max_memory = 500, additionalArguments="", \
							transferOutput=False, minMappingQuality=30, minBaseQuality=20):
		"""
		2012.5.7
			
		"""
		job = Job(namespace=workflow.namespace, name=samtoolsDepth.name, version=workflow.version)
		job.addArguments(samtools_path, bamF, outputFile, "%s"%minMappingQuality, "%s"%minBaseQuality)
		if additionalArguments:
			job.addArguments(additionalArguments)
		#it's either symlink or stage-in
		job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addCalculateDepthMeanMedianModeJob(self, workflow, executable=None, \
							inputFile=None, outputFile=None, alignmentID=None, fractionToSample=0.001, \
							noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E7,\
							parentJobLs=[], job_max_memory = 500, additionalArguments=None, \
							transferOutput=False,):
		"""
		2012.6.15 turn maxNumberOfSamplings into integer when passing it to the job
		2012.5.7
			a job to take input of samtoolsDepth
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-i", inputFile, "-o", outputFile, "-a %s"%(alignmentID), "-n %s"%(noOfLinesInHeader), "-f %s"%(fractionToSample),\
						"-w %s"%(whichColumn), '-m %d'%(maxNumberOfSamplings))
		if additionalArguments:
			job.addArguments(additionalArguments)
		#it's either symlink or stage-in
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addReadCountJob(self, workflow, VariousReadCountJava=None, genomeAnalysisTKJar=None,\
					refFastaFList=None, bamF=None, baiF=None, readCountOutputF=None,\
					parentJobLs=[], job_max_memory = 1000, additionalArguments="", \
					transferOutput=False):
		"""
		2011-11-25
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		job = Job(namespace=workflow.namespace, name=VariousReadCountJava.name, version=workflow.version)
		job.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "VariousReadCount",\
			'-R', refFastaF, '-o', readCountOutputF, "-mmq 30")
		job.addArguments("-I", bamF)
		if additionalArguments:
			job.addArguments(additionalArguments)
		job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
		self.registerFilesAsInputToJob(job, refFastaFList)
		job.output = readCountOutputF
		job.uses(readCountOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addReformatFlagstatOutputJob(self, workflow, executable=None, inputF=None, alignmentID=None, outputF=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.4.3
			input is output of "samtools flagstat"
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-a %s'%(alignmentID), '-i', inputF, '-o', outputF)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for parentJob in parentJobLs:
			if parentJob:
				workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addJobs(self, workflow, alignmentDataLs, refName2size, samtools=None, DOCWalkerJava=None, \
				ContigDOCWalkerJava=None, ContigVariousReadCountJava=None, \
				VariousReadCountJava=None,  genomeAnalysisTKJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir="/Network/Data/vervet/vervetPipeline/tmp/",\
				mkdirWrap=None, mv=None, \
				ReformatFlagstatOutput=None, PutFlagstatOutput2DB=None, PutDOCOutput2DB=None,\
				samtoolsDepth=None, \
				reduceDepthOfCoverage=None, reduceVariousReadCount=None,\
				refFastaFList=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				dataDir=None, needPerContigJob=False, skipAlignmentWithStats=False,\
				needSSHDBTunnel=0):
		"""
		2012.6.15
			replace  CalculateMedianModeFromSAMtoolsDepthOutput with CalculateMedianMeanOfInputColumn
		2012.5.7
			replace GATK's DOC walker with samtoolsDepth + CalculateMedianModeFromSAMtoolsDepthOutput
			add argument needSSHDBTunnel
		2012.4.6
			add argument skipAlignmentWithStats to skip alignments with stats
		2012.4.3
			add ReformatFlagstatOutput, PutDOCOutput2DB, PutFlagstatOutput2DB jobs
		2011-11-4
			
		"""
		sys.stderr.write("Adding jobs for %s references and %s alignments..."%(len(refName2size), len(alignmentDataLs)))
		perContigJobMaxMemory = 1000	#in MB
		refFastaF = refFastaFList[0]
		no_of_jobs = 0
		if needFastaDictJob:
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
										refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
			no_of_jobs += 1
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
			no_of_jobs += 1
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
		"""
		#2012.4.3 no more VariousReadCountJava job
		readCountOutputF = File('VariousReadCount.tsv')
		readCountOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=readCountOutputF, namespace=namespace, version=version)
		"""
		
		depthOfCoverageOutputF = File('DepthOfCoverage.tsv')
		depthOfCoverageOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=depthOfCoverageOutputF, namespace=namespace, version=version,\
							transferOutput=True)
		
		depthOfCoveragePerChrOutputF = File('DepthOfCoveragePerChr.tsv')
		depthOfCoveragePerChrOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=depthOfCoveragePerChrOutputF, namespace=namespace, version=version,\
							transferOutput=True)
		flagStatOutputF = File('FlagStat.tsv')
		flagStatOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=flagStatOutputF, namespace=namespace, version=version,\
							transferOutput=True)
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow, alignmentDataLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir, tmpDir=tmpDir)
		#alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		no_of_jobs += 2
		no_of_alns_with_depth_jobs = 0
		no_of_alns_with_flagstat_jobs = 0
		samtools_path = yh_pegasus.getAbsPathOutOfExecutable(samtools)
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			
			bamF= alignmentData.bamF
			baiF = alignmentData.baiF
			
			if skipAlignmentWithStats and alignment.median_depth is not None and alignment.mean_depth is not None and alignment.mode_depth is not None:
				pass
			else:
				#DOCOutputFnamePrefix = '%s_DOC'%(alignment.id)
				depthOutputFile = File('%s_DOC.tsv.gz'%(alignment.id))
				samtoolsDepthJob = self.addSAMtoolsDepthJob(workflow, samtoolsDepth=samtoolsDepth, samtools_path=samtools_path,\
							bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
							parentJobLs=alignmentData.jobLs, job_max_memory = 500, additionalArguments="", \
							transferOutput=False)
				meanMedianModeDepthFile = File("%s_meanMedianModeDepth.tsv"%(alignment.id))
				meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(workflow, \
							executable=workflow.CalculateMedianMeanOfInputColumn, \
							inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, alignmentID=alignment.id, fractionToSample=self.fractionToSample, \
							noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E7,\
							parentJobLs=[samtoolsDepthJob], job_max_memory = 500, additionalArguments=None, \
							transferOutput=False)
				"""
				DOCJob = self.addDepthOfCoverageJob(workflow, DOCWalkerJava=DOCWalkerJava, genomeAnalysisTKJar=genomeAnalysisTKJar,\
							refFastaFList=refFastaFList, bamF=bamF, baiF=baiF, \
							DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
							parentJobLs=alignmentData.jobLs, job_max_memory = perContigJobMaxMemory*8, transferOutput=True,\
							fractionToSample=self.fractionToSample)
				"""
				self.addInputToStatMergeJob(workflow, statMergeJob=depthOfCoverageOutputMergeJob, inputF=meanMedianModeDepthFile,\
							parentJobLs=[meanMedianModeDepthJob])
				no_of_jobs += 2
				no_of_alns_with_depth_jobs += 1
				
			if skipAlignmentWithStats and alignment.perc_reads_mapped is not None:
				pass
			else:
				oneFlagStatOutputF = File('%s_flagstat.txt'%(alignment.id))
				samtoolsFlagStatJob = self.addSamtoolsFlagstatJob(workflow, executable=samtools, inputF=bamF, outputF=oneFlagStatOutputF, \
					parentJobLs=alignmentData.jobLs, extraDependentInputLs=[baiF], transferOutput=False, \
					extraArguments=None, job_max_memory=perContigJobMaxMemory)
				
				reformatFlagStatOutputF = File('%s_flagstat.tsv'%(alignment.id))
				reformatFlagStatOutputJob = self.addReformatFlagstatOutputJob(workflow, executable=ReformatFlagstatOutput, \
									inputF=oneFlagStatOutputF, alignmentID=alignment.id, outputF=reformatFlagStatOutputF, \
									parentJobLs=[samtoolsFlagStatJob], extraDependentInputLs=[], transferOutput=False, \
									extraArguments=None, job_max_memory=20)
				self.addInputToStatMergeJob(workflow, statMergeJob=flagStatOutputMergeJob, inputF=reformatFlagStatOutputJob.output, \
							parentJobLs=[reformatFlagStatOutputJob])
				no_of_jobs += 2
				no_of_alns_with_flagstat_jobs += 1
			"""
				#2012.4.3 no more VariousReadCountJava job
			readCountOutputF = File('%s_variousReadCount.tsv'%(alignment.id))
			readCountJob = self.addReadCountJob(workflow, VariousReadCountJava=VariousReadCountJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, refFastaFList=refFastaFList, \
						bamF=bamF, baiF=baiF, readCountOutputF=readCountOutputF,\
						parentJobLs=None, job_max_memory = perContigJobMaxMemory*2, \
						transferOutput=True)
			self.addInputToStatMergeJob(workflow, statMergeJob=readCountOutputMergeJob, inputF=readCountOutputF,\
						parentJobLs=[readCountJob])
			"""
			if not needPerContigJob:	#no need for per-contig job
				continue
			
			statOutputDir = '%s'%(alignment.id)
			statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=statOutputDir, namespace=namespace, \
													version=version)
			
			"""
				#2012.4.3 no more VariousReadCountJava job
			finalVariousReadCountOutputF = File("%s_VariousReadCount.tsv"%(alignment.getReadGroup()))
			reduceVariousReadCountJob = Job(namespace=namespace, name=reduceVariousReadCount.name, version=version)
			reduceVariousReadCountJob.addArguments("-o", finalVariousReadCountOutputF)
			reduceVariousReadCountJob.uses(finalVariousReadCountOutputF, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(reduceVariousReadCountJob)
			"""
			
			"""
			finalDepthOfCoverageOutputF = File("%s_DepthOfCoverage.tsv"%(alignment.getReadGroup()))
			reduceDepthOfCoverageJob = Job(namespace=namespace, name=reduceDepthOfCoverage.name, version=version)
			reduceDepthOfCoverageJob.addArguments("-o", finalDepthOfCoverageOutputF)
			reduceDepthOfCoverageJob.uses(finalDepthOfCoverageOutputF, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(reduceDepthOfCoverageJob)
			"""
			
			for refName, refSize in refName2size.iteritems():	#this could be further improved to work on per-interval
				depthOutputFile = File(os.path.join(statOutputDir, '%s_%s_DOC.tsv.gz'%(alignment.id, refName)))
				samtoolsDepthJob = self.addSAMtoolsDepthJob(workflow, samtoolsDepth=samtoolsDepth, samtools_path=samtools_path,\
							bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
							parentJobLs=[statOutputDirJob]+alignmentData.jobLs, job_max_memory = 500, additionalArguments="", \
							transferOutput=False)
				meanMedianModeDepthFile = File(os.path.join(statOutputDir, "%s_%s_meanMedianModeDepth.tsv"%(alignment.id, refName)))
				meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(workflow, \
							executable=workflow.CalculateMedianMeanOfInputColumn, \
							inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, alignmentID="%s-%s"%(alignment.id, refName), \
							fractionToSample=self.fractionToSample, \
							noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E7,\
							parentJobLs=[samtoolsDepthJob], job_max_memory = 500, additionalArguments="-r %s"%(refName), \
							transferOutput=False)
				
				self.addInputToStatMergeJob(workflow, statMergeJob=depthOfCoveragePerChrOutputMergeJob, inputF=meanMedianModeDepthFile,\
							parentJobLs=[meanMedianModeDepthJob])
				"""
				DOCOutputFnamePrefix = os.path.join(statOutputDir, '%s_%s_DOC'%(alignment.id, refName))
				DOCJob = self.addDepthOfCoverageJob(workflow, DOCWalkerJava=ContigDOCWalkerJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar,\
						refFastaFList=refFastaFList, bamF=bamF, baiF=baiF, \
						DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
						parentJobLs=[statOutputDirJob]+alignmentData.jobLs, job_max_memory = perContigJobMaxMemory*3, additionalArguments="-L %s"%(refName),\
						transferOutput=False,\
						fractionToSample=self.fractionToSample)
				
				reduceDepthOfCoverageJob.addArguments(DOCJob.sample_interval_summary_file)
				reduceDepthOfCoverageJob.uses(DOCJob.sample_interval_summary_file, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=DOCJob, child=reduceDepthOfCoverageJob)
				"""
				
				no_of_jobs += 1
				"""
				#2012.4.3 no more VariousReadCountJava job
				readCountOutputF = File(os.path.join(statOutputDir, '%s_%s_variousReadCount.tsv'%(alignment.id, refName)))
				readCountJob = self.addReadCountJob(workflow, VariousReadCountJava=ContigVariousReadCountJava, \
							genomeAnalysisTKJar=genomeAnalysisTKJar, refFastaFList=refFastaFList, \
							bamF=bamF, baiF=baiF, readCountOutputF=readCountOutputF,\
							parentJobLs=statOutputDirJob, job_max_memory = perContigJobMaxMemory, additionalArguments="-L %s"%(refName), \
							transferOutput=False)
				
				reduceVariousReadCountJob.addArguments(readCountOutputF)
				reduceVariousReadCountJob.uses(readCountOutputF, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=readCountJob, child=reduceVariousReadCountJob)
				"""
		
		flagStat2DBLogFile = File("flagStat2DB.log")
		flagStat2DBJob = self.addPutStuffIntoDBJob(workflow, executable=PutFlagstatOutput2DB, inputFileLs=[flagStatOutputF], \
					logFile=flagStat2DBLogFile, commit=True, \
					parentJobLs=[flagStatOutputMergeJob], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
		DOC2DBLogFile = File("DOC2DB.log")
		DOC2DBJob = self.addPutStuffIntoDBJob(workflow, executable=PutDOCOutput2DB, inputFileLs=[depthOfCoverageOutputF], \
					logFile=DOC2DBLogFile, commit=True, \
					parentJobLs=[depthOfCoverageOutputMergeJob], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
		
		no_of_jobs += 2
		sys.stderr.write(" %s jobs, %s alignments with depth jobs, %s alignments with flagstat jobs.\n"%(no_of_jobs, \
																				no_of_alns_with_depth_jobs, no_of_alns_with_flagstat_jobs))
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-25
			split out of run()
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableList = []
		
		reduceDepthOfCoverage = Executable(namespace=namespace, name="ReduceDepthOfCoverage", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceDepthOfCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/ReduceDepthOfCoverage.py"), site_handler))
		#clustering is buggy for programs with long arguments.
		#reduceDepthOfCoverage.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(reduceDepthOfCoverage)
		self.reduceDepthOfCoverage = reduceDepthOfCoverage
		
		reduceVariousReadCount = Executable(namespace=namespace, name="ReduceVariousReadCount", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceVariousReadCount.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/ReduceVariousReadCount.py"), site_handler))
		#clustering is buggy for programs with long arguments.
		#reduceVariousReadCount.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(reduceVariousReadCount)
		self.reduceVariousReadCount = reduceVariousReadCount
		
		ContigDOCWalkerJava = Executable(namespace=namespace, name="ContigDOCWalkerJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		ContigDOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableList.append(ContigDOCWalkerJava)
		
		
		ContigVariousReadCountJava = Executable(namespace=namespace, name="ContigVariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		ContigVariousReadCountJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableList.append(ContigVariousReadCountJava)
		
		ReformatFlagstatOutput = Executable(namespace=namespace, name="ReformatFlagstatOutput", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		ReformatFlagstatOutput.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/ReformatFlagstatOutput.py"), site_handler))
		executableList.append(ReformatFlagstatOutput)
		
		samtoolsDepth = Executable(namespace=namespace, name="samtoolsDepth", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		samtoolsDepth.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/samtoolsDepth.sh"), site_handler))
		executableList.append(samtoolsDepth)
		
		CalculateMedianModeFromSAMtoolsDepthOutput = Executable(namespace=namespace, name="CalculateMedianModeFromSAMtoolsDepthOutput", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		CalculateMedianModeFromSAMtoolsDepthOutput.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py"), site_handler))
		executableList.append(CalculateMedianModeFromSAMtoolsDepthOutput)
		
		CalculateMedianMeanOfInputColumn = Executable(namespace=namespace, name="CalculateMedianMeanOfInputColumn", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		CalculateMedianMeanOfInputColumn.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/CalculateMedianMeanOfInputColumn"), site_handler))
		executableList.append(CalculateMedianMeanOfInputColumn)
				
		
		PutFlagstatOutput2DB = Executable(namespace=namespace, name="PutFlagstatOutput2DB", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		PutFlagstatOutput2DB.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/PutFlagstatOutput2DB.py"), site_handler))
		#PutFlagstatOutput2DB.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(PutFlagstatOutput2DB)
		self.PutFlagstatOutput2DB = PutFlagstatOutput2DB
		
		PutDOCOutput2DB = Executable(namespace=namespace, name="PutDOCOutput2DB", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		PutDOCOutput2DB.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/PutDOCOutput2DB.py"), site_handler))
		#PutDOCOutput2DB.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		self.addExecutable(PutDOCOutput2DB)
		self.PutDOCOutput2DB = PutDOCOutput2DB
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
		
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
		
		workflow = self.initiateWorkflow()
		
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		if self.needPerContigJob:
			refName2size = self.getTopNumberOfContigs(self.contigMaxRankBySize)
			#refName2size = set(['Contig149'])	#temporary when testing Contig149
			#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		else:
			refName2size = {}
		
		alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
									dataDir=self.localDataDir, aln_method_id=self.alignment_method_id)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars(workflow)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
		
		self.addJobs(workflow, alignmentDataLs, refName2size, samtools=workflow.samtools, DOCWalkerJava=workflow.DOCWalkerJava, 
				ContigDOCWalkerJava=workflow.ContigDOCWalkerJava, \
				VariousReadCountJava=workflow.VariousReadCountJava, \
				ContigVariousReadCountJava=workflow.ContigVariousReadCountJava, \
				genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
				createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
				addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar,\
				MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=self.tmpDir,\
				mkdirWrap=workflow.mkdirWrap, mv=workflow.mv, \
				ReformatFlagstatOutput=workflow.ReformatFlagstatOutput, PutFlagstatOutput2DB=workflow.PutFlagstatOutput2DB, \
				PutDOCOutput2DB=workflow.PutDOCOutput2DB,\
				samtoolsDepth=workflow.samtoolsDepth, \
				reduceDepthOfCoverage=workflow.reduceDepthOfCoverage, reduceVariousReadCount=workflow.reduceVariousReadCount,\
				
				refFastaFList=refFastaFList, \
				namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				dataDir=self.dataDir, needPerContigJob=self.needPerContigJob,\
				skipAlignmentWithStats=self.skipAlignmentWithStats,\
				needSSHDBTunnel=self.needSSHDBTunnel)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = InspectAlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
