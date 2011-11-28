#!/usr/bin/env python
"""
Examples:
	
	#2011-11-4 run on condorpool
	%s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice -o InspectTop804ContigsAlnRefSeq524Alignments.xml -I 552-661 -N 804
	
	#2011-11-5 run it on hoffman2
	%s -a 524 -j hoffman2 -l hoffman2 -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_hoffman2.xml 
		-I 552-661 -e /u/home/eeskin/polyacti/ -m /u/home/eeskin/polyacti/NetworkData/ 
		-J /u/local/apps/java/jre1.6.0_23/bin/java -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /Network/Data/vervet/db/
	
	#2011-11-5 run on uschpc (input data is on uschpc), for each top contig as well
	%s -a 524 -j uschpc -l uschpc -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_uschpc.xml
		-I 552-661 -P -e /home/cmb-03/mn/yuhuang/ -m /home/cmb-03/mn/yuhuang/tmp/
		-J /usr/usc/jdk/default/bin/java -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -D /Network/Data/vervet/db/
	
	#2011-11-25 on hoffman2's condor pool
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost -N 7559 -o InspectRefSeq524WholeAlignment.xml -C 30
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -J ~/bin/jdk/bin/java
	 
Description:
	2011-11-4
		a pegasus workflow that inspects no-of-reads-aligned, inferred insert size and etc.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])	# sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0]

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
from AlignmentToCallPipeline import AlignmentToCallPipeline
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class InspectAlignmentPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict
	option_default_dict.update({
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("contigMinRankBySize", 1, int): [1, 'M', 1, 'min rank (by size) of contigs'],\
						("tmpDir", 1, ): ["/tmp/", 'm', 1, 'for MarkDuplicates.jar, default is /tmp/ but sometimes too small'],\
						("topNumberOfContigs", 1, int): [156, 'N', 1, 'max rank (by size) of contigs'],\
						("needPerContigJob", 0, int): [0, 'P', 0, 'toggle to add DepthOfCoverage and VariousReadCount jobs for each contig.'],\
						("needFastaIndexJob", 0, int): [0, '', 0, 'toggle to add a reference index job by samtools'],\
						("needFastaDictJob", 0, int): [0, '', 0, 'toggle to add a reference dict job by picard CreateSequenceDictionary.jar'],\
						})

	def __init__(self,  **keywords):
		"""
		2011-11-4
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
	
	def addDepthOfCoverageJob(self, workflow, DOCWalkerJava=None, genomeAnalysisTKJar=None,\
							refFastaFList=None, bamF=None, baiFile=None, DOCOutputFnamePrefix=None,\
							parentJob=None, job_max_memory = 1000, additionalArguments="", \
							transferOutput=False):
		"""
		2011-11-25
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		DOCJob = Job(namespace=workflow.namespace, name=DOCWalkerJava.name, version=workflow.version)
		DOCJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "DepthOfCoverage",\
			'-R', refFastaF, '-o', DOCOutputFnamePrefix, "-pt sample", "--omitDepthOutputAtEachBase",\
			"-mmq 30", "-mbq 20")
		DOCJob.addArguments("-I", bamF)
		if additionalArguments:
			DOCJob.addArguments(additionalArguments)
		#it's either symlink or stage-in
		DOCJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		DOCJob.uses(baiFile, transfer=True, register=True, link=Link.INPUT)
		self.registerFilesAsInputToJob(DOCJob, refFastaFList)
		DOCJob.sample_summary_file = File('%s.sample_summary'%(DOCOutputFnamePrefix))
		DOCJob.sample_interval_summary_file =File('%s.sample_interval_summary'%(DOCOutputFnamePrefix)) 
		DOCJob.uses(DOCJob.sample_summary_file, transfer=transferOutput, register=True, link=Link.OUTPUT)
		DOCJob.uses(DOCJob.sample_interval_summary_file, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(DOCJob)
		yh_pegasus.setJobProperRequirement(DOCJob, job_max_memory=job_max_memory)
		if parentJob:
			workflow.depends(parent=parentJob, child=DOCJob)
		return DOCJob
	
	def addReadCountJob(self, workflow, VariousReadCountJava=None, genomeAnalysisTKJar=None,\
					refFastaFList=None, bamF=None, baiFile=None, readCountOutputF=None,\
					parentJob=None, job_max_memory = 1000, additionalArguments="", \
					transferOutput=False):
		"""
		2011-11-25
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		readCountJob = Job(namespace=workflow.namespace, name=VariousReadCountJava.name, version=workflow.version)
		readCountJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "VariousReadCount",\
			'-R', refFastaF, '-o', readCountOutputF, "-mmq 30")
		readCountJob.addArguments("-I", bamF)
		if additionalArguments:
			DOCJob.addArguments(additionalArguments)
		readCountJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		readCountJob.uses(baiFile, transfer=True, register=True, link=Link.INPUT)
		self.registerFilesAsInputToJob(readCountJob, refFastaFList)
		readCountJob.output = readCountOutputF
		readCountJob.uses(readCountOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(readCountJob)
		yh_pegasus.setJobProperRequirement(readCountJob, job_max_memory=job_max_memory)
		if parentJob:
			workflow.depends(parent=parentJob, child=readCountJob)
		return readCountJob
	
	def addJobs(self, workflow, alignmentLs, refName2size, samtools=None, DOCWalkerJava=None, \
				ContigDOCWalkerJava=None, ContigVariousReadCountJava=None, \
				VariousReadCountJava=None,  genomeAnalysisTKJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir="/Network/Data/vervet/vervetPipeline/tmp/",\
				mkdirWrap=None, mv=None, \
				
				reduceDepthOfCoverage=None, reduceVariousReadCount=None,\
				refFastaFList=None, \
				
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				dataDir=None, needPerContigJob=False):
		"""
		2011-11-4
			
		"""
		sys.stderr.write("Adding jobs for %s references and %s alignments..."%(len(refName2size), len(alignmentLs)))
		perContigJobMaxMemory = 1000	#in MB
		refFastaF = refFastaFList[0]
		if needFastaDictJob:
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
										refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
		readCountOutputF = File('VariousReadCount.tsv')
		readCountOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=readCountOutputF, namespace=namespace, version=version)
		
		depthOfCoverageOutputF = File('DepthOfCoverage.tsv')
		depthOfCoverageOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=depthOfCoverageOutputF, namespace=namespace, version=version)
				
		returnData = self.addAddRG2BamJobsAsNeeded(workflow, alignmentLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		
		for alignment in alignmentLs:
			#if alignment.id in alignmentId2RGJobDataLs:
				index_sam_job, bamF, baiFile = alignmentId2RGJobDataLs[alignment.id]
				
				DOCOutputFnamePrefix = '%s_DOC'%(alignment.id)
				DOCJob = self.addDepthOfCoverageJob(workflow, DOCWalkerJava=DOCWalkerJava, genomeAnalysisTKJar=genomeAnalysisTKJar,\
							refFastaFList=refFastaFList, bamF=bamF, baiFile=baiFile, \
							DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
							parentJob=None, job_max_memory = perContigJobMaxMemory*2, transferOutput=True)
				self.addInputToStatMergeJob(workflow, statMergeJob=depthOfCoverageOutputMergeJob, inputF=DOCJob.sample_summary_file,\
							parentJobLs=[DOCJob])
				
				readCountOutputF = File('%s_variousReadCount.tsv'%(alignment.id))
				readCountJob = self.addReadCountJob(workflow, VariousReadCountJava=VariousReadCountJava, \
							genomeAnalysisTKJar=genomeAnalysisTKJar, refFastaFList=refFastaFList, \
							bamF=bamF, baiFile=baiFile, readCountOutputF=readCountOutputF,\
							parentJob=None, job_max_memory = perContigJobMaxMemory*2, \
							transferOutput=True)
				self.addInputToStatMergeJob(workflow, statMergeJob=readCountOutputMergeJob, inputF=readCountOutputF,\
							parentJobLs=[readCountJob])
				
				if not needPerContigJob:	#no need for per-contig job
					continue
				
				statOutputDir = '%s'%(alignment.id)
				statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=statOutputDir, namespace=namespace, \
														version=version)
				
				finalVariousReadCountOutputF = File("%s_VariousReadCount.tsv"%(alignment.getReadGroup()))
				reduceVariousReadCountJob = Job(namespace=namespace, name=reduceVariousReadCount.name, version=version)
				reduceVariousReadCountJob.addArguments("-o", finalVariousReadCountOutputF)
				reduceVariousReadCountJob.uses(finalVariousReadCountOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(reduceVariousReadCountJob)
				
				finalDepthOfCoverageOutputF = File("%s_DepthOfCoverage.tsv"%(alignment.getReadGroup()))
				reduceDepthOfCoverageJob = Job(namespace=namespace, name=reduceDepthOfCoverage.name, version=version)
				reduceDepthOfCoverageJob.addArguments("-o", finalDepthOfCoverageOutputF)
				reduceDepthOfCoverageJob.uses(finalDepthOfCoverageOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(reduceDepthOfCoverageJob)
				
				for refName, refSize in refName2size.iteritems():	#this could be further improved to work on per-interval
					DOCOutputFnamePrefix = os.path.join(statOutputDir, '%s_%s_DOC'%(alignment.id, refName))
					DOCJob = self.addDepthOfCoverageJob(workflow, DOCWalkerJava=ContigDOCWalkerJava, \
							genomeAnalysisTKJar=genomeAnalysisTKJar,\
							refFastaFList=refFastaFList, bamF=bamF, baiFile=baiFile, \
							DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
							parentJob=statOutputDirJob, job_max_memory = perContigJobMaxMemory, additionalArguments="-L %s"%(refName),\
							transferOutput=False)
					
					reduceDepthOfCoverageJob.addArguments(DOCJob.sample_interval_summary_file)
					reduceDepthOfCoverageJob.uses(DOCJob.sample_interval_summary_file, transfer=True, register=True, link=Link.INPUT)
					workflow.depends(parent=DOCJob, child=reduceDepthOfCoverageJob)
					
					
					readCountOutputF = File(os.path.join(statOutputDir, '%s_%s_variousReadCount.tsv'%(alignment.id, refName)))
					readCountJob = self.addReadCountJob(workflow, VariousReadCountJava=ContigVariousReadCountJava, \
								genomeAnalysisTKJar=genomeAnalysisTKJar, refFastaFList=refFastaFList, \
								bamF=bamF, baiFile=baiFile, readCountOutputF=readCountOutputF,\
								parentJob=statOutputDirJob, job_max_memory = perContigJobMaxMemory, additionalArguments="-L %s"%(refName), \
								transferOutput=False)
					
					reduceVariousReadCountJob.addArguments(readCountOutputF)
					reduceVariousReadCountJob.uses(readCountOutputF, transfer=True, register=True, link=Link.INPUT)
					workflow.depends(parent=readCountJob, child=reduceVariousReadCountJob)
					
					if index_sam_job is not None:
						workflow.depends(parent=index_sam_job, child=DOCJob)
						workflow.depends(parent=index_sam_job, child=readCountJob)
		sys.stderr.write(".Done\n")
	
	def registerCustomExecutables(self, workflow):
		"""
		2011-11-25
			split out of run()
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		reduceDepthOfCoverage = Executable(namespace=namespace, name="ReduceDepthOfCoverage", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceDepthOfCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reduce/ReduceDepthOfCoverage.py"), site_handler))
		#clustering is buggy for programs with long arguments.
		#reduceDepthOfCoverage.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(reduceDepthOfCoverage)
		workflow.reduceDepthOfCoverage = reduceDepthOfCoverage
		
		reduceVariousReadCount = Executable(namespace=namespace, name="ReduceVariousReadCount", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceVariousReadCount.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reduce/ReduceVariousReadCount.py"), site_handler))
		#clustering is buggy for programs with long arguments.
		#reduceVariousReadCount.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(reduceVariousReadCount)
		workflow.reduceVariousReadCount = reduceVariousReadCount
		
		ContigDOCWalkerJava = Executable(namespace=namespace, name="ContigDepthOfCoverageWalker", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		ContigDOCWalkerJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#this is for one contig only, so add clusters_size
		ContigDOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(ContigDOCWalkerJava)
		workflow.ContigDOCWalkerJava = ContigDOCWalkerJava
		
		ContigVariousReadCountJava = Executable(namespace=namespace, name="ContigVariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		ContigVariousReadCountJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#this is for one contig only, so add clusters_size
		ContigVariousReadCountJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(ContigVariousReadCountJava)
		workflow.ContigVariousReadCountJava = ContigVariousReadCountJava
		
	
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
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		refName2size = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refName2size = set(['Contig149'])	#temporary when testing Contig149
		#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = refName2size.keys()
		
		alignmentLs = self.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2, dataDir=self.localDataDir)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars(workflow)
		
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		self.addJobs(workflow, alignmentLs, refName2size, samtools=workflow.samtools, DOCWalkerJava=workflow.DOCWalkerJava, 
				ContigDOCWalkerJava=workflow.ContigDOCWalkerJava, \
				VariousReadCountJava=workflow.VariousReadCountJava, \
				ContigVariousReadCountJava=workflow.ContigVariousReadCountJava, \
				genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
				createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
				addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar,\
				MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=self.tmpDir,\
				mkdirWrap=workflow.mkdirWrap, mv=workflow.mv, \
				reduceDepthOfCoverage=workflow.reduceDepthOfCoverage, reduceVariousReadCount=workflow.reduceVariousReadCount,\
				refFastaFList=refFastaFList, \
				namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				dataDir=self.dataDir, needPerContigJob=self.needPerContigJob)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = InspectAlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
