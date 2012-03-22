#!/usr/bin/env python
"""
Examples:
	# 
	%s 
	
	# 2012.3.14 
	%s -i 1-864 -o workflow/read_count_isq_1_864.xml -u yh -l condorpool -j condorpool -z uclaOffice -F readCount -c
	
Description:
	2012.3.14
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow


class ReadFileBaseCountWorkflow(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						})

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
	
	def registerCustomExecutables(self, workflow):
		"""
		2012.3.14
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		
		CountFastqReadBaseCount = Executable(namespace=namespace, name="CountFastqReadBaseCount", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		CountFastqReadBaseCount.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/CountFastqReadBaseCount.py"), \
										site_handler))
		CountFastqReadBaseCount.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CountFastqReadBaseCount)
		workflow.CountFastqReadBaseCount = CountFastqReadBaseCount
		
		
		PutReadBaseCountIntoDB = Executable(namespace=namespace, name="PutReadBaseCountIntoDB", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		PutReadBaseCountIntoDB.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/PutReadBaseCountIntoDB.py"), \
										site_handler))
		PutReadBaseCountIntoDB.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PutReadBaseCountIntoDB)
		workflow.PutReadBaseCountIntoDB = PutReadBaseCountIntoDB
		
	
	def registerISQFiles(self, workflow=None, db_vervet=None, ind_seq_id_ls=[], localDataDir='', pegasusFolderName='', \
						input_site_handler='local'):
		"""
		2012.3.14
		"""
		sys.stderr.write("Finding all ISQ-affiliated files of %s ind seq entries ..."%(len(ind_seq_id_ls)))
		returnData = PassingData(jobDataLs=[])
		counter = 0
		Table = VervetDB.IndividualSequence
		query = Table.query.filter(Table.id.in_(ind_seq_id_ls))
		for individual_sequence in query:
			if individual_sequence.individual_sequence_file_ls:	#not empty
				for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
					absPath = os.path.join(localDataDir, individual_sequence_file.path)
					inputF = File(os.path.join(pegasusFolderName, individual_sequence_file.path))
					inputF.addPFN(PFN("file://" + absPath, input_site_handler))
					inputF.absPath = absPath
					workflow.addFile(inputF)
					returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
														isqf_id=individual_sequence_file.id))
			elif individual_sequence.path:
				absPath = os.path.join(localDataDir, individual_sequence.path)
				if os.path.isfile(absPath):
					inputF = File(os.path.join(pegasusFolderName, individual_sequence.path))
					inputF.addPFN(PFN("file://" + absPath, input_site_handler))
					inputF.absPath = absPath
					workflow.addFile(inputF)
					returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
														isqf_id=None))
				else:
					sys.stderr.write("Warning: IndividualSequence.id=%s doesn't have any affiliated IndividualSequenceFile entries while its path %s is not a file.\n"%\
									(individual_sequence.id, individual_sequence.path))
		
		sys.stderr.write(" %s files registered.\n"%(len(returnData.jobDataLs)))
		return returnData
	
	def addPutReadBaseCountIntoDBJob(self, workflow, executable=None, inputFileLs=[], \
					logFile=None, commit=False, parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, **keywords):
		"""
		2012.3.14
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-v", self.drivername, "-z", self.hostname, "-d", self.dbname, \
						"-u", self.db_user, "-p", self.db_passwd, \
						"--logFilename", logFile)
		if extraArguments:
			job.addArguments(extraArguments)
		if commit:
			job.addArguments("-c")
		for inputFile in inputFileLs:
			job.addArguments(inputFile)
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(logFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = logFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	
	def addCountFastqReadBaseCountJob(self, workflow, executable=None, inputFile=None, \
								outputFile=None, isq_id=None, isqf_id=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.14
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("-i", inputFile, "-o", outputFile)
		if isq_id:
			job.addArguments("-q %s"%(isq_id))
		if isqf_id:
			job.addArguments("-f %s"%(isqf_id))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addJobs(self, workflow, inputData=None, pegasusFolderName=""):
		"""
		2012.3.14
		"""
		
		sys.stderr.write("Adding read counting jobs on %s input ..."%(len(inputData.jobDataLs)))
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		if topOutputDir:
			topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
			no_of_jobs += 1
		else:
			topOutputDirJob = None
		
		finalReduceFile = File(os.path.join(topOutputDir, 'read_base_count.tsv'))
		
		readBaseCountMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
						outputF=finalReduceFile, transferOutput=True, extraArguments=None, parentJobLs=[topOutputDirJob])
		
		logFile = os.path.join(topOutputDir, 'PutReadBaseCountIntoDB.log')
		putCountIntoDBJob = self.addPutReadBaseCountIntoDBJob(workflow, executable=workflow.PutReadBaseCountIntoDB, inputFileLs=[finalReduceFile], \
					logFile=logFile, commit=self.commit, parentJobLs=[readBaseCountMergeJob], extraDependentInputLs=[], transferOutput=True, \
					extraArguments=None, \
					job_max_memory=10)
		no_of_jobs += 2
		for jobData in inputData.jobDataLs:
			#add the read count job
			outputFile = os.path.join(topOutputDir, 'read_count_isq_%s_isqf_%s.tsv'%(jobData.isq_id, jobData.isqf_id))
			readCountJob = self.addCountFastqReadBaseCountJob(workflow, executable=workflow.CountFastqReadBaseCount, \
								inputFile=jobData.output, outputFile=outputFile, isq_id=jobData.isq_id, isqf_id=jobData.isqf_id, \
								parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
								job_max_memory=10)
			
			no_of_jobs += 1
			self.addInputToStatMergeJob(workflow, statMergeJob=readBaseCountMergeJob, \
								inputF=readCountJob.output, parentJobLs=[readCountJob])
			
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return putCountIntoDBJob
	
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
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerCustomJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		inputData = self.registerISQFiles(workflow=workflow, db_vervet=db_vervet, ind_seq_id_ls=self.ind_seq_id_ls, \
										localDataDir=self.localDataDir, pegasusFolderName=self.pegasusFolderName,\
										input_site_handler=self.input_site_handler)
		self.addJobs(workflow, inputData=inputData, pegasusFolderName=self.pegasusFolderName)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = ReadFileBaseCountWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()