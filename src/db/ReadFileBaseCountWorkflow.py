#!/usr/bin/env python
"""
Examples:
	# 2012.5.3 run on hoffman2's condorpool, need sshDBTunnel (-H1)
	%s  -i 963-1346 -o dags/ReadCount/read_count_isq_936_1346.xml -u yh --commit -z localhost
		--pegasusFolderName readcount --needSSHDBTunnel
		-l hcondor -j hcondor 
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
	# 2012.3.14 
	%s -i 1-864 -o dags/ReadCount/read_count_isq_1_864.xml -u yh -l condorpool -j condorpool -z uclaOffice
		--pegasusFolderName readCount --commit
	
Description:
	2012.3.14
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus, AbstractNGSWorkflow
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow


class ReadFileBaseCountWorkflow(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						})

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.3.14
		"""
		if workflow is None:
			workflow = self
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'mapper/CountFastqReadBaseCount.py'), \
										name='CountFastqReadBaseCount', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/PutReadBaseCountIntoDB.py'), \
										name='PutReadBaseCountIntoDB', clusterSizeMultipler=0.2)
		
	
	def registerISQFiles(self, workflow=None, db_vervet=None, ind_seq_id_ls=[], local_data_dir='', pegasusFolderName='', \
						input_site_handler='local'):
		"""
		2012.3.14
		"""
		sys.stderr.write("Finding all ISQ-affiliated files of %s ind seq entries ..."%(len(ind_seq_id_ls)))
		returnData = PassingData(jobDataLs=[])
		counter = 0
		Table = VervetDB.IndividualSequence
		query = Table.query.filter(Table.id.in_(ind_seq_id_ls))
		individual_sequence_id_set = set()
		missed_individual_sequence_id_set = set()
		for individual_sequence in query:
			if individual_sequence.individual_sequence_file_ls:	#not empty
				for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
					absPath = os.path.join(local_data_dir, individual_sequence_file.path)
					if os.path.isfile(absPath):
						inputF = File(os.path.join(pegasusFolderName, individual_sequence_file.path))
						inputF.addPFN(PFN("file://" + absPath, input_site_handler))
						inputF.absPath = absPath
						workflow.addFile(inputF)
						returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
															isqf_id=individual_sequence_file.id))
						individual_sequence_id_set.add(individual_sequence.id)
					else:
						missed_individual_sequence_id_set.add(individual_sequence.id)
						sys.stderr.write("Warning: IndividualSequenceFile.id=%s (isq-id=%s) doesn't have any affiliated IndividualSequenceFile entries while its path %s is not a file.\n"%\
									(individual_sequence_file.id, individual_sequence.id, absPath))
			elif individual_sequence.path:
				absPath = os.path.join(local_data_dir, individual_sequence.path)
				if os.path.isfile(absPath):
					inputF = File(os.path.join(pegasusFolderName, individual_sequence.path))
					inputF.addPFN(PFN("file://" + absPath, input_site_handler))
					inputF.absPath = absPath
					workflow.addFile(inputF)
					returnData.jobDataLs.append(PassingData(output=inputF, jobLs=[], isq_id=individual_sequence.id,\
														isqf_id=None))
					individual_sequence_id_set.add(individual_sequence.id)
				else:
					sys.stderr.write("Warning: IndividualSequence.id=%s doesn't have any affiliated IndividualSequenceFile entries while its path %s is not a file.\n"%\
									(individual_sequence.id, absPath))
					missed_individual_sequence_id_set.add(individual_sequence.id)
		
		sys.stderr.write(" %s files registered for %s individual_sequence entries. missed %s individual-sequence entries.\n"%\
						(len(returnData.jobDataLs), len(individual_sequence_id_set), len(missed_individual_sequence_id_set)))
		return returnData
	
	def addPutReadBaseCountIntoDBJob(self, workflow, executable=None, inputFileLs=[], \
					logFile=None, commit=False, parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=1, **keywords):
		"""
		2012.5.3
			add argument sshDBTunnel
		2012.3.14
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments("--drivername", self.drivername, "--hostname", self.hostname, "--dbname", self.dbname, \
						"--db_user", self.db_user, "--db_passwd", self.db_passwd, \
						"--logFilename", logFile)
		if extraArguments:
			job.addArguments(extraArguments)
		if commit:
			job.addArguments("--commit")
		for inputFile in inputFileLs:
			job.addArguments(inputFile)
			job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(logFile, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = logFile
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel)
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
		job.addArguments("--inputFname", inputFile, "--outputFname", outputFile)
		if isq_id:
			job.addArguments("--isq_id %s"%(isq_id))
		if isqf_id:
			job.addArguments("--isqf_id %s"%(isqf_id))
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
	
	def addJobs(self, workflow, inputData=None, pegasusFolderName="", needSSHDBTunnel=0):
		"""
		2012.3.14
		"""
		
		sys.stderr.write("Adding read counting jobs on %s input ..."%(len(inputData.jobDataLs)))
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		if topOutputDir:
			topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
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
					job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
		no_of_jobs += 2
		for jobData in inputData.jobDataLs:
			#add the read count job
			outputFile = File(os.path.join(topOutputDir, 'read_count_isq_%s_isqf_%s.tsv'%(jobData.isq_id, jobData.isqf_id)))
			readCountJob = self.addCountFastqReadBaseCountJob(workflow, executable=workflow.CountFastqReadBaseCount, \
								inputFile=jobData.output, outputFile=outputFile, isq_id=jobData.isq_id, isqf_id=jobData.isqf_id, \
								parentJobLs=jobData.jobLs + [topOutputDirJob], extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
								job_max_memory=10)
			
			no_of_jobs += 1
			self.addInputToStatMergeJob(workflow, statMergeJob=readBaseCountMergeJob, \
								inputF=readCountJob.output, parentJobLs=[readCountJob])
			
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return putCountIntoDBJob
	
	def setup_run(self):
		"""
		2013.04.07 wrap all standard pre-run() related functions into this function.
			setting up for run(), called by run()
		"""
		pdata = AbstractNGSWorkflow.setup_run(self)
		workflow = pdata.workflow
		
		db_vervet = self.db_vervet
		session = db_vervet.session
		session.begin(subtransactions=True)
		"""
		Traceback (most recent call last):
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 249, in <module>
		    instance.run()
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 232, in run
		    pdata = self.setup_run()
		  File "/u/home/eeskin/polyacti/script/vervet/src/db/ReadFileBaseCountWorkflow.py", line 217, in setup_run
		    session.begin()
		  File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/scoping.py", line 139, in do
		    return getattr(self.registry(), name)(*args, **kwargs)
		  File "/u/home/eeskin/polyacti/lib/python/sqlalchemy/orm/session.py", line 550, in begin
		    "A transaction is already begun.  Use subtransactions=True "
		sqlalchemy.exc.InvalidRequestError: A transaction is already begun.  Use subtransactions=True to allow subtransactions.

		"""
		
		inputData = self.registerISQFiles(workflow=workflow, db_vervet=db_vervet, ind_seq_id_ls=self.ind_seq_id_ls, \
										local_data_dir=self.local_data_dir, pegasusFolderName=self.pegasusFolderName,\
										input_site_handler=self.input_site_handler)
		
		registerReferenceData = self.getReferenceSequence()
		
		
		return PassingData(workflow=workflow, inputData=inputData,\
						registerReferenceData=registerReferenceData)
	def run(self):
		"""
		2011-7-11
		"""
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		inputData=pdata.inputData
		
		self.addJobs(workflow, inputData=inputData, pegasusFolderName=self.pegasusFolderName,
					needSSHDBTunnel=self.needSSHDBTunnel)
		
		self.end_run()

if __name__ == '__main__':
	main_class = ReadFileBaseCountWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()