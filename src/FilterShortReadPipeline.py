#!/usr/bin/env python
"""
Examples:
	
	# 2011-8-16 on condorpool.
	%s -o filterShortReadPipeline.xml -u yh -i 1-507 -z uclaOffice -j condorpool -l condorpool
	
	# 2011-8-16 use hoffman2 site_handler. watch additional arguments for tunnel setup
	%s -o filterShortReadPipeline.xml -u yh -i 1-8,15-130 
		-l hoffman2 -e /u/home/eeskin/polyacti
		-s polyacti@login3 -a 5432 -z dl324b-1.cmb.usc.edu -t /u/home/eeskin/polyacti/NetworkData/vervet/db
	#2012.2.15
	%s -z uclaOffice -l condorpool -c -j condorpool  -o workflow/FilterReadPipeline_isq_id_527_626.xml  -i 527-626 -u yh
	
Description:
	2011-8-16
		script that generates a read-filtering pegasus workflow.
		It uses picard_path/FilterRead.jar and AddFilteredSequences2DB.py.
	
	1. Be careful with the db connection setting as it'll be passed to the db-registration job.
		Make sure all computing nodes have access to the db.
	2. The workflow has to be run on nodes where they have direct db and db-affiliated file-storage access. 
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class FilterShortReadPipeline(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update(
						{
						('sshTunnelCredential', 0, ): ['', 's', 1, 'a ssh credential to allow machine to access db server. \
										polyacti@login3, yuhuang@hpc-login2. if empty or port is empty, no tunnel', ],\
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									no-individual_sequence_file-affiliated entries will be discarded.', ],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						})
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		else:
			self.ind_seq_id_ls = []
	
	def registerCustomExecutables(self, workflow):
		"""
		2012.1.3
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		addFilteredSequences2DB = Executable(namespace=namespace, name="AddFilteredSequences2DB", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		addFilteredSequences2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AddFilteredSequences2DB.py"), site_handler))
		workflow.addExecutable(addFilteredSequences2DB)
		workflow.addFilteredSequences2DB = addFilteredSequences2DB
		
		"""2011-8-31 replace filterReadJar
		filterShortRead = Executable(namespace=namespace, name="FilterShortRead", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		filterShortRead.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "FilterShortRead.py"), site_handler))
		workflow.addExecutable(filterShortRead)
		"""
	
	def registerCustomJars(self, workflow, ):
		"""
		2012.2.10
		"""
		site_handler = self.site_handler
		
		abs_path = os.path.join(self.picard_path, 'FilterRead.jar')
		filterReadJar = File(abs_path)
		filterReadJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(filterReadJar)
		workflow.filterReadJar = filterReadJar
		
	
	def addAddFilteredSequences2DB_job(self, workflow, executable=None, \
							inputFile=None, individual_sequence_id=None, outputDir=None, logFile=None,\
							parent_individual_sequence_file_id=None, \
							parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
							commit=0, \
							extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.2.10
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-v', self.drivername, '-z', self.hostname, '-d', self.dbname, '-k', self.schema,\
						'-u', self.db_user, '-p', self.db_passwd, \
						'-i', inputFile, '-n', str(individual_sequence_id), '-o', outputDir, \
						'-e %s'%(parent_individual_sequence_file_id))
		if commit:
			job.addArguments("-c")
		if self.port:
			job.addArguments("--port=%s"%self.port)
		if self.sshTunnelCredential:
			job.addArguments("--sshTunnelCredential=%s"%(self.sshTunnelCredential))
		
		job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
		if logFile:
			job.addArguments('-g', logFile)
			job.uses(logFile, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
			job.output = logFile
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	
	def addFilterReadJob(self, workflow, executable=None, jar=None,\
					parentJobLs=[], job_max_memory=2000, job_max_walltime = 120, \
					extraDependentInputLs=[], \
					transferOutput=False, **keywords):
		"""
		2012.2.9
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments(javaMemRequirement, '-jar', jar,\
						'M=20', 'N=50', 'a=0.3', 'l=2')
		"""
		m = "Minimum phred score. any base with phred score below this number will be turned into N")
	
		@Option(shortName = "n", doc = "final read (after all filtering) must have length>=this number")
	
		@Option(shortName = "a", doc = "any read with Ns more than this percentage will be discarded")
		
		@Option(shortName = "l", doc = "during head/tail trimming, smoothing phred score is applied. amounts to how many flanking bases used"
		"""
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def getLibrarySplitOrder2DBEntryLs(self, individual_sequence):
		"""
		2012.2.10
			generate a dictionary for this program to check if some (library, split_order)s have been filtered and recorded in db already
		"""
		library_split_order2db_entry_ls = {}
		for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
			key = (individual_sequence_file.library, individual_sequence_file.split_order)
			if key not in library_split_order2db_entry_ls:
				library_split_order2db_entry_ls[key] = []
			library_split_order2db_entry_ls[key].append(individual_sequence_file)
		return library_split_order2db_entry_ls
	
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
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerCustomJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		isq_id2LibrarySplitOrder2FileLs = db_vervet.getISQ_ID2LibrarySplitOrder2FileLs(self.ind_seq_id_ls, dataDir=self.dataDir, filtered=0)
		no_of_jobs = 0
		for ind_seq_id, LibrarySplitOrder2FileLs in isq_id2LibrarySplitOrder2FileLs.iteritems():
			parent_individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if parent_individual_sequence is not None and parent_individual_sequence.format=='fastq':
				"""
				check if the child individual_sequence already exists in db or not. if it does, what about its files?? if not, go add filtering jobs.
				"""
				
				individual_sequence = db_vervet.getIndividualSequence(individual_id=parent_individual_sequence.individual_id, \
						sequencer=parent_individual_sequence.sequencer, sequence_type=parent_individual_sequence.sequence_type,\
						sequence_format=parent_individual_sequence.format, path_to_original_sequence=None, tissue_name=None, coverage=None,\
						quality_score_format='Standard', filtered=1,\
						parent_individual_sequence_id=parent_individual_sequence.id)
				
				library_split_order2filtered_db_entry_ls = self.getLibrarySplitOrder2DBEntryLs(individual_sequence)
				
				sequenceOutputDir = os.path.join(self.dataDir, individual_sequence.path)
				sequenceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sequenceOutputDir)
				no_of_jobs += 1
				
				filteredReadOutputDir = os.path.join(os.path.basename(individual_sequence.path))
				filteredReadOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=filteredReadOutputDir)
				no_of_jobs += 1
				for key, fileObjLs in LibrarySplitOrder2FileLs.iteritems():
					if key in library_split_order2filtered_db_entry_ls:
						sys.stderr.write("Warning: this pair of filtered individual_sequence_file(s), %s, have been in db already. skip.\n"%(repr(key)))
						continue
					library, split_order = key[:2]
					
					#add filter jobs
					filterShortRead_job = self.addFilterReadJob(workflow, executable=workflow.java, jar=workflow.filterReadJar,\
						parentJobLs=[filteredReadOutputDirJob], job_max_memory=2000, job_max_walltime = 120, \
						extraDependentInputLs=[], transferOutput=False)
					no_of_jobs += 1
					for i in xrange(len(fileObjLs)):
						fileObj = fileObjLs[i]
						inputFile = self.registerOneInputFile(workflow, fileObj.path)
						outputFname = os.path.join(filteredReadOutputDir, os.path.basename(fileObj.path))
						outputFile = File(outputFname)		#take the base filename as the output filename. it'll be in scratch/.
						if i==0:	#1st mate	#also add the quality_score_format
							filterShortRead_job.addArguments('V=%s'%fileObj.db_entry.quality_score_format)
							filterShortRead_job.addArguments("I=", inputFile, 'O=', outputFile)
						elif i==1:	#2nd mate
							filterShortRead_job.addArguments("J=", inputFile, 'P=', outputFile)
						else:
							sys.stderr.write("Error: mate %s appeared in paired-end data (individualSequenceID=%s).\n"%(i+1, ind_seq_id))
							sys.exit(4)
						filterShortRead_job.uses(inputFile, transfer=True, register=True, link=Link.INPUT)
						filterShortRead_job.uses(outputFile, transfer=False, register=True, link=Link.OUTPUT)
						
						logFile = File('%s_%s.register.log'%(individual_sequence.id, fileObj.db_entry.id))
						addFilteredSequences2DB_job = self.addAddFilteredSequences2DB_job(workflow, executable=workflow.addFilteredSequences2DB, \
									inputFile=outputFile, individual_sequence_id=individual_sequence.id, outputDir=sequenceOutputDir, \
									logFile=logFile, \
									parent_individual_sequence_file_id=fileObj.db_entry.id,\
									parentJobLs=[sequenceOutputDirJob, filterShortRead_job], commit=self.commit, \
									extraDependentInputLs=[], transferOutput=True)
						no_of_jobs += 1
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = FilterShortReadPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
