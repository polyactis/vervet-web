#!/usr/bin/env python
"""
Examples:
	
	# 2011-8-16 on condorpool.
	%s -o filterShortReadPipeline.xml -u yh -i 1-507 -z uclaOffice -j condorpool -l condorpool --commit
	
	# 2011-8-16 use hoffman2 site_handler. needs ssh tunnel for db-access jobs (-H), always commit (--commit) otherwise, no records in IndividualSequence
	# make job cluster size=50 (-C 50)
	%s -o dags/FilterReads/FilterShortReadPipeline.xml -i 1-8,15-130 
		-z localhost --commit -u yh
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java
		-H -C 50
	
	#2012.2.15
	%s -z uclaOffice --commit -l condorpool -j condorpool  -o dags/FilterReadPipeline_isq_id_527_626.xml 
		-i 527-626 -u yh
	
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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus import yh_pegasus
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow

class FilterShortReadPipeline(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update(
						{
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									no-individual_sequence_file-affiliated entries will be discarded.', ],\
						})
	# 2012.6.8 skipFilteredSequenceFiles is automatically enforced by checking library_split_order2filtered_db_entry_ls.
		# 
	#('skipFilteredSequenceFiles', 0, int):[0, 'K', 0, 'skip individual_sequence_file entries whose file entry already exists'],\
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		else:
			self.ind_seq_id_ls = []
	
	def registerCustomExecutables(self, workflow):
		"""
		2012.1.3
		"""
		vervetSrcPath = self.vervetSrcPath
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "db/input/AddFilteredSequences2DB.py"), \
										name='AddFilteredSequences2DB', clusterSizeMultipler=0.5)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
										name='FilterReadJava', clusterSizeMultipler=0.6)
		#2012.7.13 don't cluster add-to-DB jobs because if one of them fail, the previous ones will be re-run.
		#workflow.java.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		
		"""2011-8-31 replace FilterReadJar
		filterShortRead = Executable(namespace=namespace, name="FilterShortRead", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		filterShortRead.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "FilterShortRead.py"), site_handler))
		workflow.addExecutable(filterShortRead)
		"""
	
	def registerCustomJars(self, workflow, ):
		"""
		2012.2.10
		"""
		
		self.registerOneJar(name="FilterReadJar", path=os.path.join(self.picard_path, 'FilterRead.jar'))
	
	def addAddFilteredSequences2DB_job(self, workflow=None, executable=None, \
							inputFile=None, individual_sequence_id=None, outputDir=None, logFile=None,\
							parent_individual_sequence_file_id=None, \
							parentJobLs=None, job_max_memory=100, walltime = 60, \
							commit=0, \
							extraDependentInputLs=None, \
							transferOutput=False, sshDBTunnel=1, **keywords):
		"""
		2013.04.05 call addGenericFile2DBJob
		2012.4.18
			add argument sshDBTunnel
		2012.2.10
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		extraArgumentList = ["--individual_sequence_id %s"%(individual_sequence_id), '--outputDir %s'%(outputDir), \
						'--parent_individual_sequence_file_id %s'%(parent_individual_sequence_file_id)]
		job = self.addGenericFile2DBJob(executable=executable, inputFile=inputFile, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", inputFileList=None, \
					data_dir=None, logFile=logFile, commit=commit,\
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, \
					job_max_memory=job_max_memory,  sshDBTunnel=sshDBTunnel, walltime=walltime,\
					key2ObjectForJob=None, objectWithDBArguments=self, **keywords)
		return job
	
	
	def addFilterReadJob(self, workflow=None, executable=None, jar=None,\
					parentJobLs=None, job_max_memory=2000, walltime = 120, \
					extraDependentInputLs=None, \
					transferOutput=False, **keywords):
		"""
		2013.04.05 use addGenericJavaJob()
		2012.2.9
			walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		
		extraArgumentList = ['M=20', 'N=50', 'a=0.3', 'l=2']
		"""
		m = "Minimum phred score. any base with phred score below this number will be turned into N")
	
		@Option(shortName = "n", doc = "final read (after all filtering) must have length>=this number")
	
		@Option(shortName = "a", doc = "any read with Ns more than this percentage will be discarded")
		
		@Option(shortName = "l", doc = "during head/tail trimming, smoothing phred score is applied. amounts to how many flanking bases used"
		"""
		
		job = self.addGenericJavaJob(executable=executable, jarFile=jar, \
					inputFile=None, inputArgumentOption=None, \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=None, outputArgumentOption=None,\
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=None, \
					extraArguments=None, extraArgumentList=extraArgumentList, extraOutputLs=None, \
					extraDependentInputLs=extraDependentInputLs, \
					no_of_cpus=None, walltime=walltime, sshDBTunnel=None, **keywords)
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
		
		db_vervet = self.db_vervet
		session = db_vervet.session
		session.begin()
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerCustomJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		isq_id2LibrarySplitOrder2FileLs = db_vervet.getISQ_ID2LibrarySplitOrder2FileLs(self.ind_seq_id_ls, data_dir=self.data_dir, \
													filtered=0, ignoreEmptyReadFile=False)	#2012.6.1 unfiltered read file shoudn't be empty
		to_work_ind_seq_id_set = set()
		parent_individual_sequence_file_id_set = set()
		for ind_seq_id, LibrarySplitOrder2FileLs in isq_id2LibrarySplitOrder2FileLs.iteritems():
			parent_individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if parent_individual_sequence is not None and parent_individual_sequence.format=='fastq':
				"""
				check if the child individual_sequence already exists in db or not. if it does, what about its files?? if not, go add filtering jobs.
				"""
				#2012.6.8
				individual_sequence = db_vervet.copyParentIndividualSequence(parent_individual_sequence=parent_individual_sequence, \
									parent_individual_sequence_id=ind_seq_id,\
									quality_score_format='Standard', filtered=1, data_dir=self.data_dir)
				"""
				# 2012.6.8 use db_vervet.copyParentIndividualSequence() instead.
				individual_sequence = db_vervet.getIndividualSequence(individual_id=parent_individual_sequence.individual_id, \
						sequencer=parent_individual_sequence.sequencer, sequence_type=parent_individual_sequence.sequence_type,\
						sequence_format=parent_individual_sequence.format, path_to_original_sequence=None, tissue_name=None, coverage=None,\
						quality_score_format='Standard', filtered=1,\
						parent_individual_sequence_id=parent_individual_sequence.id, data_dir=self.data_dir)
				"""
				library_split_order2filtered_db_entry_ls = self.getLibrarySplitOrder2DBEntryLs(individual_sequence)
				
				sequenceOutputDirJob = None
				filteredReadOutputDirJob = None
				for key, fileObjLs in LibrarySplitOrder2FileLs.iteritems():
					if key in library_split_order2filtered_db_entry_ls:
						sys.stderr.write("Warning: this pair of filtered individual_sequence_file(s), %s, parent_individual_sequence (id=%s, %s),\
			individual_sequence (id=%s, %s) are already in db. skip.\n"%\
										(repr(key), parent_individual_sequence.id, parent_individual_sequence.individual.code,\
										individual_sequence.id, individual_sequence.individual.code))
						continue
					else:
						if sequenceOutputDirJob is None:
							sequenceOutputDir = os.path.join(self.data_dir, individual_sequence.path)
							sequenceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sequenceOutputDir)
						if filteredReadOutputDirJob is None:
							filteredReadOutputDir = os.path.join(os.path.basename(individual_sequence.path))
							filteredReadOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=filteredReadOutputDir)
						
					library, split_order = key[:2]
					
					#add filter jobs
					filterShortRead_job = self.addFilterReadJob(executable=self.FilterReadJava, jar=workflow.FilterReadJar,\
						parentJobLs=[filteredReadOutputDirJob], job_max_memory=2000, walltime = 120, \
						extraDependentInputLs=None, transferOutput=False)
					for i in xrange(len(fileObjLs)):
						fileObj = fileObjLs[i]
						try:	#2012.7.2
							inputFile = self.registerOneInputFile(workflow, inputFname=fileObj.path, folderName='inputIndividualSequenceFile')
						except:
							import pdb
							pdb.set_trace()
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
						addFilteredSequences2DB_job = self.addAddFilteredSequences2DB_job(workflow, \
									executable=workflow.AddFilteredSequences2DB, \
									inputFile=outputFile, individual_sequence_id=individual_sequence.id, outputDir=sequenceOutputDir, \
									logFile=logFile, \
									parent_individual_sequence_file_id=fileObj.db_entry.id,\
									parentJobLs=[sequenceOutputDirJob, filterShortRead_job], commit=self.commit, \
									extraDependentInputLs=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
						to_work_ind_seq_id_set.add(ind_seq_id)
						parent_individual_sequence_file_id_set.add(fileObj.db_entry.id)
		sys.stderr.write("%s jobs, %s individual_sequence entries, %s parent_individual_sequence_file_id s.\n"%\
						(self.no_of_jobs, len(to_work_ind_seq_id_set), len(parent_individual_sequence_file_id_set)))
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
