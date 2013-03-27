#!/usr/bin/env python
"""
Examples:
	# 2012.2.9 convert all existing 'fastq' IndividualSequence entries into new format (IndividualSequenceFile,)
	%s -o ConvertOldIsqRecordsToNewOnes.xml -u yh -p secret -c -z uclaOffice -j condorpool -l condorpool
	
	# 2012.2.9 convert only IndividualSequence entries whose id fall into 496-508 range and their format is 'fastq'. 
	%s -o ConvertOldIsqRecordsPipeline_12Development.xml -u yh -i 496-508 -c ...
	
Description:
	2012.2.8
		convert old IndividualSequence records into new ones (IndividualSequenceFile & IndividualSequence).
		This program is only meant for running once as it modifies the IndividualSequence.path
			so that next run of this program would fetch newly split files.
	1. Be careful with the db connection setting as it'll be passed to the db-registration job.
		Make sure all computing nodes have access to the db.
	2. The workflow has to be run on nodes where they have direct db and db-affiliated file-storage access.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow
from vervet.src.UnpackAndAddIndividualSequence2DB import UnpackAndAddIndividualSequence2DB
from vervet.src import VervetDB

class ConvertOldIsqRecordsPipeline(UnpackAndAddIndividualSequence2DB):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('isq_id_ls', 0, ): ['', 'i', 1, 'coma/dash-separated list of individual_sequence.id. used to restrict the input. \
							if not given, all format=fastq isq entries are input.', ],\
						('minNoOfReads', 1, int): [5000000, '', 1, 'minimum number of reads in each split fastq file. The upper limit in each split file is 2*minNoOfReads.', ],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (table individual_sequence_file, this argument gets passed to RegisterAndMoveSplitSequenceFiles.py)'],\
						})
	
	def __init__(self,  **keywords):
		"""
		2012.2.8
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		if self.isq_id_ls:
			self.isq_id_ls = getListOutOfStr(self.isq_id_ls, data_type=int)
	
	
	def fetchIndividualSequenceFromDB(self, db_vervet, isq_id_ls=None, format='fastq'):
		"""
		2012.2.8
		"""
		sys.stderr.write("Fetching list of individual_sequence ids from db ...")
		Table = VervetDB.IndividualSequence
		query = Table.query.filter_by(format=format)
		isq_ls = []
		if isq_id_ls:
			query = query.filter(Table.id.in_(isq_id_ls))
		for row in query:
			isq_ls.append(row)
		sys.stderr.write("%s db entries.\n"%(len(isq_ls)))
		return isq_ls
	
	def parseLibraryMateIDFromFilename(self, filename):
		"""
		2012.2.8
			the filename looks like gerald_64J6AAAXX_1_1.fastq.gz
			'gerald_64J6AAAXX_1' is the library and the last '1' is the mate id.
		"""
		import re
		library_mate_id_pattern = re.compile(r'(\w+?)(_([12])|).fastq')	#the library searching is non-greedy.
		#mate_id might be missing (454 date, single end).
		library_mate_id_search_result = library_mate_id_pattern.search(filename)
		if library_mate_id_search_result:
			library = library_mate_id_search_result.group(1)
			mate_id = library_mate_id_search_result.group(3)
			if mate_id:	#it could be None because it's optional in the pattern.
				mate_id = int(mate_id)
		else:
			library = None
			mate_id = None
		return [library, mate_id]
	
	def run(self):
		"""
		2012.2.8
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, db_user=self.db_user,
					db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		session = db_vervet.session
		session.begin()
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		isq_ls = self.fetchIndividualSequenceFromDB(db_vervet, self.isq_id_ls)
		
		no_of_jobs = 1
		
		individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs([isq.id for isq in isq_ls], \
													data_dir=self.local_data_dir, checkOldPath=True)
		
		for individualSequenceID, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			individual_sequence = VervetDB.IndividualSequence.get(individualSequenceID)
			
			newISQPath = individual_sequence.constructRelativePathForIndividualSequence()
			#newISQPath = '%s_split'%(newISQPath)
			if individual_sequence.path !=newISQPath:
				individual_sequence.path = newISQPath
				session.add(individual_sequence)
				session.flush()
			
			sequenceOutputDir = os.path.join(self.data_dir, individual_sequence.path)
			sequenceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sequenceOutputDir)
			
			for filePair in FilePairLs:
				for fileRecord in filePair:
					filename = fileRecord[0]
					absPath = os.path.join(self.local_data_dir, filename)
					fastqFile = self.registerOneInputFile(workflow, absPath)
					library, mate_id = self.parseLibraryMateIDFromFilename(filename)[:2]
					if library is None:
						sys.stderr.write("Warning: can't parse library out of file %s of isq %s & skip.\n"%(filename, individualSequenceID))
						continue
					
					if mate_id:
						prefix = '%s_%s'%(library, mate_id)
					else:
						prefix = library
					
					outputFilenamePrefix = '%s_%s'%(individual_sequence.id, prefix)
					
					splitOutputDir = outputFilenamePrefix
					splitOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=splitOutputDir)
					
					splitFastQFnamePrefix = os.path.join(splitOutputDir, outputFilenamePrefix)
					
					logFile = File('%s_%s.split.log'%(individual_sequence.id, prefix))
					splitReadFileJob1 = self.addSplitReadFileJob(workflow, executable=workflow.splitReadFile, \
									inputF=fastqFile, outputFnamePrefix=splitFastQFnamePrefix, \
									outputFnamePrefixTail="", minNoOfReads=self.minNoOfReads, \
									logFile=logFile, parentJobLs=[splitOutputDirJob], \
									job_max_memory=2000, walltime = 800, \
									extraDependentInputLs=[], transferOutput=True)
					
					logFile = File('%s_%s.register.log'%(individual_sequence.id, prefix))
					registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.registerAndMoveSplitSequenceFiles, \
									inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
									individual_sequence_id=individual_sequence.id, bamFile=None, library=library, mate_id=mate_id, \
									parentJobLs=[splitReadFileJob1, sequenceOutputDirJob], job_max_memory=100, walltime = 60, \
									commit=self.commit, extraDependentInputLs=[], \
									transferOutput=True)
					
					no_of_jobs += 3
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = ConvertOldIsqRecordsPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
