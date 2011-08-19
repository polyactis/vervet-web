#!/usr/bin/env python
"""
Examples:
	
	# 2011-8-16 on condorpool
	%s -o filterShortReadPipeline.xml -u yh -i 1-8,15-130
	
	# 2011-8-16 use hoffman2 site_handler. watch additional arguments for tunnel setup
	%s -o filterShortReadPipeline.xml -u yh -i 1-8,15-130 
		-l hoffman2 -e /u/home/eeskin/polyacti
		-s polyacti@login3 -a 5432 -z dl324b-1.cmb.usc.edu -t ~/NetworkData/vervet/db
	
	
Description:
	2011-8-16
		script that generates a read-filtering pegasus workflow.
		It uses picard_path/FilterRead.jar and AddFilteredSequences2DB.py. 
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
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from Pegasus.DAX3 import *


class FilterShortReadPipeline(object):
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
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
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
		
		self.picard_path = self.picard_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
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
		
		# Create a abstract dag
		workflow = ADAG("FilterShortRead")
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		
		mkdir = Executable(namespace=namespace, name="mkdir", version=version, os=operatingSystem, arch=architecture, installed=True)
		mkdir.addPFN(PFN("file://" + '/bin/mkdir', site_handler))
		workflow.addExecutable(mkdir)
		
		filterShortRead = Executable(namespace=namespace, name="FilterShortRead", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		filterShortRead.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "FilterShortRead.py"), site_handler))
		workflow.addExecutable(filterShortRead)
		
		addFilteredSequences2DB = Executable(namespace=namespace, name="AddFilteredSequences2DB", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		addFilteredSequences2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AddFilteredSequences2DB.py"), site_handler))
		workflow.addExecutable(addFilteredSequences2DB)
		
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(java)
		
		for ind_seq_id in self.ind_seq_id_ls:
			individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if individual_sequence is not None and individual_sequence.format=='fastq':
				inputDir = os.path.join(self.dataDir, individual_sequence.path)
				
				#start to collect all files affiliated with this individual_sequence record 
				inputFnameLs = []
				if os.path.isfile(inputDir):	#it's a file already
					inputFname = inputDir
					inputFnameLs.append(inputFname)
				elif os.path.isdir(inputDir):
					files = os.listdir(inputDir)
					for fname in files:
						inputFname = os.path.join(inputDir, fname)
						inputFnameLs.append(inputFname)
				
				#create filter jobs
				
				#create this db job first.
				addFilteredSequences2DB_job = Job(namespace=namespace, name=addFilteredSequences2DB.name, version=version)
				addFilteredSequences2DB_job.addArguments('-v', self.drivername, '-z', self.hostname, '-d', self.dbname, '-k', self.schema,\
										'-u', self.db_user, '-p', self.db_passwd, \
										'-n', str(individual_sequence.id), "-c")
				#add optional parameters
				if self.port:
					addFilteredSequences2DB_job.addArguments("-o", self.port)
				if self.sshTunnelCredential:
					addFilteredSequences2DB_job.addArguments("-s", self.sshTunnelCredential)
				if self.dataDir:
					addFilteredSequences2DB_job.addArguments("-t", self.dataDir)
				
				workflow.addJob(addFilteredSequences2DB_job)
				
				filterShortRead_job_ls= []
				for inputFname in inputFnameLs:
					prefix, suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(inputFname)
					if suffix=='.fastq':	#sometimes other files get in the way
						outputFname = os.path.split(inputFname)[1]
						outputFile = File(outputFname)		#take the base filename as the output filename. it'll be in scratch/.
						filterShortRead_job = Job(namespace=namespace, name=java.name, version=version)
						filterShortRead_job.addArguments('-jar', os.path.join(self.picard_path, 'FilterRead.jar'),\
										"I=%s"%inputFname, 'V=%s'%individual_sequence.quality_score_format,\
										'O=%s'%outputFname)
						#samtools_index_job.uses(picard_output, transfer=False, register=True, link=Link.OUTPUT)	
						#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
						#samtools_index_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
						workflow.addJob(filterShortRead_job)
						filterShortRead_job_ls.append(filterShortRead_job)
						addFilteredSequences2DB_job.addArguments(outputFile)
						workflow.addDependency(parent=filterShortRead_job, child=addFilteredSequences2DB_job)
				
			
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = FilterShortReadPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
