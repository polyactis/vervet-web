#!/usr/bin/env python
"""
Examples:
	
	# 2011-8-16 on condorpool
	%s -o inspectBaseQuality.xml -u yh -i 1-8,15-130
	
	# 2011-8-16 use hoffman2 site_handler
	%s -o inspectBaseQuality.xml -u yh -i 1-8,15-130 
		-l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db
	
	
Description:
	2011-8-16
		construct a pegasus workflow to run InspectBaseQuality.py over a list of individual sequences
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


class InspectBaseQualityPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. non-fastq entries will be discarded.', ],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. If not given, use the default stored in db.'],\
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
		workflow = ADAG("InspectBaseQualityPipeline")
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
		
		inspectBaseQuality = Executable(namespace=namespace, name="InspectBaseQuality", version=version, \
								os=operatingSystem, arch=architecture, installed=True)
		inspectBaseQuality.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "InspectBaseQuality.py"), site_handler))
		workflow.addExecutable(inspectBaseQuality)
		
		#must use db_vervet.data_dir.
		# If self.dataDir differs from db_vervet.data_dir, this program (must be run on submission host) won't find files.
		individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, dataDir=db_vervet.data_dir)
		
		for ind_seq_id, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			individual_sequence = VervetDB.IndividualSequence.get(ind_seq_id)
			if individual_sequence is not None and individual_sequence.format=='fastq':
				#start to collect all files affiliated with this individual_sequence record 
				inputFnameLs = []
				for filePair in FilePairLs:
					for fileRecord in filePair:
						relativePath, format, sequence_type = fileRecord[:3]
						filename = os.path.join(self.dataDir, relativePath)
						inputFnameLs.append(filename)
				
				#create jobs
				for inputFname in inputFnameLs:
					prefix, suffix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(inputFname)
					if suffix=='.fastq':	#sometimes other files get in the way
						inspectBaseQuality_job = Job(namespace=namespace, name=inspectBaseQuality.name, version=version)
						inspectBaseQuality_job.addArguments("-i", inputFname, '-u', 'yh', '-e', '0.005', '-p', 'secret',\
													'-q', individual_sequence.quality_score_format)
						#samtools_index_job.uses(picard_output, transfer=False, register=True, link=Link.OUTPUT)	
						#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
						#samtools_index_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
						workflow.addJob(inspectBaseQuality_job)
						#workflow.addDependency(parent=picard_job, child=samtools_index_job)
				
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = InspectBaseQualityPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
