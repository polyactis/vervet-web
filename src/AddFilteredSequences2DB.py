#!/usr/bin/env python
"""
Examples:
	%s -n 130 -u yh 1.fastq.gz 2.fastq.gz
	
	%s 

Description:
	2011-8-18
		follow-up to FilterShortRead.py or picard/dist/FilterRead.jar to add filtered sequences to db
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

import random
import VervetDB
from pymodule import PassingData, ProcessOptions, utils, yh_matplotlib

class AddFilteredSequences2DB(object):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, 'o', 1, 'database port number. must be non-empty if need ssh tunnel'],\
						('sshTunnelCredential', 0, ): ['', 's', 1, 'a ssh credential to allow machine to access db server. \
										polyacti@login3, yuhuang@hpc-login2. if empty or port is empty, no tunnel', ],\
						('parent_individual_sequence_id', 1, int): [None, 'n', 1, 'The individual_sequence id of pre-filter sequences.', ],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. If not given, use the default stored in db.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path) and qsub jobs'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.inputFnameLs = inputFnameLs
	
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		if self.sshTunnelCredential and self.port:
			try:
				utils.sshTunnel(serverHostname=self.hostname, port=self.port, middleManCredential=self.sshTunnelCredential)
				self.hostname = 'localhost'	#now the tunnel is built, access is on localhost
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, hostname=self.hostname, database=self.dbname, schema=self.schema,\
						username=self.db_user, password=self.db_passwd,	port=self.port)
		db_vervet.setup(create_tables=False)
		session = db_vervet.session
		session.begin()	#no transaction for input node as there is no data insertion
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		parent_individual_sequence = VervetDB.IndividualSequence.get(self.parent_individual_sequence_id)
		
		if not parent_individual_sequence:
			sys.stderr.write("parent individual sequence %s doesn't exist in db.\n")
			sys.exit(4)
		
		individual_sequence = db_vervet.copyParentIndividualSequence(parent_individual_sequence=parent_individual_sequence)
		
		outputDir = os.path.join(self.dataDir, individual_sequence.path)
		if not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		for inputFname in self.inputFnameLs:
			commandline = 'mv %s %s/'%(inputFname, outputDir)
			if self.commit:	#qsub only when db transaction will be committed.
				return_data = utils.runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				if return_data.stderr_content:
					sys.stderr.write("Error in moving file. Exit now.\n")
					sys.exit(3)
				
		if self.commit:
			session.commit()
		else:
			session.rollback()
			
			

if __name__ == '__main__':
	main_class = AddFilteredSequences2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
