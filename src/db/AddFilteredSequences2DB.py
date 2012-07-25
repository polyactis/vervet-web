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
from pymodule import PassingData, ProcessOptions, utils, yh_matplotlib
from pymodule.AbstractDBInteractingClass import AbstractDBInteractingClass
from mapper.RegisterAndMoveSplitSequenceFiles import RegisterAndMoveSplitSequenceFiles
from vervet.src import VervetDB	#have to import it from vervet.src, not directly cuz that's how it's imported in RegisterAndMoveSplitSequenceFiles


class AddFilteredSequences2DB(RegisterAndMoveSplitSequenceFiles):
	__doc__ = __doc__
	option_default_dict = AbstractDBInteractingClass.option_default_dict.copy()
	option_default_dict.update({
						('inputFname', 1, ): ['', 'i', 1, 'the filtered fastq files', ],\
						('individual_sequence_id', 1, int): [None, 'n', 1, 'The individual_sequence id of the input sequence file.', ],\
						('parent_individual_sequence_file_id', 1, int): [None, 'e', 1, 'ID of the parent of this filtered individual_sequence_file' ],\
						('outputDir', 1, ): ['', 'o', 1, 'output folder to which split files from inputDir will be moved', ],\
						('logFilename', 0, ): [None, 'g', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow'],\
						})

	def __init__(self, **keywords):
		"""
		2011-7-11
		"""
		RegisterAndMoveSplitSequenceFiles.__init__(self, **keywords)
	
	
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
		
		individual_sequence = VervetDB.IndividualSequence.get(self.individual_sequence_id)
		
		if not os.path.isdir(self.outputDir):
			os.makedirs(self.outputDir)
		
		parent_individual_sequence_file = VervetDB.IndividualSequenceFile.get(self.parent_individual_sequence_file_id)
		if not parent_individual_sequence_file:
			sys.stderr.write("parent individual_sequence_file %s doesn't exist in db.\n"%(self.parent_individual_sequence_file_id))
			sys.exit(4)
		
		db_entry = db_vervet.copyParentIndividualSequenceFile(parent_individual_sequence_file=parent_individual_sequence_file,\
									individual_sequence_id=self.individual_sequence_id, quality_score_format='Standard', filtered=1)
		if db_entry:
			#move the file
			inputDir, filename = os.path.split(self.inputFname)
			exitCode = db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=filename, \
													inputDir=inputDir, outputDir=self.outputDir, \
								relativeOutputDir=individual_sequence.path, shellCommand='cp -rL', \
								srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
								constructRelativePathFunction=None)
			"""
			#2012.7.13 old way
			exitCode = self.moveNewISQFileIntoDBStorage(session, individual_sequence_file=db_entry, filename=filename, inputDir=inputDir, \
										outputDir=self.outputDir, relativeOutputDir=individual_sequence.path,\
										shellCommand='cp', srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs)
			"""
			if exitCode!=0:
				sys.stderr.write("Error: moveNewISQFileIntoDBStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				#delete all recorded target files
				self.rmGivenFiles(filenameLs=self.dstFilenameLs)
				sys.exit(exitCode)
		else:
			sys.stderr.write("Error: IndividualSequenceFile db entry is None.\n")
			sys.exit(3)
		
		if self.logFilename:
			outf = open(self.logFilename, 'w')
			outf.write("file %s was added into db.\n"%(self.inputFname))
			outf.close()
		
		if self.commit:
			try:
				session.commit()
				#delete all source files
				self.rmGivenFiles(filenameLs=self.srcFilenameLs)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				#delete all recorded target files
				self.rmGivenFiles(filenameLs=self.dstFilenameLs)
				sys.exit(3)
		else:
			session.rollback()
			#delete all target files
			self.rmGivenFiles(filenameLs=self.dstFilenameLs)
			

if __name__ == '__main__':
	main_class = AddFilteredSequences2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
