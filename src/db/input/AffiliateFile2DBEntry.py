#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i  OneLibAlignment/6246_depth.tsv.gz
		--logFilename  OneLibAlignment/6246_depth_2db.log
		--db_entry_id 6246 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--drivername postgresql --hostname localhost --dbname vervetdb --db_user yh
		--schema public --tableClassName IndividualAlignment
		--filePathColumnName path_to_depth_file
		--fileSizeColumnName depth_file_size
		--outputFileRelativePath individual_alignment/6246_1517_VEB1010_GA_vs_3488_by_method6_realigned1_reduced0_p4475_m227_depth.tsv.gz

Description:
	2013.08.08 this program would affiliate a file to an existing database table entry.
		store its path in db,  calculate file size, etc.
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AffiliateFile2DBEntry(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						#('inputDir', 1, ): ['', 'i', 1, 'input folder that contains split fastq files', ],\
						('db_entry_id', 1, int):[None, '', 1, 'ID of the db entry with which the input file is to be affiliated with'],\
						('tableClassName', 1, ):[None, '', 1, 'table class name'],\
						('filePathColumnName', 1, ):[None, '', 1, 'column name to store file path'],\
						('fileSizeColumnName', 0, ):[None, '', 1, 'column name to store file size'],\
						('outputFileRelativePath', 1, ):[None, '', 1, 'the relative path of the output file, to be attached to self.data_dir'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		session.begin()
		
		if not self.data_dir:
			self.data_dir = self.db_vervet.data_dir
		data_dir = self.data_dir
		
		inputFileRealPath = os.path.realpath(self.inputFname)
		logMessage = "Affiliating file %s to db table %s, entry id=%s .\n"%(self.inputFname, self.tableClassName,\
																		self.db_entry_id)
		
		if os.path.isfile(self.inputFname):
			TableClass = getattr(VervetDB, self.tableClassName, None)
			db_entry = TableClass.get(self.db_entry_id)
			
			existing_db_entry_file_path = getattr(db_entry, self.filePathColumnName, None)
			if db_entry and existing_db_entry_file_path and os.path.isfile(os.path.join(data_dir, existing_db_entry_file_path)):
				sys.stderr.write("ERROR (and non-zero exit): another file %s is already in db.\n"%\
								(existing_db_entry_file_path))
				self.sessionRollback(session)
				self.cleanUpAndExitOnFailure(exitCode=3)
			
			try:
				#make sure these files are stored in self.dstFilenameLs and self.srcFilenameLs
				#copy further files if there are

				logMessage = self.db_vervet.copyFileWithAnotherFilePrefix(inputFname=self.inputFname, \
												outputFileRelativePath=self.outputFileRelativePath, \
												outputDir=self.data_dir,\
												logMessage=logMessage, srcFilenameLs=self.srcFilenameLs, \
												dstFilenameLs=self.dstFilenameLs)
				setattr(db_entry, self.filePathColumnName, self.outputFileRelativePath)
				session.add(db_entry)
				session.flush()
				self.db_vervet.updateDBEntryPathFileSize(db_entry=db_entry, data_dir=self.data_dir, absPath=None, \
								file_path_column_name=self.filePathColumnName, \
								file_size_column_name=self.fileSizeColumnName)
				
				## 2012.7.17 commented out because md5sum is calculated above
				#db_vervet.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=data_dir)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				self.sessionRollback(session)
				self.cleanUpAndExitOnFailure(exitCode=5)
		else:
			logMessage += "%s doesn't exist.\n"%(inputFileRealPath)
		self.outputLogMessage(logMessage)
		
		if self.commit:
			try:
				session.flush()
				session.commit()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				self.cleanUpAndExitOnFailure(exitCode=3)
		else:
			#delete all target files but exit gracefully (exit 0)
			self.sessionRollback(session)
			self.cleanUpAndExitOnFailure(exitCode=0)

if __name__ == '__main__':
	main_class = AffiliateFile2DBEntry
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()