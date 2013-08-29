#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s --db_user yh --schema public -i  folderReduce/63_alignments_CAE5_depth_GADAOut_minSegLength1000.tsv.gz
		--methodShortName PopSabaeusCoverageOn3488 --alignmentIDList 6115-6177 --chromosome CAE5 --format tsv
		--data_dir /u/home/p/polyacti/NetworkData/vervet/db/ --commit
		--logFilename  folderLog/AddAlignmentDepthIntervalFile2DB_chr_CAE5.log

Description:
	2013.08.29 add a new AlignmentDepthIntervalFile into db
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy, numpy, re
from pymodule import ProcessOptions, PassingData, utils, NextGenSeq
from pymodule import MatrixFile
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from AddAlignmentDepthIntervalMethod2DB import AddAlignmentDepthIntervalMethod2DB

parentClass = AbstractVervetMapper

class AddAlignmentDepthIntervalFile2DB(parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(parentClass.option_default_dict)
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('chromosome', 0, ): ['', '', 1, 'which chromosome is the data from', ],\
						('alignmentIDList', 0, ): ['', '', 1, 'coma/dash-separated list of alignment IDs, used to verify/create AlignmentDepthIntervalMethod', ],\
						('methodShortName', 1, ):[None, 's', 1, 'column short_name of AlignmentDepthIntervalMethod table, \
		will be created if not present in db.'],\
						('format', 1, ):['tsv', '', 1, 'format for AlignmentDepthIntervalFile entry'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.alignmentIDList = utils.getListOutOfStr(list_in_str=self.alignmentIDList, data_type=int)
		
		self.characterPattern = re.compile(r'[a-zA-Z]')

	def parseInputFile(self, inputFname=None, **keywords):
		"""
		2013.08.23
		"""
		reader = MatrixFile(inputFname=inputFname)
		no_of_intervals = 0
		interval_value_ls = []
		chromosome_size = 0
		min_interval_value = None
		max_interval_value = None
		for row in reader:
			if row[0][0]=='#' or self.characterPattern.search(row[0]):	#2013.08.28 skip comments and headers
				continue
			span = int(row[2])
			value = float(row[3])
			no_of_intervals += 1
			chromosome_size += span
			interval_value_ls.append(value)
			if min_interval_value is None or value <min_interval_value:
				min_interval_value = value
			if max_interval_value is None or value >max_interval_value:
				max_interval_value = value
		return PassingData(no_of_intervals=no_of_intervals, chromosome_size=chromosome_size, \
						mean_interval_value=numpy.mean(interval_value_ls),\
						median_interval_value=numpy.median(interval_value_ls),\
						min_interval_value=min_interval_value,\
						max_interval_value=max_interval_value)
	
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
		
		realPath = os.path.realpath(self.inputFname)
		logMessage = "Handling file %s ...\n"%(self.inputFname)
			
		alignmentList = self.db_vervet.getAlignmentsFromAlignmentIDList(self.alignmentIDList)
		
		method = self.db_vervet.getAlignmentDepthIntervalMethod(short_name=self.methodShortName, description=None, ref_ind_seq_id=None, \
					individualAlignmentLs=alignmentList, parent_db_entry=None, parent_id=None, \
					no_of_alignments=None, no_of_intervals=None, \
					sum_median_depth=None, sum_mean_depth=None,\
					data_dir=self.data_dir)
		self.checkIfAlignmentListMatchMethodDBEntry(individualAlignmentLs=alignmentList, methodDBEntry=method, session=session)
		
		inputFileData = self.parseInputFile(inputFname=self.inputFname)
		logMessage += "chromosome_size=%s, no_of_intervals=%s.\n"%(inputFileData.chromosome_size, inputFileData.no_of_intervals )
		
		db_entry = self.db_vervet.getAlignmentDepthIntervalFile(alignment_depth_interval_method=method, alignment_depth_interval_method_id=None,  \
					path=None, file_size=None, \
					chromosome=self.chromosome, chromosome_size=inputFileData.chromosome_size, \
					no_of_chromosomes=1, no_of_intervals=inputFileData.no_of_intervals,\
					format=self.format,\
					mean_interval_value=inputFileData.mean_interval_value, median_interval_value=inputFileData.median_interval_value, \
					min_interval_value=inputFileData.min_interval_value, max_interval_value=inputFileData.max_interval_value,\
					md5sum=None, original_path=None, data_dir=self.data_dir)
		if db_entry.id and db_entry.path:
			isPathInDB = self.db_vervet.isPathInDBAffiliatedStorage(relativePath=db_entry.path, data_dir=self.data_dir)
			if isPathInDB==-1:
				sys.stderr.write("Error while updating AlignmentDepthIntervalFile.path with the new path, %s.\n"%(db_entry.path))
				self.cleanUpAndExitOnFailure(exitCode=isPathInDB)
			elif isPathInDB==1:	#successful exit, entry already in db
				sys.stderr.write("Warning: file %s is already in db.\n"%\
									(db_entry.path))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=0)
			else:	#not in db affiliated storage, keep going.
				#to overwrite an old db entry
				db_entry.chromosome_size = inputFileData.chromosome_size
				db_entry.no_of_intervals = inputFileData.no_of_intervals
				db_entry.mean_interval_value=inputFileData.mean_interval_value
				db_entry.median_interval_value=inputFileData.median_interval_value
				db_entry.min_interval_value=inputFileData.min_interval_value
				db_entry.max_interval_value=inputFileData.max_interval_value
				session.add(db_entry)
				session.flush()
		
		#move the file and update the db_entry's path as well
		inputFileBasename = os.path.basename(self.inputFname)
		relativePath = db_entry.constructRelativePath(sourceFilename=inputFileBasename)
		exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=inputFileBasename, \
								inputDir=os.path.split(self.inputFname)[0], dstFilename=os.path.join(self.data_dir, relativePath), \
								relativeOutputDir=None, shellCommand='cp -rL', \
								srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
								constructRelativePathFunction=db_entry.constructRelativePath)
		
		if exitCode!=0:
			sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
			session.rollback()
			self.cleanUpAndExitOnFailure(exitCode=exitCode)
		
		self.db_vervet.updateDBEntryPathFileSize(db_entry=db_entry, data_dir=self.data_dir)
		
		#logMessage += " is empty (no loci) or not VCF file.\n"
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
			session.rollback()
			#delete all target files but exit gracefully (exit 0)
			self.cleanUpAndExitOnFailure(exitCode=0)
	


if __name__ == '__main__':
	main_class = AddAlignmentDepthIntervalFile2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()