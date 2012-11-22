#!/usr/bin/env python
"""
	%s  -v postgresql -z 149.142.212.14 -d vervetdb -k public -u yh 
		-i 3185_gerald_D14CGACXX_7_GCCAAT -o /Network/Data/vervet/db/individual_sequence/3185_6059_2002099_GA_0_0
		-t individual_sequence/3185_6059_2002099_GA_0_0 -l gerald_D14CGACXX_7_GCCAAT -n 3185 -f fastq
		-c
		-m 1 -a  gerald_D14CGACXX_7_GCCAAT.bam  -g  3185_gerald_D14CGACXX_7_GCCAAT_1.register.log

Description:
	2012.1.27
		program to register and move split output by picard's SplitReadFile.jar to db storage.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, re
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper

from vervet.src import VervetDB

class RegisterAndMoveSplitSequenceFiles(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('inputDir', 1, ): ['', 'i', 1, 'input folder that contains split fastq files', ],\
						('outputDir', 1, ): ['', 'o', 1, 'output folder to which split files from inputDir will be moved', ],\
						('relativeOutputDir', 1, ): ['', 't', 1, 'path of the output folder relative to db.dataDir. it should form the last part of outputDir.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("individual_sequence_id", 1, int): [None, 'n', 1, 'individual_sequence.id that is affiliated to the split files in inputDir'],\
						("bamFilePath", 0, ): [None, 'a', 1, 'the original bam file from which the files in inputDir came'],\
						('library', 1,): [None, 'l', 1, 'library name for files in inputDir', ],\
						('mate_id', 0, ): [None, 'm', 1, '1: first end; 2: 2nd end. of paired-end or mate-paired libraries'],\
						("sequence_format", 1, ): ["fastq", 'f', 1, 'fasta, fastq, etc.'],\
						('logFilename', 0, ): [None, 'g', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow'],\
						('mate_id_associated_with_bam', 0, int):[0, 'w', 0, 'toggle to make mate_id associated with bamFilePath as well'],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=None, **keywords)	#connectDB(), and setup srcFilenameLs and dstFilenameLs

	
	def addBamFileToDB(self, db_vervet, bamFilePath, library=None, individual_sequence_id=None, mate_id=None):
		"""
		2012.4.30
			add mate_id
		2012.1.27
			1. run md5sum
			2. check if it already exists in db
				if not, add it into db
				if yes, exit the program. no more work.
		"""
		sys.stderr.write("Record the bamfile into db ...")
		md5sum = utils.get_md5sum(bamFilePath)
		db_entry = VervetDB.IndividualSequenceFileRaw.query.filter_by(md5sum=md5sum).first()
		if db_entry:
			sys.stderr.write("Warning: another file %s with the identical md5sum %s (library=%s) as this file %s is already in db.\n"%\
								(db_entry.path, md5sum, library, bamFilePath))
			#sys.exit(3)
		else:
			db_entry = db_vervet.getIndividualSequenceFileRaw(individual_sequence_id, library=library, md5sum=md5sum, \
											path=bamFilePath, mate_id=mate_id)
		realpath = os.path.realpath(bamFilePath)
		if realpath!=db_entry.path:
			db_entry.path = realpath
			db_vervet.session.add(db_entry)
			db_vervet.session.flush()
		mate_id2split_order_ls = {}
		for individual_sequence_file in db_entry.individual_sequence_file_ls:
			mate_id = individual_sequence_file.mate_id
			if mate_id is None:	#sequence entries without mate_id are just from one mate.
				mate_id = 1
			if mate_id not in mate_id2split_order_ls:
				mate_id2split_order_ls[mate_id] = []
			mate_id2split_order_ls[mate_id].append(individual_sequence_file.split_order)
		if len(mate_id2split_order_ls)>2:
			sys.stderr.write("Error: db sequence files spawned from bam file %s (md5sum=%s) form %s(>2) mates.\
				Unless this bam file contains reads from >2 mates, reads from this bam files should be stored in db already.\n"%\
				(db_entry.path, md5sum, len(mate_id2split_order_ls)))
			sys.exit(4)
		return db_entry
	
	def parseSplitOrderOutOfFilename(self, filename, library, mate_id=None):
		"""
		2012.2.9
			mate_id is optional
		2012.1.27
			filename might look like gerald_81LL0ABXX_4_TTAGGC_2_1.fastq.gz.
				library_(mate_id)_(split_order).fastq.gz
			library is gerald_81LL0ABXX_4_TTAGGC.
			mate_id is 2.
			split_order is 1.
		"""
		if mate_id:
			prefix = '%s_%s'%(library, mate_id)
		else:
			prefix = library
		split_order_pattern = re.compile(r'%s_(\d+).fastq'%(prefix))
		split_order_search_result = split_order_pattern.search(filename)
		if split_order_search_result:
			split_order = split_order_search_result.group(1)
		else:
			split_order = None
		return split_order
	
	def moveNewISQFileIntoDBStorage(self, session, individual_sequence_file=None, filename=None, inputDir=None, outputDir=None, \
								relativeOutputDir=None, shellCommand='cp', srcFilenameLs=None, dstFilenameLs=None):
		"""
		2012.7.13 superceded by VervetDB.moveFileIntoDBAffiliatedStorage()
		2012.7.4
			add srcFilename and dstFilename into given arguments (srcFilenameLs, dstFilenameLs) for later undo
		2012.6.8
			return non-zero if failure in move or destinaion file already exists
		2012.2.10
			this function moves a file to a db-affiliated storage path
		"""
		newfilename = '%s_%s'%(individual_sequence_file.id, filename)
		newPath = os.path.join(relativeOutputDir, newfilename)
		if individual_sequence_file.path!=newPath:
			individual_sequence_file.path = newPath
			session.add(individual_sequence_file)
			session.flush()
		
		srcFilename = os.path.join(inputDir, filename)
		dstFilename = os.path.join(outputDir, newfilename)
		if os.path.isfile(dstFilename):
			sys.stderr.write("Error: destination %s already exits.\n"%(dstFilename))
			exitCode = 2
		else:
			#move the file
			commandline = '%s %s %s'%(shellCommand, srcFilename, dstFilename)
			return_data = utils.runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			if srcFilenameLs is not None:
				srcFilenameLs.append(srcFilename)
			if dstFilenameLs is not None:
				dstFilenameLs.append(dstFilename)
			
			if return_data.stderr_content:
				#something wrong. abort
				sys.stderr.write("commandline %s failed: %s\n"%(commandline, return_data.stderr_content))
				#remove the db entry
				session.delete(individual_sequence_file)
				session.flush()
				exitCode = 3
			else:
				exitCode = 0
		return exitCode
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		relativePathIndex = self.outputDir.find(self.relativeOutputDir)
		noOfCharsInRelativeOutputDir = len(self.relativeOutputDir)
		if self.outputDir[relativePathIndex:relativePathIndex+noOfCharsInRelativeOutputDir]!=self.relativeOutputDir:
			sys.stderr.write('Error: relativeOutputDir %s is not the last part of outputDir %s.\n'%\
							(self.relativeOutputDir, self.outputDir))
			sys.exit(4)
		
		db_vervet = self.db_vervet
		session = db_vervet.session
		session.begin()
		
		if self.bamFilePath:
			if self.mate_id_associated_with_bam:
				mate_id_passed = self.mate_id
			else:
				mate_id_passed = None
			file_raw_db_entry = self.addBamFileToDB(db_vervet, self.bamFilePath, library=self.library, \
										individual_sequence_id=self.individual_sequence_id, mate_id=mate_id_passed)
		else:
			file_raw_db_entry = None
		
		counter = 0
		real_counter = 0
		for filename in os.listdir(self.inputDir):
			split_order = self.parseSplitOrderOutOfFilename(filename, self.library, self.mate_id)
			counter += 1
			if split_order:
				#save db entry
				if file_raw_db_entry:
					individual_sequence_file_raw_id = file_raw_db_entry.id
				else:
					individual_sequence_file_raw_id = None
				db_entry = db_vervet.getIndividualSequenceFile(self.individual_sequence_id, library=self.library, mate_id=self.mate_id, \
										split_order=int(split_order), format=self.sequence_format,\
										filtered=0, parent_individual_sequence_file_id=None, \
										individual_sequence_file_raw_id=individual_sequence_file_raw_id)
				
				#move the file
				exitCode = db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=filename, \
													inputDir=self.inputDir, outputDir=self.outputDir, \
								relativeOutputDir=self.relativeOutputDir, shellCommand='cp -rL', \
								srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
								constructRelativePathFunction=None)
				"""
				exitCode = self.moveNewISQFileIntoDBStorage(session, individual_sequence_file=db_entry, filename=filename, inputDir=self.inputDir, \
										outputDir=self.outputDir, relativeOutputDir=self.relativeOutputDir,\
										shellCommand='cp', srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs)
				"""
				if exitCode!=0:
					sys.stderr.write("Error: moveNewISQFileIntoDBStorage() exits with %s code.\n"%(exitCode))
					session.rollback()
					#delete all recorded target files
					self.cleanUpAndExitOnFailure(exitCode=exitCode)
				real_counter += 1
		
		if self.logFilename:
			outf = open(self.logFilename, 'w')
			outf.write("%s files processed, %s of them added into db.\n"%(counter, real_counter))
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
				#delete all target files.
				self.cleanUpAndExitOnFailure(exitCode=3)
		else:
			session.rollback()
			#delete all target files
			self.cleanUpAndExitOnFailure(exitCode=0)
		
if __name__ == '__main__':
	main_class = RegisterAndMoveSplitSequenceFiles
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()