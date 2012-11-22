#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i  OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_dupMarked.bam
		--logFilename  OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_2db.log
		--individual_alignment_id 2278 --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--drivername postgresql --hostname localhost --dbname vervetdb --db_user yh
		--schema public OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_dupMarked.metric

Description:
	2012.9.21
		Add the alignment file into database
		1. register an alignment entry in db
			a. if individual_alignment_id is not None
			b. if parent_individual_alignment_id is not None:
			c. construct alignment using all other arguments
		2. copy the file (& bai file if it exists and other files in the commandline arguments) over
		3. write the log if instructed so (for workflow purpose)
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, NextGenSeq, Genome
from pymodule.VCFFile import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AddAlignmentFile2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						#('inputDir', 1, ): ['', 'i', 1, 'input folder that contains split fastq files', ],\
						('individual_alignment_id', 0, int):[None, '', 1, 'fetch the db individual_alignment based on this ID'],\
						('individual_sequence_id', 0, int):[None, '', 1, 'used to construct individual_alignment'],\
						('ref_sequence_id', 0, int):[None, '', 1, 'used to construct individual_alignment'],\
						('alignment_method_id', 0, int):[None, '', 1, 'used to construct individual_alignment'],\
						('parent_individual_alignment_id', 0, int):[None, '', 1, 'the parent ID of individual_alignment.\n\
	if given, an individual_alignment db entry will be created as a copy of this one.'],\
						('mask_genotype_method_id', 0, int):[None, '', 1, 'for alignments coming out of base quality recalibration'],\
						('individual_sequence_file_raw_id', 0, int):[None, '', 1, 'for library specific alignment'],\
						('format', 0, ):[None, 'f', 1, 'format for GenotypeFile entry'],\
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
		if not self.dataDir:
			self.dataDir = self.db_vervet.data_dir
		dataDir = self.dataDir
		
		realPath = os.path.realpath(self.inputFname)
		logMessage = "Adding file %s to db .\n"%(self.inputFname)
		
		if os.path.isfile(realPath):
			if self.individual_alignment_id:
				individual_alignment = VervetDB.IndividualAlignment.get(self.individual_alignment_id)
			elif self.parent_individual_alignment_id:
				individual_alignment = self.db_vervet.copyParentIndividualAlignment(parent_individual_alignment_id=self.parent_individual_alignment_id,\
																	mask_genotype_method_id=self.mask_genotype_method_id,\
																	dataDir=self.dataDir)
			else:
				#alignment for this library of the individual_sequence
				individual_sequence = VervetDB.IndividualSequence.get(self.individual_sequence_id)
				
				individual_alignment = self.db_vervet.getAlignment(individual_sequence_id=self.individual_sequence_id,\
										path_to_original_alignment=None, sequencer=individual_sequence.sequencer,\
										sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
										ref_individual_sequence_id=self.ref_sequence_id, \
										alignment_method_id=self.alignment_method_id, alignment_format=self.format,\
										individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
										dataDir=dataDir, \
										mask_genotype_method_id=self.mask_genotype_method_id, \
										parent_individual_alignment_id=self.parent_individual_alignment_id,\
										individual_sequence_file_raw_id=self.individual_sequence_file_raw_id)
			needSessionFlush = False
			if not individual_alignment.path:
				individual_alignment.path = individual_alignment.constructRelativePath()
				needSessionFlush = True
			
			if self.mask_genotype_method_id and \
					individual_alignment.mask_genotype_method_id!=self.mask_genotype_method_id:
				individual_alignment.mask_genotype_method_id = self.mask_genotype_method_id
				needSessionFlush = True
			if self.individual_sequence_file_raw_id and \
					individual_alignment.individual_sequence_file_raw_id != self.individual_sequence_file_raw_id:
				individual_alignment.individual_sequence_file_raw_id = self.individual_sequence_file_raw_id
				needSessionFlush = True
			
			if needSessionFlush:
				session.add(individual_alignment)
				session.flush()
			
			try:
				md5sum = utils.get_md5sum(realPath)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				self.cleanUpAndExitOnFailure(exitCode=4)
			
			db_entry = VervetDB.IndividualAlignment.query.filter_by(md5sum=md5sum).first()
			if db_entry and db_entry.id!=individual_alignment.id and db_entry.path and os.path.isfile(os.path.join(dataDir, db_entry.path)):
				sys.stderr.write("Warning: another file %s with the identical md5sum %s as this file %s, is already in db.\n"%\
								(db_entry.path, md5sum, realPath))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=3)
			
			
			if individual_alignment.md5sum is None or individual_alignment.md5sum!=md5sum:
				individual_alignment.md5sum = md5sum
				session.add(individual_alignment)
				session.flush()
			
			#move the file and update the db_entry's path as well
			exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=individual_alignment, filename=os.path.basename(realPath), \
									inputDir=os.path.split(realPath)[0], dstFilename=os.path.join(self.dataDir, individual_alignment.path), \
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=individual_alignment.constructRelativePath)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			try:
				#make sure these files are stored in self.dstFilenameLs and self.srcFilenameLs
				#copy further files if there are
				if self.inputFnameLs:
					for inputFname in self.inputFnameLs:
						logMessage = self.db_vervet.copyFileWithAnotherFilePrefix(inputFname=inputFname, \
												filenameWithPrefix=individual_alignment.path, \
												outputDir=self.dataDir,\
												logMessage=logMessage, srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs)
				
				self.db_vervet.updateDBEntryPathFileSize(db_entry=individual_alignment, data_dir=dataDir)
				
				## 2012.7.17 commented out because md5sum is calculated above
				#db_vervet.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=dataDir)
				#copy the bai index file if it exists
				baiFilename = '%s.bai'%(realPath)
				if os.path.isfile(baiFilename):
					srcFilename = baiFilename
					dstFilename = os.path.join(self.dataDir, '%s.bai'%(individual_alignment.path))
					utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
					logMessage += "bai file %s has been copied to %s.\n"%(srcFilename, dstFilename)
					self.srcFilenameLs.append(srcFilename)
					self.dstFilenameLs.append(dstFilename)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=5)
		else:
			logMessage += "%s doesn't exist.\n"%(realPath)
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
	main_class = AddAlignmentFile2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()