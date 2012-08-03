#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i  folder/Contig456.filter_by_vcftools.recode.vcf.gz  -s 323VRCSKNevisTrioCallerMAC10MAF.05 -f VCF
		-c -v postgresql -z uclaOffice -d vervetdb -u yh  -k public

Description:
	2012.5.2
		Add locus from one VCF file into database. 
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


class AddVCFFile2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						#('inputDir', 1, ): ['', 'i', 1, 'input folder that contains split fastq files', ],\
						('genotypeMethodShortName', 1, ):[None, 's', 1, 'column short_name of GenotypeMethod table, \
		will be created if not present in db.'],\
						('checkEmptyVCFByReading', 0, int):[0, 'E', 0, 'toggle to check if a vcf file is empty by reading its content'],\
						('format', 1, ):[None, 'f', 1, 'format for GenotypeFile entry'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
	
	def getAlignmentLsFromVCF(self, db_vervet=None, vcfFile=None):
		"""
		2012.7.13
		"""
		sys.stderr.write("Parsing a list of alignments out of vcfFile sample id list ...")
		alignmentList = []
		for read_group in vcfFile.getSampleIDList():	#sample_id_ls is not good because its index 0 is "ref"
			individualAlignment = db_vervet.parseAlignmentReadGroup(read_group).individualAlignment
			alignmentList.append(individualAlignment)
		sys.stderr.write(" %s alignment IDs.\n"%(len(alignmentList)))
		return alignmentList
	
	def getNoOfIndividualsFromVCFFile(self):
		"""
		2012.7.13
		"""
	
	def getNoOfLociFromVCFFile(self, vcfFile=None):
		"""
		2012.7.13
		
		"""
		sys.stderr.write("Calculating no. of loci ... ")
		no_of_loci = 0
		for vcfRecord in vcfFile.parseIter():
			no_of_loci += 1
		sys.stderr.write("%s loci in this file.\n"%(no_of_loci))
		return no_of_loci
	
	def checkIfAlignmentListMatchDB(self, individualAlignmentLs=[], genotypeMethod=None, session=None):
		"""
		2012.7.18
		"""
		#make sure genotypeMethod.individual_alignment_ls is identical to individualAlignmentLs
		alignmentIDSetInFile = set([alignment.id for alignment in individualAlignmentLs])
		alignmentIDSetInGenotypeMethod = set([alignment.id for alignment in genotypeMethod.individual_alignment_ls])
		if alignmentIDSetInFile!=alignmentIDSetInGenotypeMethod:
			sys.stderr.write("ERROR: alignmentIDSetInFile (%s) doesn't match alignmentIDSetInFile (%s).\n"%\
							(repr(alignmentIDSetInFile), repr(alignmentIDSetInGenotypeMethod)))
			if session:
				session.rollback()
			#delete all target files if there is any
			self.cleanUpAndExitOnFailure(exitCode=2)
		
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
		logMessage = "file %s.\n"%(self.inputFname)
		if NextGenSeq.isFileNameVCF(realPath, includeIndelVCF=False) and \
				not NextGenSeq.isVCFFileEmpty(realPath, checkContent=self.checkEmptyVCFByReading):
			vcfFile = VCFFile(inputFname=self.inputFname)
			
			individualAlignmentLs = self.getAlignmentLsFromVCF(db_vervet=self.db_vervet, vcfFile=vcfFile)
			
			genotypeMethod = self.db_vervet.getGenotypeMethod(short_name=self.genotypeMethodShortName, \
															individualAlignmentLs=individualAlignmentLs,\
															no_of_individuals=len(individualAlignmentLs), no_of_loci=None,\
															dataDir=self.dataDir)
			self.checkIfAlignmentListMatchDB(individualAlignmentLs, genotypeMethod, session)
			
			no_of_loci = self.getNoOfLociFromVCFFile(vcfFile)
			if no_of_loci>0:	#file with zero loci could have identical md5sum
				try:
					md5sum = utils.get_md5sum(realPath)
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					self.cleanUpAndExitOnFailure(exitCode=4)
			else:
				md5sum = None
			db_entry = VervetDB.GenotypeFile.query.filter_by(md5sum=md5sum).first()
			if db_entry:
				sys.stderr.write("Warning: another file %s with the identical md5sum %s as this file %s is already in db.\n"%\
									(db_entry.path, md5sum, realPath))
				session.rollback()
				#2012.8.3 when the jobs are clustered into one merged job and it failed halfway
				# and retried elsewhere, the redundancy check should not exit with non-zero. otherwise the merged job would fail again. 
				self.cleanUpAndExitOnFailure(exitCode=0)
			no_of_individuals = len(individualAlignmentLs)
			chromosome = Genome.getChrFromFname(self.inputFname)
			genotypeFile = self.db_vervet.getGenotypeFile(genotype_method=genotypeMethod,\
										chromosome=chromosome, format=self.format, path=None, file_size=None, md5sum=md5sum,\
										original_path=realPath, no_of_individuals=no_of_individuals, no_of_loci=no_of_loci,\
										dataDir=self.dataDir)
			
			#move the file and update the db_entry's path as well
			exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=genotypeFile, filename=os.path.basename(self.inputFname), \
									inputDir=os.path.split(self.inputFname)[0], outputDir=os.path.join(self.dataDir, genotypeMethod.path), \
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=genotypeFile.constructRelativePath)
			#copy the tbi (tabix) index file if it exists
			tbiFilename = '%s.tbi'%(realPath)
			if os.path.isfile(tbiFilename):
				srcFilename = tbiFilename
				dstFilename = os.path.join(self.dataDir, '%s.tbi'%(genotypeFile.path))
				utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
				logMessage += "tbi file %s has been copied to %s.\n"%(srcFilename, dstFilename)
			## 2012.7.17 commented out because md5sum is calcualted above
			#db_vervet.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=dataDir)
			# #2012.7.17 record the size of db_entry.path (folder or file)
			self.db_vervet.updateDBEntryPathFileSize(db_entry=genotypeFile, data_dir=self.dataDir)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			
			vcfFile.close()
			logMessage += "%s individuals, %s loci, md5sum=%s.\n"%(no_of_individuals, no_of_loci, md5sum)
		else:
			logMessage += " is empty (no loci).\n"
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
	main_class = AddVCFFile2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()