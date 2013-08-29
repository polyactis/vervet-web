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

import copy
from pymodule import ProcessOptions, PassingData, utils, NextGenSeq
from pymodule import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AddVCFFile2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetMapper.option_default_dict)
	#option_default_dict.pop(('inputFname', 0, ))
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
		2012.8.30
			add chromosome info as well
		2012.7.13
		
		"""
		sys.stderr.write("Calculating no. of loci & no. of chromosomes ... ")
		no_of_loci = 0
		chromosome2noOfLoci = {}
		for vcfRecord in vcfFile.parseIter():
			no_of_loci += 1
			if vcfRecord.chr not in chromosome2noOfLoci:
				chromosome2noOfLoci[vcfRecord.chr] =0
			chromosome2noOfLoci[vcfRecord.chr] += 1
		sys.stderr.write("%s loci & %s chromosomes in this file.\n"%(no_of_loci, len(chromosome2noOfLoci)))
		return PassingData(no_of_loci=no_of_loci, chromosome2noOfLoci=chromosome2noOfLoci)
	
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
		logMessage = "file %s.\n"%(self.inputFname)
		if NextGenSeq.isFileNameVCF(realPath, includeIndelVCF=True) and \
				not NextGenSeq.isVCFFileEmpty(realPath, checkContent=self.checkEmptyVCFByReading):
			vcfFile = VCFFile(inputFname=self.inputFname)
			
			individualAlignmentLs = self.getAlignmentLsFromVCF(db_vervet=self.db_vervet, vcfFile=vcfFile)
			
			genotypeMethod = self.db_vervet.getGenotypeMethod(short_name=self.genotypeMethodShortName, \
															individualAlignmentLs=individualAlignmentLs,\
															no_of_individuals=len(individualAlignmentLs), no_of_loci=None,\
															data_dir=self.data_dir)
			self.checkIfAlignmentListMatchMethodDBEntry(individualAlignmentLs, genotypeMethod, session)
			
			pdata = self.getNoOfLociFromVCFFile(vcfFile)
			chromosome2noOfLoci = pdata.chromosome2noOfLoci
			no_of_loci = pdata.no_of_loci
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
			"""
			db_entry = VervetDB.GenotypeFile.query.filter_by(md5sum=md5sum).first()
			if db_entry:
				sys.stderr.write("Warning: another file %s with the identical md5sum %s as this file %s is already in db.\n"%\
									(db_entry.path, md5sum, realPath))
				session.rollback()
				#2012.8.3 when the jobs are clustered into one merged job and it failed halfway
				# and retried elsewhere, the redundancy check should not exit with non-zero. otherwise the merged job would fail again. 
				self.cleanUpAndExitOnFailure(exitCode=0)
			"""
			no_of_individuals = len(individualAlignmentLs)
			no_of_chromosomes = len(chromosome2noOfLoci)
			if no_of_chromosomes == 1:	#2012.8.30 use 1st chromosome
				chromosome = chromosome2noOfLoci.keys()[0]
			else:
				chromosome = None
			genotypeFile = self.db_vervet.getGenotypeFile(genotype_method=genotypeMethod,\
										chromosome=chromosome, format=self.format, path=None, file_size=None, md5sum=md5sum,\
										original_path=realPath, no_of_individuals=no_of_individuals, no_of_loci=no_of_loci,\
										data_dir=self.data_dir, no_of_chromosomes=no_of_chromosomes)
			if genotypeFile.id and genotypeFile.path:
				isPathInDB = self.db_vervet.isPathInDBAffiliatedStorage(relativePath=genotypeFile.path, data_dir=self.data_dir)
				if isPathInDB==-1:
					sys.stderr.write("Error while updating genotypeFile.path with the new path, %s.\n"%(genotypeFile.path))
					self.cleanUpAndExitOnFailure(exitCode=isPathInDB)
				elif isPathInDB==1:	#successful exit, entry already in db
					sys.stderr.write("Warning: file %s is already in db.\n"%\
										(genotypeFile.path))
					session.rollback()
					self.cleanUpAndExitOnFailure(exitCode=0)
				else:	#not in db affiliated storage, keep going.
					pass
			#move the file and update the db_entry's path as well
			inputFileBasename = os.path.basename(self.inputFname)
			relativePath = genotypeFile.constructRelativePath(sourceFilename=inputFileBasename)
			exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=genotypeFile, filename=inputFileBasename, \
									inputDir=os.path.split(self.inputFname)[0], dstFilename=os.path.join(self.data_dir, relativePath), \
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=genotypeFile.constructRelativePath)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			
			#copy the tbi (tabix) index file if it exists
			tbiFilename = '%s.tbi'%(realPath)
			if os.path.isfile(tbiFilename):
				srcFilename = tbiFilename
				dstFilename = os.path.join(self.data_dir, '%s.tbi'%(genotypeFile.path))
				utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
				logMessage += "tbi file %s has been copied to %s.\n"%(srcFilename, dstFilename)
			## 2012.7.17 commented out because md5sum is calcualted above
			#db_vervet.updateDBEntryMD5SUM(db_entry=genotypeFile, data_dir=data_dir)
			# #2012.7.17 record the size of db_entry.path (folder or file)
			self.db_vervet.updateDBEntryPathFileSize(db_entry=genotypeFile, data_dir=self.data_dir)
			
			vcfFile.close()
			logMessage += "%s individuals, %s loci, md5sum=%s.\n"%(no_of_individuals, no_of_loci, md5sum)
		else:
			logMessage += " is empty (no loci) or not VCF file.\n"
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