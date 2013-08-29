#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice -s 323VRCSKNevisTrioCallerMAF0.1
		-i Contig103.filter_by_vcftools.recode.vcf.gz

Description:
	2012.7.17 
		Add a new genotype method (based on one VCF file) into database. 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, NextGenSeq, Genome
from pymodule import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
from AddVCFFile2DB import AddVCFFile2DB

class AddGenotypeMethod2DB(AddVCFFile2DB):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('inputFname', 1, ): ['', 'i', 1, 'the VCF file from whose header the alignments&no. of individuals will be parsed.', ],\
						('genotypeMethodShortName', 1, ):[None, 's', 1, 'column short_name of GenotypeMethod table,\
			will be created if not present in db.'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AddVCFFile2DB.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
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
		
		
		vcfFile = VCFFile(inputFname=self.inputFname)
		
		individualAlignmentLs = self.getAlignmentLsFromVCF(db_vervet=self.db_vervet, vcfFile=vcfFile)
		
		genotypeMethod = self.db_vervet.getGenotypeMethod(short_name=self.genotypeMethodShortName, \
												individualAlignmentLs=individualAlignmentLs,\
												no_of_individuals=len(individualAlignmentLs), no_of_loci=None,\
												data_dir=self.data_dir)
		self.checkIfAlignmentListMatchMethodDBEntry(individualAlignmentLs, genotypeMethod, session)
		self.outputLogMessage(logMessage="genotypeMethod %s (%s) added into db.\n"%\
							(genotypeMethod.id, self.genotypeMethodShortName))
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
			#delete all target files
			self.cleanUpAndExitOnFailure(exitCode=0)

if __name__ == '__main__':
	main_class = AddGenotypeMethod2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()