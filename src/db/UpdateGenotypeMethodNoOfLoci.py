#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice -s 323VRCSKNevisTrioCallerMAF0.1
		-i Contig103.filter_by_vcftools.recode.vcf.gz

Description:
	2012.7.17 
		Update the number of loci in one genotype method by summing the no_of_loci of its GenotypeFile entries.
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

class UpdateGenotypeMethodNoOfLoci(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('genotypeMethodID', 0, int): [None, 'i', 1, 'GenotypeMethod.id, used to fetch db entry, non-zero exit if not present in db', ],\
						('genotypeMethodShortName', 0, ):[None, 's', 1, 'column short_name of GenotypeMethod table,\
			non-zero exit if not present in db.'],\
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
		
		if self.genotypeMethodID:
			genotypeMethod = VervetDB.GenotypeMethod.get(self.genotypeMethodID)
		elif self.genotypeMethodShortName:
			genotypeMethod = VervetDB.GenotypeMethod.query.filter_by(short_name=self.genotypeMethodShortName).first()
		else:
			sys.stderr.write("ERROR: Both genotypeMethodID (%s) and genotypeMethodShortName (%s) are null.\n"%\
							(self.genotypeMethodID, self.genotypeMethodShortName))
			sys.exit(2)
		
		if genotypeMethod is None:
			sys.stderr.write("ERROR: genotypeMethod with genotypeMethodID (%s) or genotypeMethodShortName (%s) doesn't exist in db.\n"%\
							(self.genotypeMethodID, self.genotypeMethodShortName))
			sys.exit(4)
		logMessage = "genotypeMethod %s (%s) has %s loci.\n"%(genotypeMethod.id, self.genotypeMethodShortName, genotypeMethod.no_of_loci)
		self.db_vervet.updateGenotypeMethodNoOfLoci(db_entry=genotypeMethod)
		
		logMessage += "genotypeMethod %s (%s) updated with %s loci.\n"%\
							(genotypeMethod.id, self.genotypeMethodShortName, genotypeMethod.no_of_loci)
		self.outputLogMessage(logMessage=logMessage)
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
	main_class = UpdateGenotypeMethodNoOfLoci
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
