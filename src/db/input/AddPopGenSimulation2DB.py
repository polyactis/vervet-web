#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i  OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_dupMarked.bam
		--logFilename  OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_2db.log
		--individual_alignment_id 2278 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--drivername postgresql --hostname localhost --dbname vervetdb --db_user yh
		--schema public OneLibAlignment/2278_634_vs_524_by_2_r4043_sequence_628C2AAXX_6_dupMarked.metric

Description:
	2013.08.04 add population genetics simulation files (PolymorphismTableFile) to db
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AddPopGenSimulation2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						#('inputDir', 1, ): ['', 'i', 1, 'input folder that contains split fastq files', ],\
						('r', 1, float):[None, '', 1, 'recombination rate per base per generation used in simulation'],\
						('rho', 1, float):[None, '', 1, 'population recombination rate used in simulation'],\
						('mu', 1, float):[None, '', 1, 'mutation rate per base per generation used in simulation'],\
						('theta', 1, float):[None, '', 1, 'population mutation rate used in simulation'],\
						('n0', 0, int):[None, '', 1, 'initial population size in simulation'],\
						('is_selection', 0, int):[None, '', 1, 'population selection parameter (s) used in simulation'],\
						('selection_parameters', 0, ):[None, '', 1, 'selection parameters used in simulation'],\
						('replicate_index', 0, int):[0, '', 1, 'which replicate instance for the same simulation type'],\
						('no_of_populations', 0, int):[1, '', 1, 'number of populations in the simulation'],\
						('no_of_chromosomes', 1, int):[None, '', 1, 'number of haplotypes in each simulated population'],\
						('chromosome_length', 1, int):[None, '', 1, 'length of simulated chromosome'],\
						('sample_size', 1, int):[None, '', 1, 'output sample size in the simulation'],\
						('no_of_polymorphic_loci', 0, int):[None, '', 1, 'number of polymorhpic loci in simulation output'],\
						('population_size_parameters', 0, ):[None, '', 1, 'population size parameters for the population genetics simulation'],\
						('parent_pop_gen_simulation_type_id', 0, ):[None, '', 1, 'if this pop genetics simulation takes input from another pop gen simulation, \n\
	what ID is its pop gene simulation type?'],\
						('simulation_programs', 0, ):[None, '', 1, 'which population genetics simulation program(s) used'],\
						('commandline', 0, ):[None, '', 1, 'commandline for population genetics simulation'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		 pop_gen_simulation_type_id=None, replicate_index=None, \
						no_of_populations=None,\
						no_of_chromosomes=None, chromosome_length=None, sample_size=None, \
						no_of_polymorphic_loci=None,
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def run(self):
		"""
		2013.08.04
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
		logMessage = "Adding file %s to db .\n"%(self.inputFname)
		
		if os.path.isfile(inputFileRealPath):
			popGenSimulationType = self.db_vervet.getPopGenSimulationType( short_name=None, r=self.r, rho=self.rho, \
							mu=self.mu, theta=self.theta, n0=self.n0, is_selection=self.is_selection,\
							selection_parameters=self.selection_parameters, indel=None, indel_parameters=None, \
							population_size_parameters=self.population_size_parameters, \
							parent_pop_gen_simulation_type_id=self.parent_pop_gen_simulation_type_id)
			popGenSimulation = self.db_vervet.getPopGenSimulation(pop_gen_simulation_type_id=popGenSimulationType.id, \
						replicate_index=self.replicate_index, no_of_populations=self.no_of_populations,\
						no_of_chromosomes=self.no_of_chromosomes, chromosome_length=self.chromosome_length, \
						sample_size=self.sample_size, \
						no_of_polymorphic_loci=self.no_of_polymorphic_loci, programs=self.simulation_programs,\
						original_path=inputFileRealPath, data_dir=self.data_dir)
			try:
				md5sum = utils.get_md5sum(inputFileRealPath)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				self.cleanUpAndExitOnFailure(exitCode=4)
			
			db_entry = VervetDB.PopGenSimulation.query.filter_by(md5sum=md5sum).first()
			if db_entry and db_entry.id!=popGenSimulation.id and db_entry.path and os.path.isfile(os.path.join(data_dir, db_entry.path)):
				sys.stderr.write("Warning: another file %s with the identical md5sum %s as this file %s, is already in db.\n"%\
								(db_entry.path, md5sum, inputFileRealPath))
				self.sessionRollback(session)
				self.cleanUpAndExitOnFailure(exitCode=3)
			
			
			if popGenSimulation.md5sum is None or popGenSimulation.md5sum!=md5sum:
				popGenSimulation.md5sum = md5sum
				session.add(popGenSimulation)
				session.flush()
			try:
				#move the file and update the db_entry's path as well
				exitCode = self.db_vervet.moveFileIntoDBAffiliatedStorage(db_entry=popGenSimulation, \
							filename=os.path.basename(inputFileRealPath), \
							inputDir=os.path.split(inputFileRealPath)[0], dstFilename=os.path.join(self.data_dir, popGenSimulation.path), \
							relativeOutputDir=None, shellCommand='cp -rL', \
							srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
							constructRelativePathFunction=popGenSimulation.constructRelativePath)
			except:
				sys.stderr.write('Except in copying %s to db-storage with except info: %s\n'%(inputFileRealPath, repr(sys.exc_info())))
				import traceback
				traceback.print_exc()
				self.sessionRollback(session)
				self.cleanUpAndExitOnFailure(exitCode=5)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with code=%s.\n"%(exitCode))
				self.sessionRollback(session)
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			try:
				#make sure these files are stored in self.dstFilenameLs and self.srcFilenameLs
				#copy further files if there are
				if self.inputFnameLs:
					for inputFname in self.inputFnameLs:
						if inputFname!=self.inputFname:	#2013.3.18 make sure it has not been copied.
							logMessage = self.db_vervet.copyFileWithAnotherFilePrefix(inputFname=inputFname, \
												filenameWithPrefix=popGenSimulation.path, \
												outputDir=self.data_dir,\
												logMessage=logMessage, srcFilenameLs=self.srcFilenameLs, \
												dstFilenameLs=self.dstFilenameLs)
				
				self.db_vervet.updateDBEntryPathFileSize(db_entry=popGenSimulation, data_dir=data_dir)
				
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
	main_class = AddPopGenSimulation2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()