#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  -i  folderMap/926_635_tantalus_GA_vs_524.fq.gz  --individual_alignment_id 926
		--format fastq --minDP 4 --maxDP=16 --minBaseQ=20 --minMapQ 30 --minRMSMapQ 10
		--minDistanceToIndel 5 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--logFilename  folderMap/926_635_tantalus_GA_vs_524_alignment926_2DB.log
		--drivername postgresql --hostname localhost --dbname vervetdb
		--db_user yh --db_passwd secret --schema public --commit
	

Description:
	2013.2.8 add alignment consensus sequence file (fastq format, input of the fq2psmcfa program from PSMC) to db
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Bio import SeqIO
from pymodule import ProcessOptions, PassingData, utils, NextGenSeq
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AddIndividualAlignmentConsensusSequence2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetMapper.option_default_dict)
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('individual_alignment_id', 1, int):[None, '', 1, 'alignment id'],\
						('minDP', 1, int):[None, '', 1, 'minimum read depth at that locus'],\
						('maxDP', 1, int):[None, '', 1, 'maximum read depth at that locus'],\
						('minBaseQ', 1, int):[20, '', 1, 'inferred consensus base quality '],\
						('minMapQ', 1, int):[30, '', 1, 'read alignment mapping quality'],\
						('minRMSMapQ', 1, int):[10, '', 1, 'root mean squared mapping quality of reads covering the locus'],\
						('minDistanceToIndel', 1, int):[5, '', 1, 'min distance to predicted short insertions or deletions'],\
						('format', 1, ):['fastq', '', 1, 'format of the input file'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.inputFname= os.path.realpath(self.inputFname)

	def add2DB(self, db=None, individual_alignment_id=None, inputFname=None, format=None, minDP=None, maxDP=None, minBaseQ=None, minMapQ=None,\
			minRMSMapQ=None, minDistanceToIndel=None, comment=None, data_dir=None, commit=0):
		"""
		2012.11.13
		"""
		session = db.session
		session.begin()
		
		no_of_chromosomes = 0
		no_of_bases = 0
		
		sys.stderr.write("Counting #chromosomes, #bases ...")
		inf = utils.openGzipFile(inputFname)
		for seq_record in SeqIO.parse(inf, 'fastq'):
			no_of_chromosomes += 1
			no_of_bases += len(seq_record)
		inf.close()
		sys.stderr.write("%s chromosomes, %s bases\n"%(no_of_chromosomes, no_of_bases))
		
		#2012.11.13 check if it's in db already
		db_entry = db.checkIndividualAlignmentConsensusSequence(individual_alignment_id=individual_alignment_id, minDP=minDP, \
									maxDP=maxDP, minBaseQ=minBaseQ, minMapQ=minMapQ,\
									minRMSMapQ=minRMSMapQ, minDistanceToIndel=minDistanceToIndel, no_of_chromosomes=no_of_chromosomes)
		if db_entry:
			sys.stderr.write("Warning: IndividualAlignmentConsensusSequence of (individual_alignment_id=%s, minDP %s, maxDP %s, etc.) already in db with id=%s.\n"%\
							(individual_alignment_id, minDP, maxDP, db_entry.id))
			sys.exit(3)
		else:
			db_entry = db.getIndividualAlignmentConsensusSequence(individual_alignment_id=individual_alignment_id, format=format, \
									minDP=minDP, maxDP=maxDP, minBaseQ=minBaseQ, \
									minMapQ=minMapQ, minRMSMapQ=minRMSMapQ, minDistanceToIndel=minDistanceToIndel, \
									no_of_chromosomes=no_of_chromosomes,no_of_bases=no_of_bases, \
									original_path=os.path.abspath(inputFname), data_dir=data_dir)
		
		if commit:
			inputFileBasename = os.path.basename(inputFname)
			#moveFileIntoDBAffiliatedStorage() will also set db_entry.path
			exitCode = db.moveFileIntoDBAffiliatedStorage(db_entry=db_entry, filename=inputFileBasename, \
									inputDir=os.path.split(inputFname)[0], \
									outputDir=data_dir,\
									relativeOutputDir=None, shellCommand='cp -rL', \
									srcFilenameLs=self.srcFilenameLs, dstFilenameLs=self.dstFilenameLs,\
									constructRelativePathFunction=db_entry.constructRelativePath, data_dir=data_dir)
			
			if exitCode!=0:
				sys.stderr.write("Error: moveFileIntoDBAffiliatedStorage() exits with %s code.\n"%(exitCode))
				session.rollback()
				self.cleanUpAndExitOnFailure(exitCode=exitCode)
			
			session.flush()
			session.commit()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
			
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#add the extracted association result into db
		self.add2DB(db=self.db_vervet, individual_alignment_id=self.individual_alignment_id, inputFname=self.inputFname, \
				format=self.format, minDP=self.minDP, maxDP=self.maxDP, minBaseQ=self.minBaseQ, minMapQ=self.minMapQ, \
				minRMSMapQ=self.minRMSMapQ, minDistanceToIndel=self.minDistanceToIndel, comment=None, \
				data_dir=self.data_dir, commit=self.commit)
		
		self.outputLogMessage("submission done.\n")
		

if __name__ == '__main__':
	main_class = AddIndividualAlignmentConsensusSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()