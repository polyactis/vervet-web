#!/usr/bin/env python
"""
Examples:
	# run the program on the hoffman2 cluster. It'll submit jobs in the end.
	%s -i ~/NetworkData/vervet/VRC/ -t /u/home/eeskintmp/polyacti/NetworkData/vervet/db/ 
		-u yh -j ~/qjob -m ~/script/vervet/data/VRC_sequencing_20110802.tsv -c -z dl324b-1.cmb.usc.edu
	
Description:
	2011-8-2
		inputDir (recursively go through), tsv file (map bam filename to monkey ID), jobFileDir (where job files will be stored)
		
		this program
			1. queries the db if specific individual for that monkey is already in db, or not. create a new entry in db if not
			2. write a qsub job file (convertBamToFastqAndGzip.sh, + mkdir + move the fastq to new directory)
			3. submit the job
		
		This program has to be run on the submission host, (currently only hoffman2) because it invokes qsub.
		option "-c" not only commits the db transaction and also qsub all jobs.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, re, csv
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData
from pymodule.utils import runLocalCommand, getColName2IndexFromHeader

class UnpackAndAddIndividualSequence2DB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('inputDir', 1, ): ['', 'i', 1, 'the input directory which ', ],\
						('bamFname2MonkeyIDMapFname', 1, ): ['', 'm', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('jobFileDir', 1, ): ['', 'j', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/vervet/bin/GenomeAnalysisTK", '', 1, 'GATK folder containing its jar binaries'],\
						("vervet_path", 1, ): ["%s/script/vervet/src", '', 1, 'folder containing programs in vervet repositary'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 1, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. If not given, use the default stored in db.'],\
						("sequencer", 1, ): ["GA", '', 1, 'choices: 454, GA, Sanger'],\
						("sequence_type", 1, ): ["PE", '', 1, 'choices: BAC, genome, scaffold, PE, SR, ...'],\
						("sequence_format", 1, ): ["fastq", '', 1, 'fasta, fastq, etc.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path) and qsub jobs'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-8-3
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.samtools_path = self.samtools_path%self.home_path
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervet_path = self.vervet_path%self.home_path
	
	def getBamFname2MonkeyID(self, bamFname2MonkeyIDMapFname, ):
		"""
		2011-8-3
			
		"""
		sys.stderr.write("Getting BamFname2MonkeyID dictionary ...")
		bamFname2MonkeyID = {}
		reader = csv.reader(open(bamFname2MonkeyIDMapFname), delimiter='\t')
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkeyIDIndex = col_name2index.get("Library")
		bamFnameIndex = col_name2index.get("Bam Path")
		monkeyIDPattern = re.compile(r'\w+-(\w+)-\d+-\w+')
		for row in reader:
			monkeyID = row[monkeyIDIndex]
			pa_search = monkeyIDPattern.search(monkeyID)
			if pa_search:
				monkeyID = pa_search.group(1)
			else:
				sys.stderr.write("Warning: could not parse monkey ID from %s. Ignore.\n"%(monkeyID))
				continue
			bamFname = row[bamFnameIndex]
			bamBaseFname = os.path.split(bamFname)[1]
			bamFname2MonkeyID[bamBaseFname] = monkeyID
		sys.stderr.write("%s entries.\n"%(len(bamFname2MonkeyID)))
		return bamFname2MonkeyID
	
	def getAllBamFiles(self, inputDir, bamFnameLs=[]):
		"""
		2011-8-3
			recursively going through the directory to get all bam files
			
		"""
		
		for inputFname in os.listdir(inputDir):
			#get the absolute path
			inputFname = os.path.join(inputDir, inputFname)
			if os.path.isfile(inputFname):
				prefix, suffix = os.path.splitext(inputFname)
				if suffix=='.bam' or suffix=='.sam':
					bamFnameLs.append(inputFname)
			elif os.path.isdir(inputFname):
				self.getAllBamFiles(inputFname, bamFnameLs)
		
	def addMonkeySequence(self, db_vervet, monkeyID, sequencer='GA', sequence_type='short-read', sequence_format='fastq'):
		"""
		2011-8-3
		"""
		individual = db_vervet.getIndividual(code=monkeyID)
		individual_sequence = db_vervet.getIndividualSequence(individual_id=individual.id, sequencer=sequencer, \
															sequence_type=sequence_type,\
						sequence_format=sequence_format, path_to_original_sequence=None, tissue_name=None, coverage=None,\
						subFolder='individual_sequence')
		
		return individual_sequence
	
	def writeQsubJob(self, jobFname, bamFname, outputDir, vervet_path):
		"""
		2011-8-3
		"""
		bamBaseFname = os.path.split(bamFname)[1]
		bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
		outputFnamePrefix = os.path.join(outputDir, bamBaseFnamePrefix)
		commandline="%s/convertBamToFastqAndGzip.sh %s %s"%(vervet_path, bamFname, outputFnamePrefix)
		
		jobF = open(jobFname, 'w')
		jobF.write("""#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_NAME.joblog.$JOB_ID
#$ -j y
#$ -l time=48:00:00
#$ -l highp
#$ -r n

mkdir %s	#create the output directory
echo %s
%s"""%(outputDir, commandline, commandline))
		jobF.close()
	
	def run(self):
		"""
		2011-8-3
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		import VervetDB
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		session = db_vervet.session
		session.begin()
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		bamFname2MonkeyID = self.getBamFname2MonkeyID(self.bamFname2MonkeyIDMapFname)
		bamFnameLs = []
		self.getAllBamFiles(self.inputDir, bamFnameLs)
		sys.stderr.write("%s total bam files.\n"%(len(bamFnameLs)))
		
		for bamFname in bamFnameLs:
			bamBaseFname = os.path.split(bamFname)[1]
			if bamBaseFname not in bamFname2MonkeyID:
				sys.stderr.write("%s doesn't have monkeyID affiliated with.\n"%(bamFname))
				continue
			monkeyID = bamFname2MonkeyID.get(bamBaseFname)
			individual_sequence = self.addMonkeySequence(db_vervet, monkeyID, sequencer=self.sequencer, sequence_type=self.sequence_type, \
														sequence_format=self.sequence_format)
			if individual_sequence.path is None:
				individual = individual_sequence.individual
				individual_sequence.path = db_vervet.constructRelativePathForIndividualSequence(individual_id=individual.id, \
								individual_sequence_id=individual_sequence.id, individual_code=individual.code,\
								sequencer=self.sequencer, tissue=None)
				session.add(individual_sequence)
				session.flush()
			jobFname = os.path.join(self.jobFileDir, 'job%s.bam2fastq.sh'%(monkeyID))
			self.writeQsubJob(jobFname, bamFname, os.path.join(self.dataDir, individual_sequence.path), self.vervet_path)
			commandline = 'qsub %s'%(jobFname)
			if self.commit:	#qsub only when db transaction will be committed.
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
		if self.commit:
			session.commit()
		else:
			session.rollback()
		
if __name__ == '__main__':
	main_class = UnpackAndAddIndividualSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
