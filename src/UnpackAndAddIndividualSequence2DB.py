#!/usr/bin/env python
"""
Examples:
	# run the program on crocea and output a local condor workflow
	%s -i ~/NetworkData/vervet/VRC/ -t /u/home/eeskintmp/polyacti/NetworkData/vervet/db/ 
		-u yh -m ~/script/vervet/data/VRC_sequencing_20110802.tsv -c -z dl324b-1.cmb.usc.edu -o /tmp/condorpool.xml

	#2011-8-26	generate a list of all bam file physical paths through find. (doing this because they are not located on crocea)
	find /u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/ -name *.bam  > /u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/bamFileList.txt
	# run the program on the crocea and output a hoffman2 workflow. (with db commit)
	%s -i ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/bamFileList.txt
		-e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh
		-m ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
		-z dl324b-1.cmb.usc.edu -l hoffman2
		-o unpackAndAdd12_2007Monkeys2DB_hoffman2.xml
		-c 
	
	# 2011-8-28 output a workflow to run on local condorpool, no db commit (because records are already in db)
	~/script/vervet/src/UnpackAndAddIndividualSequence2DB.py -i /Network/Data/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/
		-u yh
		-m ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
		-z dl324b-1.cmb.usc.edu -l condorpool -o unpackAndAdd12_2007Monkeys2DB_condor.xml
	
Description:
	2011-8-2
		input:
			input (if directory, recursively go and find all bam files; if file, each line is path to bam),
			tsv file (map bam filename to monkey ID),
		
		this program
			1. queries db if individual record for that monkey is already in db, or not. create a new entry in db if not
			2. queries db if individual_sequence record for that monkey is already in db. create one new entry if not.
				update individual_sequence.path if it's none
			3. if commit, the db action would be committed
			4. outputs the whole unpack workflow to an xml output file.
		
		This program has to be run on the pegasus submission host.
		option "-c" commits the db transaction.
	
	The bamFname2MonkeyIDMapFname contains a map from monkey ID to the bam filename (only base filename is used, relative directory is not used).
	Example ("Library" and "Bam Path" are required):
		FlowCell	Lane	Index Sequence	Library	Common Name	Bam Path	MD5
		64J6AAAXX	1		VCAC-2007002-1-lib1	African Green Monkey	gerald_64J6AAAXX_1.bam	gerald_64J6AAAXX_1.bam.md5
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

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
from Pegasus.DAX3 import *

class UnpackAndAddIndividualSequence2DB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number. must be non-empty if need ssh tunnel'],\
						('input', 1, ): ['', 'i', 1, 'if it is a directory, take every *.bam file in it. If it is a file, every line should be a path to actual bam file. ', ],\
						('bamFname2MonkeyIDMapFname', 1, ): ['', 'm', 1, 'a tsv version of WUSTL xls file detailing what monkey is in which bam file.', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/vervet/bin/GenomeAnalysisTK", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. If not given, use the default stored in db.'],\
						("sequencer", 1, ): ["GA", '', 1, 'choices: 454, GA, Sanger'],\
						("sequence_type", 1, ): ["PE", '', 1, 'choices: BAC, genome, scaffold, PE, SR, ...'],\
						("sequence_format", 1, ): ["fastq", '', 1, 'fasta, fastq, etc.'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path and individual_sequence records if inexistent)'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	#('jobFileDir', 0, ): ['', 'j', 1, 'folder to contain qsub scripts', ],\

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
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	def getBamBaseFname2MonkeyID(self, bamFname2MonkeyIDMapFname, ):
		"""
		2011-8-3
			
		"""
		sys.stderr.write("Getting bamBaseFname2MonkeyID dictionary ...")
		bamBaseFname2MonkeyID = {}
		reader = csv.reader(open(bamFname2MonkeyIDMapFname), delimiter='\t')
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkeyIDIndex = col_name2index.get("Library")
		bamFnameIndex = col_name2index.get("Bam Path")
		monkeyIDPattern = re.compile(r'\w+-(\w+)-\d+-\w+')	# i.e. VCAC-2007002-1-lib1
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
			bamBaseFname2MonkeyID[bamBaseFname] = monkeyID
		sys.stderr.write("%s entries.\n"%(len(bamBaseFname2MonkeyID)))
		return bamBaseFname2MonkeyID
	
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
		
	def addMonkeySequence(self, db_vervet, monkeyID, sequencer='GA', sequence_type='PE', sequence_format='fastq',\
						path_to_original_sequence=None):
		"""
		2011-8-3
		"""
		individual = db_vervet.getIndividual(code=monkeyID)
		individual_sequence = db_vervet.getIndividualSequence(individual_id=individual.id, sequencer=sequencer, \
						sequence_type=sequence_type, sequence_format=sequence_format, \
						path_to_original_sequence=path_to_original_sequence, tissue_name=None, coverage=None,\
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
#$ -l time=23:00:00
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
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema,\
					port=self.port)
		db_vervet.setup(create_tables=False)
		session = db_vervet.session
		session.begin()
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		# Create a abstract dag
		workflow = ADAG("UnpackAndAddIndividualSequence2DB")
		site_handler = self.site_handler
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		
		#add the MergeSamFiles.jar file into workflow
		convertBamToFastqAndGzip = Executable(namespace=namespace, name="convertBamToFastqAndGzip", version=version, os=operatingSystem, \
									arch=architecture, installed=True)
		convertBamToFastqAndGzip.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, 'convertBamToFastqAndGzip.sh'), site_handler))
		workflow.addExecutable(convertBamToFastqAndGzip)
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, arch=architecture, \
							installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mkdirWrap.sh"), site_handler))
		workflow.addExecutable(mkdirWrap)
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, \
							installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
		
		
		bamBaseFname2MonkeyID = self.getBamBaseFname2MonkeyID(self.bamFname2MonkeyIDMapFname)
		bamFnameLs = []
		if os.path.isdir(self.input):
			inputDir = self.input
			self.getAllBamFiles(inputDir, bamFnameLs)
		elif os.path.isfile(self.input):
			inf = open(self.input)
			for line in inf:
				bamFnameLs.append(line.strip())
		else:
			sys.stderr.write("%s is neither a folder nor a file.\n"%(self.input))
			sys.exit(4)
			
		sys.stderr.write("%s total bam files.\n"%(len(bamFnameLs)))
		
		for bamFname in bamFnameLs:
			bamBaseFname = os.path.split(bamFname)[1]
			if bamBaseFname not in bamBaseFname2MonkeyID:
				sys.stderr.write("%s doesn't have monkeyID affiliated with.\n"%(bamFname))
				continue
			monkeyID = bamBaseFname2MonkeyID.get(bamBaseFname)
			individual_sequence = self.addMonkeySequence(db_vervet, monkeyID, sequencer=self.sequencer, sequence_type=self.sequence_type, \
										sequence_format=self.sequence_format, path_to_original_sequence=bamFname)
			if individual_sequence.path is None:
				individual = individual_sequence.individual
				individual_sequence.path = db_vervet.constructRelativePathForIndividualSequence(individual_id=individual.id, \
								individual_sequence_id=individual_sequence.id, individual_code=individual.code,\
								sequencer=self.sequencer, tissue=None)
				individual_sequence.original_path = bamFname
				session.add(individual_sequence)
				session.flush()
			
			outputDir = os.path.join(self.dataDir, individual_sequence.path)
			
			# Add a mkdir job for the individual_sequence.path directory.
			mkCallDirJob = Job(namespace=namespace, name=mkdirWrap.name, version=version)
			mkCallDirJob.addArguments(outputDir)
			workflow.addJob(mkCallDirJob)
			
			bamBaseFname = os.path.split(bamFname)[1]
			bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
			outputFnamePrefix = os.path.join(outputDir, bamBaseFnamePrefix)
			
			convertBamToFastqAndGzip_job = Job(namespace=namespace, name=convertBamToFastqAndGzip.name, version=version)
			convertBamToFastqAndGzip_job.addArguments(bamFname, outputFnamePrefix)
			job_max_memory = 2000	#in MB
			convertBamToFastqAndGzip_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))	#same as max_memory
			# GLOBUS 5.04 's RSL says it should be max_memory.
			convertBamToFastqAndGzip_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			job_max_walltime = 1380	#23 hours in minutes (max time allowed on hoffman2 is 24 hours)
			convertBamToFastqAndGzip_job.addProfile(Profile(Namespace.GLOBUS, key="maxwalltime", value="%s"%job_max_walltime))
			# GLOBUS 5.04 's RSL says it should be max_wall_time
			
			workflow.addJob(convertBamToFastqAndGzip_job)
			workflow.addDependency(parent=mkCallDirJob, child=convertBamToFastqAndGzip_job)
			
			"""
			
			jobFname = os.path.join(self.jobFileDir, 'job%s.bam2fastq.sh'%(monkeyID))
			self.writeQsubJob(jobFname, bamFname, os.path.join(self.dataDir, individual_sequence.path), self.vervet_path)
			commandline = 'qsub %s'%(jobFname)
			if self.commit:	#qsub only when db transaction will be committed.
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			"""
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		
		if self.commit:
			session.commit()
		else:
			session.rollback()
		
		
if __name__ == '__main__':
	main_class = UnpackAndAddIndividualSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
