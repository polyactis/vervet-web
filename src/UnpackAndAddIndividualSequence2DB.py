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
	%s -i /Network/Data/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/
		-u yh
		-m ~/mnt/hoffman2/u/home/eeskintmp/polyacti/NetworkData/vervet/raw_sequence/xfer.genome.wustl.edu/gxfer3/46019922623327/Vervet_12_4X_README.tsv
		-z dl324b-1.cmb.usc.edu -j condorpool -l condorpool -o unpackAndAdd12_2007Monkeys2DB_condor.xml
	
	# 2012.4.30 run on hcondor, to import McGill 1X data (-y2), (-e) is not necessary because it's running on hoffman2 and can recognize home folder.
	#. -H means it needs sshTunnel for db-interacting jobs.
	%s -i ~/NetworkData/vervet/raw_sequence/McGill96_1X/ -z localhost -u yh -j hcondor -l hcondor -c
		-o workflow/unpackMcGill96_1X.xml -y2 -H 
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -e /u/home/eeskin/polyacti
	
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
	
	1. Be careful with the db connection setting as it'll be passed to the db-registration job.
		Make sure all computing nodes have access to the db.
	2. The workflow has to be run on nodes where they have direct db and db-affiliated file-storage access.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, re, csv
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from pymodule.utils import runLocalCommand, getColName2IndexFromHeader
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class UnpackAndAddIndividualSequence2DB(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('input', 1, ): ['', 'i', 1, 'if it is a folder, take all .bam/.sam/.fastq files recursively. If it is a file, every line should be a path to the input file.', ],\
						('bamFname2MonkeyIDMapFname', 0, ): ['', 'm', 1, 'a tsv version of WUSTL xls file detailing what monkey is in which bam file.', ],\
						('minNoOfReads', 1, int): [1000000, '', 1, 'minimum number of reads in each split fastq file. The upper limit in each split file is 2*minNoOfReads.', ],\
						("sequencer", 1, ): ["GA", '', 1, 'choices: 454, GA, Sanger'],\
						("sequence_type", 1, ): ["PE", '', 1, 'choices: BAC, genome, scaffold, PE, SR, ...'],\
						("sequence_format", 1, ): ["fastq", 'f', 1, 'fasta, fastq, etc.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path and individual_sequence records if inexistent)'],\
						('inputType', 1, int): [1, 'y', 1, 'input type. 1: bam files from WUSTL. 2: fastq files from McGill (paired ends are split.).', ],\
						})
	#('jobFileDir', 0, ): ['', 'j', 1, 'folder to contain qsub scripts', ],\
	
	def __init__(self,  **keywords):
		"""
		2011-8-3
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
	
	def getBamBaseFname2MonkeyID(self, inputFname, ):
		"""
		2011-8-3
			the input looks like this:
			#	FlowCell	Lane	Index Sequence	Library	Common Name	Bam Path	MD5
			1	64J6AAAXX	1	VCAC-2007002-1-lib1	African	Green	Monkey	/gscmnt/sata755/production/csf_111215677/gerald_64J6AAAXX_1.bam	/gscmnt/sata755/production/csf_111215677/gerald_64J6AAAXX_1.bam.md5
			2	64J6AAAXX	2	VCAC-2007006-1-lib1	African	Green	Monkey	/gscmnt/sata751/production/csf_111215675/gerald_64J6AAAXX_2.bam	/gscmnt/sata751/production/csf_111215675/gerald_64J6AAAXX_2.bam.md5
		"""
		sys.stderr.write("Getting bamBaseFname2MonkeyID dictionary ...")
		bamBaseFname2MonkeyID = {}
		reader = csv.reader(open(inputFname), delimiter='\t')
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkeyIDIndex = col_name2index.get("Library")
		bamFnameIndex = col_name2index.get("Bam Path")
		if bamFnameIndex is None:	#2012.2.9
			bamFnameIndex = col_name2index.get("BAM Path")
		#monkeyIDPattern = re.compile(r'\w+-(\w+)-\d+-\w+')	# i.e. VCAC-2007002-1-lib1
		monkeyIDPattern = re.compile(r'\w+-(\w+)-\w+-\w+')	# 2012.5.29 i.e. VCAC-VGA00006-AGM0075-lib1 ,
		# VCAC-VZC1014-AGM0055-lib1, VCAC-1996031-VRV0265-lib2a, VCAC-VKD7-361-VKD7-361-lib1 (VKD7 is to be taken),
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
						path_to_original_sequence=None, dataDir=None):
		"""
		2012.4.30
			add argument dataDir
		2011-8-3
		"""
		individual = db_vervet.getIndividual(code=monkeyID)
		individual_sequence = db_vervet.getIndividualSequence(individual_id=individual.id, sequencer=sequencer, \
						sequence_type=sequence_type, sequence_format=sequence_format, \
						path_to_original_sequence=path_to_original_sequence, tissue_name=None, coverage=None,\
						subFolder='individual_sequence', dataDir=dataDir)
		
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
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.1.3
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		convertBamToFastqAndGzip = Executable(namespace=namespace, name="convertBamToFastqAndGzip", version=version, os=operatingSystem, \
									arch=architecture, installed=True)
		convertBamToFastqAndGzip.addPFN(PFN("file://" + os.path.join(vervetSrcPath, 'shell/convertBamToFastqAndGzip.sh'), site_handler))
		#convertBamToFastqAndGzip.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(convertBamToFastqAndGzip)
		workflow.convertBamToFastqAndGzip = convertBamToFastqAndGzip
		
		splitReadFile = Executable(namespace=namespace, name="splitReadFile", version=version, os=operatingSystem, \
								arch=architecture, installed=True)
		splitReadFile.addPFN(PFN("file://" + os.path.join(vervetSrcPath, 'shell/SplitReadFileWrapper.sh'), site_handler))
		#splitReadFile.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(splitReadFile)
		workflow.splitReadFile = splitReadFile
		
		registerAndMoveSplitSequenceFiles  = Executable(namespace=namespace, name="registerAndMoveSplitSequenceFiles", version=version, os=operatingSystem, \
								arch=architecture, installed=True)
		registerAndMoveSplitSequenceFiles.addPFN(PFN("file://" + os.path.join(vervetSrcPath, 'mapper/RegisterAndMoveSplitSequenceFiles.py'), site_handler))
		#registerAndMoveSplitSequenceFiles.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(registerAndMoveSplitSequenceFiles)
		workflow.registerAndMoveSplitSequenceFiles = registerAndMoveSplitSequenceFiles
		
	def addConvertBamToFastqAndGzipJob(self, workflow=None, executable=None, \
							inputF=None, outputFnamePrefix=None, \
							parentJobLs=[], job_max_memory=2000, job_max_walltime = 800, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.1.3
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			The executable should be convertBamToFastqAndGzip.
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments(inputF, outputFnamePrefix)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		output1 = "%s_1.fastq.gz"%(outputFnamePrefix)
		output2 = "%s_2.fastq.gz"%(outputFnamePrefix)
		
		job.uses(output1, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
		job.uses(output2, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
		
		job.output1 = output1
		job.output2 = output2
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addSplitReadFileJob(self, workflow, executable=None, \
							inputF=None, outputFnamePrefix=None, outputFnamePrefixTail="",\
							minNoOfReads=5000000, logFile=None, parentJobLs=[], job_max_memory=2000, job_max_walltime = 800, \
							extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.2.9
			argument outputFnamePrefixTail is now useless.
		2012.1.24
			executable is shell/SplitReadFileWrapper.sh
			
			which calls "wc -l" to count the number of reads beforehand to derive a proper minNoOfReads (to avoid files with too few reads).
			
			run SplitReadFile and generate the output directly into the db-affiliated folders
			
			a log file is generated and registered for transfer (so that pegasus won't skip it)
		
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments(self.javaPath, repr(job_max_memory), workflow.SplitReadFileJar, inputF, outputFnamePrefix, \
						repr(minNoOfReads))
		
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		if logFile:
			job.addArguments(logFile)
			job.uses(logFile, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
			job.output = logFile
		
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addRegisterAndMoveSplitFileJob(self, workflow=None, executable=None, \
							inputDir=None, outputDir=None, relativeOutputDir=None, logFile=None,\
							individual_sequence_id=None, bamFile=None, library=None, mate_id=None, \
							parentJobLs=[], job_max_memory=100, job_max_walltime = 60, \
							commit=0, sequence_format='fastq',\
							extraDependentInputLs=[], extraArguments=None, \
							transferOutput=False, sshDBTunnel=1, **keywords):
		"""
		2012.4.30
			add argument extraArguments, sshDBTunnel
		2012.1.24
			job_max_walltime is in minutes (max time allowed on hoffman2 is 24 hours).
			
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-v', self.drivername, '-z', self.hostname, '-d', self.dbname,\
						'-k', self.schema, '-u', self.db_user, '-p', self.db_passwd, \
						'-i', inputDir, '-o', outputDir, '-t', relativeOutputDir,'-l', library, '-n %s'%individual_sequence_id, \
						'-f', sequence_format)
		if commit:
			job.addArguments("-c")
		if mate_id:
			job.addArguments('-m %s'%(mate_id))
		if bamFile:
			job.addArguments('-a', bamFile)
			job.uses(bamFile, transfer=True, register=True, link=Link.INPUT)
		if logFile:
			job.addArguments('-g', logFile)
			job.uses(logFile, transfer=transferOutput, register=transferOutput, link=Link.OUTPUT)
			job.output = logFile
		if extraArguments:
			job.addArguments(extraArguments)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=job_max_walltime, sshDBTunnel=sshDBTunnel)
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def getInputFnameLsFromInput(self, input, suffixSet=set(['.fastq']), fakeSuffix='.gz'):
		"""
		2012.4.30
			this function supercedes self.getAllBamFiles() and it's more flexible.
		"""
		inputFnameLs = []
		if os.path.isdir(input):
			self.getFilesWithSuffixFromFolderRecursive(inputFolder=input, suffixSet=suffixSet, fakeSuffix=fakeSuffix, \
												inputFnameLs=inputFnameLs)
		elif os.path.isfile(input):
			inf = open(input)
			for line in inf:
				inputFnameLs.append(line.strip())
			del inf
		else:
			sys.stderr.write("%s is neither a folder nor a file.\n"%(input))
			sys.exit(4)
		return inputFnameLs
	
	def getMonkeyID2FastqObjectLs(self, fastqFnameLs=None,):
		"""
		2012.4.30
			each fastq file looks like 7_Index-11.2006013_R1.fastq.gz, 7_Index-11.2006013_R2.fastq.gz,
				7_Index-10.2005045_replacement_R1.fastq.gz, 7_Index-10.2005045_replacement_R2.fastq.gz
			
			The data from McGill is dated 2012.4.27.
			Each monkey is sequenced at 1X. There are 96 of them. Each library seems to contain 2 monkeys.
				But each monkey has only 1 library.
			
			8_Index_23.2008126_R1.fastq.gz
			8_Index_23.2008126_R2.fastq.gz
			8_Index_23.2009017_R1.fastq.gz
			8_Index_23.2009017_R2.fastq.gz

		"""
		sys.stderr.write("Passing monkeyID2FastqObjectLs from %s files ..."%(len(fastqFnameLs)))
		monkeyID2FastqObjectLs = {}
		import re, random
		monkeyIDPattern = re.compile(r'(?P<library>[-\w]+)\.(?P<monkeyID>\d+)((_replacement)|(_pool)|())_R(?P<mateID>\d).fastq.gz')
		counter = 0
		real_counter = 0
		library2UniqueLibrary = {}	#McGill's library ID , 7_Index-11, is not unique enough.
		for fastqFname in fastqFnameLs:
			counter += 1
			monkeyIDSearchResult = monkeyIDPattern.search(fastqFname)
			if monkeyIDSearchResult:
				real_counter += 1
				library = monkeyIDSearchResult.group('library')
				monkeyID = monkeyIDSearchResult.group('monkeyID')
				mateID = monkeyIDSearchResult.group('mateID')
				#concoct a unique library ID
				if library not in library2UniqueLibrary:
					library2UniqueLibrary[library] = '%s_%s'%(library, repr(random.random())[2:])
				uniqueLibrary = library2UniqueLibrary[library]
				fastqObject = PassingData(library=uniqueLibrary, monkeyID=monkeyID, mateID=mateID, absPath=fastqFname)
				if monkeyID not in monkeyID2FastqObjectLs:
					monkeyID2FastqObjectLs[monkeyID] = []
				monkeyID2FastqObjectLs[monkeyID].append(fastqObject)
			else:
				sys.stderr.write("Error: can't parse monkeyID, library, mateID out of %s.\n"%fastqFname)
				sys.exit(4)
		sys.stderr.write(" %s monkeys and %s files in the dictionary.\n"%(len(monkeyID2FastqObjectLs), real_counter))
		return monkeyID2FastqObjectLs
		
	
	def addJobsToProcessMcGillData(self, workflow, db_vervet=None, bamFname2MonkeyIDMapFname=None, input=None, dataDir=None, \
			minNoOfReads=None, commit=None,\
			sequencer=None, sequence_type=None, sequence_format=None):
		"""
		2012.4.30
		"""
		fastqFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.fastq']), fakeSuffix='.gz')
		monkeyID2FastqObjectLs = self.getMonkeyID2FastqObjectLs(fastqFnameLs)
		
		no_of_jobs = 0
		for monkeyID, fastqObjectLs in monkeyID2FastqObjectLs.iteritems():
			individual_sequence = self.addMonkeySequence(db_vervet, monkeyID, sequencer=sequencer, sequence_type=sequence_type, \
										sequence_format=sequence_format, dataDir=dataDir)
			
			sequenceOutputDir = os.path.join(dataDir, individual_sequence.path)
			sequenceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sequenceOutputDir)
			
			splitOutputDir = '%s'%(individual_sequence.id)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=splitOutputDir)
			
			no_of_jobs += 2
			
			
			for fastqObject in fastqObjectLs:
				library = fastqObject.library
				mateID = fastqObject.mateID
				fastqPath = fastqObject.absPath
				fastqFile = self.registerOneInputFile(workflow, fastqPath)
			
				splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mateID))
				logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mateID))
				splitReadFileJob1 = self.addSplitReadFileJob(workflow, executable=workflow.splitReadFile, \
								inputF=fastqFile, outputFnamePrefix=splitFastQFnamePrefix, \
								outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
								logFile=logFile, parentJobLs=[splitOutputDirJob], \
								job_max_memory=2000, job_max_walltime = 800, \
								extraDependentInputLs=[], transferOutput=True)
				
				logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mateID))
				# 2012.4.30 add '-w' to registerAndMoveSplitSequenceFiles so that it will be used to distinguish IndividualSequenceFileRaw
				registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.registerAndMoveSplitSequenceFiles, \
								inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
								individual_sequence_id=individual_sequence.id, bamFile=fastqFile, library=library, mate_id=mateID, \
								parentJobLs=[splitReadFileJob1, sequenceOutputDirJob], job_max_memory=100, job_max_walltime = 60, \
								commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], extraArguments='-w', \
								transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
				no_of_jobs += 2
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		
	
	def addJobsToProcessWUSTLData(self, workflow, db_vervet=None, bamFname2MonkeyIDMapFname=None, input=None, dataDir=None, \
			minNoOfReads=None, commit=None,\
			sequencer=None, sequence_type=None, sequence_format=None):
		"""
		2012.4.30
		"""
		bamBaseFname2MonkeyID = self.getBamBaseFname2MonkeyID(bamFname2MonkeyIDMapFname)
		bamFnameLs = self.getInputFnameLsFromInput(input, suffixSet=set(['.bam', '.sam']), fakeSuffix='.gz')
		
		sys.stderr.write("%s total bam files.\n"%(len(bamFnameLs)))
		
		sam2fastqOutputDir = 'sam2fastq'
		sam2fastqOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sam2fastqOutputDir)
		no_of_jobs = 1
		for bamFname in bamFnameLs:
			bamBaseFname = os.path.split(bamFname)[1]
			if bamBaseFname not in bamBaseFname2MonkeyID:
				sys.stderr.write("%s doesn't have monkeyID affiliated with.\n"%(bamFname))
				continue
			monkeyID = bamBaseFname2MonkeyID.get(bamBaseFname)
			individual_sequence = self.addMonkeySequence(db_vervet, monkeyID, sequencer=sequencer, sequence_type=sequence_type, \
										sequence_format=sequence_format, dataDir=dataDir)
			#2012.2.10 stop passing path_to_original_sequence=bamFname to self.addMonkeySequence()
			
			"""
			#2012.2.10 temporary, during transition from old records to new ones.
			newISQPath = individual_sequence.constructRelativePathForIndividualSequence()
			newISQPath = '%s_split'%(newISQPath)
			if individual_sequence.path is None or individual_sequence.path !=newISQPath:
				individual_sequence.path = newISQPath
				session.add(individual_sequence)
				session.flush()
			"""
			
			sequenceOutputDir = os.path.join(dataDir, individual_sequence.path)
			sequenceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=sequenceOutputDir)
			
			bamInputF = yh_pegasus.registerFile(workflow, bamFname)
			
			bamBaseFname = os.path.split(bamFname)[1]
			bamBaseFnamePrefix = os.path.splitext(bamBaseFname)[0]
			library = bamBaseFnamePrefix
			
			outputFnamePrefix = os.path.join(sam2fastqOutputDir, '%s_%s'%(individual_sequence.id, library))
			
			convertBamToFastqAndGzip_job = self.addConvertBamToFastqAndGzipJob(workflow, executable=workflow.convertBamToFastqAndGzip, \
							inputF=bamInputF, outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[sam2fastqOutputDirJob], job_max_memory=2000, job_max_walltime = 800, \
							extraDependentInputLs=[], \
							transferOutput=False)
			
			splitOutputDir = '%s_%s'%(individual_sequence.id, library)
			#same directory containing split files from both mates is fine as RegisterAndMoveSplitSequenceFiles could pick up.
			splitOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=splitOutputDir)
			
			
			mate_id = 1
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mate_id))
			logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
			splitReadFileJob1 = self.addSplitReadFileJob(workflow, executable=workflow.splitReadFile, \
							inputF=convertBamToFastqAndGzip_job.output1, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
							job_max_memory=2000, job_max_walltime = 800, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.registerAndMoveSplitSequenceFiles, \
							inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
							individual_sequence_id=individual_sequence.id, bamFile=bamInputF, library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob1, sequenceOutputDirJob], job_max_memory=100, job_max_walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			#handle the 2nd end
			mate_id = 2
			splitFastQFnamePrefix = os.path.join(splitOutputDir, '%s_%s_%s'%(individual_sequence.id, library, mate_id))
			logFile = File('%s_%s_%s.split.log'%(individual_sequence.id, library, mate_id))
			splitReadFileJob2 = self.addSplitReadFileJob(workflow, executable=workflow.splitReadFile, \
							inputF=convertBamToFastqAndGzip_job.output2, outputFnamePrefix=splitFastQFnamePrefix, \
							outputFnamePrefixTail="", minNoOfReads=minNoOfReads, \
							logFile=logFile, parentJobLs=[convertBamToFastqAndGzip_job, splitOutputDirJob], \
							job_max_memory=2000, job_max_walltime = 800, \
							extraDependentInputLs=[], transferOutput=True)
			
			logFile = File('%s_%s_%s.register.log'%(individual_sequence.id, library, mate_id))
			registerJob1 = self.addRegisterAndMoveSplitFileJob(workflow, executable=workflow.registerAndMoveSplitSequenceFiles, \
							inputDir=splitOutputDir, outputDir=sequenceOutputDir, relativeOutputDir=individual_sequence.path, logFile=logFile,\
							individual_sequence_id=individual_sequence.id, bamFile=bamInputF, library=library, mate_id=mate_id, \
							parentJobLs=[splitReadFileJob2, sequenceOutputDirJob], job_max_memory=100, job_max_walltime = 60, \
							commit=commit, sequence_format=sequence_format, extraDependentInputLs=[], \
							transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
			no_of_jobs += 5
			"""
			
			jobFname = os.path.join(self.jobFileDir, 'job%s.bam2fastq.sh'%(monkeyID))
			self.writeQsubJob(jobFname, bamFname, os.path.join(self.dataDir, individual_sequence.path), self.vervet_path)
			commandline = 'qsub %s'%(jobFname)
			if self.commit:	#qsub only when db transaction will be committed.
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			"""
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
	
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
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables(workflow=workflow)
		
		if self.inputType==1:
			self.addJobsToProcessWUSTLData(workflow, db_vervet=db_vervet, bamFname2MonkeyIDMapFname=self.bamFname2MonkeyIDMapFname, input=self.input, \
					dataDir=self.dataDir, minNoOfReads=self.minNoOfReads, commit=self.commit,\
					sequencer=self.sequencer, sequence_type=self.sequence_type, sequence_format=self.sequence_format)
		elif self.inputType==2:
			self.addJobsToProcessMcGillData(workflow, db_vervet=db_vervet, bamFname2MonkeyIDMapFname=None, input=self.input, \
					dataDir=self.dataDir, minNoOfReads=self.minNoOfReads, commit=self.commit,\
					sequencer=self.sequencer, sequence_type=self.sequence_type, sequence_format=self.sequence_format)
		else:
			sys.stderr.write("error: inputType %s not supported.\n"%(self.inputType))
			sys.exit(3)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		if self.commit:
			session.commit()
		else:
			session.rollback()
		
if __name__ == '__main__':
	main_class = UnpackAndAddIndividualSequence2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
