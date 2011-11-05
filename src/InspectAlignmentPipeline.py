#!/usr/bin/env python
"""
Examples:
	
	#2011-11-4 run on condorpool
	%s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice -o InspectAlnRefSeq524Alignments.xml -I 552-661
	
Description:
	2011-11-4
		a pegasus workflow that inspects no-of-reads-aligned, inferred insert size and etc.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])	# sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0]

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from AlignmentToCallPipeline import AlignmentToCallPipeline

class InspectAlignmentPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("samtools_path", 1, ): ["%s/bin/samtools", '', 1, 'samtools binary'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						("topNumberOfContigs", 1, int): [156, 'N', 1, 'number of contigs'],\
						("needFastaIndexJob", 0, int): [0, '', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, '', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-11-4
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
		if self.ind_aln_id_ls:
			self.ind_aln_id_ls = getListOutOfStr(self.ind_aln_id_ls, data_type=int)
		
		self.samtools_path = self.samtools_path%self.home_path
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	def addJobs(self, workflow, alignmentLs, refName2size, samtools=None, DOCWalkerJava=None, \
				VariousReadCountJava=None,  genomeAnalysisTKJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mkdirWrap=None, mv=None, \
				reduceDepthOfCoverage=None, reduceVariousReadCount=None,\
				refFastaF=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				dataDir=None):
		"""
		2011-11-4
			
		"""
		sys.stderr.write("Adding jobs for %s references and %s alignments..."%(len(refName2size), len(alignmentLs)))
		job_max_memory = 1000	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		vcf_job_max_memory = 1000
		if needFastaDictJob:	# the .dict file is required for GATK
			refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
			refFastaDictF = File(refFastaDictFname)
			#not os.path.isfile(refFastaDictFname) or 
			fastaDictJob = Job(namespace=namespace, name=createSequenceDictionaryJava.name, version=version)
			fastaDictJob.addArguments('-jar', createSequenceDictionaryJar, \
					'REFERENCE=', refFastaF, 'OUTPUT=', refFastaDictF)
			fastaDictJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaDictJob.uses(refFastaDictF, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaDictJob)
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
			refFastaIndexF = File(refFastaIndexFname)
			#not os.path.isfile(refFastaIndexFname)
			fastaIndexJob = Job(namespace=namespace, name=samtools.name, version=version)
			fastaIndexJob.addArguments("faidx", refFastaF)
			fastaIndexJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaIndexJob.uses(refFastaIndexFname, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaIndexJob)
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
		
		
		returnData = self.addAddRG2BamJobsAsNeeded(workflow, alignmentLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		
		for alignment in alignmentLs:
			if alignment.id in alignmentId2RGJobDataLs:
				statOutputDir = '%s'%(alignment.id)
				statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=statOutputDir, namespace=namespace, \
														version=version)
				finalVariousReadCountOutputF = File("%s_VariousReadCount.tsv"%(alignment.getReadGroup()))
				reduceVariousReadCountJob = Job(namespace=namespace, name=reduceVariousReadCount.name, version=version)
				reduceVariousReadCountJob.addArguments("-o", finalVariousReadCountOutputF)
				reduceVariousReadCountJob.uses(finalVariousReadCountOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(reduceVariousReadCountJob)
				
				finalDepthOfCoverageOutputF = File("%s_DepthOfCoverage.tsv"%(alignment.getReadGroup()))
				reduceDepthOfCoverageJob = Job(namespace=namespace, name=reduceDepthOfCoverage.name, version=version)
				reduceDepthOfCoverageJob.addArguments("-o", finalDepthOfCoverageOutputF)
				reduceDepthOfCoverageJob.uses(finalDepthOfCoverageOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(reduceDepthOfCoverageJob)
				
				for refName, refSize in refName2size.iteritems():
				
					index_sam_job, bamF, bai_output = alignmentId2RGJobDataLs[alignment.id]
					DOCOutputFnamePrefix = os.path.join(statOutputDir, '%s_%s_DOC'%(alignment.id, refName))
					DOCOutputF = File('%s.sample_interval_summary'%(DOCOutputFnamePrefix))
					DOCJob = Job(namespace=namespace, name=DOCWalkerJava.name, version=version)
					DOCJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "DepthOfCoverage",\
						'-R', refFastaF, '-o', DOCOutputFnamePrefix, "-pt sample", "--omitDepthOutputAtEachBase",\
						"-L", refName, "-mmq 30")
					DOCJob.addArguments("-I", bamF)
					workflow.addJob(DOCJob)
					
					readCountOutputF = File(os.path.join(statOutputDir, '%s_%s_variousReadCount.tsv'%(alignment.id, refName)))
					readCountJob = Job(namespace=namespace, name=VariousReadCountJava.name, version=version)
					readCountJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "VariousReadCount",\
						'-R', refFastaF, '-o', readCountOutputF, \
						"-L", refName, "-mmq 30")
					readCountJob.addArguments("-I", bamF)
					workflow.addJob(readCountJob)
					
					#it's either symlink or stage-in
					DOCJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
					DOCJob.uses(bai_output, transfer=True, register=True, link=Link.INPUT)
					DOCJob.uses(DOCOutputF, transfer=False, register=True, link=Link.OUTPUT)
					reduceDepthOfCoverageJob.addArguments(DOCOutputF)
					reduceDepthOfCoverageJob.uses(DOCOutputF, transfer=False, register=True, link=Link.INPUT)
					
					readCountJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
					readCountJob.uses(bai_output, transfer=True, register=True, link=Link.INPUT)
					readCountJob.uses(readCountOutputF, transfer=False, register=True, link=Link.OUTPUT)
					reduceVariousReadCountJob.addArguments(readCountOutputF)
					reduceVariousReadCountJob.uses(readCountOutputF, transfer=False, register=True, link=Link.INPUT)
					
					workflow.depends(parent=statOutputDirJob, child=DOCJob)
					workflow.depends(parent=DOCJob, child=reduceDepthOfCoverageJob)
					workflow.depends(parent=statOutputDirJob, child=readCountJob)
					workflow.depends(parent=readCountJob, child=reduceVariousReadCountJob)
					if index_sam_job is not None:
						workflow.depends(parent=index_sam_job, child=DOCJob)
						workflow.depends(parent=index_sam_job, child=readCountJob)
		sys.stderr.write(".Done\n")
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		refName2size = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refName2size = set(['Contig149'])	#temporary when testing Contig149
		#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = refName2size.keys()
		
		alignmentLs = self.getAlignments(self.aln_ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2, dataDir=self.localDataDir)
		
		refSequence = VervetDB.IndividualSequence.get(self.aln_ref_ind_seq_id)
		
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		if self.needFastaDictJob and self.needFastaIndexJob:
			refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path
			# Add it into replica only when needed.
			refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
			workflow.addFile(refFastaF)
			# If it's not needed, assume the index is done and all relevant files are in absolute path.
			# and no replica transfer
		else:
			refFastaF = File(refFastaFname)
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = ADAG(workflowName)
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		abs_path = os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar')
		genomeAnalysisTKJar = File(abs_path)
		genomeAnalysisTKJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(genomeAnalysisTKJar)
		
		abs_path = os.path.join(self.picard_path, 'CreateSequenceDictionary.jar')
		createSequenceDictionaryJar = File(abs_path)
		createSequenceDictionaryJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(createSequenceDictionaryJar)
		
		abs_path = os.path.join(self.picard_path, 'BuildBamIndex.jar')
		BuildBamIndexFilesJar = File(abs_path)
		BuildBamIndexFilesJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(BuildBamIndexFilesJar)
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroupsAndCleanSQHeader.jar')
		addOrReplaceReadGroupsAndCleanSQHeaderJar = File(abs_path)
		addOrReplaceReadGroupsAndCleanSQHeaderJar.addPFN(PFN("file://" + abs_path, \
												site_handler))
		workflow.addFile(addOrReplaceReadGroupsAndCleanSQHeaderJar)
		
		abs_path = os.path.join(self.picard_path, 'AddOrReplaceReadGroups.jar')
		addOrReplaceReadGroupsJar = File(abs_path)
		addOrReplaceReadGroupsJar.addPFN(PFN("file://" + abs_path, \
											site_handler))
		workflow.addFile(addOrReplaceReadGroupsJar)
		
		"""
		#2011-9-2
		self.outputSeqCoverage(self.seqCoverageFname)
		seqCoverageF = File(os.path.basename(self.seqCoverageFname))
		seqCoverageF.addPFN(PFN("file://" + os.path.abspath(self.seqCoverageFname), \
											self.input_site_handler))
		workflow.addFile(seqCoverageF)
		"""
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		clusters_size = 5
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mkdirWrap)
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		mv.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mv)
		
		
		"""
		# 2011-9-1 not used anymore. too slow. replace with addOrReplaceReadGroupsAndCleanSQHeaderJar.
		addRGAndCleanSQHeaderAlignment = Executable(namespace=namespace, name="AddRGAndCleanSQHeaderAlignment", version=version, \
												os=operatingSystem, arch=architecture, installed=True)
		addRGAndCleanSQHeaderAlignment.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AddRGAndCleanSQHeaderAlignment.py"), site_handler))
		workflow.addExecutable(addRGAndCleanSQHeaderAlignment)
		"""
		
		samtools = Executable(namespace=namespace, name="samtools", version=version, os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
		
		java = Executable(namespace=namespace, name="java", version=version, os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(java)
		
		createSequenceDictionaryJava = Executable(namespace=namespace, name="createSequenceDictionaryJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		createSequenceDictionaryJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(createSequenceDictionaryJava)
		
		
		DOCWalkerJava = Executable(namespace=namespace, name="DepthOfCoverageWalker", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		DOCWalkerJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		DOCWalkerJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(DOCWalkerJava)
		
		VariousReadCountJava = Executable(namespace=namespace, name="VariousReadCountJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		VariousReadCountJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		VariousReadCountJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(VariousReadCountJava)
		
		addOrReplaceReadGroupsJava = Executable(namespace=namespace, name="addOrReplaceReadGroupsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		addOrReplaceReadGroupsJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(addOrReplaceReadGroupsJava)
		
		BuildBamIndexFilesJava = Executable(namespace=namespace, name="BuildBamIndexFilesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		BuildBamIndexFilesJava.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(BuildBamIndexFilesJava)
		
		reduceDepthOfCoverage = Executable(namespace=namespace, name="ReduceDepthOfCoverage", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceDepthOfCoverage.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reduce/ReduceDepthOfCoverage.py"), site_handler))
		#reduceDepthOfCoverage.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(reduceDepthOfCoverage)
		
		reduceVariousReadCount = Executable(namespace=namespace, name="ReduceVariousReadCount", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		reduceVariousReadCount.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reduce/ReduceVariousReadCount.py"), site_handler))
		#reduceVariousReadCount.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(reduceVariousReadCount)
		
		if self.site_handler=='condorpool':
			bamListF_site_handler = 'condorpool'
		else:	#every other site requires transfer from local
			bamListF_site_handler = 'local'
		#bamListF = self.outputListOfBam(alignmentLs, self.dataDir, self.bamListFname, workflow=workflow, site_handler=bamListF_site_handler)
		
		
		self.addJobs(workflow, alignmentLs, refName2size, samtools=samtools, DOCWalkerJava=DOCWalkerJava, \
				VariousReadCountJava=VariousReadCountJava,  genomeAnalysisTKJar=genomeAnalysisTKJar, \
				createSequenceDictionaryJava=createSequenceDictionaryJava, createSequenceDictionaryJar=createSequenceDictionaryJar, \
				addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
				BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar,\
				mkdirWrap=mkdirWrap, mv=mv, \
				reduceDepthOfCoverage=reduceDepthOfCoverage, reduceVariousReadCount=reduceVariousReadCount,\
				refFastaF=refFastaF, \
				namespace=namespace, version=version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				dataDir=self.dataDir)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = InspectAlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
