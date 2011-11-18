#!/usr/bin/env python
"""
Examples:
	%s
	
	dirPrefix=AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_single_sample_condor_20111105T0143/556_
	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o 4HighCovVRC_isq_15_18_vs_524_top804Contigs_gatk_vs_samtools_overlap_stat.xml -z uclaOffice -u yh -k genome
	
Description:
	2011-11-7 pegasus workflow that compares overlap between two vcf files (mapper/CheckTwoVCFOverlap.py),
			calculate mismatch rate, pi statistics based on the intersection
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

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *

class CheckTwoVCFOverlapPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [524, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('vcf1Dir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", 'S', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
		
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = ADAG(workflowName)
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		clusters_size = 40
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mkdirWrap)
		
		
		checkTwoVCFOverlap = Executable(namespace=namespace, name="CheckTwoVCFOverlap", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		checkTwoVCFOverlap.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CheckTwoVCFOverlap.py"), site_handler))
		checkTwoVCFOverlap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(checkTwoVCFOverlap)
		
		mergeSameHeaderTablesIntoOne = Executable(namespace=namespace, name="MergeSameHeaderTablesIntoOne", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		mergeSameHeaderTablesIntoOne.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/MergeSameHeaderTablesIntoOne.py"), site_handler))
		#long arguments will happen so no clustering
		#mergeSameHeaderTablesIntoOne.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mergeSameHeaderTablesIntoOne)
		
		statOutputDir = "statDir"
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=statOutputDir, namespace=namespace, version=version)
		
		"""
		inputFLs = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler)
		
		refSequence = VervetDB.IndividualSequence.get(self.aln_ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path.
		refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
		workflow.addFile(refFastaF)
		"""
		import re
		chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		mergeSameHeaderTablesIntoOneJob = Job(namespace=namespace, name=mergeSameHeaderTablesIntoOne.name, version=version)
		finalOutputF = File('overlapStat.tsv')
		mergeSameHeaderTablesIntoOneJob.addArguments('-o', finalOutputF)
		mergeSameHeaderTablesIntoOneJob.uses(finalOutputF, transfer=True, register=True, link=Link.OUTPUT)
		workflow.addJob(mergeSameHeaderTablesIntoOneJob)
		
		counter = 0
		for inputFname in os.listdir(self.vcf1Dir):
			gatkVCFAbsPath = os.path.join(os.path.abspath(self.vcf1Dir), inputFname)
			samtoolsVCFAbsPath = os.path.join(os.path.abspath(self.vcf2Dir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(gatkVCFAbsPath):
				if  not NextGenSeq.isVCFFileEmpty(samtoolsVCFAbsPath):	#make sure the samtools vcf exists
					chr = chr_pattern.search(inputFname).group(1)
					chr_size = chr2size.get(chr)
					if chr_size is None:
						sys.stderr.write("size for chr %s is unknown. set it to 1000.\n"%(chr))
						chr_size = 1000
					checkTwoVCFOverlapJob = Job(namespace=namespace, name=checkTwoVCFOverlap.name, version=version)
					
					vcf1 = File(os.path.join('vcf1', inputFname))	#relative path
					vcf1.addPFN(PFN("file://" + gatkVCFAbsPath, input_site_handler))
					workflow.addFile(vcf1)
					
					vcf2 = File(os.path.join('vcf2', inputFname))	#relative path
					vcf2.addPFN(PFN("file://" + samtoolsVCFAbsPath, input_site_handler))
					workflow.addFile(vcf2)
					
					outputF = File(os.path.join(statOutputDir, '%s_overlap.tsv'%(os.path.splitext(inputFname)[0])))
					
					checkTwoVCFOverlapJob.addArguments("-i", vcf1, "-j", vcf2, "-c", chr, "-l %s"%(chr_size), '-o', outputF)
					checkTwoVCFOverlapJob.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
					checkTwoVCFOverlapJob.uses(vcf2, transfer=True, register=True, link=Link.INPUT)
					checkTwoVCFOverlapJob.uses(outputF, transfer=False, register=True, link=Link.OUTPUT)
					
					
					job_max_memory=1000	#in MB. 
					checkTwoVCFOverlapJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
					checkTwoVCFOverlapJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
					workflow.addJob(checkTwoVCFOverlapJob)
					
					mergeSameHeaderTablesIntoOneJob.addArguments(outputF)
					mergeSameHeaderTablesIntoOneJob.uses(outputF, transfer=False, register=True, link=Link.INPUT)
					workflow.depends(parent=checkTwoVCFOverlapJob, child=mergeSameHeaderTablesIntoOneJob)
					counter += 1
		sys.stderr.write("%s jobs.\n"%(counter+1))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = CheckTwoVCFOverlapPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
