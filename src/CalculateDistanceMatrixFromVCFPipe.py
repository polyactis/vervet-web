#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -i 8GenomeVsTop156Contigs_GATK_all_bases/call/ -n 8 -s 1 -m 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -i 8GenomeVsTop156Contigs_GATK_all_bases/call/ -n 8 -s 1 -m 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
Description:
	2011-10-14
		a program which generates a pegasus workflow dag (xml file) to 
		
		1. filter VCF file
		2. calculate pairwise distance matrix for each vcf file
	
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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *


class CalculateDistanceMatrixFromVCFPipe(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('inputDir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", 'S', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('numberOfReadGroups', 1, int): [None, 'n', 1, 'number of read groups/genomes in the inputFname', ],\
						('min_MAF', 1, float): [0.0, 'M', 1, 'minimum MAF for SNP filter', ],\
						('max_NA_rate', 1, float): [0.4, 'm', 1, 'maximum NA rate for SNP filter', ],\
						('convertHetero2NA', 1, int):[0, 'c', 1, 'convertHetero2NA mode. 0: no conversion, 1: convert to NA.'],\
						('defaultCoverage', 1, float): [9, 'f', 1, 'default coverage when coverage is not available for a read group'],\
						('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\
						("site_type", 1, int): [2, 's', 1, '1: all genome sites, 2: variants only'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.inputDir = os.path.abspath(self.inputDir)
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	
	def registerAllInputFiles(self, workflow, inputDir, input_site_handler=None):
		"""
		2011-9-29
			vcf files only
		"""
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		inputFLs = []
		fnameLs = os.listdir(inputDir)
		for fname in fnameLs:
			if fname[-6:]!='vcf.gz' and fname[-3:]!='vcf':	#ignore non-vcf files
				continue
			inputFname = os.path.join(inputDir, fname)
			
			inputF = File(os.path.basename(inputFname))
			inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
			workflow.addFile(inputF)
			inputFLs.append(inputF)
		sys.stderr.write("%s files.\n"%(len(inputFLs)))
		return inputFLs
		
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
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
		clusters_size = 20
		
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
		
		genotypeCallByCoverage = Executable(namespace=namespace, name="GenotypeCallByCoverage", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		workflow.addExecutable(genotypeCallByCoverage)
		
		calcula = Executable(namespace=namespace, name="CalculatePairwiseDistanceOutOfSNPXStrainMatrix", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		calcula.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py"), site_handler))
		workflow.addExecutable(calcula)
		
		callOutputDir = "call"
		callOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=callOutputDir, namespace=namespace, version=version)
		matrixDir = "pairwiseDistMatrix"
		matrixDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=matrixDir, namespace=namespace, version=version)
		
		inputFLs = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler)
		
		refSequence = VervetDB.IndividualSequence.get(self.aln_ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaF = File(os.path.basename(refFastaFname))	#use relative path, otherwise, it'll go to absolute path.
		refFastaF.addPFN(PFN("file://" + refFastaFname, self.input_site_handler))
		workflow.addFile(refFastaF)
		
		
		"""
		#2011-9-2
		self.outputSeqCoverage(self.seqCoverageFname)
		seqCoverageF = File(os.path.basename(self.seqCoverageFname))
		seqCoverageF.addPFN(PFN("file://" + os.path.abspath(self.seqCoverageFname), \
											self.input_site_handler))
		workflow.addFile(seqCoverageF)
		"""
		seqCoverageF = None
		
		for inputF in inputFLs:
			#add the cover filter (or filter+call) job after index is done
			genotypeCallByCoverage_job = Job(namespace=namespace, name=genotypeCallByCoverage.name, version=version)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(inputF.name))
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job.addArguments("-i", inputF, "-n", str(self.numberOfReadGroups), \
					"-o", genotypeCallOutput, '-e', refFastaF, '-y', str(self.genotypeCallerType), \
					'-s', repr(self.site_type), "-f", str(self.defaultCoverage))
			if seqCoverageF:
				genotypeCallByCoverage_job.addArguments("-q", seqCoverageF)
				genotypeCallByCoverage_job.uses(seqCoverageF, transfer=True, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(genotypeCallOutput, transfer=False, register=True, link=Link.OUTPUT)
			if self.site_type==1:	#all sites require additional memory
				job_max_memory = 5000	#in MB
			else:	#variants only, less memory
				job_max_memory=1000	#in MB. 
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
			genotypeCallByCoverage_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			workflow.addJob(genotypeCallByCoverage_job)
			workflow.depends(parent=callOutputDirJob, child=genotypeCallByCoverage_job)
			
			
			#add the pairwise distance matrix job after filter is done
			calcula_job = Job(namespace=namespace, name=calcula.name, version=version)
			
			calculaOutputFname =os.path.join(matrixDir, '%s.pairwiseDist.convertHetero2NA%s.minMAF%.2f.maxNA%.2f.tsv'%(inputF.name, \
								self.convertHetero2NA, self.min_MAF, self.max_NA_rate))
			calculaOutput = File(calculaOutputFname)
			calcula_job.addArguments("-i", genotypeCallOutput, "-n", str(self.min_MAF), \
						"-o", calculaOutput, '-m', repr(self.max_NA_rate), '-c', str(self.convertHetero2NA))
			calcula_job.uses(genotypeCallOutput, transfer=False, register=False, link=Link.INPUT)
			calcula_job.uses(calculaOutput, transfer=True, register=False, link=Link.OUTPUT)
			
			workflow.addJob(calcula_job)
			workflow.depends(parent=genotypeCallByCoverage_job, child=calcula_job)
			workflow.depends(parent=matrixDirJob, child=calcula_job)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = CalculateDistanceMatrixFromVCFPipe
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
