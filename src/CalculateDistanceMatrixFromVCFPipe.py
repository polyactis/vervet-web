#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -i 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -i 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1)
	%s -a 524 -i Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5/
		-C 1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -u yh -z localhost
		-o workflow/PairwiseDistance_Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5.xml

Description:
	2011-10-14
		a program which generates a pegasus workflow dag (xml file) to 
		
		1. filter VCF file
		2. calculate pairwise distance matrix for each vcf file
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow


class CalculateDistanceMatrixFromVCFPipe(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('min_MAF', 1, float): [0.0, 'M', 1, 'minimum MAF for SNP filter', ],\
						('max_NA_rate', 1, float): [0.4, 'N', 1, 'maximum NA rate for SNP filter', ],\
						('convertHetero2NA', 1, int):[0, 'c', 1, 'convertHetero2NA mode. 0: no conversion, 1: convert to NA.'],\
						('hetHalfMatchDistance', 1, float): [0.5, 'q', 1, 'distance between two half-matched genotypes. AG vs A or AG vs AC', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		
		self.inputDir = os.path.abspath(self.inputDir)
	
	def registerCustomExecutables(self, workflow):
		"""
		2011-11-28
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		
		AggregateAndHClusterDistanceMatrix = Executable(namespace=namespace, name="AggregateAndHClusterDistanceMatrix", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AggregateAndHClusterDistanceMatrix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "AggregateAndHClusterDistanceMatrix.py"), \
													site_handler))
		AggregateAndHClusterDistanceMatrix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(AggregateAndHClusterDistanceMatrix)
		workflow.AggregateAndHClusterDistanceMatrix = AggregateAndHClusterDistanceMatrix
	
	
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
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		callOutputDir = "call"
		callOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=callOutputDir)
		matrixDir = "pairwiseDistMatrix"
		matrixDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=matrixDir)
		
		aggregateDistanceMatrixOutputF = File('aggregateDistanceMatrix.tsv')
		figureFnamePrefix = 'aggregateDistanceMatrix'
		aggregateAndHClusterDistanceMatrixJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.AggregateAndHClusterDistanceMatrix, \
									outputF=aggregateDistanceMatrixOutputF, \
									parentJobLs=[], \
									extraDependentInputLs=[], transferOutput=True, extraArguments="-f %s"%(figureFnamePrefix))
		aggregateAndHClusterDistanceMatrixJob.uses('%s.png'%(figureFnamePrefix), transfer=True, register=True, link=Link.OUTPUT)
		aggregateAndHClusterDistanceMatrixJob.uses('%s.svg'%(figureFnamePrefix), transfer=True, register=True, link=Link.OUTPUT)
		aggregateAndHClusterDistanceMatrixJob.uses('%s_PCA.tsv'%(figureFnamePrefix), transfer=True, register=True, link=Link.OUTPUT)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler,
								checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName=self.pegasusFolderName)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		
		"""
		#2011-9-2
		self.outputSeqCoverage(self.seqCoverageFname)
		seqCoverageF = File(os.path.basename(self.seqCoverageFname))
		seqCoverageF.addPFN(PFN("file://" + os.path.abspath(self.seqCoverageFname), \
											self.input_site_handler))
		workflow.addFile(seqCoverageF)
		"""
		seqCoverageF = None
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			#add the cover filter (or filter+call) job after index is done
			inputFBaseName = os.path.basename(inputF.name)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(inputFBaseName))
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job = self.addVCF2MatrixJob(workflow, executable=workflow.genotypeCallByCoverage, \
															inputVCF=inputF, outputFile=genotypeCallOutput, \
						refFastaF=None, run_type=3, numberOfReadGroups=10, \
						parentJobLs=[callOutputDirJob]+jobData.jobLs, extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000)
			
			calculaOutputFname =os.path.join(matrixDir, '%s.pairwiseDist.convertHetero2NA%s.minMAF%.2f.maxNA%.2f.tsv'%(inputFBaseName, \
								self.convertHetero2NA, self.min_MAF, self.max_NA_rate))
			calculaOutput = File(calculaOutputFname)
			calcula_job = self.addCalculatePairwiseDistanceFromSNPXStrainMatrixJob(workflow, executable=workflow.calcula, \
															inputFile=genotypeCallOutput, outputFile=calculaOutput, \
						min_MAF=self.min_MAF, max_NA_rate=self.max_NA_rate, convertHetero2NA=self.convertHetero2NA, \
						hetHalfMatchDistance=self.hetHalfMatchDistance,\
						parentJobLs=[genotypeCallByCoverage_job, matrixDirJob], extraDependentInputLs=[], transferOutput=True, \
						extraArguments=None, job_max_memory=2000)
			
			self.addInputToStatMergeJob(workflow, statMergeJob=aggregateAndHClusterDistanceMatrixJob, inputF=calculaOutput, \
							parentJobLs=[calcula_job])
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = CalculateDistanceMatrixFromVCFPipe
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
