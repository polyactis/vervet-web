#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	%s -a 524 -I Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5/
		-C 1 -H -m1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -u yh -z localhost
		-o workflow/PairwiseDistance/PairwiseDistance_Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5.xml
		#-U 0 -Z 3000

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
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from vervet.src import VervetDB, AbstractVervetWorkflow

class CalculateDistanceMatrixFromVCFPipe(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('min_MAF', 1, float): [0.0, '', 1, 'minimum MAF for SNP filter', ],\
						('max_NA_rate', 1, float): [1, '', 1, 'maximum NA rate for SNP filter', ],\
						('convertHetero2NA', 1, int):[0, 'c', 1, 'convertHetero2NA mode. 0: no conversion, 1: convert to NA.'],\
						('hetHalfMatchDistance', 1, float): [0.5, 'q', 1, 'distance between two half-matched genotypes. AG vs A or AG vs AC', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 3000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		
		self.inputDir = os.path.abspath(self.inputDir)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		AggregateAndHClusterDistanceMatrix = Executable(namespace=namespace, name="AggregateAndHClusterDistanceMatrix", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AggregateAndHClusterDistanceMatrix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "reducer/AggregateAndHClusterDistanceMatrix.py"), \
													site_handler))
		executableClusterSizeMultiplierList.append((AggregateAndHClusterDistanceMatrix, 0))
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		parentPreReduceData = AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix, passingData=passingData, \
							transferOutput=transferOutput, **keywords)
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		callOutputDir = "call"
		callOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=callOutputDir)
		passingData.callOutputDirJob = callOutputDirJob
		
		matrixDir = "pairwiseDistMatrix"
		matrixDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=matrixDir)
		passingData.matrixDirJob = matrixDirJob
		
		reduceOutputDirJob = passingData.reduceOutputDirJob
		#2012.10.9 reduceOutputDirJob was added to passingData during AbstractVCFWorkflow.preReduce()
		
		#reduceOutputDir = "aggregateData"
		#reduceOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=reduceOutputDir)
		#passingData.reduceOutputDirJob = reduceOutputDirJob
		
		figureFnamePrefix = os.path.join(reduceOutputDirJob.output, 'aggregateDistanceMatrix')
		aggregateDistanceMatrixOutputF = File('%s.tsv'%(figureFnamePrefix))
		PCAFile = File('%s_PCA.tsv'%(figureFnamePrefix))
		aggregateAndHClusterDistanceMatrixJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.AggregateAndHClusterDistanceMatrix, \
									outputF=aggregateDistanceMatrixOutputF, \
									parentJobLs=[reduceOutputDirJob],extraOutputLs=[PCAFile, File('%s.png'%(figureFnamePrefix)), \
																				File('%s.svg'%(figureFnamePrefix))], \
									extraDependentInputLs=[], transferOutput=True, extraArguments="-f %s"%(figureFnamePrefix))
		returnData.aggregateAndHClusterDistanceMatrixJob = aggregateAndHClusterDistanceMatrixJob
		
		#2012.9.5 add the job to append meta info (country, sex, latitude, etc. of each monkey)
		outputF = File('%s_withMetaInfo.tsv'%(figureFnamePrefix))
		appendInfo2PCAOutputJob = self.addGenericDBJob(executable=self.AppendInfo2SmartPCAOutput, inputFile=PCAFile, \
				outputFile=outputF, \
				parentJobLs=[aggregateAndHClusterDistanceMatrixJob], extraDependentInputLs=None, \
				extraOutputLs=None,\
				transferOutput=True, \
				extraArgumentList=None, extraArguments=None, sshDBTunnel=self.needSSHDBTunnel, \
				key2ObjectForJob=None, job_max_memory=2000)
		
		
		return returnData
	
	def mapEachInterval(self, workflow=None, \
					VCFFile=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.9.22
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		intervalFnamePrefix = passingData.intervalFnamePrefix
		jobData = passingData.jobData
		callOutputDirJob = passingData.callOutputDirJob
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		
		genotypeCallOutputFname = os.path.join(callOutputDirJob.output, '%s.call'%(intervalFnamePrefix))
		genotypeCallOutput = File(genotypeCallOutputFname)
		genotypeCallByCoverage_job = self.addVCF2MatrixJob(workflow, executable=self.GenotypeCallByCoverage, \
											inputVCF=VCFFile, outputFile=genotypeCallOutput, \
					refFastaF=None, run_type=3, numberOfReadGroups=10, minDepth=self.minDepth,\
					parentJobLs=[callOutputDirJob, splitVCFJob]+jobData.jobLs, extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000)
		
		matrixDirJob = passingData.matrixDirJob
		calculaOutputFname =os.path.join(matrixDirJob.output, '%s.pairwiseDist.convertHetero2NA%s.minMAF%.2f.maxNA%.2f.tsv'%(intervalFnamePrefix, \
							self.convertHetero2NA, self.min_MAF, self.max_NA_rate))
		calculaOutput = File(calculaOutputFname)
		calculaJob = self.addCalculatePairwiseDistanceFromSNPXStrainMatrixJob(workflow, \
										executable=self.CalculatePairwiseDistanceOutOfSNPXStrainMatrix, \
										inputFile=genotypeCallOutput, outputFile=calculaOutput, \
					min_MAF=self.min_MAF, max_NA_rate=self.max_NA_rate, convertHetero2NA=self.convertHetero2NA, \
					hetHalfMatchDistance=self.hetHalfMatchDistance,\
					parentJobLs=[genotypeCallByCoverage_job, matrixDirJob], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000)
		returnData.jobDataLs.append(PassingData(jobLs=[calculaJob], file=calculaJob.output, \
											fileList=[calculaJob.output]))
		returnData.calculaJob = calculaJob
		return returnData
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		for jobData in mapEachIntervalData.jobDataLs:
			calculaJob = jobData.jobLs[0]
			self.addInputToStatMergeJob(workflow, statMergeJob=preReduceReturnData.aggregateAndHClusterDistanceMatrixJob, \
						inputF=calculaJob.output, \
						parentJobLs=[calculaJob])
		return returnData
	
	"""
	2011.9-28
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = self.db_vervet
		
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
		
		topOutputDir = "aggregateData"
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		figureFnamePrefix = os.path.join(topOutputDir, 'aggregateDistanceMatrix')
		aggregateDistanceMatrixOutputF = File('%s.tsv'%(figureFnamePrefix))
		aggregateAndHClusterDistanceMatrixJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.AggregateAndHClusterDistanceMatrix, \
									outputF=aggregateDistanceMatrixOutputF, \
									parentJobLs=[topOutputDirJob], \
									extraDependentInputLs=[], transferOutput=True, extraArguments="-f %s"%(figureFnamePrefix))
		aggregateAndHClusterDistanceMatrixJob.uses('%s.png'%(figureFnamePrefix), transfer=True, register=True, link=Link.OUTPUT)
		aggregateAndHClusterDistanceMatrixJob.uses('%s.svg'%(figureFnamePrefix), transfer=True, register=True, link=Link.OUTPUT)
		PCAFile = File('%s_PCA.tsv'%(figureFnamePrefix))
		aggregateAndHClusterDistanceMatrixJob.uses(PCAFile, transfer=True, register=True, link=Link.OUTPUT)
		
		#2012.9.5 add the job to append meta info (country, sex, latitude, etc. of each monkey)
		outputF = File('%s_withMetaInfo.tsv'%(figureFnamePrefix))
		appendInfo2PCAOutputJob = self.addGenericJob(executable=self.AppendInfo2SmartPCAOutput, inputFile=PCAFile, \
				outputFile=outputF, \
				parentJobLs=[aggregateAndHClusterDistanceMatrixJob, topOutputDirJob], extraDependentInputLs=None, \
				extraOutputLs=None,\
				transferOutput=True, \
				extraArgumentList=None, extraArguments=None, key2ObjectForJob=None, job_max_memory=2000)
		self.addDBArgumentsToOneJob(job=appendInfo2PCAOutputJob, objectWithDBArguments=self)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler,
								checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName=self.pegasusFolderName)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		
		seqCoverageF = None
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			#add the cover filter (or filter+call) job after index is done
			inputFBaseName = os.path.basename(inputF.name)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(inputFBaseName))
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job = self.addVCF2MatrixJob(workflow, executable=workflow.GenotypeCallByCoverage, \
															inputVCF=inputF, outputFile=genotypeCallOutput, \
						refFastaF=None, run_type=3, numberOfReadGroups=10, minDepth=self.minDepth,\
						parentJobLs=[callOutputDirJob]+jobData.jobLs, extraDependentInputLs=[], transferOutput=False, \
						extraArguments=" --minDepth %s "%(self.minDepth), job_max_memory=2000)
			
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
		
	"""

	
if __name__ == '__main__':
	main_class = CalculateDistanceMatrixFromVCFPipe
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
