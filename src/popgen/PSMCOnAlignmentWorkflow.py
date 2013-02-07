#!/usr/bin/env python
"""
Examples:
	%s
	
	#2013.2.4
	%s --ind_aln_id_ls 929,1508,1510,926,1577,1507,1534,1524,1512
		-u yh -a 524  -z localhost -o  dags/PSMCOnAlignment/PSMCOnAlignmentOnSubspecies.xml
		-j hcondor -l hcondor
		--clusters_size 1 -e /u/home/eeskin/polyacti/
		--dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2
		# -Y
		#--minContigID 96 --maxContigID 100	#useless right now
	
Description:
	2013.1.25 run psmc (Li, Durbin 2011) on alignments
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from pymodule import VCFFile
from vervet.src import VervetDB, AbstractVervetAlignmentWorkflow

class PSMCOnAlignmentWorkflow(AbstractVervetAlignmentWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetAlignmentWorkflow.option_default_dict)
	option_default_dict.update({
					("psmcFolderPath", 1, ): ["%s/script/psmc", '', 1, \
							'path to the folder that contains the source/binaries of PSMC (Li, Durbin 2011)'],\
					('patternOfPSMCTimeIntervals', 1, ):["4+25*2+4+6", '', 1, 'pattern of parameters, -p of psmc '],\
					('maxNoOfIterations', 1, int):[25, '', 1, 'number of iterations for the EM algorithm, -N of psmc'],\
					('initThetaRhoRatio', 1, float):[4, '', 1, 'initial theta/rho ratio, -r of psmc'],\
					('max2N0CoalescentTime', 1, float):[25, '', 1, 'maximum 2N0 coalescent time, -t of psmc '],\
					
					('absMutationRatePerNucleotide', 1, float):[5e-09, '', 1, 'absolute mutation rate per nucleotide per generation, -u FLOAT of psmc_plot'],\
					('maxNoOfGenerations', 1, int):[0, '', 1, 'maximum generations, 0 for auto , -X of psmc_plot'],\
					('minNoOfGenerations', 1, int):[10000, '', 1, 'minimum number of generations, 0 for auto , -x of psmc_plot'],\
					('noOfYearsPerGeneration', 1, int):[5, '', 1, 'number of years per generation, -g of psmc_plot.pl'],\
					('maxPopulationSize', 1, int):[30, '', 1, 'maximum popsize in terms of 10^4, 0 for auto, -Y of psmc_plot.pl'],\
					
					('bootstrapTrunkSize', 1, int):[50000, '', 1, 'number of 100bp fragments for each bootstrap psmc input. option for splitfa'],\
					('noOfBootstraps', 1, int):[100, '', 1, 'number of bootstraps'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2013.1.25
		"""
		
		AbstractVervetAlignmentWorkflow.__init__(self, **keywords)
		
		if hasattr(self, 'psmcFolderPath'):
			self.psmcFolderPath =  self.insertHomePath(self.psmcFolderPath, self.home_path)
		else:
			self.psmcFolderPath = None
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		2013.1.28
		"""
		AbstractVervetAlignmentWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		
		
		psmc = Executable(namespace=namespace, name="psmc", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		psmc.addPFN(PFN("file://" + os.path.join(self.psmcFolderPath, "psmc"), site_handler))
		executableClusterSizeMultiplierList.append((psmc, 1))

		fq2psmcfa = Executable(namespace=namespace, name="fq2psmcfa", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		fq2psmcfa.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "popgen/fq2psmcfa.sh"), site_handler))
		executableClusterSizeMultiplierList.append((fq2psmcfa, 1))
		
		psmc2history = Executable(namespace=namespace, name="psmc2history", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		psmc2history.addPFN(PFN("file://" + os.path.join(self.psmcFolderPath, "utils/psmc2history.pl"), site_handler))
		executableClusterSizeMultiplierList.append((psmc2history, 1))
		
		history2ms = Executable(namespace=namespace, name="history2ms", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		history2ms.addPFN(PFN("file://" + os.path.join(self.psmcFolderPath, "utils/history2ms.pl"), site_handler))
		executableClusterSizeMultiplierList.append((history2ms, 1))
		
		GenerateMSCommandFromPSMCOutput = Executable(namespace=namespace, name="GenerateMSCommandFromPSMCOutput", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		GenerateMSCommandFromPSMCOutput.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "popgen/GenerateMSCommandFromPSMCOutput.sh"), site_handler))
		executableClusterSizeMultiplierList.append((GenerateMSCommandFromPSMCOutput, 1))
		
		psmc_plot = Executable(namespace=namespace, name="psmc_plot", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		psmc_plot.addPFN(PFN("file://" + os.path.join(self.psmcFolderPath, "utils/psmc_plot.pl"), site_handler))
		executableClusterSizeMultiplierList.append((psmc_plot, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def mapEachAlignment(self, workflow=None, alignmentData=None, passingData=None, transferOutput=True, **keywords):
		"""
		2013.1.25
			run fq2psmcfa, psmc, ms-simulation-commandline, plot
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
		"""					

		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		refFastaFile = passingData.refFastaFList[0]
		
		alignmentData = passingData.alignmentData
		
		alignment = alignmentData.alignment
		parentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		
		bamFnamePrefix = alignment.getReadGroup()
		
		minDP = max(1, alignment.median_depth/2)
		maxDP = alignment.median_depth*2
		consensusSequenceFile = File(os.path.join(topOutputDirJob.output, '%s.fq.gz'%(passingData.bamFnamePrefix)))
		ExtractConsensusSequenceFromAlignmentJob = self.addGenericJob(executable=self.ExtractConsensusSequenceFromAlignment, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=parentJobLs+[topOutputDirJob], \
					extraDependentInputLs=[refFastaFile, bamF, baiF], extraOutputLs=[consensusSequenceFile], transferOutput=True, \
					extraArguments=None, extraArgumentList=[refFastaFile, bamF, consensusSequenceFile, '%s %s'%(minDP, maxDP)], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=3000)	#max_walltime is in minutes, 50 hours
		
		psmcInputFile = File(os.path.join(topOutputDirJob.output, '%s.psmcfa'%(passingData.bamFnamePrefix)))
		fq2psmcfaJob = self.addGenericJob(executable=self.fq2psmcfa, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=[ExtractConsensusSequenceFromAlignmentJob, topOutputDirJob], \
					extraDependentInputLs=[consensusSequenceFile], extraOutputLs=[psmcInputFile], transferOutput=False, \
					extraArguments=None, extraArgumentList=[self.psmcFolderPath, consensusSequenceFile, psmcInputFile], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=360)
		passingData.fq2psmcfaJob = fq2psmcfaJob
		
		psmcOutputFile = File(os.path.join(topOutputDirJob.output, '%s.psmc'%(passingData.bamFnamePrefix)))
		psmcJob = self.addGenericJob(executable=self.psmc, \
					inputArgumentOption=None, \
					outputFile=psmcOutputFile, outputArgumentOption="-o", inputFileList=[psmcInputFile], \
					parentJob=None, parentJobLs=[fq2psmcfaJob, topOutputDirJob], \
					extraDependentInputLs=[], extraOutputLs=[], transferOutput=True, \
					extraArguments='-N%s -t%s -r%s -p %s'%(self.maxNoOfIterations, self.max2N0CoalescentTime, \
									self.initThetaRhoRatio, self.patternOfPSMCTimeIntervals), \
					extraArgumentList=[], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=600)	#10 hours
		
		passingData.wholeGenomePSMCJob = psmcJob
		
		msCommandFile = File(os.path.join(topOutputDirJob.output, '%s.ms_command.sh'%(passingData.bamFnamePrefix)))
		psmcOutput2MSCommandJob = self.addGenericJob(executable=self.GenerateMSCommandFromPSMCOutput, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=[psmcJob, topOutputDirJob], \
					extraDependentInputLs=[psmcOutputFile], extraOutputLs=[msCommandFile], transferOutput=True, \
					extraArguments='', extraArgumentList=[self.psmcFolderPath, psmcOutputFile, msCommandFile], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=None)
		
		plotOutputFnamePrefix = os.path.join(topOutputDirJob.output, '%s'%(passingData.bamFnamePrefix))
		self.addPSMCPlotJob(inputFile=psmcOutputFile, plotOutputFnamePrefix=plotOutputFnamePrefix, \
						absMutationRatePerNucleotide=self.absMutationRatePerNucleotide, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[psmcJob, topOutputDirJob], extraDependentInputLs=[psmcOutputFile], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		
		#2013.2.5 twice the minNoOfGenerations	#estimates close to 10^4 is not accurate
		self.addPSMCPlotJob(inputFile=psmcOutputFile, plotOutputFnamePrefix=plotOutputFnamePrefix, \
						absMutationRatePerNucleotide=self.absMutationRatePerNucleotide, \
						minNoOfGenerations=self.minNoOfGenerations*2, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[psmcJob, topOutputDirJob], extraDependentInputLs=[psmcOutputFile], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		
		return returnData
	
	
	def addPSMCPlotJob(self, inputFile=None, \
					plotOutputFnamePrefix=None, \
					absMutationRatePerNucleotide=None,\
					minNoOfGenerations=None, maxNoOfGenerations=None, maxPopulationSize=None, noOfYearsPerGeneration=None,\
					parentJob=None, parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=1000,  \
					no_of_cpus=None, max_walltime=None, **keywords):
		"""
		2013.2.5
				
		
		perl psmc_plot.pl -u 0.00000001 -x 1000 -X 2000000 -Y 8 -R -T \"Chinese Wolf\" -g 1
			-P \"left top\" AllBootstrapData AllBootstrapData.psmc
		
		Usage:   psmc_plot.pl [options] <out.prefix> <in.psmc>

		Options: -u FLOAT	absolute mutation rate per nucleotide  per generation [2.5e-08]
			-s INT	skip used in data preparation [100]
			-X FLOAT	maximum generations, 0 for auto [0]
			-x FLOAT	minimum generations, 0 for auto [10000]
			-Y FLOAT	maximum popsize, 0 for auto [0]
			-m INT	minimum number of iteration [5]
			-n INT	take n-th iteration (suppress GOF) [20]
			-M titles	multiline mode [null]
			-f STR	font for title, labels and tics [Helvetica,16]
			-g INT	number of years per generation [25]
			-w INT	line width [4]
			-P STR	position of the keys [right top]
			-T STR	figure title [null]
			-N FLOAT	false negative rate [0]
			-S		no scaling
			-L		show the last bin
			-p		convert to PDF (with epstopdf)
			-R		do not remove temporary files
			-G		plot grid
		
		"""
		plotOutputFnamePrefix = "%s_u%s_x%s_X%s_Y%s_g%s"%(plotOutputFnamePrefix, absMutationRatePerNucleotide, minNoOfGenerations, \
									maxNoOfGenerations, \
									maxPopulationSize, noOfYearsPerGeneration)
		plotOutputFile = File('%s.eps'%(plotOutputFnamePrefix))
		plotPDFOutputFile  = File('%s.pdf'%(plotOutputFnamePrefix))
		psmc_plotJob = self.addGenericJob(executable=self.psmc_plot, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=None, \
					parentJobLs=parentJobLs, \
					extraDependentInputLs=[inputFile], extraOutputLs=[plotOutputFile, plotPDFOutputFile], transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=["-u %s -x %s -X %s -Y %s -g %s -p -R"%\
								(absMutationRatePerNucleotide, minNoOfGenerations, maxNoOfGenerations, \
								maxPopulationSize, noOfYearsPerGeneration),\
											plotOutputFnamePrefix, inputFile], \
					job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=no_of_cpus, max_walltime=max_walltime)
		
		plotPNGOutputFile = File('%s.png'%(plotOutputFnamePrefix))
		self.addConvertImageJob(inputFile=plotOutputFile, \
					outputFile=plotPNGOutputFile, density=300, \
					resizeDimension=None, \
					parentJobLs=[psmc_plotJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, job_max_memory=500)
		return psmc_plotJob
	
	def mapReduceOneAlignment(self, workflow=None, alignmentData=None, passingData=None, \
						chrIDSet=None, chr2IntervalDataLs=None, chr2VCFFile=None, \
						outputDirPrefix=None, transferOutput=False, **keywords):
		"""
		2013.1.25
			run splitfa, 
			map:
				bootstrap psmc runs
			reduce: concatenate all bootstrap-psmc output + the normal psmc output
				plot the result
				
		"""
		fq2psmcfaJob = passingData.fq2psmcfaJob
		
		bootstrapOutputDir = "%sBootstrapPSMC_%s"%(outputDirPrefix, passingData.bamFnamePrefix)
		bootstrapOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=bootstrapOutputDir)
		
		splitPSMCInputFile = File(os.path.join(bootstrapOutputDirJob.output, \
									'%s.split%s.psmcfa'%(passingData.bamFnamePrefix, self.bootstrapTrunkSize)))
		splitfaJob = self.addGenericJob(executable=self.splitfa, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=[fq2psmcfaJob, bootstrapOutputDirJob], \
					extraDependentInputLs=[fq2psmcfaJob.output], extraOutputLs=[splitPSMCInputFile], transferOutput=True, \
					extraArgumentList=[self.psmcFolderPath, fq2psmcfaJob.output, splitPSMCInputFile, "%s"%self.bootstrapTrunkSize], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=None)
		
		mergeBootstrapPSMCOutputFile = File(os.path.join(bootstrapOutputDirJob.output, \
									'%s.wholeGenome_%sBootstrap.psmc'%(passingData.bamFnamePrefix, self.noOfBootstraps))) 
		mergeBootstrapPSMCOutputJob = self.addStatMergeJob(statMergeProgram=self.MergeSameHeaderTablesIntoOne, \
				outputF=mergeBootstrapPSMCOutputFile, transferOutput=True, extraArguments='--noHeader', parentJobLs=[bootstrapOutputDirJob])
		#first add the whole-genome psmc output
		wholeGenomePSMCJob = passingData.wholeGenomePSMCJob 
		self.addInputToStatMergeJob(workflow, statMergeJob=mergeBootstrapPSMCOutputJob, \
								inputF=wholeGenomePSMCJob.output, parentJobLs=[wholeGenomePSMCJob])
		
		for i in xrange(self.noOfBootstraps):
			psmcOutputFile = File(os.path.join(bootstrapOutputDirJob.output, '%s.%s.psmc'%(passingData.bamFnamePrefix,i,)))
			psmcJob = self.addGenericJob(executable=self.psmc, \
					inputArgumentOption=None, \
					outputFile=psmcOutputFile, outputArgumentOption="-o", inputFileList=[splitPSMCInputFile], \
					parentJob=None, parentJobLs=[splitfaJob, bootstrapOutputDirJob], \
					extraDependentInputLs=[], extraOutputLs=[], transferOutput=False, \
					extraArguments='-N%s -t%s -r%s -p %s -b'%(self.maxNoOfIterations, self.max2N0CoalescentTime, \
									self.initThetaRhoRatio, self.patternOfPSMCTimeIntervals), \
					extraArgumentList=[], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=600)
			self.addInputToStatMergeJob(workflow, statMergeJob=mergeBootstrapPSMCOutputJob, \
								inputF=psmcJob.output, parentJobLs=[psmcJob])
		
		
		plotOutputFnamePrefix = os.path.join(bootstrapOutputDirJob.output, '%s.wholeGenome_%sBootstrap'%\
											(passingData.bamFnamePrefix, self.noOfBootstraps))
		
		psmc_plotJob= self.addPSMCPlotJob(inputFile=mergeBootstrapPSMCOutputJob.output, plotOutputFnamePrefix=plotOutputFnamePrefix, \
						absMutationRatePerNucleotide=self.absMutationRatePerNucleotide, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[mergeBootstrapPSMCOutputJob, bootstrapOutputDirJob], extraDependentInputLs=[psmcOutputFile], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		#2013.2.5 twice the minNoOfGenerations	#estimates close to 10^4 is not accurate
		psmc_plotJob= self.addPSMCPlotJob(inputFile=mergeBootstrapPSMCOutputJob.output, plotOutputFnamePrefix=plotOutputFnamePrefix, \
						absMutationRatePerNucleotide=self.absMutationRatePerNucleotide, \
						minNoOfGenerations=self.minNoOfGenerations*2, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[mergeBootstrapPSMCOutputJob, bootstrapOutputDirJob], extraDependentInputLs=[psmcOutputFile], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		return
	
if __name__ == '__main__':
	main_class = PSMCOnAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()