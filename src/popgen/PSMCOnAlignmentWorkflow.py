#!/usr/bin/env python
"""
Examples:
	#2013.2.4
	%s --ind_aln_id_ls 929,1508,1510,926,1577,1507,1534,1524,1512
		-u yh -a 524  -z localhost -o  dags/PSMCOnAlignment/PSMCOnAlignmentOnSubspecies.xml
		-j hcondor -l hcondor
		--clusters_size 25 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 --needSSHDBTunnel --commit
		# -Y
		#--minContigID 96 --maxContigID 100	#useless right now
	
	#2013.2.12 more refined pattern: 4*2+40*1+5*2
	%s --ind_aln_id_ls 929,1508,1510,926,1577,1507,1534 --patternOfPSMCTimeIntervals 4*2+40*1+5*2
		-u yh -a 524  -z localhost -o  dags/PSMCOnAlignment/PSMCOn7AlignmentOnSubspeciesRefinedPattern.xml
		-j hcondor -l hcondor
		--clusters_size 25 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 --needSSHDBTunnel --commit
	
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
					('minBaseQ', 1, int):[20, '', 1, 'inferred consensus base quality '],\
					('minMapQ', 1, int):[30, '', 1, 'read alignment mapping quality'],\
					('minRMSMapQ', 1, int):[10, '', 1, 'root mean squared mapping quality of reads covering the locus'],\
					('minDistanceToIndel', 1, int):[5, '', 1, 'min distance to predicted short insertions or deletions'],\
					
					('patternOfPSMCTimeIntervals', 1, ):["4+25*1+4+6", '', 1, 'pattern of parameters, -p of psmc.\n\
	The first (most recent) parameter spans the first 4 atomic time intervals;\n\
	each of the next 25 parameters spans 1 intervals;\n\
	the 27th spans 4 intervals and the last parameter spans the last 6 time intervals.'],\
					('maxNoOfIterations', 1, int):[25, '', 1, 'number of iterations for the EM algorithm, -N of psmc'],\
					('initThetaRhoRatio', 1, float):[4, '', 1, 'initial theta/rho ratio, -r of psmc'],\
					('max2N0CoalescentTime', 1, float):[15, '', 1, 'maximum 2N0 coalescent time, -t of psmc '],\
					
					('mutationRatePerNucleotidePerGeneration', 1, float):[5e-09, '', 1, 'absolute mutation rate per nucleotide per generation, -u FLOAT of psmc_plot'],\
					('maxNoOfGenerations', 1, int):[0, '', 1, 'maximum generations, 0 for auto , -X of psmc_plot'],\
					('minNoOfGenerations', 1, int):[10000, '', 1, 'minimum number of generations, 0 for auto , -x of psmc_plot'],\
					('noOfYearsPerGeneration', 1, int):[5, '', 1, 'number of years per generation, -g of psmc_plot.pl'],\
					('maxPopulationSize', 1, int):[30, '', 1, 'maximum popsize in terms of 10^4, 0 for auto, -Y of psmc_plot.pl'],\
					
					('chromosomeLengthToSimulate', 1, int):[20000000, '', 1, 'passed to history2ms.pl. ms/msHOT would segment-fault if it is too long, i.e. 300Mb, or even 30Mb occasionally'],\
					('recombinationHotSpotRateMultipler', 1, float):[10, '', 1, 'passed to history2ms.pl. recombination rate in hotspots are FLOAT times larger. used in ms simulation. '],\
					('divergenceTimeToSimulate', 1, int):[0, '', 1, "no idea of what it is, passed to history2ms.pl"],\
					
					('bootstrapTrunkSize', 1, int):[50000, '', 1, 'number of 100bp fragments for each bootstrap psmc input. option for splitfa'],\
					('noOfBootstraps', 1, int):[100, '', 1, 'number of bootstraps for each whole-alignment psmc job'],\
					('noOfSimulations', 1, int):[50, '', 1, 'number of simulations via ms for each whole-alignment psmc job'],\
					('msPath', 0, ): ["%s/script/lh3_foreign/msHOT-lite/msHOT-lite", '', 1, 'path to the ms or msHOT, msHOT-lite program, '],\
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
	
		if hasattr(self, 'msPath'):
			self.msPath =  self.insertHomePath(self.msPath, self.home_path)
		else:
			self.msPath = None
	
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
		

		psmc_plot = Executable(namespace=namespace, name="psmc_plot", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		psmc_plot.addPFN(PFN("file://" + os.path.join(self.psmcFolderPath, "utils/psmc_plot.pl"), site_handler))
		executableClusterSizeMultiplierList.append((psmc_plot, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'popgen/GenerateMSCommandFromPSMCOutput.sh'), \
										name='GenerateMSCommandFromPSMCOutput', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'popgen/ReplaceMSPathInMSCommandFile.py'), \
										name='ReplaceMSPathInMSCommandFile', clusterSizeMultipler=2)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/AddIndividualAlignmentConsensusSequence2DB.py'), \
										name='AddIndividualAlignmentConsensusSequence2DB', clusterSizeMultipler=0.05)
		#2013.2.10 reduce clustering size for runShellCommand because ms takes a while
		self.setOrChangeExecutableClusterSize(executable=self.runShellCommand, clusterSizeMultipler=0.1)
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2013.2.10
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		return returnData
	
	def mapEachAlignment(self, workflow=None, alignmentData=None, passingData=None, transferOutput=True, **keywords):
		"""
		2013.1.25
			run fq2psmcfa, psmc, ms-simulation-commandline, plot
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
		"""
		returnData = PassingData(jobDataLs=[])
		
		topOutputDirJob = passingData.topOutputDirJob
		plotOutputDirJob = passingData.plotOutputDirJob
		
		refFastaFile = passingData.refFastaFList[0]
		
		alignmentData = passingData.alignmentData
		
		alignment = alignmentData.alignment
		alignmentParentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		
		minDP = int(max(1, alignment.median_depth/2))
		maxDP = int(alignment.median_depth*2)
		individual_alignment_consensus_sequence = self.db_vervet.checkIndividualAlignmentConsensusSequence(individual_alignment_id=alignment.id, \
									minDP=minDP, \
									maxDP=maxDP, minBaseQ=self.minBaseQ, minMapQ=self.minMapQ,\
									minRMSMapQ=self.minRMSMapQ, minDistanceToIndel=self.minDistanceToIndel)
		if individual_alignment_consensus_sequence:
			consensusSequenceFname = individual_alignment_consensus_sequence.getFileAbsPath(oldDataDir=self.db_vervet.data_dir, \
												newDataDir=self.data_dir)
			consensusSequenceFile = self.registerOneInputFile(inputFname=consensusSequenceFname, \
										folderName=self.pegasusFolderName)
			extractConsensusSequenceFromAlignmentJob = PassingData(output=consensusSequenceFile)
		else:
			consensusSequenceFile = File(os.path.join(topOutputDirJob.output, '%s.fq.gz'%(passingData.bamFnamePrefix)))
			extractConsensusSequenceFromAlignmentJob = self.addGenericJob(executable=self.ExtractConsensusSequenceFromAlignment, \
						inputArgumentOption=None, \
						outputArgumentOption=None, inputFileList=[], \
						parentJob=None, parentJobLs=alignmentParentJobLs+[topOutputDirJob], \
						extraDependentInputLs=[refFastaFile, bamF, baiF], extraOutputLs=[consensusSequenceFile], transferOutput=True, \
						extraArguments=None, extraArgumentList=[refFastaFile, bamF, consensusSequenceFile, \
										'%s %s %s %s %s %s'%(minDP, maxDP, self.minBaseQ, self.minMapQ, self.minRMSMapQ, self.minDistanceToIndel)], \
						job_max_memory=1000,  sshDBTunnel=None, \
						key2ObjectForJob=None, no_of_cpus=None, max_walltime=3000)	#max_walltime is in minutes, 50 hours
			
			logFile = File(os.path.join(topOutputDirJob.output, '%s_alignment%s_2DB.log'%(passingData.bamFnamePrefix, alignment.id)))
			self.addGenericFile2DBJob(executable=self.AddIndividualAlignmentConsensusSequence2DB, \
									inputFile=extractConsensusSequenceFromAlignmentJob.output, \
									inputArgumentOption="-i", \
						data_dir=self.data_dir, logFile=logFile, commit=self.commit,\
						parentJobLs=[extractConsensusSequenceFromAlignmentJob], extraDependentInputLs=None, \
						extraOutputLs=None, transferOutput=True, \
						extraArgumentList=["--individual_alignment_id %s"%(alignment.id), " --format fastq", "--minDP %s"%(minDP), \
										"--maxDP=%s"%(maxDP), "--minBaseQ=%s"%(self.minBaseQ), \
										"--minMapQ %s"%(self.minMapQ), "--minRMSMapQ %s"%(self.minRMSMapQ), \
										"--minDistanceToIndel %s"%(self.minDistanceToIndel)], \
						job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, \
						key2ObjectForJob=None)
		
		#added self.minBaseQ to this fq2psmcfa job 
		psmcInputFile = File(os.path.join(topOutputDirJob.output, '%s.psmcfa'%(passingData.bamFnamePrefix)))
		fq2psmcfaJob = self.addGenericJob(executable=self.fq2psmcfa, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=[extractConsensusSequenceFromAlignmentJob, topOutputDirJob], \
					extraDependentInputLs=[extractConsensusSequenceFromAlignmentJob.output], extraOutputLs=[psmcInputFile], transferOutput=False, \
					extraArguments=None, extraArgumentList=[self.psmcFolderPath, \
											extractConsensusSequenceFromAlignmentJob.output, psmcInputFile, "%s"%self.minBaseQ], \
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
		
		psmcJob.alignmentData = alignmentData	#to be accessed in reduce()
		
		passingData.wholeGenomePSMCJob = psmcJob
		returnData.wholeGenomePSMCJob = psmcJob	#to be accessed in reduce()
		
		msCommandFile = File(os.path.join(topOutputDirJob.output, '%s.ms_command.sh'%(passingData.bamFnamePrefix)))
		psmcOutput2MSCommandJob = self.addGenericJob(executable=self.GenerateMSCommandFromPSMCOutput, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=[], \
					parentJob=None, parentJobLs=[psmcJob, topOutputDirJob], \
					extraDependentInputLs=[psmcJob.output], extraOutputLs=[msCommandFile], transferOutput=True, \
					extraArguments='', extraArgumentList=[self.psmcFolderPath, psmcJob.output, msCommandFile, \
						"-L %s -u %s -g %s -R %s -d %s"%\
						(self.chromosomeLengthToSimulate, self.mutationRatePerNucleotidePerGeneration, self.noOfYearsPerGeneration, self.recombinationHotSpotRateMultipler, self.divergenceTimeToSimulate)], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=None)
		
		self.simulateOneAlignmentPSMCResult(alignmentData=alignmentData, passingData=passingData, \
								psmcOutput2MSCommandJob=psmcOutput2MSCommandJob, outputDirPrefix=passingData.outputDirPrefix, \
								transferOutput=transferOutput)
		
		plotOutputFnamePrefix = os.path.join(plotOutputDirJob.output, '%s'%(passingData.bamFnamePrefix))
		self.addPSMCPlotJob(inputFileList=[psmcJob.output], \
						figureLegendList=[self.getPSMCPlotLabelForOneAlignment(alignment)], \
						figureTitle = 'depth%s-%s'%(minDP, maxDP),\
						plotOutputFnamePrefix=plotOutputFnamePrefix, \
						mutationRatePerNucleotidePerGeneration=self.mutationRatePerNucleotidePerGeneration, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[psmcJob, plotOutputDirJob], extraDependentInputLs=[psmcJob.output], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		
		return returnData
	
	def getPSMCPlotLabelForOneAlignment(self, individual_alignment=None):
		"""
		2013.2.12
		"""
		return '%s_%s%dMD'%(individual_alignment.individual_sequence.individual.code, individual_alignment.individual_sequence.individual.sex, individual_alignment.median_depth)
	
	def addPSMCPlotJob(self, inputFileList=None, figureLegendList=None, figureTitle=None, \
					plotOutputFnamePrefix=None, \
					mutationRatePerNucleotidePerGeneration=None,\
					minNoOfGenerations=None, maxNoOfGenerations=None, maxPopulationSize=None, noOfYearsPerGeneration=None,\
					parentJob=None, parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=1000,  \
					no_of_cpus=None, max_walltime=None, **keywords):
		"""
		2013.2.10, add/change argument inputFileList, figureLegendList
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
		plotOutputFnamePrefix = "%s_u%s_x%s_X%s_Y%s_g%s"%(plotOutputFnamePrefix, mutationRatePerNucleotidePerGeneration, minNoOfGenerations, \
									maxNoOfGenerations, \
									maxPopulationSize, noOfYearsPerGeneration)
		plotOutputFile = File('%s.eps'%(plotOutputFnamePrefix))
		#plotPDFOutputFile  = File('%s.pdf'%(plotOutputFnamePrefix))
		if figureLegendList:
			figureLegendArgument = "-M %s"%(",".join(figureLegendList))
		else:
			figureLegendArgument = ""
		
		if figureTitle:
			figureTitleArgument = "-T %s"%(figureTitle)
		else:
			figureTitleArgument = ""
		psmc_plotJob = self.addGenericJob(executable=self.psmc_plot, \
					inputArgumentOption=None, \
					outputArgumentOption=None, inputFileList=inputFileList, \
					parentJobLs=parentJobLs, \
					extraDependentInputLs=None, extraOutputLs=[plotOutputFile], transferOutput=transferOutput, \
					extraArguments=extraArguments, frontArgumentList=["-u %s -x %s -X %s -Y %s -g %s %s %s -R"%\
								(mutationRatePerNucleotidePerGeneration, minNoOfGenerations, maxNoOfGenerations, \
									maxPopulationSize, noOfYearsPerGeneration, figureLegendArgument, figureTitleArgument),\
							plotOutputFnamePrefix], \
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
		alignmentData = passingData.alignmentData
		plotOutputDirJob = passingData.plotOutputDirJob
		
		alignment = alignmentData.alignment
		bootstrapOutputDir = "%sBootstrapPSMC_%s"%(outputDirPrefix, passingData.bamFnamePrefix)
		bootstrapOutputDirJob = self.addMkDirJob(outputDir=bootstrapOutputDir)
		
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
				outputF=mergeBootstrapPSMCOutputFile, transferOutput=True, extraArguments='--noHeader --exitNonZeroIfAnyInputFileInexistent', parentJobLs=[bootstrapOutputDirJob])
		#first add the whole-genome psmc output
		wholeGenomePSMCJob = passingData.wholeGenomePSMCJob 
		self.addInputToStatMergeJob(workflow, statMergeJob=mergeBootstrapPSMCOutputJob, \
								inputF=wholeGenomePSMCJob.output, parentJobLs=[wholeGenomePSMCJob])
		
		for i in xrange(self.noOfBootstraps):
			psmcOutputFile = File(os.path.join(bootstrapOutputDirJob.output, '%s.%s.psmc'%(passingData.bamFnamePrefix,i,)))
			psmcJob = self.addGenericJob(executable=self.psmc, \
					inputArgumentOption=None, \
					outputFile=psmcOutputFile, outputArgumentOption="-o", inputFileList=[splitfaJob.output], \
					parentJob=None, parentJobLs=[splitfaJob, bootstrapOutputDirJob], \
					extraDependentInputLs=[], extraOutputLs=[], transferOutput=False, \
					extraArguments='-N%s -t%s -r%s -p %s -b'%(self.maxNoOfIterations, self.max2N0CoalescentTime, \
									self.initThetaRhoRatio, self.patternOfPSMCTimeIntervals), \
					extraArgumentList=[], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=600)	#-b is for bootstrapped input
			self.addInputToStatMergeJob(workflow, statMergeJob=mergeBootstrapPSMCOutputJob, \
								inputF=psmcJob.output, parentJobLs=[psmcJob])
		
		
		plotOutputFnamePrefix = os.path.join(plotOutputDirJob.output, '%s.wholeGenome_%sBootstrap'%\
											(passingData.bamFnamePrefix, self.noOfBootstraps))
		psmc_plotJob= self.addPSMCPlotJob(inputFileList=[mergeBootstrapPSMCOutputJob.output], \
						figureLegendList=[self.getPSMCPlotLabelForOneAlignment(alignment)],\
						figureTitle = '%sBootstraps'%(self.noOfBootstraps),\
						plotOutputFnamePrefix=plotOutputFnamePrefix, \
						mutationRatePerNucleotidePerGeneration=self.mutationRatePerNucleotidePerGeneration, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[mergeBootstrapPSMCOutputJob, plotOutputDirJob], \
						extraDependentInputLs=[], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		return
	
	def simulateOneAlignmentPSMCResult(self, alignmentData=None, passingData=None, psmcOutput2MSCommandJob=None,\
						outputDirPrefix=None, transferOutput=False, **keywords):
		"""
		2013.2.10 
			take the ms commandlines from psmcOutput2MSCommandJob
			for i in xrange(noOfSimulations):
				run noOfSimulations ms (msHOT) jobs
				turn ms output into psmc input
				run psmc	(output in bootstrap-like output)
			run psmc_plot on all of them (assuming bootstrap mode)
		"""
		fq2psmcfaJob = passingData.fq2psmcfaJob
		alignmentData = passingData.alignmentData
		plotOutputDirJob = passingData.plotOutputDirJob
		
		alignment = alignmentData.alignment
		simulateOutputDir = "%sSimulatePSMC_%s"%(outputDirPrefix, passingData.bamFnamePrefix)
		simulateOutputDirJob = self.addMkDirJob(outputDir=simulateOutputDir)
		
		mergeSimulatePSMCOutputFile = File(os.path.join(simulateOutputDirJob.output, \
									'%s.wholeGenome_%sSimulations.psmc'%(passingData.bamFnamePrefix, self.noOfSimulations))) 
		mergeSimulatePSMCOutputJob = self.addStatMergeJob(statMergeProgram=self.MergeSameHeaderTablesIntoOne, \
				outputF=mergeSimulatePSMCOutputFile, transferOutput=True, extraArguments='--noHeader --exitNonZeroIfAnyInputFileInexistent', parentJobLs=[simulateOutputDirJob])
		#first add the whole-genome psmc output
		wholeGenomePSMCJob = passingData.wholeGenomePSMCJob 
		self.addInputToStatMergeJob(self, statMergeJob=mergeSimulatePSMCOutputJob, \
								inputF=wholeGenomePSMCJob.output, parentJobLs=[wholeGenomePSMCJob])
		
		newMSCommandFile = File(os.path.join(simulateOutputDirJob.output, \
									'%s.new_ms_command.sh'%(passingData.bamFnamePrefix)))
		replaceMSPathJob = self.addAbstractMapperLikeJob(executable=self.ReplaceMSPathInMSCommandFile, \
					inputF=psmcOutput2MSCommandJob.output, outputF=newMSCommandFile, extraOutputLs=None,\
					parentJobLs=[psmcOutput2MSCommandJob, simulateOutputDirJob], transferOutput=True, job_max_memory=2000,\
					extraArguments=None, extraArgumentList=["--oldMSPath msHOT-lite", "--msPath %s"%(self.msPath),], \
					extraDependentInputLs=None, \
					sshDBTunnel=None, **keywords)
		for i in xrange(self.noOfSimulations):

			msOutputFile = File(os.path.join(simulateOutputDirJob.output, \
									'%s.sim%s_msOutput.txt.gz'%(passingData.bamFnamePrefix, i)))
			msJob = self.addGenericJob(executable=self.runShellCommand, inputFile=replaceMSPathJob.output, inputArgumentOption=None, \
					outputFile=msOutputFile, outputArgumentOption=None, \
					parentJob=None, parentJobLs=[replaceMSPathJob], extraDependentInputLs=None, transferOutput=False, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=["1"], job_max_memory=2000,  \
					no_of_cpus=None, max_walltime=None)	#the extra "1" argument enables gzipping
			
			#output a fasta file, every 100bp -> N/T/K
			msOutputFastQFile = File(os.path.join(simulateOutputDirJob.output, \
									'%s.sim%s_msOutput.fq.gz'%(passingData.bamFnamePrefix, i)))
			msOutput2FastQJob = self.addGenericJob(executable=self.ConvertMSOutput2FASTQ, inputFile=msJob.output, \
					outputFile=msOutputFastQFile,\
					parentJob=None, parentJobLs=[msJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
					frontArgumentList=None, extraArguments="--inputFileFormat 2", extraArgumentList=None, job_max_memory=2000,  \
					no_of_cpus=None, max_walltime=None)
			
			#added self.minBaseQ to this fq2psmcfa job 
			psmcInputFile = File(os.path.join(simulateOutputDirJob.output, '%s.sim%s.psmcfa'%(passingData.bamFnamePrefix,i)))
			fq2psmcfaJob = self.addGenericJob(executable=self.fq2psmcfa, \
						inputArgumentOption=None, \
						outputArgumentOption=None, inputFileList=[], \
						parentJob=None, parentJobLs=[msOutput2FastQJob, simulateOutputDirJob], \
						extraDependentInputLs=[msOutput2FastQJob.output], extraOutputLs=[psmcInputFile], transferOutput=False, \
						extraArguments=None, extraArgumentList=[self.psmcFolderPath, \
												msOutput2FastQJob.output, psmcInputFile, "%s"%self.minBaseQ], \
						job_max_memory=1000,  sshDBTunnel=None, \
						key2ObjectForJob=None, no_of_cpus=None, max_walltime=360)
			
			splitPSMCInputFile = File(os.path.join(simulateOutputDirJob.output, \
									'%s.sim%s.bootstrap.psmcfa'%(passingData.bamFnamePrefix, i)))
			splitfaJob = self.addGenericJob(executable=self.splitfa, \
						inputArgumentOption=None, \
						outputArgumentOption=None, inputFileList=[], \
						parentJob=None, parentJobLs=[fq2psmcfaJob, simulateOutputDirJob], \
						extraDependentInputLs=[fq2psmcfaJob.output], extraOutputLs=[splitPSMCInputFile], transferOutput=False, \
						extraArgumentList=[self.psmcFolderPath, fq2psmcfaJob.output, splitPSMCInputFile, "%s"%(self.bootstrapTrunkSize*3)], \
						job_max_memory=1000,  sshDBTunnel=None, \
						key2ObjectForJob=None, no_of_cpus=None, max_walltime=None)
			
			psmcOutputFile = File(os.path.join(simulateOutputDirJob.output, '%s.sim%s.psmc'%(passingData.bamFnamePrefix,i,)))
			psmcJob = self.addGenericJob(executable=self.psmc, \
					inputArgumentOption=None, \
					outputFile=psmcOutputFile, outputArgumentOption="-o", inputFileList=[splitfaJob.output], \
					parentJob=None, parentJobLs=[splitfaJob, simulateOutputDirJob], \
					extraDependentInputLs=[], extraOutputLs=[], transferOutput=False, \
					extraArguments='-N%s -t%s -r%s -p %s -b'%(self.maxNoOfIterations, self.max2N0CoalescentTime, \
									self.initThetaRhoRatio, self.patternOfPSMCTimeIntervals), \
					extraArgumentList=[], \
					job_max_memory=1000,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=None, max_walltime=600)
			self.addInputToStatMergeJob(self, statMergeJob=mergeSimulatePSMCOutputJob, \
								inputF=psmcJob.output, parentJobLs=[psmcJob])
		
		
		plotOutputFnamePrefix = os.path.join(plotOutputDirJob.output, '%s.%sMSSimulations'%\
											(passingData.bamFnamePrefix, self.noOfSimulations))
		psmc_plotJob= self.addPSMCPlotJob(inputFileList=[mergeSimulatePSMCOutputJob.output], \
						figureLegendList=[self.getPSMCPlotLabelForOneAlignment(alignment)],\
						figureTitle = '%sSimulations'%(self.noOfSimulations),\
						plotOutputFnamePrefix=plotOutputFnamePrefix, \
						mutationRatePerNucleotidePerGeneration=self.mutationRatePerNucleotidePerGeneration, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=[mergeSimulatePSMCOutputJob, plotOutputDirJob], \
						extraDependentInputLs=[], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		return
	
	def reduce(self, workflow=None, passingData=None, reduceAfterEachAlignmentDataLs=None,
			transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData()
		returnData.jobDataLs = []
		
		plotOutputDir = "%sMultiLinePSMCPlot"%(passingData.outputDirPrefix)
		plotOutputDirJob = self.addMkDirJob(outputDir=plotOutputDir)
		
		
		psmc_plot_inputFileList = []
		psmc_plot_parentJobList = []
		figureLegendList = []
		for mapEachAlignmentData in passingData.mapEachAlignmentDataLs:
			wholeGenomePSMCJob = mapEachAlignmentData.wholeGenomePSMCJob
			psmc_plot_inputFileList.append(wholeGenomePSMCJob.output)
			psmc_plot_parentJobList.append(wholeGenomePSMCJob)
			alignment = wholeGenomePSMCJob.alignmentData.alignment
			legend = '%s_M%s'%(alignment.individual_sequence.individual.code, alignment.median_depth)
			figureLegendList.append(legend)
			
		plotOutputFnamePrefix = os.path.join(plotOutputDirJob.output, 'MultiLine_%sAlignments'%\
											(len(psmc_plot_inputFileList)))
		#no space within arguments
		psmc_plotJob= self.addPSMCPlotJob(inputFileList=psmc_plot_inputFileList, \
						figureLegendList=figureLegendList,\
						figureTitle = '%sAlignments'%(len(psmc_plot_inputFileList)),\
						plotOutputFnamePrefix=plotOutputFnamePrefix, \
						mutationRatePerNucleotidePerGeneration=self.mutationRatePerNucleotidePerGeneration, \
						minNoOfGenerations=self.minNoOfGenerations, maxNoOfGenerations=self.maxNoOfGenerations, \
						maxPopulationSize=self.maxPopulationSize, noOfYearsPerGeneration=self.noOfYearsPerGeneration, \
						parentJobLs=psmc_plot_parentJobList + [plotOutputDirJob], \
						extraDependentInputLs=[], \
						transferOutput=True, job_max_memory=2000, no_of_cpus=None, max_walltime=None)
		return returnData
	
if __name__ == '__main__':
	main_class = PSMCOnAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()