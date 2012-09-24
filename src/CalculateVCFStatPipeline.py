#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s  -a 524 -I ./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/gatk/
		-o workflow/VCFStat/VCFStat_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_gatk.xml
		-l condorpool -j condorpool  -u yh -z uclaOffice
	
	#different window size for pi and TiTv
	%s  -a 524 -I ./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/gatk/
		-o workflow/VCFStat/VCFStat_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_gatk_window20Mb.xml
		-l condorpool -j condorpool  -u yh -z uclaOffice
		-w 20000000
	
	#2012.5.11 run on hoffman2 condor, need ssh tunnel for the PlotVCFtoolsStat (-H)
	%s ...
		-j hcondor -l hcondor -u yh -z localhost  -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-H 
	
	# 2012.8.1 on hoffman2 condor, LD window size=800Kb (-L), minChrLengthForPlot=4million (-P), minChrSiz for computing =1mil (-n)
	# need ssh tunnel for the PlotVCFtoolsStat (-H)
	# try "-n 0" to get rid of Ti/Tv, pi, snpDensity jobs
	# try "-L 0" to get rid of LD jobs (which take long time)
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	# try "-s 0.01" to increase data sampling rate
	%s -a 524 -I FilterGenotypeMethod5_ByMethod7Sites_NoDepthFilter_MaxSNPMissing0.5.2012.7.30T1806/method_5_vcftoolsFilter/
		-o workflow/VCFStat/VCFStat_Method8_L800000P4000000n1000000.xml
		-j hcondor -l hcondor -u yh -z localhost  -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-L 800000 -P 4000000 -n 1000000 -H 
		#-V 90 -x 100 -C 1 -s 0.01
		
	

Description:
	2011-10-6
		a program which generates a pegasus workflow dag (xml file) that for each vcf file,
		generates these statistics and plots and more 
		
		pi (by two ways) => theta, snp density, AAF over position, LD (r^2) over distance, TsTv, heterozygosity, more
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow

class CalculateVCFStatPipeline(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('windowSize', 1, int): [200000, 'w', 1, "window size for TiTv, pi, snpDensity calculation by vcftools.\
			set it to 0 if you do not need these stats.", ],\
						('LDWindowSize', 0, int): [0, 'L', 1, 'window size for LD calculation by vcftools, set it to 0 to skip LD jobs.', ],\
						('minChrLengthForPlot', 1, int): [2000000, 'P', 1, 'minimum chromosome size for a chromosome to be included in plot', ],\
						('minChrSize', 1, int): [500000, 'n', 1, 'minimum chromosome size for any computing ', ],\
						('vcf1Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('includeIndelVCF', 1, int): [0, '', 1, 'toggle this to include indel VCF, filename with "indel"'],\
						('samplingRate', 1, float): [0.001, 's', 1, "sampling rate for plotting the locus-specific stats.\
			The LD plot's sampling rate is samplingRate*samplingRate*10 "],\
						})

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		
		self.vcf1Dir = os.path.abspath(self.vcf1Dir)
	
	def addChrLengthAppendingJob(self, workflow=None, AddChromosomeLengthToTSVFile=None, inputF=None, outputF=None, \
							divideByLength=True, chrLength=1000, divideStartingColumn=2, \
							chromosome=None, addChrName=False, \
							parentJobLs=None, extraDependentInputLs=None, \
							extraArguments=None, transferOutput=False, job_max_memory=2000, **keywords):
		"""
		2012.7.31
			rewrite it using addGenericJob()
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = ["-l %s"%chrLength, "-s %s"%divideStartingColumn]
		extraOutputLs = []
		key2ObjectForJob = {}
		if chromosome:
			extraArgumentList.append("-c %s"%(chromosome))
		if addChrName:
			extraArgumentList.append("-a")
		if divideByLength:
			extraArgumentList.append("-v")
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=AddChromosomeLengthToTSVFile, inputFile=inputF, outputFile=outputF, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job

	
	
	def addTallyAACFromVCFJob(self, workflow, TallyAACFromVCF=None, genomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputF=None, parentJobLs=None, \
							job_max_memory=1000, extraDependentInputLs=[],\
							input_site_handler=None, transferOutput=False, **keywords):
		"""
		2012.7.25 why *.vcf.tbi file is needed? get rid of it. 
		"""
		# Add a mkdir job for any directory.
		TallyAACFromVCFJob = Job(namespace=workflow.namespace, name=TallyAACFromVCF.name, version=workflow.version)
		refFastaF = refFastaFList[0]
		TallyAACFromVCFJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", genomeAnalysisTKJar, "-R", refFastaF, 
						"-T TallyAACFromVCF", "--variant:VCF", inputVCFF, "-o", outputF)
		for refFastaFile in refFastaFList:
			TallyAACFromVCFJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		TallyAACFromVCFJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)
		"""	# 2012.7.25
		vcfTabixF = File('%s.tbi'%(inputVCFF.name))	#relative path
		vcfTabixF.absPath = '%s.tbi'%inputVCFF.absPath
		vcfTabixF.addPFN(PFN("file://" + vcfTabixF.absPath, input_site_handler))
		workflow.addFile(vcfTabixF)
		TallyAACFromVCFJob.uses(vcfTabixF, transfer=True, register=True, link=Link.INPUT)
		"""
		for input in extraDependentInputLs:
			TallyAACFromVCFJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		TallyAACFromVCFJob.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		TallyAACFromVCFJob.output = outputF
		workflow.addJob(TallyAACFromVCFJob)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=parentJob, child=TallyAACFromVCFJob)
		yh_pegasus.setJobProperRequirement(TallyAACFromVCFJob, job_max_memory=job_max_memory)
		return TallyAACFromVCFJob
	
	def addCountHomoHetInOneVCFJob(self, workflow=None, executable=None, inputF=None, outputF=None, \
							chrLength=1000, chromosome=None, \
							parentJobLs=None, extraDependentInputLs=None, \
							extraArguments=None, transferOutput=False, job_max_memory=2000, **keywords):
		"""
		2012.8.1 encapsulate the code from run() to here
		
			countHomoHetInOneVCFJob = Job(namespace=workflow.namespace, name=workflow.CountHomoHetInOneVCF.name, version=workflow.version)
			outputF = File(os.path.join(statOutputDir, "%s.homoHetCountPerSample.tsv"%(os.path.splitext(inputFname)[0])))
			countHomoHetInOneVCFJob.addArguments("-i", vcf1, "-c", chr, "-l %s"%(chr_size), '-o', outputF)
			countHomoHetInOneVCFJob.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
			countHomoHetInOneVCFJob.uses(outputF, transfer=False, register=True, link=Link.OUTPUT)
			countHomoHetInOneVCFJob.output = outputF
			yh_pegasus.setJobProperRequirement(countHomoHetInOneVCFJob, job_max_memory=1000)
			workflow.addJob(countHomoHetInOneVCFJob)
			workflow.depends(parent=statOutputDirJob, child=countHomoHetInOneVCFJob)
			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = ["-l %s"%chrLength]
		extraOutputLs = []
		key2ObjectForJob = {}
		if chromosome:
			extraArgumentList.append("-c %s"%(chromosome))
		if extraArguments:
			extraArgumentList.append(extraArguments)
		job= self.addGenericJob(executable=executable, inputFile=inputF, outputFile=outputF, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if not workflow:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#executableList = []
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		CountHomoHetInOneVCF = Executable(namespace=namespace, name="CountHomoHetInOneVCF", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CountHomoHetInOneVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CountHomoHetInOneVCF.py"), site_handler))
		executableClusterSizeMultiplierList.append((CountHomoHetInOneVCF, 1))
		
		AddChromosomeLengthToTSVFile = Executable(namespace=namespace, name="AddChromosomeLengthToTSVFile", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		AddChromosomeLengthToTSVFile.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/AddChromosomeLengthToTSVFile.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddChromosomeLengthToTSVFile, 1))
		
		AddHetFractionToVCFtoolsHWE = Executable(namespace=namespace, name="AddHetFractionToVCFtoolsHWE", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		AddHetFractionToVCFtoolsHWE.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/AddHetFractionToVCFtoolsHWE.py"), site_handler))
		executableClusterSizeMultiplierList.append((AddHetFractionToVCFtoolsHWE, 1))
		
		TallyAACFromVCF = Executable(namespace=namespace, name="TallyAACFromVCF", version=version, os=operatingSystem, arch=architecture, installed=True)
		TallyAACFromVCF.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((TallyAACFromVCF, 0.2))
		
		#2012.8.29 vcftoolsWrapper just for LD caculation 
		import copy
		vcftoolsWrapperLD = copy.deepcopy(self.vcftoolsWrapper)
		vcftoolsWrapperLD.name = 'vcftoolsWrapperLD'
		vcftoolsWrapperLD.clearProfiles()	#get rid of existing profiles
		executableClusterSizeMultiplierList.append((vcftoolsWrapperLD, 0))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	
	def addStatCalculationJobs(self, workflow=None, inputData=None, refFastaFList=None, chr2size=None,\
						windowSize=None, minChrLengthForPlot=None, minChrSize=None, LDWindowSize=None, \
						transferOutput=True, outputDirPrefix="", samplingRate=0.001, minSiteGap=10000):
		"""
		2012.8.13
			if windowSize=0, no jobs for Ti/Tv, snpDensity, pi
		2012.8.2
			these arguments won't be used: minMAC=None, minMAF=None, maxSNPMissingRate=None, 
		"""
		if workflow is None:
			workflow = self
		
		no_of_jobs = 0
		statOutputDir = "%sstat"%(outputDirPrefix)
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=statOutputDir, \
												namespace=workflow.namespace, version=workflow.version)
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir, \
												namespace=workflow.namespace, version=workflow.version)
		no_of_jobs += 2
		input_site_handler = self.input_site_handler
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		homoHetCountFinalOutputF = File(os.path.join(statOutputDir, 'homoHetCountPerSamplePerContig.tsv'))
		homoHetCountMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=homoHetCountFinalOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[homoHetCountMergeJob ], \
											fileList=[homoHetCountFinalOutputF]))
		no_of_jobs += 1
		
		homoHetCountReduceOutputF = File(os.path.join(statOutputDir, 'homoHetCountPerSample.tsv'))
		homoHetCountReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=homoHetCountReduceOutputF, \
							extraArguments='-k 0 -v 7,3,2,4,5', transferOutput=False)
					#reduce by sampleId and aggregate only certain columns
					# 2012.8.6 ReduceMatrixBySumSameKeyColsAndThenDivide is similar to ReduceMatrixByChosenColumn.
					# it first does ReduceMatrixByChosenColumn and then divides the 1st chosen column by the 2nd chosen column.
					# that's why 7 (NoOfHet) and 3 (NoOfTotal) are the 1st two in extraArguments.
		returnData.jobDataLs.append(PassingData(jobLs=[homoHetCountReduceJob], \
											fileList=[homoHetCountReduceOutputF ]))
		outputFile = File( os.path.join(plotOutputDir, 'HeterozygoteFractionPerIndividual_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[homoHetCountReduceOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="NoOfHet_by_NoOfTotal", whichColumnPlotLabel="HetFraction", \
					logWhichColumn=False, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, homoHetCountReduceJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2
		
		siteGapOutputF = File(os.path.join(statOutputDir, 'siteGap.tsv'))
		siteGapMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteGapOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteGapMergeJob ], fileList=[siteGapOutputF]))
		outputFile = File( os.path.join(plotOutputDir, 'siteGapHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[siteGapOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="distanceToNextSite", whichColumnPlotLabel="log_distanceToNextSite", \
					logWhichColumn=True, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=samplingRate,\
					parentJobLs=[plotOutputDirJob, siteGapMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2
		
		#2012.8.13 a merging job for the selected large site gap
		largeSiteGapOutputF = File(os.path.join(statOutputDir, 'largeSiteGapMin%s.tsv'%(minSiteGap)))
		largeSiteGapMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=largeSiteGapOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[largeSiteGapMergeJob ], fileList=[largeSiteGapOutputF]))
		outputFile = File( os.path.join(plotOutputDir, 'largeSiteGap%sHist.png'%(minSiteGap)))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[largeSiteGapOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="distanceToNextSite", whichColumnPlotLabel="log_distanceToNextSite", \
					logWhichColumn=True, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, largeSiteGapMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2
		
		perIndividualHetReduceFile = File(os.path.join(statOutputDir, 'perIndividualHetReduce.tsv'))
		perIndividualHetReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=perIndividualHetReduceFile, extraArguments='-k 0 -v 1-3', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[perIndividualHetReduceJob], \
											fileList=[perIndividualHetReduceFile]))
		no_of_jobs += 1
		
		if windowSize>0:
			TiTvFinalOutputF = File(os.path.join(statOutputDir, 'TiTvWindowSize%s.tsv'%(windowSize)))
			TiTvMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=TiTvFinalOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[TiTvMergeJob], \
												fileList=[TiTvFinalOutputF ]))
			outputFnamePrefix = os.path.join(plotOutputDir, 'TiTvWindowSize%s_Plot'%(windowSize))
			# whichColumnPlotLabel and posColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
			self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[TiTvFinalOutputF], \
								outputFnamePrefix=outputFnamePrefix, \
								whichColumn=None, whichColumnLabel="Ts/Tv", whichColumnPlotLabel="Ts/Tv", need_svg=False, \
								logWhichColumn=False, valueForNonPositiveYValue=-1, \
								posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHROM", \
								minChrLength=minChrLengthForPlot, posColumnLabel="BinStart", minNoOfTotal=50,\
								figureDPI=100, ylim_type=2, samplingRate=1,\
								parentJobLs=[TiTvMergeJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
			no_of_jobs += 2
			
			TiTvSummaryReduceOutputF = File(os.path.join(statOutputDir, 'TiTv.summary.tsv'))
			TiTvSummaryReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=TiTvSummaryReduceOutputF, extraArguments='-k 0 -v 1', transferOutput=False)
			
			returnData.jobDataLs.append(PassingData(jobLs=[TiTvSummaryReduceJob], fileList=[TiTvSummaryReduceOutputF]))
			no_of_jobs += 1
		
		if windowSize>0:
			windowedPiFinalOutputF = File(os.path.join(statOutputDir, 'PiWindowSize%s.tsv'%(windowSize)))
			windowedPiMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=windowedPiFinalOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[windowedPiMergeJob, ], \
												fileList=[windowedPiFinalOutputF]))
			no_of_jobs += 1
		
		if windowSize>0:
			snpDensityOutputF = File(os.path.join(statOutputDir, 'SNPDensityByWindowSize%s.tsv'%(windowSize)))
			snpDensityMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=snpDensityOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[snpDensityMergeJob], \
												fileList=[ snpDensityOutputF]))
			outputFnamePrefix = os.path.join(plotOutputDir, 'SNPDensityWindowSize%s_Plot'%(windowSize))
			self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[snpDensityOutputF], \
								outputFnamePrefix=outputFnamePrefix, \
								whichColumn=None, whichColumnLabel="SNPS/KB", whichColumnPlotLabel="SNPS/KB", need_svg=False, \
								logWhichColumn=False, valueForNonPositiveYValue=-1, \
								posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHROM", \
								minChrLength=minChrLengthForPlot, posColumnLabel="BIN_START", minNoOfTotal=50,\
								figureDPI=100, ylim_type=2, samplingRate=1,\
								parentJobLs=[snpDensityMergeJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
			no_of_jobs += 2
		
		hweMergeFile = File(os.path.join(statOutputDir, 'hweMerge.tsv'))
		hweMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=hweMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[hweMergeJob ], \
											fileList=[hweMergeFile]))
		outputFnamePrefix = os.path.join(plotOutputDir, 'HWEPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[hweMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="P", whichColumnPlotLabel="-logHWEpvalue", need_svg=False, \
							logWhichColumn=True, valueForNonPositiveYValue=-1, \
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHR", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[hweMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		outputFnamePrefix = os.path.join(plotOutputDir, 'hetFraction')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[hweMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="hetFraction", whichColumnPlotLabel="hetFraction", need_svg=False, \
							logWhichColumn=False, valueForNonPositiveYValue=-1, \
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHR", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[hweMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		no_of_jobs += 3
		
		#2012.8.2
		siteMeanDepthMergeFile = File(os.path.join(statOutputDir, 'siteMeanDepthMerge.tsv'))
		siteMeanDepthMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteMeanDepthMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteMeanDepthMergeJob], fileList=[siteMeanDepthMergeFile]))
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteMeanDepth')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteMeanDepthMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="MEAN_DEPTH", whichColumnPlotLabel="meanDepth", need_svg=False, \
							logWhichColumn=False, valueForNonPositiveYValue=-1, \
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHROM", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[siteMeanDepthMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteDepthVar')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteMeanDepthMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="VAR_DEPTH", whichColumnPlotLabel="depthVariation", need_svg=False, \
							logWhichColumn=False, valueForNonPositiveYValue=-1, \
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHROM", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[siteMeanDepthMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		no_of_jobs += 3
		
		#reduce missingness per individual  per contig data to missingness per individual.
		imissingReduceFile = File(os.path.join(statOutputDir, 'imissingReduce.tsv'))
		imissingReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
						outputF=imissingReduceFile, extraArguments='-k 0 -v 1-4', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[imissingReduceJob ], \
											fileList=[imissingReduceFile]))
		no_of_jobs += 1
		
		#merge missingness per position
		lmissingMergeFile = File(os.path.join(statOutputDir, 'lmissingMerge.tsv'))
		lmissingMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
											outputF=lmissingMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[lmissingMergeJob], \
											fileList=[lmissingMergeFile]))
		outputFnamePrefix = os.path.join(plotOutputDir, 'lmissingPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[lmissingMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="F_MISS", whichColumnPlotLabel="missingFrequency", need_svg=False, \
							logWhichColumn=False, valueForNonPositiveYValue=-1,\
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHR", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[lmissingMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		no_of_jobs += 2
		
		#LD r2
		if LDWindowSize>0:
			LDMergeFile = File(os.path.join(statOutputDir, 'LDMerge.tsv'))
			LDMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=LDMergeFile, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[LDMergeJob], \
												fileList=[LDMergeFile]))
			outputFile = File( os.path.join(plotOutputDir, 'LDPlot%s.png'%(LDWindowSize)))
			#square the sampling rate
			LDPlotJob = self.addPlotLDJob(executable=workflow.PlotLD, inputFileList=None, outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="R^2", whichColumnPlotLabel="r2", \
						logWhichColumn=False,\
						xColumnPlotLabel="distance", chrLengthColumnHeader=None, chrColumnHeader="CHR", \
						minChrLength=minChrLengthForPlot, xColumnHeader="POS1", pos2ColumnHeader="POS2", minNoOfTotal=50,\
						figureDPI=100, ylim_type=2, samplingRate=samplingRate*samplingRate*10,\
						parentJobLs=[LDMergeJob, plotOutputDirJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput)
			no_of_jobs += 2
		else:
			LDMergeJob = None
			LDPlotJob = None
		
		#2012.8.2
		siteQualityMergeFile = File(os.path.join(statOutputDir, 'siteQualityMerge.tsv'))
		siteQualityMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteQualityMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteQualityMergeJob], fileList=[siteQualityMergeFile]))
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteQualityPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteQualityMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnLabel="QUAL", whichColumnPlotLabel="siteQuality", need_svg=False, \
							logWhichColumn=False,\
							posColumnPlotLabel="position", chrLengthColumnLabel="chrLength", chrColumnLabel="CHROM", \
							minChrLength=minChrLengthForPlot, posColumnLabel="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=samplingRate, logCount=True,\
							parentJobLs=[siteQualityMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, sshDBTunnel=self.needSSHDBTunnel)
		no_of_jobs += 2
		
		AACTallyReduceOutputF = File(os.path.join(statOutputDir, 'AAC_tally.tsv'))
		AACTallyReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=AACTallyReduceOutputF, extraArguments='-k 0 -v 1', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[AACTallyReduceJob], \
											fileList=[AACTallyReduceOutputF]))
		outputFile = File( os.path.join(plotOutputDir, 'NoOfLoci_vs_AAC.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.PlotXYAsBarChart, inputFileList=[AACTallyReduceOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="NumberOfLoci", whichColumnPlotLabel="NoOfLoci", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="AAC", xColumnPlotLabel="AlternativeAlleleCount", \
					minNoOfTotal=2,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, AACTallyReduceJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2	
		
		
		counter = 0
		no_of_vcf_files = 0
		no_of_vcf_non_empty_files = 0
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			tbi_F = getattr(jobData, 'tbi_F', None)
			inputJobLs = jobData.jobLs
			inputFname = os.path.basename(inputF.name)
			chromosome = self.getChrFromFname(inputFname)
			commonPrefix = inputFname.split('.')[0]
			vcf1 = inputF
			chr = chromosome
			no_of_vcf_non_empty_files += 1
			chr_size = chr2size.get(chr)
			if chr_size is None:
				sys.stderr.write("size for chr %s is unknown. set it to 1000.\n"%(chr))
				chr_size = 1000
			if chr_size<minChrSize:
				continue
			
			outputF = File(os.path.join(statOutputDir, "%s.homoHetCountPerSample.tsv"%(commonPrefix)))
			countHomoHetInOneVCFJob = self.addCountHomoHetInOneVCFJob(executable=workflow.CountHomoHetInOneVCF, inputF=vcf1, \
										outputF=outputF, \
							chrLength=chr_size, chromosome=chr, \
							parentJobLs=[statOutputDirJob]+inputJobLs, extraDependentInputLs=None, \
							extraArguments=None, transferOutput=False, job_max_memory=2000)
			self.addInputToStatMergeJob(workflow, statMergeJob=homoHetCountMergeJob, inputF=outputF, \
						parentJobLs=[countHomoHetInOneVCFJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=homoHetCountReduceJob, inputF=countHomoHetInOneVCFJob.output, \
						parentJobLs=[countHomoHetInOneVCFJob])
			no_of_jobs += 1
			
			outputF = File(os.path.join(statOutputDir, "%s.siteGap.tsv"%(commonPrefix)))
			outputVCFSiteGapJob = self.addCountHomoHetInOneVCFJob(executable=workflow.OutputVCFSiteGap, inputF=vcf1, \
								outputF=outputF, \
								chrLength=chr_size, chromosome=chr, \
								parentJobLs=[statOutputDirJob]+inputJobLs, extraDependentInputLs=None, \
								extraArguments=None, transferOutput=False, job_max_memory=2000)
			self.addInputToStatMergeJob(workflow, statMergeJob=siteGapMergeJob, inputF=outputF, \
						parentJobLs=[outputVCFSiteGapJob])
			no_of_jobs += 1
			
			outputFile = File(os.path.join(statOutputDir, "%s.laregSiteGapMin%s.tsv"%(commonPrefix, minSiteGap)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			selectSiteGapJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
								inputFileList=outputVCFSiteGapJob.outputLs, \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="distanceToNextSite", \
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						samplingRate=1,\
						parentJobLs=[statOutputDirJob, outputVCFSiteGapJob], \
						extraDependentInputLs=None, \
						extraArguments="-V %s"%(minSiteGap), transferOutput=False,  job_max_memory=2000)
			self.addInputToStatMergeJob(workflow, statMergeJob=largeSiteGapMergeJob, inputF=outputFile, \
						parentJobLs=[selectSiteGapJob])
			no_of_jobs += 1
			
			### vcftools job
			outputFnamePrefix = os.path.join(statOutputDir, "%s.vcftools"%(commonPrefix))
			perIndividualHetJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, \
													outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs, snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						perIndividualDepth=False, perIndividualHeterozygosity=True, \
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None)
			no_of_jobs += 1
			
			if windowSize>0:
				windowPIJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[statOutputDirJob]+inputJobLs, snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
							perIndividualDepth=False, perIndividualHeterozygosity=False, \
							perSiteHWE=False, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
							LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=windowSize, perSitePI=False, \
							SNPDensityWindowSize=None, calculateMissingNess=False,\
							extraDependentInputLs=None, transferOutput=False, \
							outputFormat=None)
				TsTvJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[statOutputDirJob]+inputJobLs, \
							TsTvWindowSize=windowSize, \
							extraDependentInputLs=None, transferOutput=False, \
							outputFormat=None)
				snpDensityJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[statOutputDirJob]+inputJobLs, snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
							perIndividualDepth=False, perIndividualHeterozygosity=False, \
							perSiteHWE=False, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
							LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=None, perSitePI=False, \
							SNPDensityWindowSize=windowSize, calculateMissingNess=False,\
							extraDependentInputLs=None, transferOutput=False, \
							outputFormat=None)
				no_of_jobs += 3
			perSiteHWEJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs, snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						perIndividualDepth=False, perIndividualHeterozygosity=False, \
						perSiteHWE=True, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
						LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=None, perSitePI=False, \
						SNPDensityWindowSize=None, calculateMissingNess=False,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None)
			missingNessJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						calculateMissingNess=True,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None, job_max_memory=1000)
			meanDepthJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs,\
						getSiteMeanDepth=True,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None, job_max_memory=1000)
			siteQualityJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs,\
						getSiteQuality=True,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None, job_max_memory=1000)
			no_of_jobs += 4
			if LDWindowSize>0 and chr_size>=minChrLengthForPlot:
				LDJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapperLD, inputVCFF=vcf1, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[statOutputDirJob]+inputJobLs, snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						genoLD=True, minLDr2=0, LDWindowByNoOfSites=None,\
						LDWindowByBP=LDWindowSize, \
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None, job_max_memory=1000)
				no_of_jobs += 1
			
			if windowSize>0:
				addChrLengthToTsTvFileJob = self.addChrLengthAppendingJob(workflow, \
								AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=TsTvJob.TsTvFile, outputF=File('%s.divByLength'%(TsTvJob.TsTvFile.name)), divideByLength=False, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[TsTvJob], \
								chromosome=chr, addChrName=True)
				windowPIFile = windowPIJob.windowPIFile
				addChrLengthToPiFileJob = self.addChrLengthAppendingJob(workflow, \
								AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=windowPIFile, outputF=File('%s.divByLength'%(windowPIFile.name)), divideByLength=False, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[windowPIJob], \
								chromosome=chr, addChrName=False)
				addChrLengthToSNPDensityFileJob = self.addChrLengthAppendingJob(workflow, \
								AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=snpDensityJob.snpDensityFile, outputF=File('%s.divByLength'%(snpDensityJob.snpDensityFile.name)), 
								divideByLength=False, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[snpDensityJob], \
								chromosome=chr, addChrName=False)
			addChrLengthToHWEJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddHetFractionToVCFtoolsHWE, \
							inputF=perSiteHWEJob.hweFile, outputF=File('%s.withChrLength'%(perSiteHWEJob.hweFile.name)), divideByLength=False, \
							chrLength=chr_size, divideStartingColumn=2, parentJobLs=[perSiteHWEJob], \
							chromosome=chr, addChrName=False)
			addChrLengthToLMissingJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
							inputF=missingNessJob.lmissFile, outputF=File('%s.withChrLength'%(missingNessJob.lmissFile.name)), divideByLength=False, \
							chrLength=chr_size, parentJobLs=[missingNessJob], \
							chromosome=chr, addChrName=False)
			addChrLengthToMeanDepthJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
							inputF=meanDepthJob.ldpethMeanFile, outputF=File('%s.withChrLength'%(meanDepthJob.ldpethMeanFile.name)), divideByLength=False, \
							chrLength=chr_size, parentJobLs=[meanDepthJob], \
							chromosome=chr, addChrName=False)
			addChrLengthToSiteQualityJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
							inputF=siteQualityJob.lqualFile, outputF=File('%s.withChrLength'%(siteQualityJob.lqualFile.name)), divideByLength=False, \
							chrLength=chr_size, parentJobLs=[siteQualityJob], \
							chromosome=chr, addChrName=False)
			no_of_jobs += 7
		
		
			
			self.addInputToStatMergeJob(workflow, statMergeJob=perIndividualHetReduceJob, inputF=perIndividualHetJob.output, \
						parentJobLs=[perIndividualHetJob])
			if windowSize>0:
				self.addInputToStatMergeJob(workflow, statMergeJob=TiTvMergeJob, inputF=addChrLengthToTsTvFileJob.output, \
							parentJobLs=[addChrLengthToTsTvFileJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=windowedPiMergeJob, inputF=addChrLengthToPiFileJob.output, \
							parentJobLs=[addChrLengthToPiFileJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=snpDensityMergeJob, inputF=addChrLengthToSNPDensityFileJob.output, \
							parentJobLs=[addChrLengthToSNPDensityFileJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=TiTvSummaryReduceJob, inputF=TsTvJob.TsTvSummaryFile, \
							parentJobLs=[TsTvJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=hweMergeJob, inputF=addChrLengthToHWEJob.output, \
						parentJobLs=[addChrLengthToHWEJob])
			if LDWindowSize>0 and chr_size>=minChrLengthForPlot:
				self.addInputToStatMergeJob(workflow, statMergeJob=LDMergeJob, inputF=LDJob.output, \
						parentJobLs=[LDJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=LDPlotJob, inputF=LDJob.output, \
						parentJobLs=[LDJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=imissingReduceJob, inputF=missingNessJob.imissFile, \
						parentJobLs=[missingNessJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=lmissingMergeJob, inputF=addChrLengthToLMissingJob.output, \
						parentJobLs=[addChrLengthToLMissingJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=siteMeanDepthMergeJob, inputF=addChrLengthToMeanDepthJob.output, \
						parentJobLs=[addChrLengthToMeanDepthJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=siteQualityMergeJob, inputF=addChrLengthToSiteQualityJob.output, \
						parentJobLs=[addChrLengthToSiteQualityJob])
			
			### TallyAACFromVCF
			# 2012.8.1 it needs the tabix index file
			outputF = File(os.path.join(statOutputDir, "%s_AAC_tally.tsv"%(commonPrefix)))
			TallyAACFromVCFJob = self.addTallyAACFromVCFJob(workflow, TallyAACFromVCF=workflow.TallyAACFromVCF, \
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						refFastaFList=refFastaFList, inputVCFF=vcf1, outputF=outputF, parentJobLs=[statOutputDirJob]+inputJobLs, \
						job_max_memory=5000, extraDependentInputLs=[tbi_F],\
						input_site_handler=input_site_handler, transferOutput=False)
			no_of_jobs += 1
			self.addInputToStatMergeJob(workflow, statMergeJob=AACTallyReduceJob, inputF=TallyAACFromVCFJob.output, \
						parentJobLs=[TallyAACFromVCFJob])
		
		sys.stderr.write("%s jobs from %s non-empty vcf files (%s total files).\n"%(no_of_jobs, \
																	no_of_vcf_non_empty_files, no_of_vcf_files))
		
		#2012.7.31 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
						outputDirPrefix="%sVCFStat"%(outputDirPrefix))
		return newReturnData
	
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		# Create a abstract dag
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		
		#have to be in front of the db_vervet connection code. Otherwise schema "genome" wont' be default path and its visible will not be visible.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
											sequence_type_id=9)
		
		#2011.11.16 initiate vervet db connection after genome db connection
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
						input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir)
		inputData = self.registerAllInputFiles(workflow, self.vcf1Dir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName="%s_%s"%(self.pegasusFolderName, vcf1Name),\
											maxContigID=self.maxContigID, \
											minContigID=self.minContigID)
		
		returnData = self.addStatCalculationJobs(workflow=workflow, inputData=inputData, refFastaFList=refFastaFList, \
									chr2size=chr2size, windowSize=self.windowSize, minChrLengthForPlot=self.minChrLengthForPlot, \
									minChrSize=self.minChrSize, LDWindowSize=self.LDWindowSize, outputDirPrefix="",\
									samplingRate=self.samplingRate, minSiteGap=30000)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
	
if __name__ == '__main__':
	main_class = CalculateVCFStatPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
