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
	# try "-s 0.01" to increase data sampling rate.
	# --intervalSize refers to the number of SNP loci, not chromosome distance. 
	m=64; country=135; pop=VRC; maxCID=100; popSS=50; n=10;
	m=64; country=144,148; pop=StKittsNevis; maxCID=100; popSS=0; n=10;	#popSS=0 means take all samples from that population
	maxIBDSharing=0.15;
	%s -a 524 -I ~/NetworkData/vervet/db/genotype_file/method_$m/
		-o workflow/VCFStat/VCFStat_Method$m\_L800000P4000000n1000000_sample$popSS\_$pop\_maxContigID$maxCID\_maxIBDSharing$maxIBDSharing.xml
		-j hcondor -l hcondor -u yh -z localhost
		-D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-L 800000 -P 4000000 -n 1000000 -H --intervalSize 2000
		--intervalOverlapSize 500   -x $maxCID -C 3 -s 0.01  --pop_sampleSize $popSS --pop_country_id_ls $country
		--popHeader=$pop
		--plinkIBDCheckOutputFname PlinkIBDCheck/PlinkIBDCheck_Method38_W50Z20R0.5.2012.9.13T102232/ibdCheckIBDCheck/LDPrunedMerged_ibdCheck.tsv
		--maxIBDSharing=0.15
		#-V 90 
		
	

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from vervet.src.AbstractVervetWorkflow import AbstractVervetWorkflow

class CalculateVCFStatPipeline(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
#	option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
				("pop_sampleSize", 0, int): [None, '', 1, 'if specified (>1), a sampling of individuals from population2 is used instead. otherwise, all.'],\
				('plinkIBDCheckOutputFname', 0, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relateness.\n\
	at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
	The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
				('maxIBDSharing', 1, float): [0.1, '', 1, 'This argument caps the maximum IBD sharing among any pair within the sampled.'],\
				
				("pop_site_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
				("pop_country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				("pop_tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				('popHeader', 1, ): ['', '', 1, 'used to identify population 1', ],\
				
						('windowSize', 1, int): [200000, 'w', 1, "window size for TiTv, pi, snpDensity calculation by vcftools.\
			set it to 0 if you do not need these stats.", ],\
						('LDWindowSize', 0, int): [0, 'L', 1, 'window size for LD calculation by vcftools, set it to 0 to skip LD jobs.', ],\
						('minChrLengthForPlot', 1, int): [2000000, 'P', 1, 'minimum chromosome size for a chromosome to be included in plot', ],\
						('minChrSize', 1, int): [500000, 'n', 1, 'minimum chromosome size for any computing ', ],\
						('minSiteGap', 1, int): [10000, '', 1, 'minimum site gap to be taken as large gap', ],\
						('samplingRate', 1, float): [0.001, 's', 1, "sampling rate for plotting the locus-specific stats.\
			The LD plot's sampling rate is samplingRate*samplingRate*10 "],\
						})
#						('vcf1Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		
		#self.vcf1Dir = os.path.abspath(self.vcf1Dir)
	
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

	
	
	def addTallyAACFromVCFJob(self, workflow, TallyAACFromVCF=None, GenomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputF=None, parentJobLs=None, \
							job_max_memory=1000, extraDependentInputLs=[],\
							input_site_handler=None, transferOutput=False, **keywords):
		"""
		2012.7.25 why *.vcf.tbi file is needed? get rid of it. 
		"""
		# Add a mkdir job for any directory.
		TallyAACFromVCFJob = Job(namespace=workflow.namespace, name=TallyAACFromVCF.name, version=workflow.version)
		refFastaF = refFastaFList[0]
		TallyAACFromVCFJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", GenomeAnalysisTKJar, "-R", refFastaF, 
						"-T TallyAACFromVCF", "--variant:VCF", inputVCFF, "-o", outputF)
		for refFastaFile in refFastaFList:
			TallyAACFromVCFJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		TallyAACFromVCFJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)
		self.addJobUse(TallyAACFromVCFJob, file=GenomeAnalysisTKJar, transfer=True, register=True, link=Link.INPUT)
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
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow)
		
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
	
		if self.plinkIBDCheckOutputFname:
			self.plinkIBDCheckOutputFile = self.registerOneInputFile(workflow=workflow, inputFname=self.plinkIBDCheckOutputFname, input_site_handler=None, \
									folderName=self.pegasusFolderName)
		else:
			self.plinkIBDCheckOutputFile = None
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		
		returnData = AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix,\
								passingData=passingData, transferOutput=transferOutput, **keywords)
		
		#ExtractSamplesFromVCF for the 1st population
		extractPopSampleIDJob = self.addExtractSampleIDJob(workflow=workflow, outputDirPrefix=outputDirPrefix, \
							passingData=passingData, transferOutput=transferOutput, \
							pop_tax_id_ls_str=self.pop_tax_id_ls, pop_site_id_ls_str=self.pop_site_id_ls, \
							pop_country_id_ls_str=self.pop_country_id_ls, popHeader=self.popHeader,\
							pop_sampleSize=self.pop_sampleSize, returnData=returnData)
		passingData.extractPopSampleIDJob = extractPopSampleIDJob
		
		statOutputDir = "%sStat"%(outputDirPrefix)
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=statOutputDir)
		passingData.statOutputDirJob = statOutputDirJob
		
		plotOutputDir = "%sPlot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		passingData.plotOutputDirJob = plotOutputDirJob
		
		input_site_handler = self.input_site_handler
		
		homoHetCountFinalOutputF = File(os.path.join(statOutputDir, 'homoHetCountPerSamplePerContig.tsv'))
		homoHetCountMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=homoHetCountFinalOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[homoHetCountMergeJob ], \
											fileList=[homoHetCountFinalOutputF]))
		passingData.homoHetCountMergeJob = homoHetCountMergeJob
		
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
		passingData.homoHetCountReduceJob = homoHetCountReduceJob
		
		outputFile = File( os.path.join(plotOutputDir, 'HeterozygoteFractionPerIndividual_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[homoHetCountReduceOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="NoOfHet_by_NoOfTotal", whichColumnPlotLabel="HetFraction", \
					logY=False, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, homoHetCountReduceJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		siteGapOutputF = File(os.path.join(statOutputDir, 'siteGap.tsv'))
		siteGapMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteGapOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteGapMergeJob ], fileList=[siteGapOutputF]))
		passingData.siteGapMergeJob = siteGapMergeJob
		
		outputFile = File( os.path.join(plotOutputDir, 'siteGapHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[siteGapOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="distanceToNextSite", whichColumnPlotLabel="log_distanceToNextSite", \
					logY=1, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=self.samplingRate,\
					parentJobLs=[plotOutputDirJob, siteGapMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		#2012.8.13 a merging job for the selected large site gap
		largeSiteGapOutputF = File(os.path.join(statOutputDir, 'largeSiteGapMin%s.tsv'%(self.minSiteGap)))
		largeSiteGapMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=largeSiteGapOutputF, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[largeSiteGapMergeJob ], fileList=[largeSiteGapOutputF]))
		passingData.largeSiteGapMergeJob = largeSiteGapMergeJob
		
		outputFile = File( os.path.join(plotOutputDir, 'largeSiteGap%sHist.png'%(self.minSiteGap)))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[largeSiteGapOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="distanceToNextSite", whichColumnPlotLabel="log_distanceToNextSite", \
					logY=1, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, largeSiteGapMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		
		perIndividualHetReduceFile = File(os.path.join(statOutputDir, 'perIndividualHetReduce.tsv'))
		perIndividualHetReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=perIndividualHetReduceFile, extraArguments='-k 0 -v 1-3', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[perIndividualHetReduceJob], \
											fileList=[perIndividualHetReduceFile]))
		passingData.perIndividualHetReduceJob = perIndividualHetReduceJob
		
		if self.windowSize>0:
			TiTvFinalOutputF = File(os.path.join(statOutputDir, 'TiTvWindowSize%s.tsv'%(self.windowSize)))
			TiTvMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=TiTvFinalOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[TiTvMergeJob], \
												fileList=[TiTvFinalOutputF ]))
			passingData.TiTvMergeJob = TiTvMergeJob
			
			outputFnamePrefix = os.path.join(plotOutputDir, 'TiTvWindowSize%s_Plot'%(self.windowSize))
			# whichColumnPlotLabel and xColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
			self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[TiTvFinalOutputF], \
								outputFnamePrefix=outputFnamePrefix, \
								whichColumn=None, whichColumnHeader="Ts/Tv", whichColumnPlotLabel="Ts/Tv", need_svg=False, \
								logY=False, valueForNonPositiveYValue=-1, \
								xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHROM", \
								minChrLength=self.minChrLengthForPlot, xColumnHeader="BinStart", minNoOfTotal=50,\
								figureDPI=100, ylim_type=2, samplingRate=1,\
								parentJobLs=[TiTvMergeJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
			TiTvSummaryReduceOutputF = File(os.path.join(statOutputDir, 'TiTv.summary.tsv'))
			TiTvSummaryReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=TiTvSummaryReduceOutputF, extraArguments='-k 0 -v 1', transferOutput=False)
			
			returnData.jobDataLs.append(PassingData(jobLs=[TiTvSummaryReduceJob], fileList=[TiTvSummaryReduceOutputF]))
			passingData.TiTvSummaryReduceJob = TiTvSummaryReduceJob
		
		if self.windowSize>0:
			windowedPiFinalOutputF = File(os.path.join(statOutputDir, 'PiWindowSize%s.tsv'%(self.windowSize)))
			windowedPiMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=windowedPiFinalOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[windowedPiMergeJob, ], \
												fileList=[windowedPiFinalOutputF]))
			passingData.windowedPiMergeJob = windowedPiMergeJob
		
		
		if self.windowSize>0:
			snpDensityOutputF = File(os.path.join(statOutputDir, 'SNPDensityByWindowSize%s.tsv'%(self.windowSize)))
			snpDensityMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=snpDensityOutputF, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[snpDensityMergeJob], \
												fileList=[ snpDensityOutputF]))
			passingData.snpDensityMergeJob = snpDensityMergeJob
			
			outputFnamePrefix = os.path.join(plotOutputDir, 'SNPDensityWindowSize%s_Plot'%(self.windowSize))
			self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[snpDensityOutputF], \
								outputFnamePrefix=outputFnamePrefix, \
								whichColumn=None, whichColumnHeader="SNPS/KB", whichColumnPlotLabel="SNPS/KB", need_svg=False, \
								logY=False, valueForNonPositiveYValue=-1, \
								xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHROM", \
								minChrLength=self.minChrLengthForPlot, xColumnHeader="BIN_START", minNoOfTotal=50,\
								figureDPI=100, ylim_type=2, samplingRate=1,\
								parentJobLs=[snpDensityMergeJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
			
		
		hweMergeFile = File(os.path.join(statOutputDir, 'hweMerge.tsv'))
		hweMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=hweMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[hweMergeJob ], \
											fileList=[hweMergeFile]))
		passingData.hweMergeJob = hweMergeJob
		
		outputFnamePrefix = os.path.join(plotOutputDir, 'HWEPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[hweMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="P", whichColumnPlotLabel="-logHWEpvalue", need_svg=False, \
							logY=2, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHR", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[hweMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		outputFnamePrefix = os.path.join(plotOutputDir, 'hetFraction')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[hweMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="hetFraction", whichColumnPlotLabel="hetFraction", need_svg=False, \
							logY=False, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHR", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[hweMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		
		#2012.8.2
		siteMeanDepthMergeFile = File(os.path.join(statOutputDir, 'siteMeanDepthMerge.tsv'))
		siteMeanDepthMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteMeanDepthMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteMeanDepthMergeJob], fileList=[siteMeanDepthMergeFile]))
		passingData.siteMeanDepthMergeJob = siteMeanDepthMergeJob
		
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteMeanDepth')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteMeanDepthMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="MEAN_DEPTH", whichColumnPlotLabel="meanDepth", need_svg=False, \
							logY=False, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHROM", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[siteMeanDepthMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteDepthVar')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteMeanDepthMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="VAR_DEPTH", whichColumnPlotLabel="depthVariation", need_svg=False, \
							logY=False, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHROM", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[siteMeanDepthMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		
		#reduce missingness per individual  per contig data to missingness per individual.
		imissingReduceFile = File(os.path.join(statOutputDir, 'imissingReduce.tsv'))
		imissingReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
						outputF=imissingReduceFile, extraArguments='-k 0 -v 1-4', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[imissingReduceJob ], \
											fileList=[imissingReduceFile]))
		passingData.imissingReduceJob = imissingReduceJob
		
		
		#merge missingness per position
		lmissingMergeFile = File(os.path.join(statOutputDir, 'lmissingMerge.tsv'))
		lmissingMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
											outputF=lmissingMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[lmissingMergeJob], \
											fileList=[lmissingMergeFile]))
		passingData.lmissingMergeJob = lmissingMergeJob
		
		outputFnamePrefix = os.path.join(plotOutputDir, 'lmissingPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[lmissingMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="F_MISS", whichColumnPlotLabel="missingFrequency", need_svg=False, \
							logY=False, valueForNonPositiveYValue=-1,\
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHR", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[lmissingMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		
		#LD r2
		if self.LDWindowSize>0:
			LDMergeFile = File(os.path.join(statOutputDir, 'LDMerge.tsv'))
			LDMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=LDMergeFile, transferOutput=False)
			returnData.jobDataLs.append(PassingData(jobLs=[LDMergeJob], \
												fileList=[LDMergeFile]))
			passingData.LDMergeJob = LDMergeJob
			
			windowStep = int(self.LDWindowSize/30)
			candidateLDDistList = range(0, self.LDWindowSize+1, windowStep)[1:]
			for maxDist in candidateLDDistList:
				for movingAverageType in [1,2,3]:
					#square the sampling rate
					#more sampling for shorter distance, because shorter distance is additional filter after random sampling
					if movingAverageType==3:
						whichColumnPlotLabel="r2Above0.8"
					elif movingAverageType==2:
						whichColumnPlotLabel="mean_r2"
					else:
						whichColumnPlotLabel="median_r2"
					
					LDSamplingRate = min(1, self.samplingRate*self.samplingRate*self.LDWindowSize/float(maxDist))
					outputFile = File( os.path.join(plotOutputDir, 'LDPlotMaxDist%s_movingAverageType%s_Sampling%.6f.png'%(maxDist, \
																								movingAverageType, LDSamplingRate)))
					passingData.LDPlotJob = self.addPlotLDJob(executable=workflow.PlotLD, inputFileList=[LDMergeJob.output], outputFile=outputFile, \
								whichColumn=None, whichColumnHeader="R^2", whichColumnPlotLabel=whichColumnPlotLabel, \
								logY=False,\
								xColumnPlotLabel="distance", chrLengthColumnHeader=None, chrColumnHeader="CHR", \
								minChrLength=self.minChrLengthForPlot, xColumnHeader="POS1", pos2ColumnHeader="POS2", minNoOfTotal=50,\
								figureDPI=100, ylim_type=2, samplingRate=LDSamplingRate,\
								minDist=None, maxDist=maxDist, movingAverageType=movingAverageType,\
								parentJobLs=[LDMergeJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=True)
		else:
			passingData.LDMergeJob = None
			passingData.LDPlotJob = None
		
		#2012.8.2
		siteQualityMergeFile = File(os.path.join(statOutputDir, 'siteQualityMerge.tsv'))
		siteQualityMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=siteQualityMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[siteQualityMergeJob], fileList=[siteQualityMergeFile]))
		passingData.siteQualityMergeJob = siteQualityMergeJob
		
		outputFnamePrefix = os.path.join(plotOutputDir, 'siteQualityPlot')
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[siteQualityMergeFile], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader="QUAL", whichColumnPlotLabel="siteQuality", need_svg=False, \
							logY=False,\
							xColumnPlotLabel="position", chrLengthColumnHeader="chrLength", chrColumnHeader="CHROM", \
							minChrLength=self.minChrLengthForPlot, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=self.samplingRate, logCount=True,\
							parentJobLs=[siteQualityMergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		
		
		AACTallyReduceOutputF = File(os.path.join(statOutputDir, 'AAC_tally.tsv'))
		AACTallyReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=AACTallyReduceOutputF, extraArguments='-k 0 -v 1', transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[AACTallyReduceJob], \
											fileList=[AACTallyReduceOutputF]))
		passingData.AACTallyReduceJob = AACTallyReduceJob
		
		outputFile = File( os.path.join(plotOutputDir, 'NoOfLoci_vs_AAC.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.PlotYAsBar, inputFileList=[AACTallyReduceOutputF], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="NumberOfLoci", whichColumnPlotLabel="NoOfLoci", \
					logY=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="AAC", xColumnPlotLabel="AlternativeAlleleCount", \
					minNoOfTotal=2,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, AACTallyReduceJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		return returnData
	
	def mapEachInterval(self, workflow=None, VCFFile=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.10.3
			#. extract individuals from individual population into separate VCF, each VCF with re-calculated AC/AF.
				#. ExtractSamplesFromVCF.py
				#. run gatk's selectVariants to update AC/AF & DP, (run two this-kind jobs) 
				#. arguments: tax-ID, site-ID, country-ID
			#. JuxtaposeAlleleFrequencyFromMultiVCFInput.py
				#. output it in the format of Draw2DHistogramOfMatrix.py. header & 3 columns: AF_1, AF_2, count. 
				#. arguments: two input file,
		"""
		if workflow is None:
			workflow = self
		returnData = AbstractVervetWorkflow.mapEachInterval(self, workflow=workflow, VCFFile=VCFFile, \
												passingData=passingData, transferOutput=transferOutput, **keywords)
		
		topOutputDirJob = passingData.topOutputDirJob
		topOutputDir = topOutputDirJob.output
		
		intervalFnamePrefix = passingData.intervalFnamePrefix
		jobData = passingData.jobData
		
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		chr = passingData.chromosome
		
		chr_size = self.chr2size.get(chr)
		if chr_size is None:
			sys.stderr.write("size for chr %s is unknown. set it to 1000.\n"%(chr))
			chr_size = 1000
		if chr_size<self.minChrSize:
			return returnData
		commonPrefix = intervalFnamePrefix
		
		#1st population
		extractPopSampleIDJob = passingData.extractPopSampleIDJob
		outputVCF = File(os.path.join(topOutputDirJob.output, '%s_pop%s.vcf'%(intervalFnamePrefix, self.popHeader)))
		#selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		popVCFConvertJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=self.SelectVariantsJava, \
				GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, inputF=VCFFile, outputF=outputVCF, \
				refFastaFList=self.refFastaFList, sampleIDKeepFile=extractPopSampleIDJob.output,\
				parentJobLs=[topOutputDirJob, splitVCFJob, extractPopSampleIDJob]+jobData.jobLs, \
				extraDependentInputLs=[], transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000)
		
		outputF = File(os.path.join(topOutputDirJob.output, "%s.homoHetCountPerSample.tsv"%(commonPrefix)))
		countHomoHetInOneVCFJob = self.addCountHomoHetInOneVCFJob(executable=workflow.CountHomoHetInOneVCF, inputF=popVCFConvertJob.output, \
									outputF=outputF, \
						chrLength=chr_size, chromosome=chr, \
						parentJobLs=[topOutputDirJob, popVCFConvertJob], extraDependentInputLs=None, \
						extraArguments=None, transferOutput=False, job_max_memory=2000)
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.homoHetCountMergeJob, inputF=outputF, \
					parentJobLs=[countHomoHetInOneVCFJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.homoHetCountReduceJob, inputF=countHomoHetInOneVCFJob.output, \
					parentJobLs=[countHomoHetInOneVCFJob])
		
		
		outputF = File(os.path.join(topOutputDirJob.output, "%s.siteGap.tsv"%(commonPrefix)))
		outputVCFSiteGapJob = self.addCountHomoHetInOneVCFJob(executable=workflow.OutputVCFSiteGap, inputF=popVCFConvertJob.output, \
							outputF=outputF, \
							chrLength=chr_size, chromosome=chr, \
							parentJobLs=[popVCFConvertJob], extraDependentInputLs=None, \
							extraArguments=None, transferOutput=False, job_max_memory=2000)
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.siteGapMergeJob, inputF=outputF, \
					parentJobLs=[outputVCFSiteGapJob])
		
		outputFile = File(os.path.join(topOutputDirJob.output, "%s.laregSiteGapMin%s.tsv"%(commonPrefix, self.minSiteGap)))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		selectSiteGapJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
							inputFileList=outputVCFSiteGapJob.outputLs, \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="distanceToNextSite", \
					logY=False, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					samplingRate=1,\
					parentJobLs=[topOutputDirJob, outputVCFSiteGapJob], \
					extraDependentInputLs=None, \
					extraArguments="-V %s"%(self.minSiteGap), transferOutput=False,  job_max_memory=2000)
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.largeSiteGapMergeJob, inputF=outputFile, \
					parentJobLs=[selectSiteGapJob])
		
		
		### vcftools job
		outputFnamePrefix = os.path.join(topOutputDirJob.output, "%s.vcftools"%(commonPrefix))
		perIndividualHetJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, \
												outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob], snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
					perIndividualDepth=False, perIndividualHeterozygosity=True, \
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None)
		
		if self.windowSize>0:
			windowPIJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, \
													outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[topOutputDirJob]+[popVCFConvertJob], snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						perIndividualDepth=False, perIndividualHeterozygosity=False, \
						perSiteHWE=False, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
						LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=self.windowSize, perSitePI=False, \
						SNPDensityWindowSize=None, calculateMissingNess=False,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None)
			TsTvJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[topOutputDirJob]+[popVCFConvertJob], \
						TsTvWindowSize=self.windowSize, \
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None)
			snpDensityJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[topOutputDirJob]+[popVCFConvertJob], snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
						perIndividualDepth=False, perIndividualHeterozygosity=False, \
						perSiteHWE=False, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
						LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=None, perSitePI=False, \
						SNPDensityWindowSize=self.windowSize, calculateMissingNess=False,\
						extraDependentInputLs=None, transferOutput=False, \
						outputFormat=None)
		perSiteHWEJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob], snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
					perIndividualDepth=False, perIndividualHeterozygosity=False, \
					perSiteHWE=True, haploLD=False, genoLD=False, LDWindowByNoOfSites=None,\
					LDWindowByBP=None, TsTvWindowSize=None, piWindowSize=None, perSitePI=False, \
					SNPDensityWindowSize=None, calculateMissingNess=False,\
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None)
		missingNessJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob], minMAC=None, minMAF=None, maxSNPMissingRate=None,\
					calculateMissingNess=True,\
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None, job_max_memory=1000)
		meanDepthJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob],\
					getSiteMeanDepth=True,\
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None, job_max_memory=1000)
		siteQualityJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob],\
					getSiteQuality=True,\
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None, job_max_memory=1000)
		if self.LDWindowSize>0 and chr_size>=self.minChrLengthForPlot:
			LDJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapperLD, inputVCFF=popVCFConvertJob.output, outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob], snpMisMatchStatFile=None, minMAC=None, minMAF=None, maxSNPMissingRate=None,\
					genoLD=True, minLDr2=0, LDWindowByNoOfSites=None,\
					LDWindowByBP=self.LDWindowSize, \
					extraDependentInputLs=None, transferOutput=False, \
					outputFormat=None, job_max_memory=1000)
		
		if self.windowSize>0:
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
		
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.perIndividualHetReduceJob, inputF=perIndividualHetJob.output, \
					parentJobLs=[perIndividualHetJob])
		if self.windowSize>0:
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.TiTvMergeJob, inputF=addChrLengthToTsTvFileJob.output, \
						parentJobLs=[addChrLengthToTsTvFileJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.windowedPiMergeJob, inputF=addChrLengthToPiFileJob.output, \
						parentJobLs=[addChrLengthToPiFileJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.snpDensityMergeJob, inputF=addChrLengthToSNPDensityFileJob.output, \
						parentJobLs=[addChrLengthToSNPDensityFileJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.TiTvSummaryReduceJob, inputF=TsTvJob.TsTvSummaryFile, \
						parentJobLs=[TsTvJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.hweMergeJob, inputF=addChrLengthToHWEJob.output, \
					parentJobLs=[addChrLengthToHWEJob])
		if self.LDWindowSize>0 and chr_size>=self.minChrLengthForPlot:
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.LDMergeJob, inputF=LDJob.output, \
					parentJobLs=[LDJob])
			#2012.10.25 stop the LDPlotJob on every input
			#self.addInputToStatMergeJob(workflow, statMergeJob=passingData.LDPlotJob, inputF=LDJob.output, \
			#		parentJobLs=[LDJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.imissingReduceJob, inputF=missingNessJob.imissFile, \
					parentJobLs=[missingNessJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.lmissingMergeJob, inputF=addChrLengthToLMissingJob.output, \
					parentJobLs=[addChrLengthToLMissingJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.siteMeanDepthMergeJob, inputF=addChrLengthToMeanDepthJob.output, \
					parentJobLs=[addChrLengthToMeanDepthJob])
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.siteQualityMergeJob, inputF=addChrLengthToSiteQualityJob.output, \
					parentJobLs=[addChrLengthToSiteQualityJob])
		
		### TallyAACFromVCF
		# 2012.8.1 it needs the tabix index file
		outputF = File(os.path.join(topOutputDirJob.output, "%s_AAC_tally.tsv"%(commonPrefix)))
		TallyAACFromVCFJob = self.addTallyAACFromVCFJob(workflow, TallyAACFromVCF=workflow.TallyAACFromVCF, \
					GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
					refFastaFList=self.refFastaFList, inputVCFF=popVCFConvertJob.output, outputF=outputF, \
					parentJobLs=[topOutputDirJob]+[popVCFConvertJob], \
					job_max_memory=5000, extraDependentInputLs=[],\
					transferOutput=False)
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.AACTallyReduceJob, inputF=TallyAACFromVCFJob.output, \
					parentJobLs=[TallyAACFromVCFJob])
		
		return returnData

	def reduce(self, workflow=None, passingData=None, reduceEachChromosomeDataLs=None,\
			transferOutput=True, **keywords):
		"""
		2012.10.3
			#. reduce all previous jobs (passingData.reduceEachChromosomeDataLs) into one final output
			#. merge all the output and run Draw2DHistogramOfMatrix.py
		
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def addStatCalculationJobs(self, workflow=None, inputData=None, refFastaFList=None, chr2size=None,\
						windowSize=None, minChrLengthForPlot=None, minChrSize=None, LDWindowSize=None, \
						transferOutput=True, outputDirPrefix="", samplingRate=0.001, minSiteGap=10000):
		"""
		2012.10.25
			call self.addAllJobs() instead
		2012.8.13
			if windowSize=0, no jobs for Ti/Tv, snpDensity, pi
		2012.8.2
			these arguments won't be used: minMAC=None, minMAF=None, maxSNPMissingRate=None, 
		"""
		if workflow is None:
			workflow = self
		
		return self.addAllJobs(workflow=workflow, inputVCFData=inputData, \
					chr2IntervalDataLs=None, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				refFastaFList=refFastaFList,\
				needFastaIndexJob=getattr(self, 'needFastaIndexJob',False), needFastaDictJob=getattr(self, 'needFastaDictJob', False), \
				data_dir=self.data_dir, no_of_gatk_threads = 1,\
				intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
				outputDirPrefix=outputDirPrefix, transferOutput=transferOutput)
	
	def connectDB(self):
		"""
		2012.9.24
			establish db connection for all derivative classes
		"""
		
		#have to be in front of the db_vervet connection code. Otherwise schema "genome" wont' be default path and its visible will not be visible.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, db_user=self.db_user,
						db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		self.chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
											sequence_type_id=9)
		
		AbstractVervetWorkflow.connectDB(self)
	
	
if __name__ == '__main__':
	main_class = CalculateVCFStatPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
