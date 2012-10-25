#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 ...  -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -I ... -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	# set pop1ss/pop2ss to 0 (or negative) if u want to take all samples 
	m=64; country1=135; country2=151; pop1=VRC; pop2=Gambia; maxCID=100; pop1ss=10; pop2ss=10; n=10;
	%s  -I /u/home/eeskin/polyacti/NetworkData/vervet/db/genotype_file/method_$m/ -E -H 
		--pop1_country_id_ls $country1 --pop2_country_id_ls $country2 --pop1Header=$pop1 --pop2Header=$pop2
		-o workflow/popgen/CompareAlleleFrequencyMethod$m\_sample$pop1ss\_$pop1\_vs_sample$pop2ss\_$pop2\_maxContigID$maxCID.xml
		-j hcondor -l hcondor -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -z localhost -C 20 -x $maxCID -a 524 --pop1_sampleSize $pop1ss
		--pop2_sampleSize $pop2ss
		--plinkIBDCheckOutputFname PlinkIBDCheck/PlinkIBDCheck_Method38_W50Z20R0.5.2012.9.13T102232/ibdCheckIBDCheck/LDPrunedMerged_ibdCheck.tsv
		--maxIBDSharing 0.1 --minAbsDelta 0.2;
		
		#-U 0 -Z 3000

Description:
	2012-10-03 workflow that investigates allele sharing between two different populations in VCF folder.
		
	#. each population is defined using tax-ID, site & country IDs
	#. split the VCF
	  #. extract individuals from individual population into separate VCF, each VCF with re-calculated AC/AF. ExtractSamplesFromVCF.py => run gatk's selectVariants to update AC/AF & DP, (run two this-kind jobs) 
	       #. arguments: tax-ID, site-ID, country-ID
	  #. JuxtaposeAlleleFrequencyFromMultiVCFInput.py VariantDiscovery.compareAlternativeAlleleFrequencyOfTwoCallSets(), output it in the format of Draw2DHistogramOfMatrix.py. header & 3 columns: AF_1, AF_2, count. 
	       #. arguments: two input file, 
	#. merge all the output and run Draw2DHistogramOfMatrix.py
	#. focus on SNPs whose AAF in either population is [0.1,0.9]
	     #. calculate how much around y=x, how much not
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from vervet.src.AbstractVervetWorkflow import AbstractVervetWorkflow

class CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
				("pop1_sampleSize", 0, int): [None, '', 1, 'if specified (>1), a sampling of individuals from population1 is used instead. otherwise, all.'],\
				("pop2_sampleSize", 0, int): [None, '', 1, 'if specified (>1), a sampling of individuals from population2 is used instead. otherwise, all.'],\
				('plinkIBDCheckOutputFname', 0, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relateness.\n\
	at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
	The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
				('maxIBDSharing', 1, float): [0.1, '', 1, 'This argument caps the maximum IBD sharing among any pair within the sampled.'],\
				('minAbsDelta', 1, float): [0.2, '', 1, 'minimum of abs(y-x) for a sample to be declared an outlier'],\
				
				("pop1_site_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
				("pop1_country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				("pop1_tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				('pop1Header', 1, ): ['', '', 1, 'used to identify population 1', ],\
				
				("pop2_site_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
				("pop2_country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				("pop2_tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
				('pop2Header', 1, ): ['', '', 1, 'used to identify population 2', ],\
			})
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		
		"""
		listArgumentName_data_type_ls = [('pop1_site_id_ls', int), ("pop1_country_id_ls", int), \
								("pop1_tax_id_ls", int), \
								('pop2_site_id_ls', int), ("pop2_country_id_ls", int), \
								("pop2_tax_id_ls", int),]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
		"""
	
	def addTranslateIDInIBDCheckResultJob(self, workflow=None, plinkIBDCheckOutputFile=None, pop_country_id_ls_str=None, \
										pop_site_id_ls_str=None, popHeader=None,\
										readGroupFile=None, parentJobLs=None, sampleIDHeader='sampleID',\
										transferOutput=False):
		"""
		2012.10.15
		"""
		job = None
		if plinkIBDCheckOutputFile:
			pop_country_id_set = set(getListOutOfStr(pop_country_id_ls_str))
			pop_site_id_set = set(getListOutOfStr(pop_site_id_ls_str))
			if 447 in pop_site_id_set or 135 in pop_country_id_set:	#either site = VRC or country = USA
				commonPrefix = os.path.splitext(plinkIBDCheckOutputFile.name)[0]
				outputFile = File('%s_%s_withReadGroup.tsv'%(commonPrefix, popHeader))
				extraArgumentList = [" --readGroupFname", readGroupFile, "--readGroupHeader %s"%(sampleIDHeader), \
									'--replaceColumnHeaderLs IID1,IID2']
				job = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.ReplaceIndividualIDInMatrixFileWithReadGroup, \
							inputFileList=None, inputFile=plinkIBDCheckOutputFile, outputFile=outputFile, \
							outputFnamePrefix=None, whichColumn=None, whichColumnHeader="sampleID", \
							minNoOfTotal=0,\
							samplingRate=1.0, \
							parentJob=None, parentJobLs=parentJobLs, \
							extraDependentInputLs=[readGroupFile], \
							extraArgumentList=extraArgumentList, transferOutput=transferOutput,  job_max_memory=1000,\
							sshDBTunnel=self.needSSHDBTunnel, \
							objectWithDBArguments=self,)
		
		return job
	
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
		#mergeSameHeaderTablesIntoOne is used here on per chromosome basis, so allow clustering
		executableClusterSizeMultiplierList.append((workflow.mergeSameHeaderTablesIntoOne, 1))
		
		
		ReplaceIndividualIDInMatrixFileWithReadGroup = Executable(namespace=namespace, name="ReplaceIndividualIDInMatrixFileWithReadGroup", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		ReplaceIndividualIDInMatrixFileWithReadGroup.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/ReplaceIndividualIDInMatrixFileWithReadGroup.py"), site_handler))
		executableClusterSizeMultiplierList.append((ReplaceIndividualIDInMatrixFileWithReadGroup, 0.5))
		
		SelectRowsWithinCoverageRange = Executable(namespace=namespace, name="SelectRowsWithinCoverageRange", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		SelectRowsWithinCoverageRange.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/SelectRowsWithinCoverageRange.py"), site_handler))
		executableClusterSizeMultiplierList.append((SelectRowsWithinCoverageRange, 0.5))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		if self.plinkIBDCheckOutputFname:
			self.plinkIBDCheckOutputFile = self.registerOneInputFile(workflow=workflow, inputFname=self.plinkIBDCheckOutputFname, input_site_handler=None, \
									folderName=self.pegasusFolderName)
		else:
			self.plinkIBDCheckOutputFile = None
	
	def addExtractSampleIDJob(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=False, \
							pop_tax_id_ls_str=None, pop_site_id_ls_str=None, pop_country_id_ls_str=None, popHeader=None,\
							pop_sampleSize=None, returnData=None, **keywords):
		"""
		2012.10.15
		"""
		if workflow is None:
			workflow = self
		#use a random VCF file as input
		jobData = passingData.chr2jobDataLs.values()[0][0]
		VCFFile = jobData.vcfFile
		inputFBaseName = os.path.basename(VCFFile.name)
		commonPrefix = inputFBaseName.split('.')[0]
		reduceOutputDirJob = passingData.reduceOutputDirJob
		#ExtractSamplesFromVCF for the 1st population
		outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID.tsv'%(commonPrefix, popHeader)))
		extraArgumentList = ['--outputFormat 2']
		if pop_tax_id_ls_str:
			extraArgumentList.append("--tax_id_ls %s"%(pop_tax_id_ls_str))
		if pop_site_id_ls_str:
			extraArgumentList.append("--site_id_ls %s"%(pop_site_id_ls_str))
		if pop_country_id_ls_str:
			extraArgumentList.append("--country_id_ls %s"%(pop_country_id_ls_str))
		
		extractPopSampleIDJob = self.addGenericDBJob(workflow=workflow, executable=self.ExtractSamplesFromVCF, inputFile=VCFFile, \
					outputFile=outputFile, inputFileList=None, \
					parentJobLs=[reduceOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=2000, \
					sshDBTunnel=self.needSSHDBTunnel, \
					key2ObjectForJob=None)
		returnData.jobDataLs.append(PassingData(jobLs=[extractPopSampleIDJob], file=extractPopSampleIDJob.output, \
											fileList=[extractPopSampleIDJob.output]))
		if pop_sampleSize and pop_sampleSize>1:
			sampleIDHeader='sampleID'	#determined by extractPopSampleIDJob
			#. SelectRowsWithinCoverageRange
			minMedianDepth = 2
			maxMedianDepth = 15
			extraArguments = " --minMedianDepth %s --maxMedianDepth %s "%(minMedianDepth, maxMedianDepth)
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID_depth%s_%s.tsv'%(commonPrefix, popHeader,\
																			minMedianDepth, maxMedianDepth)))
			selectRowsWithinCoverageRangeJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.SelectRowsWithinCoverageRange, \
						inputFileList=None, inputFile=extractPopSampleIDJob.output, outputFile=outputFile, \
						outputFnamePrefix=None, whichColumn=None, whichColumnHeader=sampleIDHeader, \
						logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=-1, \
						minNoOfTotal=0,\
						samplingRate=1.0, \
						parentJobLs=[extractPopSampleIDJob], \
						extraDependentInputLs=None, \
						extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=100,\
						sshDBTunnel=self.needSSHDBTunnel, \
						objectWithDBArguments=self,)
			returnData.jobDataLs.append(PassingData(jobLs=[selectRowsWithinCoverageRangeJob], file=selectRowsWithinCoverageRangeJob.output, \
											fileList=[selectRowsWithinCoverageRangeJob.output]))
			#. optional, a ReplaceIndividualIDInMatrixFileWithReadGroup job (for VRC) on the IBD check result
			translateIDInIBDResultJob = self.addTranslateIDInIBDCheckResultJob(workflow=workflow, plinkIBDCheckOutputFile=self.plinkIBDCheckOutputFile, \
										pop_country_id_ls_str=pop_country_id_ls_str, \
										pop_site_id_ls_str=pop_site_id_ls_str, popHeader=popHeader,\
										readGroupFile=selectRowsWithinCoverageRangeJob.output, parentJobLs=[selectRowsWithinCoverageRangeJob], \
										sampleIDHeader=sampleIDHeader, transferOutput=transferOutput)
			#. SampleRows job
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleSize%s.tsv'%(commonPrefix, popHeader, pop_sampleSize)))
			extraArgumentList = [" --sampleSize %s "%(pop_sampleSize), "--maxIBDSharing %s"%(self.maxIBDSharing)]
			if translateIDInIBDResultJob:
				extraArgumentList.extend(["--plinkIBDCheckOutputFname", translateIDInIBDResultJob.output])
				extraDependentInputLs = [translateIDInIBDResultJob.output]
				returnData.jobDataLs.append(PassingData(jobLs=[translateIDInIBDResultJob], file=translateIDInIBDResultJob.output, \
											fileList=[translateIDInIBDResultJob.output]))
			else:
				extraDependentInputLs = None
			sampleRowsJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.SampleRows, \
						inputFileList=None, inputFile=selectRowsWithinCoverageRangeJob.output, outputFile=outputFile, \
						outputFnamePrefix=None, whichColumn=None, whichColumnHeader=sampleIDHeader, \
						logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=-1, \
						minNoOfTotal=0, \
						samplingRate=1.0, \
						parentJob=translateIDInIBDResultJob, parentJobLs=[selectRowsWithinCoverageRangeJob], \
						extraDependentInputLs=extraDependentInputLs, \
						extraArgumentList=extraArgumentList, transferOutput=transferOutput,  job_max_memory=1000)
			
			returnData.jobDataLs.append(PassingData(jobLs=[sampleRowsJob], file=sampleRowsJob.output, \
											fileList=[sampleRowsJob.output]))
			#rename the extractPopSampleIDJob
			extractPopSampleIDJob = sampleRowsJob
		return extractPopSampleIDJob
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix,\
								passingData=passingData, transferOutput=transferOutput, **keywords)
		
		frequencyDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
										outputDir="%sFrequency"%(outputDirPrefix))
		passingData.frequencyDirJob = frequencyDirJob
		
		#ExtractSamplesFromVCF for the 1st population
		extractPop1SampleIDJob = self.addExtractSampleIDJob(workflow=workflow, outputDirPrefix=outputDirPrefix, \
							passingData=passingData, transferOutput=transferOutput, \
							pop_tax_id_ls_str=self.pop1_tax_id_ls, pop_site_id_ls_str=self.pop1_site_id_ls, \
							pop_country_id_ls_str=self.pop1_country_id_ls, popHeader=self.pop1Header,\
							pop_sampleSize=self.pop1_sampleSize, returnData=returnData)
		passingData.extractPop1SampleIDJob = extractPop1SampleIDJob
		
		#ExtractSamplesFromVCF for the 2nd population
		extractPop2SampleIDJob = self.addExtractSampleIDJob(workflow=workflow, outputDirPrefix=outputDirPrefix, \
							passingData=passingData, transferOutput=transferOutput, \
							pop_tax_id_ls_str=self.pop2_tax_id_ls, pop_site_id_ls_str=self.pop2_site_id_ls, \
							pop_country_id_ls_str=self.pop2_country_id_ls, popHeader=self.pop2Header,\
							pop_sampleSize=self.pop2_sampleSize, returnData=returnData)
		passingData.extractPop2SampleIDJob = extractPop2SampleIDJob
		return returnData
	
	def mapEachInterval(self, workflow=None, \
					VCFFile=None, passingData=None, transferOutput=False, **keywords):
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
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []

		topOutputDirJob = passingData.topOutputDirJob
		topOutputDir = topOutputDirJob.output
		
		intervalFnamePrefix = passingData.intervalFnamePrefix
		jobData = passingData.jobData
		
		frequencyDirJob = passingData.frequencyDirJob
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		chr = passingData.chromosome
		
		
		#1st population
		extractPop1SampleIDJob = passingData.extractPop1SampleIDJob
		outputVCF = File(os.path.join(topOutputDirJob.output, '%s_pop%s.vcf'%(intervalFnamePrefix, self.pop1Header)))
		#selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		pop1VCFConvertJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=self.SelectVariantsJava, \
				genomeAnalysisTKJar=self.genomeAnalysisTKJar, inputF=VCFFile, outputF=outputVCF, \
				refFastaFList=self.refFastaFList, sampleIDKeepFile=extractPop1SampleIDJob.output,\
				parentJobLs=[topOutputDirJob, splitVCFJob, extractPop1SampleIDJob]+jobData.jobLs, \
				extraDependentInputLs=[], transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000)
		
		#2nd population
		extractPop2SampleIDJob = passingData.extractPop2SampleIDJob
		outputVCF = File(os.path.join(topOutputDirJob.output, '%s_pop%s.vcf'%(intervalFnamePrefix, self.pop2Header)))
		#selectVariants would generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		pop2VCFConvertJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=self.SelectVariantsJava, \
				genomeAnalysisTKJar=self.genomeAnalysisTKJar, inputF=VCFFile, outputF=outputVCF, \
				refFastaFList=self.refFastaFList, sampleIDKeepFile=extractPop2SampleIDJob.output,\
				parentJobLs=[topOutputDirJob, splitVCFJob, extractPop2SampleIDJob]+jobData.jobLs, \
				extraDependentInputLs=[], transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000)
		
		
		
		#add the JuxtaposeAlleleFrequencyFromMultiVCFInput job
		outputFile = File(os.path.join(frequencyDirJob.output, '%s_AF_pop%s_vs_pop%s.tsv'%(intervalFnamePrefix, \
															self.pop1Header, self.pop2Header)))
		extraArguments = "--inputHeaderLs %s,%s"%(self.pop1Header, self.pop2Header)
		juxtaposeAFJob = self.addGenericJob(workflow=workflow, executable=self.JuxtaposeAlleleFrequencyFromMultiVCFInput, \
										inputFile=None, \
					outputFile=outputFile, inputFileList=[pop1VCFConvertJob.output, pop2VCFConvertJob.output], \
					parentJobLs=[frequencyDirJob, pop1VCFConvertJob, pop2VCFConvertJob], extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=None, job_max_memory=2000, \
					sshDBTunnel=None, \
					key2ObjectForJob=None)
		
		returnData.jobDataLs.append(PassingData(jobLs=[juxtaposeAFJob], file=juxtaposeAFJob.output, \
											fileList=[juxtaposeAFJob.output]))
		returnData.juxtaposeAFJob = juxtaposeAFJob
		return returnData
		
	def reduceEachChromosome(self, workflow=None, chromosome=None, passingData=None, mapEachVCFDataLs=None,\
					reduceEachVCFDataLs=None, \
					transferOutput=True, \
					**keywords):
		"""
		2012.10.3
			#. merge all VCF-level reduce job (from one chromosome) output (passingData.reduceEachVCFDataLs) into one first
				#taking the input jobs of each reduceEachVCFData as input of this per-chromosome reduce job.
			#. don't use passingData.mapEachVCFDataLsLs, cuz it's empty.
			
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		reduceOutputDirJob = passingData.reduceOutputDirJob
		chr = passingData.chromosome
		
		fnamePrefix = os.path.join(reduceOutputDirJob.output, '%s_frequency_juxtapose'%(chr))
		outputFile = File('%s.tsv'%(fnamePrefix))
		reduceEachChromosomeJob = self.addStatMergeJob(workflow, \
									statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
									outputF=outputFile, \
									parentJobLs=[reduceOutputDirJob],extraOutputLs=[], \
									extraDependentInputLs=[], transferOutput=transferOutput,)
		#2012.10.7 don't add it to returnData.jobDataLs unless it needs to be gzipped and transferred out
		#returnData.jobDataLs.append(PassingData(jobLs=[reduceEachChromosomeJob], file=reduceEachChromosomeJob.output, \
		#									fileList=[reduceEachChromosomeJob.output]))
		returnData.reduceEachChromosomeJob = reduceEachChromosomeJob
		
		for reduceEachVCFData in reduceEachVCFDataLs:
			for mapEachIntervalData in reduceEachVCFData.mapEachIntervalDataLs:
				juxtaposeAFJob = mapEachIntervalData.juxtaposeAFJob
				self.addInputToStatMergeJob(workflow, statMergeJob=reduceEachChromosomeJob, \
						parentJobLs=[juxtaposeAFJob])
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
		
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		fnamePrefix = os.path.join(reduceOutputDirJob.output, 'frequency_juxtapose_%s_vs_%s'%(self.pop1Header, self.pop2Header))
		outputFile = File('%s.tsv'%(fnamePrefix))
		reduceJob = self.addStatMergeJob(workflow, \
									statMergeProgram=self.MergeSameHeaderTablesIntoOne, \
									outputF=outputFile, \
									parentJobLs=[reduceOutputDirJob],extraOutputLs=[], \
									extraDependentInputLs=[], transferOutput=transferOutput,)
		returnData.jobDataLs.append(PassingData(jobLs=[reduceJob], file=reduceJob.output, \
											fileList=[reduceJob.output]))
		
		for reduceEachChromosomeData in reduceEachChromosomeDataLs:
			parentJob = reduceEachChromosomeData.reduceEachChromosomeJob
			self.addInputToStatMergeJob(workflow, statMergeJob=reduceJob, \
						parentJobLs=[parentJob])
		
		#add a Draw2DHistogramOfMatrix.py job
		outputFile = File('%s.png'%(fnamePrefix))
		drawJob = self.addDraw2DHistogramOfMatrixJob(workflow=workflow, executable=self.Draw2DHistogramOfMatrix, \
											inputFileList=None, inputFile=reduceJob.output, outputFile=outputFile, \
				outputFnamePrefix=None, whichColumn=None, whichColumnHeader=self.pop1Header, whichColumnPlotLabel=None, \
				positiveLog=False, valueForNonPositiveYValue=-1, \
				missingDataNotation='NA',\
				xColumnHeader=self.pop2Header, xColumnPlotLabel=None, \
				minNoOfTotal=100,\
				figureDPI=300, formatString='.', samplingRate=1, need_svg=False, \
				zColumnHeader=None, logX=False, logY=False, logZ=False,\
				parentJobLs=[reduceJob], \
				extraDependentInputLs=None, \
				extraArgumentList=None, extraArguments=None, transferOutput=True,  job_max_memory=2000)
		returnData.drawJob = drawJob
		
		#2012.10.15 add a EstimateOutliersIn2DData job
		extraArgumentList= ['--minAbsDelta %s'%(self.minAbsDelta)]
		outputFile = File('%s_outlierStat_minAbsDelta%s.tsv'%(fnamePrefix, self.minAbsDelta))
		estimateOutlierJob = self.addAbstractPlotJob(workflow=workflow, executable=self.EstimateOutliersIn2DData, \
					inputFileList=None, inputFile=reduceJob.output, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=self.pop1Header, whichColumnPlotLabel=None, \
					logWhichColumn=False, positiveLog=False, valueForNonPositiveYValue=-1, \
					missingDataNotation='NA',\
					xColumnHeader=self.pop2Header, xColumnPlotLabel=None, \
					minNoOfTotal=0,\
					samplingRate=1, \
					parentJob=reduceJob, parentJobLs=None, \
					extraDependentInputLs=None, \
					extraArgumentList=extraArgumentList, extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		
		returnData.jobDataLs.append(PassingData(jobLs=[estimateOutlierJob], file=estimateOutlierJob.output, \
											fileList=[estimateOutlierJob.output]))
		returnData.estimateOutlierJob = estimateOutlierJob
		
		return returnData

if __name__ == '__main__':
	main_class = CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()