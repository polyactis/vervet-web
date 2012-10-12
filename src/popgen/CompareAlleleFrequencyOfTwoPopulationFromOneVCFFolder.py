#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 ...  -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -I ... -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	%s -a 524 -I ~/NetworkData/vervet/db/genotype_file/method_63/ -E -H --pop1_country_id_ls 144 --pop2_country_id_ls 148
		--pop1Header=StKitts --pop2Header=Nevis 
		-C 1 -H -m1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -u yh -z localhost
		-o workflow/popgen/CompareAlleleFrequencyMethod63_StKitts_vs_Nevis_MaxContigID2.xml
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
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
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
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix,\
								passingData=passingData, transferOutput=transferOutput, **keywords)
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		frequencyDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, \
										outputDir="%sFrequency"%(outputDirPrefix))
		passingData.frequencyDirJob = frequencyDirJob
		
		#use a random VCF file as input
		jobData = passingData.chr2jobDataLs.values()[0][0]
		VCFFile = jobData.vcfFile
		inputFBaseName = os.path.basename(VCFFile.name)
		commonPrefix = inputFBaseName.split('.')[0]
		reduceOutputDirJob = passingData.reduceOutputDirJob
		#ExtractSamplesFromVCF for the 1st population
		outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID.tsv'%(commonPrefix, self.pop1Header)))
		extraArgumentList = ['--outputFormat 2']
		if self.pop1_tax_id_ls:
			extraArgumentList.append("--tax_id_ls %s"%(self.pop1_tax_id_ls))
		if self.pop1_site_id_ls:
			extraArgumentList.append("--site_id_ls %s"%(self.pop1_site_id_ls))
		if self.pop1_country_id_ls:
			extraArgumentList.append("--country_id_ls %s"%(self.pop1_country_id_ls))
		
		extractPop1SamplesJob = self.addGenericDBJob(workflow=workflow, executable=self.ExtractSamplesFromVCF, inputFile=VCFFile, \
					outputFile=outputFile, inputFileList=None, \
					parentJobLs=[reduceOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=2000, \
					sshDBTunnel=self.needSSHDBTunnel, \
					key2ObjectForJob=None)
		passingData.extractPop1SamplesJob = extractPop1SamplesJob
		returnData.jobDataLs.append(PassingData(jobLs=[extractPop1SamplesJob], file=extractPop1SamplesJob.output, \
											fileList=[extractPop1SamplesJob.output]))
		
		#ExtractSamplesFromVCF for the 2nd population
		outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID.tsv'%(commonPrefix, self.pop2Header)))
		extraArgumentList = ['--outputFormat 2']
		if self.pop2_tax_id_ls:
			extraArgumentList.append("--tax_id_ls %s"%(self.pop2_tax_id_ls))
		if self.pop2_site_id_ls:
			extraArgumentList.append("--site_id_ls %s"%(self.pop2_site_id_ls))
		if self.pop2_country_id_ls:
			extraArgumentList.append("--country_id_ls %s"%(self.pop2_country_id_ls))
		
		extractPop2SamplesJob = self.addGenericDBJob(workflow=workflow, executable=self.ExtractSamplesFromVCF, inputFile=VCFFile, \
					outputFile=outputFile, inputFileList=None, \
					parentJobLs=[reduceOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=2000, \
					sshDBTunnel=self.needSSHDBTunnel, \
					key2ObjectForJob=None)
		passingData.extractPop2SamplesJob = extractPop2SamplesJob
		returnData.jobDataLs.append(PassingData(jobLs=[extractPop2SamplesJob], file=extractPop2SamplesJob.output, \
											fileList=[extractPop2SamplesJob.output]))
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
		extractPop1SamplesJob = passingData.extractPop1SamplesJob
		outputVCF = File(os.path.join(topOutputDirJob.output, '%s_pop%s.vcf'%(intervalFnamePrefix, self.pop1Header)))
		subsetPop1Job = self.addVCFSubsetJob(workflow, executable=workflow.vcfSubset, vcfSubsetPath=workflow.vcfSubsetPath, \
					sampleIDFile=extractPop1SamplesJob.output,\
					inputVCF=VCFFile, outputF=outputVCF, \
					parentJobLs=[topOutputDirJob, splitVCFJob, extractPop1SamplesJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=None)
		
		#selectVariants would generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		VCF4OutputF = File(os.path.join(topOutputDir, '%s_pop%s.niceformat.vcf'%(intervalFnamePrefix, self.pop1Header)))
		pop1VCFConvertJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=self.SelectVariantsJava, \
				genomeAnalysisTKJar=self.genomeAnalysisTKJar, inputF=subsetPop1Job.output, outputF=VCF4OutputF, \
				refFastaFList=self.refFastaFList, parentJobLs=[subsetPop1Job], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, interval=chr)
		
		#2nd population
		extractPop2SamplesJob = passingData.extractPop2SamplesJob
		outputVCF = File(os.path.join(topOutputDirJob.output, '%s_pop%s.vcf'%(intervalFnamePrefix, self.pop2Header)))
		subsetPop2Job = self.addVCFSubsetJob(workflow, executable=workflow.vcfSubset, vcfSubsetPath=workflow.vcfSubsetPath, \
					sampleIDFile=extractPop2SamplesJob.output,\
					inputVCF=VCFFile, outputF=outputVCF, \
					parentJobLs=[topOutputDirJob, splitVCFJob, extractPop2SamplesJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=None)
		
		#selectVariants would generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		VCF4OutputF = File(os.path.join(topOutputDir, '%s_pop%s.niceformat.vcf'%(intervalFnamePrefix, self.pop2Header)))
		pop2VCFConvertJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=self.SelectVariantsJava, \
				genomeAnalysisTKJar=self.genomeAnalysisTKJar, inputF=subsetPop2Job.output, outputF=VCF4OutputF, \
				refFastaFList=self.refFastaFList, parentJobLs=[subsetPop2Job], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, interval=chr)
		
		
		
		#add the JuxtaposeAlleleFrequencyFromMultiVCFInput job
		outputFile = File(os.path.join(frequencyDirJob.output, '%s_AF_pop1_vs_pop2.tsv'%(intervalFnamePrefix)))
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
		
		fnamePrefix = os.path.join(reduceOutputDirJob.output, 'frequency_juxtapose')
		outputFile = File('%s.tsv'%(fnamePrefix))
		reduceJob = self.addStatMergeJob(workflow, \
									statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
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
		
		return returnData

if __name__ == '__main__':
	main_class = CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()