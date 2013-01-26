#!/usr/bin/env python
"""
Examples:
	%s
	
	# 2012.9.18
	%s  -L ~/NetworkData/vervet/db/genotype_file/method_41 -i 633,634,635,636,637,638 
		-a 524 -o workflow/HaplotypeScore/HaplotypeScore_ISQ633_638_vsMethod41.xml -l hcondor
		-j hcondor -z localhost -u yh
		-e /u/home/eeskin/polyacti
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-C 5 -H -J ~/bin/jdk/bin/java
	
Description:
	#2012.9.18 workflow to annotate variants with HaplotypeScore,
		http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_HaplotypeScore.html.
	
	Input: a VCF folder, list of alignment IDs (or use site_id, ref_ind_seq_id to filter )
	
	for each whole-genome alignment
		partition a whole-genome alignment into several interval alignment (5Mb of each chr/contig)
			each sub-alignment (GATK -L argument) is paired with the VCF from the same contig/chromosome.
			run GATK annotate-variants job
		draw histogram & GW-plot of HaplotypeScore for each alignment.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
from Pegasus.DAX3 import *
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule import VCFFile
#from AlignmentToCallPipeline import AlignmentToCallPipeline
#from AbstractVervetWorkflow import AbstractVervetWorkflow
from vervet.src import VervetDB, AbstractVervetAlignmentAndVCFWorkflow

parentClass = AbstractVervetAlignmentAndVCFWorkflow
class HaplotypeScoreWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		outputDirPrefix = passingData.outputDirPrefix
		
		#pass it along
		passingData.annotationName = 'HaplotypeScore'
		
		statOutputDir = "%sstat"%(outputDirPrefix)
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=statOutputDir)
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		self.no_of_jobs += 2
		
		returnData.plotOutputDirJob = plotOutputDirJob
		returnData.statOutputDirJob = statOutputDirJob
		return returnData
	
	def reduceBeforeEachAlignment(self, workflow=None, passingData=None, preReduceReturnData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
			add a merge variant annotation job, GW plot job
			
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		outputDirPrefix = passingData.outputDirPrefix
		
		statOutputDirJob = preReduceReturnData.statOutputDirJob
		plotOutputDirJob = preReduceReturnData.plotOutputDirJob
		
		mergeOutputF = File(os.path.join(statOutputDirJob.output, '%s_%s.tsv'%(passingData.bamFnamePrefix, passingData.annotationName)))
		mergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=mergeOutputF, transferOutput=transferOutput, parentJobLs=[statOutputDirJob],)
		returnData.jobDataLs.append(PassingData(jobLs=[mergeJob ], file=mergeJob.output, fileList=[mergeJob.output], mergeJob=mergeJob))
		self.no_of_jobs += 1
		
		outputFnamePrefix = os.path.join(plotOutputDirJob.output, '%s_%s_Plot'%(passingData.bamFnamePrefix, passingData.annotationName))
		# whichColumnPlotLabel and xColumnPlotLabel should not contain spaces or ( or ). because they will disrupt shell commandline
		self.addPlotVCFtoolsStatJob(executable=workflow.PlotVCFtoolsStat, inputFileList=[mergeOutputF], \
							outputFnamePrefix=outputFnamePrefix, \
							whichColumn=None, whichColumnHeader=passingData.annotationName, whichColumnPlotLabel=passingData.annotationName, \
							need_svg=False, \
							logY=False, valueForNonPositiveYValue=-1, \
							xColumnPlotLabel="position", chrLengthColumnHeader=None, chrColumnHeader="CHROM", \
							minChrLength=None, xColumnHeader="POS", minNoOfTotal=50,\
							figureDPI=100, ylim_type=2, samplingRate=0.01,\
							parentJobLs=[mergeJob, plotOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, sshDBTunnel=self.needSSHDBTunnel)
		self.no_of_jobs += 1
		
		outputFile = File( os.path.join(plotOutputDirJob.output, '%s_%s_Hist.png'%(passingData.bamFnamePrefix, passingData.annotationName)))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[mergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader=passingData.annotationName, whichColumnPlotLabel=passingData.annotationName, \
					logY=None, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=0.1,\
					parentJobLs=[plotOutputDirJob, mergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		self.no_of_jobs += 1
		
		return returnData
	
	def mapEachInterval(self, workflow=None, alignmentData=None, intervalData=None,\
			VCFFile=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		alignment = alignmentData.alignment
		parentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		bamFnamePrefix = passingData.bamFnamePrefix
		
		
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFnameSignature = intervalData.intervalFnameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
		
		annotationName = passingData.annotationName
		outputFile = File(os.path.join(topOutputDirJob.output, '%s_%s.%s.vcf'%(bamFnamePrefix, overlapFilenameSignature, annotationName)))
		variantAnnotatorJob = self.addGATKVariantAnnotatorJob(workflow, executable=workflow.annotateVariantJava, \
								genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, bamFile=bamF, \
								VCFFile=VCFFile, annotationName=annotationName, interval=bcftoolsInterval, outputFile=outputFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob]+parentJobLs, 
								extraDependentInputLs=[baiF, VCFFile.tbi_F], \
								transferOutput=False, \
								extraArguments=None, job_max_memory=4000)
		
		outputFile = File(os.path.join(topOutputDirJob.output, '%s_%s.%s.tsv'%(bamFnamePrefix, overlapFilenameSignature, annotationName)))
		extractInfoJob = self.addGenericJob(workflow=workflow, executable=workflow.ExtractInfoFromVCF, inputFile=variantAnnotatorJob.output, \
						inputArgumentOption="-i", \
						outputFile=outputFile, outputArgumentOption="-o", \
						parentJobLs=[variantAnnotatorJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=False, \
						extraArguments="-k %s"%(annotationName), extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
						key2ObjectForJob=None)
		
		returnData.jobDataLs.append(PassingData(jobLs=[variantAnnotatorJob, extractInfoJob], file=variantAnnotatorJob.output, \
											fileList=[variantAnnotatorJob.output, extractInfoJob.output]))
		returnData.variantAnnotatorJob=variantAnnotatorJob
		returnData.extractInfoJob=extractInfoJob
		#add the sub-alignment to the alignment merge job
		self.no_of_jobs += 2
		return returnData
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, \
					reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		"""
		for jobData in reduceBeforeEachAlignmentData.jobDataLs:
			self.addInputToStatMergeJob(workflow, statMergeJob=jobData.mergeJob, inputF=mapEachIntervalData.extractInfoJob.output, \
							parentJobLs=[mapEachIntervalData.extractInfoJob])
	
	
	def addGATKVariantAnnotatorJob(self, workflow, executable=None, genomeAnalysisTKJar=None, bamFile=None, \
								VCFFile=None, annotationName='HaplotypeScore', interval=None, outputFile=None, \
								refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
								extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, **keywords):
		"""
		2011-12-4
			inputFile is a bam file.
			outputFile is a VCF file.
		
			java -Xmx2g -jar  ~/script/gatk/dist/GenomeAnalysisTK.jar
				-R ~/NetworkData/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
				-T VariantAnnotator  -I ~/NetworkData/vervet/db/individual_alignment/925_634_vs_524_by_2.bam 
				-A HaplotypeScore -o /tmp/Aln925_HaplotypeScore.vcf
				--variant ~/NetworkData/vervet/db/genotype_file/method_44_47544_VCF_998chromosomes_998VCFIntoOne.vcf.gz
			
		"""
		#GATK job
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
		#javaMemRequirement = "-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*50/100, job_max_memory, \
		#																			MaxPermSize*50/100, MaxPermSize)
		refFastaFile = refFastaFList[0]
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', genomeAnalysisTKJar, "-T", "VariantAnnotator",\
						"-I", bamFile, "-o", outputFile, \
						"--annotation", annotationName, self.defaultGATKArguments, \
						" --variant", VCFFile, '-L', interval, "-R", refFastaFile]
		#bamFile and outputFile have to be put after -jar and java memory requirement.
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([bamFile, VCFFile] + refFastaFList)
		
		job= self.addGenericJob(executable=executable, inputFile=None, inputArgumentOption="-I",\
							outputFile=None, outputArgumentOption="-o", \
							parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
							extraOutputLs=[outputFile],\
							transferOutput=transferOutput, \
							extraArgumentList=extraArgumentList, \
							job_max_memory=memRequirementData.memRequirement, **keywords)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		2011-11-28
		"""
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		
		annotateVariantJava = Executable(namespace=namespace, name="annotateVariantJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		annotateVariantJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((annotateVariantJava, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	


if __name__ == '__main__':
	main_class = HaplotypeScoreWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()