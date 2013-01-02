#!/usr/bin/env python
"""
Examples:
	# 2012.9.21 run base quality recalibration on VRC alignments (-S 447), individual_sequence_id from 639-642 (-i ...)
	# filtered sequences (-Q 1), alignment method 2 (-G 2)
	# -N 1000 (top 1000 contigs)
	#  -Z 10000000 (10 million bp for each interval) -B 30000 (30kb overlap between intervals),
	%s -L ~/NetworkData/vervet/db/genotype_file/method_17/ -i 639-642
		-S 447 -u yh -z localhost -Q1 -G2 -a 524 -o workflow/BaseQualityRecalibration/BaseQualityRecalibration_VRC447_vsMethod17.xml
		-l hcondor -j hcondor -z localhost -u yh -N 1000  -Z 10000000 -B 30000
		-e /u/home/eeskin/polyacti
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-H -J ~/bin/jdk/bin/java --mask_genotype_method_id 17

	# 2012.9.18
	%s  -L ~/NetworkData/vervet/db/genotype_file/method_41 -i 633,634,635,636,637,638 
		-a 524 -o workflow/BaseQualityRecalibration/BaseQualityRecalibration_ISQ633_638_vsMethod41.xml -l hcondor
		-j hcondor -z localhost -u yh -Z 10000000 -B 30000
		-e /u/home/eeskin/polyacti
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-C 5 -H -J ~/bin/jdk/bin/java --mask_genotype_method_id 41
		
Description:
	#2012.9.21  recalibration needs ~500K reads to get accurate estimate. So adjust the -Z according to coverage.  
	Input: a VCF folder, list of alignment IDs (or use site_id, ref_ind_seq_id to filter )
	
	for each whole-genome alignment
	select a whole-genome alignment into several interval alignment (5Mb of each chr/contig)
			using samtools view -h ... (AlignmentToCallPipeline.addSelectAndSplitBamJobs())
		each sub-alignment is paired with the VCF from the same contig/chromosome.
		run CountCovariates job (skeleton in ShortRead2AlignmentPipeline.py)
		TableRecalibration job
	Merge all sub-alignment into one
	add alignment file into db (input: parent alignment ID, alignment file path)
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
from Pegasus.DAX3 import *
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.VCFFile import VCFFile
#from AlignmentToCallPipeline import AlignmentToCallPipeline
#from AbstractVervetWorkflow import AbstractVervetWorkflow
from AbstractAlignmentAndVCFWorkflow import AbstractAlignmentAndVCFWorkflow

parentClass = AbstractAlignmentAndVCFWorkflow
class AlignmentReadBaseQualityRecalibrationWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(parentClass.partitionWorkflowOptionDict.copy())
	option_default_dict.update({
				('mask_genotype_method_id', 1, int):[None, '', 1, 'which genotype method is used to mask out polymorphic sites for recalibration'],\
							})
	option_default_dict[('intervalSize', 1, int)][0] = 10000000
	"""
	option_default_dict.update({
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						('inputDir', 0, ): [None, 'i', 1, 'folder containing vcf or vcf.gz files that will be used by CountCovariates to exclude variant sites', ],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'To filter alignments. None: whatever; 0: unfiltered sequences, 1: filtered sequences'],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						('intervalOverlapSize', 1, int): [3000, '', 1, 'extension of an interval on each side. overlap size is actually 2X of this number.\
							Its value should be bigger than maximum read length, to cover reads that are partially aligned within one interval.', ],\
						('intervalSize', 1, int): [5000000, '', 1, 'size for adjacent intervals from one contig/chromosome', ],\
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						('run_type', 1, int): [1, 'y', 1, '', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	"""
	
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
		#AlignmentToCallPipeline.__init__(self, **keywords)
		#self.inputDir = os.path.abspath(self.inputDir)
	
	def mapEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
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
		
		
		"""
		#2012.9.21 perhaps a downsampling job
		outputFname = os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, overlapFilenameSignature))
		outputFile = File(outputFname)
		selectAlignmentJob, bamIndexJob1 = self.addSelectAlignmentJob(executable=workflow.samtools, inputFile=bamF, \
				outputFile=outputFile, region=overlapInterval, parentJobLs=[topOutputDirJob] + parentJobLs, \
				extraDependentInputLs=[baiF], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=True)
		"""
		
		"""
		#2012.9.21 count covariates job is moved to map()
		recalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.csv'%(bamFnamePrefix, chromosome)))
		countCovariatesJob = self.addGATKCountCovariatesJob(workflow, executable=workflow.countCovariatesJava, \
								genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, inputFile=bamF, \
								VCFFile=VCFFile, interval=chromosome, outputFile=recalFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob]+parentJobLs, 
								extraDependentInputLs=[baiF, VCFFile.tbi_F], \
								transferOutput=False, \
								extraArguments=None, job_max_memory=4000)
		
		self.no_of_jobs += 1
		returnData.countCovariatesJob = countCovariatesJob
		returnData.jobDataLs.append(PassingData(jobLs=[countCovariatesJob], file=countCovariatesJob.recalFile, \
											fileList=[countCovariatesJob.recalFile]))
		"""
		
		return returnData
	
	def mapEachInterval(self, workflow=None, alignmentData=None, intervalData=None,\
							VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None,\
							mapEachChromosomeData=None, transferOutput=False, **keywords):
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
		
		"""
		outputFname = os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, overlapFilenameSignature))
		outputFile = File(outputFname)
		selectAlignmentJob, bamIndexJob1 = self.addSelectAlignmentJob(executable=workflow.samtools, inputFile=bamF, \
				outputFile=outputFile, region=overlapInterval, parentJobLs=[topOutputDirJob] + parentJobLs, \
				extraDependentInputLs=[baiF], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=True)
		"""
		
		recalFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.csv'%(bamFnamePrefix, overlapFilenameSignature)))
		countCovariatesJob = self.addGATKCountCovariatesJob(workflow, executable=workflow.countCovariatesJava, \
								genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, inputFile=bamF, \
								VCFFile=VCFFile, interval=overlapInterval, outputFile=recalFile, \
								refFastaFList=passingData.refFastaFList, parentJobLs=[topOutputDirJob]+parentJobLs, 
								extraDependentInputLs=[baiF, VCFFile.tbi_F], \
								transferOutput=False, \
								extraArguments=None, job_max_memory=4000)
		"""
		countCovariatesJob = mapEachChromosomeData.countCovariatesJob
		"""
		
		recalBAMFile = File(os.path.join(topOutputDirJob.output, '%s_%s.recal_data.bam'%(bamFnamePrefix, overlapFilenameSignature)))
		tableRecalibrationJob, bamIndexJob2 = self.addGATKTableRecalibrationJob(workflow, executable=workflow.tableRecalibrationJava, \
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, inputFile=bamF, \
							recalFile=countCovariatesJob.recalFile, interval=overlapInterval, outputFile=recalBAMFile, \
							refFastaFList=passingData.refFastaFList, parentJobLs=[countCovariatesJob] + parentJobLs, \
							extraDependentInputLs=[baiF, VCFFile.tbi_F], transferOutput=False, \
							extraArguments=None, job_max_memory=3000, needBAMIndexJob=True)
		
		nonOverlapBamFile = File(os.path.join(topOutputDirJob.output, '%s_%s.bam'%(bamFnamePrefix, intervalFnameSignature)))
		
		selectAlignmentJob, bamIndexJob3 = self.addSelectAlignmentJob(executable=workflow.samtools, \
															inputFile=tableRecalibrationJob.output, \
				outputFile=nonOverlapBamFile, region=mpileupInterval, parentJobLs=[tableRecalibrationJob, bamIndexJob2], \
				extraDependentInputLs=[bamIndexJob2.output], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, needBAMIndexJob=False)	#the next mergeBamJob doesn't need bai files.
		
		passingData.AlignmentJobAndOutputLs.append([selectAlignmentJob, selectAlignmentJob.output])
		#add the sub-alignment to the alignment merge job
		self.no_of_jobs += 5
		return returnData
	

	def reduceAfterEachAlignment(self, workflow=None, passingData=None, transferOutput=False, dataDir=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		if workflow is None:
			workflow = self
		getattr(passingData, 'mapEach', [])
		AlignmentJobAndOutputLs = getattr(passingData, 'AlignmentJobAndOutputLs', [])
		bamFnamePrefix = passingData.bamFnamePrefix
		topOutputDirJob = passingData.topOutputDirJob
		individual_alignment = passingData.individual_alignment
		if len(AlignmentJobAndOutputLs)>0:	#2012.3.29	merge alignment output only when there is something to merge!
			mergedBamFile = File(os.path.join(topOutputDirJob.output, '%s_recal.bam'%(bamFnamePrefix)))
			alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
					outputBamFile=mergedBamFile, \
					samtools=workflow.samtools, java=workflow.java, \
					mergeSamFilesJava=workflow.mergeSamFilesJava, mergeSamFilesJar=workflow.mergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					mv=workflow.mv, parentJobLs=[topOutputDirJob], \
					transferOutput=False)
			self.no_of_jobs += 1
			#2012.9.19 add/copy the alignment file to db-affliated storage
			#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
			logFile = File(os.path.join(topOutputDirJob.output, '%s_2db.log'%(bamFnamePrefix)))
			alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
								inputFile=alignmentMergeJob.output, otherInputFileList=[],\
								parent_individual_alignment_id=individual_alignment.id, \
								mask_genotype_method_id=self.mask_genotype_method_id,\
								logFile=logFile, dataDir=dataDir, \
								parentJobLs=[alignmentMergeJob, bamIndexJob], \
								extraDependentInputLs=[bamIndexJob.output], \
								extraArguments=None, transferOutput=transferOutput, \
								job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
			self.no_of_jobs += 1
			returnData.jobDataLs.append(PassingData(jobLs=[alignment2DBJob], file=alignment2DBJob.logFile, \
											fileList=[alignment2DBJob.logFile]))
		return returnData

	
	
	def addGATKCountCovariatesJob(self, workflow, executable=None, genomeAnalysisTKJar=None, inputFile=None, \
								VCFFile=None, interval=None, outputFile=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, **keywords):
		"""
		2011-12-4
			inputFile is a bam file.
			outputFile is recalFile (csv file).
		
			java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAn alysisTK.jar -l INFO
				-R ../../../../NCBI/hs_genome.fasta -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam
				-T CountCovariates   -cov ReadGroupCovariate    -cov QualityScoreCovariate  
				-cov CycleCovariate    -cov DinucCovariate
				-recalFile  454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
				-B:mask,VCF 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.GATK.vcf
		"""
		#GATK job
		#MaxPermSize= min(35000, max(1024, job_max_memory*7/9))
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
		#javaMemRequirement = "-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*50/100, job_max_memory, \
		#																			MaxPermSize*50/100, MaxPermSize)
		refFastaFile = refFastaFList[0]
		extraArgumentList = [memRequirementData.memRequirementInStr, '-jar', genomeAnalysisTKJar, "-T", "CountCovariates",\
						"-cov ReadGroupCovariate", "-cov QualityScoreCovariate", \
						"-cov CycleCovariate", "-cov DinucCovariate", "-I", inputFile, "-R", refFastaFile,\
						"-L", interval, self.defaultGATKArguments, \
						"-recalFile", outputFile,\
						"-knownSites", VCFFile]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputFile, VCFFile] + refFastaFList)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[outputFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, \
						job_max_memory=memRequirementData.memRequirement, **keywords)
		
		job.recalFile = outputFile
		return job
	
	def addGATKTableRecalibrationJob(self, workflow, executable=None, genomeAnalysisTKJar=None, inputFile=None, \
								recalFile=None, interval=None, outputFile=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, needBAMIndexJob=False, **keywords):
		"""
		2012.7.27
			java -jar ~/script/vervet/bin/GenomeAnalysisTK-1.0.4705/GenomeAnalysisTK.jar -l INFO 
				-R ../../../../NCBI/hs_genome.fasta -I 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.bam 
				-T TableRecalibration  --out 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal.bam 
				-recalFile 454_vs_hg19.3eQTL.minPerBaseAS0.4.minMapQ125.score2.recal_data.csv
		"""
		#GATK job
		memRequirementData = self.getJVMMemRequirment(job_max_memory=job_max_memory)
		#MaxPermSize= min(35000, max(1024, job_max_memory*7/9))
		javaMemRequirement = memRequirementData.memRequirementInStr
		job_max_memory = memRequirementData.memRequirement
		
		#"-Xms%sm -Xmx%sm -XX:PermSize=%sm -XX:MaxPermSize=%sm"%(job_max_memory*50/100, job_max_memory, \
		#																			MaxPermSize*50/100, MaxPermSize)
		refFastaFile = refFastaFList[0]
		extraArgumentList = [javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "TableRecalibration",\
						"-I", inputFile, "-R", refFastaFile,\
						"-L", interval, self.defaultGATKArguments,\
						"-recalFile", recalFile, "--out", outputFile]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if extraDependentInputLs is None:
			extraDependentInputLs=[]
		extraDependentInputLs.extend([inputFile, recalFile] + refFastaFList)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[outputFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		if needBAMIndexJob:
			# add the index job on the bam file
			bamIndexJob = self.addBAMIndexJob(BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, \
						BuildBamIndexFilesJar=self.BuildBamIndexFilesJar, \
						inputBamF=job.output, parentJobLs=[job], \
						transferOutput=transferOutput, job_max_memory=job_max_memory)
		else:
			bamIndexJob = None
		return job, bamIndexJob
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		countCovariatesJava = Executable(namespace=namespace, name="countCovariatesJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		countCovariatesJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((countCovariatesJava, 1))
		
		tableRecalibrationJava = Executable(namespace=namespace, name="tableRecalibrationJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		tableRecalibrationJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((tableRecalibrationJava, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	

if __name__ == '__main__':
	main_class = AlignmentReadBaseQualityRecalibrationWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()