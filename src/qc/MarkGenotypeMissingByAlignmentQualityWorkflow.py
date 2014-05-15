#!/usr/bin/env python
"""
Examples:
	#
	%s  ...
	
	#2013.12.04  intervalSize means no of loci
	%s -I ~/NetworkData/vervet/db/genotype_file/method_225/
		-H -C 10 -j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		--genotypeMethodID 225
		--ref_ind_seq_id 3488 --intervalSize 80000
		-o dags/MarkGenotypeMissing/MarkGenotypeMissing_Method225_Ref3488.xml
		--db_user yh -z localhost
		#--notToUseDBToInferVCFNoOfLoci
Description:
	2013.12.04
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import Executable, File, PFN
from pymodule import ProcessOptions, PassingData
from pymodule.pegasus import yh_pegasus
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow
from vervet.src.db import VervetDB

parentClass = AbstractVervetWorkflow

class MarkGenotypeMissingByAlignmentQualityWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
					('genotypeMethodID', 1, int): [None, '', 1, 'from which genotype method the VCF files are from', ],\
					('minMediumCoverageThreshold', 1, int): [4, '', 1, 'minimum coverage for median or higher coverage individuals', ],\
					})
	
	#2012.9.25 no overlap and make the interval (no of loci) a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 70000
	option_default_dict[('max_walltime', 1, int)][0] = 1300	#under 23 hours
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		
		self.needSplitChrIntervalData = False
		parentClass.__init__(self, **keywords)
		self.needSplitChrIntervalData = False
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.06.14
			move topOutputDirJob from addAllJobs to here. 
		2012.9.17
		"""
		
		if workflow is None:
			workflow = self
		returnData = parentClass.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix, \
										passingData=passingData, transferOutput=transferOutput, **keywords)
		
		#add a stat merge job, merge of all 
		#outputFile = File(os.path.join(self.reduceStatDirJob.output, 'mergedGenotypeMissingStat.tsv.gz'))
		#self.mergeGenotypeMissingStatJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
		#							outputF=outputFile, \
		#							parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
		#							extraDependentInputLs=None, transferOutput=True)
		
		#add a stat merge job, how many missing genotypes (samples) per locus
		outputFile = File(os.path.join(self.reduceStatDirJob.output, 'missingPerLocusStat.tsv.gz'))
		self.mergeMissingPerLocusStatJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
									outputF=outputFile, \
									parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
									extraDependentInputLs=None, transferOutput=True)
		
		#add a stat merge job, count the missing genotypes by reason,
		# it's a reducer of reducer. that's why "--keyColumnLs 0 --valueColumnLs 1"
		
		outputFile = File(os.path.join(self.reduceStatDirJob.output, 'genotypeMissingReasonStat.tsv.gz'))
		self.reduceGenotypeMissingStatByReasonJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByChosenColumn, \
									outputF=outputFile, extraArgumentList=["--keyColumnLs 0 --valueColumnLs 1"],\
									parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
									extraDependentInputLs=None, transferOutput=True)
		
		#histogram for 
		outputFile = File( os.path.join(self.plotDirJob.output, 'genotypeMissingFractionPerLocusHistogram.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(executable=workflow.DrawHistogram, inputFileList=[self.mergeMissingPerLocusStatJob.output], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="occurrence_byFixedValue", whichColumnPlotLabel="missingFractionPerLocus", \
					xScaleLog=0, yScaleLog=1, \
					logCount=False, logY=False, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,legendType=1, \
					parentJobLs=[self.plotDirJob, self.mergeMissingPerLocusStatJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=8000)	#lots of input data,
		
		return returnData
	
	def mapEachInterval(self, workflow=None, VCFJobData=None, chromosome=None,intervalData=None,\
					mapEachChromosomeData=None, passingData=None, transferOutput=False, \
					**keywords):
		"""
		2012.9.22
			argument VCFJobData looks like PassingData(file=splitVCFFile, vcfFile=splitVCFFile, fileLs=[splitVCFFile], \
																		job=splitVCFJob, jobLs=[splitVCFJob], tbi_F=None)
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		passingData.intervalFileBasenamePrefix
		passingData.splitVCFFile
		passingData.unitNumber
		"""
		passingData.fileBasenamePrefix
		## 2013.06.19 structures available from passingData, specific to the interval
		passingData.splitVCFFile = splitVCFFile
		#VCFJobData.file is same as passingData.splitVCFFile 
		passingData.unitNumber = unitNumber
		passingData.intervalFileBasenamePrefix = '%s_%s_splitVCF_u%s'%(chromosome, commonPrefix, unitNumber)
		passingData.noOfIndividuals = jobData.file.noOfIndividuals
		passingData.span = self.intervalSize + self.intervalOverlapSize*2 	#2013.06.19 for memory/walltime gauging
		"""
		
		
		# 2013.12.04 VCF combine
		#a combine VCF horizontally job
		#2011-9-22 union of all GATK intervals for one contig
		combineCallOutputFname = os.path.join(self.mapDirJob.output, '%s.vcf.gz'%passingData.intervalFileBasenamePrefix)
		combineCallOutputF = File(combineCallOutputFname)
		returnData.combineCallJob = self.addGATKCombineVariantsJob(executable=self.CombineVariantsJavaInReduce, \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					outputFile=combineCallOutputF, \
					genotypeMergeOptions='UNSORTED', \
					parentJobLs=[self.mapDirJob], \
					extraArguments=None, extraArgumentList=None, extraDependentInputLs=None,\
					transferOutput=False, job_max_memory=10000, walltime=180)
		
		#add a job to calculate missing fraction per locus in medium and higher coverage individuals
		outputFile = File(os.path.join(self.reduceStatDirJob.output, '%s_missingFractionPerLocus.tsv.gz'%(passingData.intervalFileBasenamePrefix)))
		missingFractionPerLocusJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
									outputF=outputFile, extraArgumentList=["--keyColumnLs 1", "--valueColumnLs 5", "--operatorType 3",\
												"--fixedValueDenominator %s"%(self.noOfMediumCoverageIndividuals)],\
									parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
									extraDependentInputLs=None, transferOutput=False)
		
		returnData.markGenotypeMissingJobLs = []
		#for each sample from the genotype method
		for alignmentData in self.alignmentDataLs:
			if not alignmentData.newAlignment:
				continue
			individual_alignment = alignmentData.alignment
			sampleID = individual_alignment.read_group
			
			# select that sample into VCF
			outputFile = File(os.path.join(self.mapDirJob.output, '%s_sample_%s.vcf.gz'%(passingData.intervalFileBasenamePrefix, sampleID)))
			selectSampleJob = self.addSelectVariantsJob(\
					inputF=passingData.splitVCFFile, outputF=outputFile, \
					interval=None,\
					refFastaFList=self.registerReferenceData.refFastaFList, sampleIDKeepFile=None, \
					snpIDKeepFile=None, sampleIDExcludeFile=None, \
					parentJobLs=[self.mapDirJob] + VCFJobData.jobLs, extraDependentInputLs=None, transferOutput=False, \
					extraArguments="--sample_name %s"%(alignmentData.alignment.read_group), \
					extraArgumentList=None, job_max_memory=2000, walltime=None)
			
			# mask genotype missing
			# given single-sample VCF file, alignment file, median-depth
			# output two files, one is a VCF , one is a missing genotype stat file
			genotypeMissingStatFile = File(os.path.join(self.mapDirJob.output, '%s_sample_%s_genotypeMissingStat.tsv.gz'%(passingData.intervalFileBasenamePrefix,\
																	sampleID)))
			outputFile = File(os.path.join(self.mapDirJob.output, '%s_sample_%s_markedMissing.vcf.gz'%(passingData.intervalFileBasenamePrefix,\
																	sampleID)))
			alignmentFile = alignmentData.newAlignment.file
			if individual_alignment.individual_sequence.individual.target_coverage>=self.minMediumCoverageThreshold:
				alignmentDepthFold=2
			else:
				alignmentDepthFold=10
			markGenotypeMissingJob = self.addGenericJob(executable=self.MarkGenotypeMissingByAlignmentQuality, \
							inputFile=selectSampleJob.output, inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", \
					parentJob=None, parentJobLs=[self.mapDirJob, selectSampleJob], \
					extraDependentInputLs=alignmentData.newAlignment.fileLs, extraOutputLs=[genotypeMissingStatFile], \
					extraArgumentList=["--alignmentFilename", alignmentFile, "--missingStatFname", genotypeMissingStatFile, \
									"--alignmentMedianDepth %s"%(alignmentData.newAlignment.median_depth), \
									"--alignmentDepthFold %s"%(alignmentDepthFold), \
									"--minMapQGoodRead 2", "--minFractionOfGoodRead 0.9",\
									"--sampleID %s"%(sampleID)], \
					transferOutput=False, sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, objectWithDBGenomeArguments=None,\
					no_of_cpus=None, job_max_memory=2000, walltime=180, \
					max_walltime=None)
			markGenotypeMissingJob.genotypeMissingStatFile = genotypeMissingStatFile
			returnData.markGenotypeMissingJobLs.append(markGenotypeMissingJob)
			
			if individual_alignment.individual_sequence.individual.target_coverage>=self.minMediumCoverageThreshold:
				#missing fraction only from medium or high coverage individuals
				self.addInputToStatMergeJob(statMergeJob=missingFractionPerLocusJob, parentJobLs=[markGenotypeMissingJob],\
									inputF=markGenotypeMissingJob.genotypeMissingStatFile)
			
			#add sample to VCF union job 
			self.addInputToStatMergeJob(statMergeJob=returnData.combineCallJob, inputF=markGenotypeMissingJob.output,\
								parentJobLs=[markGenotypeMissingJob], \
								extraDependentInputLs=[],\
								inputArgumentOption="--variant")
			"""
			# combine VCF horizontally
			self.hierarchicalCombineIntervalJobIntoOneChromosome(unionJob=returnData.combineCallJob.callUnionJob, \
					refFastaFList=refFastaFList, fileBasenamePrefix="subUnion",\
					intervalJobLs=chromosome2jobData[chromosome].genotypingJobLs, \
					maxNoOfIntervalJobsInOneUnion=250,\
					outputDirJob=combineChromosomeJobData.callOutputDirJob,
					transferOutput=False, job_max_memory=None, walltime=120, needBGzipAndTabixJob=False)
			"""
		
		#add job to remove loci with too many missing values
		maxMissingFraction = 0.5
		outputFile = File(os.path.join(self.mapDirJob.output, '%s_lociWithMissingFractionAbove%sRemoved.vcf.gz'%\
									(passingData.intervalFileBasenamePrefix, maxMissingFraction)))
		
		removeHighMissingLocusJob = self.addGenericJob(executable=self.FilterLocusBasedOnLocusStatFile, \
					inputArgumentOption="-i", inputFile=returnData.combineCallJob.output, \
					outputArgumentOption="-o", outputFile=outputFile, \
					parentJobLs=[self.mapDirJob, returnData.combineCallJob, missingFractionPerLocusJob], \
					extraDependentInputLs=[missingFractionPerLocusJob.output], extraOutputLs=None, \
					extraArgumentList=["--statFname", missingFractionPerLocusJob.output, \
									"--maxValue %s"%(maxMissingFraction), "--runType 1"], \
					transferOutput=False, sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, objectWithDBGenomeArguments=None,\
					no_of_cpus=None, job_max_memory=2000, walltime=180, \
					max_walltime=None)
		removeHighMissingLocusJob.intervalData = intervalData
		
		returnData.missingFractionPerLocusJob = missingFractionPerLocusJob
		returnData.removeHighMissingLocusJob = removeHighMissingLocusJob
		
		return returnData
	
	
	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		2013.12.04
		
		"""
		#check how often one locus has missing genotype
		outputFile = File(os.path.join(self.reduceStatDirJob.output, '%s_missingPerLocusStat.tsv.gz'%(passingData.fileBasenamePrefix)))
		mergeMissingFractionFromOneVCFJob = self.addStatMergeJob(statMergeProgram=self.mergeSameHeaderTablesIntoOne, \
									outputF=outputFile, extraArgumentList=None,\
									parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
									extraDependentInputLs=None, transferOutput=False)
		self.addInputToStatMergeJob(statMergeJob=self.mergeMissingPerLocusStatJob, \
						parentJobLs=[mergeMissingFractionFromOneVCFJob])
		
		#add a stat merge job, count the missing genotypes by reason
		outputFile = File(os.path.join(self.reduceStatDirJob.output, '%s_genotypeMissingReasonStat.tsv.gz'%(passingData.fileBasenamePrefix)))
		reduceGenotypeMissingStatByReasonOneVCFJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByChosenColumn, \
									outputF=outputFile, extraArgumentList=["--keyColumnLs 6 --valueColumnLs 5"],\
									parentJobLs=[self.reduceStatDirJob],extraOutputLs=None, \
									extraDependentInputLs=None, transferOutput=False)
		self.addInputToStatMergeJob(statMergeJob=self.reduceGenotypeMissingStatByReasonJob, \
						parentJobLs=[reduceGenotypeMissingStatByReasonOneVCFJob])
		
		
		#2013.12.04 get the VCF union
		combineCallOutputFname = os.path.join(self.reduceEachInputDirJob.folder, '%s.vcf'%chromosome)
		combineCallOutputF = File(combineCallOutputFname)
		combineCallJob = self.addGATKCombineVariantsJob(executable=self.CombineVariantsJavaInReduce, \
							refFastaFList=self.registerReferenceData.refFastaFList, outputFile=combineCallOutputF, \
							genotypeMergeOptions='UNSORTED', \
			parentJobLs=[self.reduceEachInputDirJob], \
			extraArguments=None, extraArgumentList=['--assumeIdenticalSamples'], extraDependentInputLs=None,\
			transferOutput=False, job_max_memory=8000, walltime=180)
		
		gatkGzipUnionOutputF = File("%s.gz"%combineCallOutputFname)
		bgzip_tabix_combineCallOutputF_job = self.addBGZIP_tabix_Job(\
				parentJob=combineCallJob, inputF=combineCallJob.output, outputF=gatkGzipUnionOutputF, \
				transferOutput=True)
		
		intervalJobLs = []
		for mapEachIntervalReturnData in mapEachIntervalDataLs:
			intervalJobLs.append(mapEachIntervalReturnData.removeHighMissingLocusJob)
			#add sample to VCF union job 
			self.addInputToStatMergeJob(statMergeJob=combineCallJob, inputF=mapEachIntervalReturnData.removeHighMissingLocusJob.output,\
								parentJobLs=[mapEachIntervalReturnData.removeHighMissingLocusJob], \
								extraDependentInputLs=[],\
								inputArgumentOption="--variant")
			self.addInputToStatMergeJob(statMergeJob=mergeMissingFractionFromOneVCFJob, \
						parentJobLs=[mapEachIntervalReturnData.missingFractionPerLocusJob])
			for markGenotypeMissingJob in mapEachIntervalReturnData.markGenotypeMissingJobLs:
				
				self.addInputToStatMergeJob(statMergeJob=reduceGenotypeMissingStatByReasonOneVCFJob, \
						inputF=markGenotypeMissingJob.genotypeMissingStatFile, parentJobLs=[markGenotypeMissingJob])
				#do not keep this all-stat file
				#self.addInputToStatMergeJob(statMergeJob=self.mergeGenotypeMissingStatJob, \
				#		inputF=markGenotypeMissingJob.genotypeMissingStatFile, parentJobLs=[markGenotypeMissingJob])
				
		
		
		"""
		#2014.1.8 one-layer merge is enough. hierarchical later. 
		#combine intervals hierarchically, to avoid open file number limit on cluster
		self.hierarchicalCombineIntervalJobIntoOneChromosome(unionJob=combineCallJob, \
					refFastaFList=self.registerReferenceData.refFastaFList, fileBasenamePrefix="subUnion",\
					intervalJobLs=intervalJobLs, \
					maxNoOfIntervalJobsInOneUnion = 250,\
					outputDirJob=self.reduceEachInputDirJob,
					transferOutput=False, job_max_memory=None, walltime=120, needBGzipAndTabixJob=False)
		"""
	
	def setup_run(self, **keywords):
		"""
		2013.07.08
			
		"""
		self.needSplitChrIntervalData = False
		pdata = parentClass.setup_run(self)
		
		self.genotypeMethod = VervetDB.GenotypeMethod.get(self.genotypeMethodID)
		
		#filter out alignment whose sequence coverage is outside the min,max range
		individual_alignment_ls = []
		for individual_alignment in self.genotypeMethod.individual_alignment_ls:
			if self.sequence_max_coverage is not None and individual_alignment.individual_sequence.coverage>self.sequence_max_coverage:
				continue
			if self.sequence_min_coverage is not None and individual_alignment.individual_sequence.coverage<self.sequence_min_coverage:
				continue
			individual_alignment_ls.append(individual_alignment)
		
		#register every alignment file so that it could be re-used
		self.alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow=self, alignmentLs=individual_alignment_ls, \
													data_dir=self.data_dir, checkFileExistence=False)
		#get corresponding new alignment as the input VCF is lifted over from old reference
		real_counter = 0
		self.noOfMediumCoverageIndividuals = 0
		for alignmentData in self.alignmentDataLs:
			newAlignment = VervetDB.IndividualAlignment.query.filter_by(ind_seq_id=alignmentData.alignment.ind_seq_id).\
					filter_by(ref_ind_seq_id=self.ref_ind_seq_id).\
					filter_by(alignment_method_id=6).first()
			"""
			self.db_vervet.checkIndividualAlignment(individual_sequence_id=alignmentData.alignment.ind_seq_id, \
					ref_individual_sequence_id=self.ref_ind_seq_id, \
					alignment_method_id=6, \
					data_dir=self.data_dir, \
					local_realigned=None, reduce_reads=None)
			"""
			if newAlignment:
				bamFile = self.registerOneInputFile(inputFname=os.path.join(self.data_dir, newAlignment.path), folderName="newAlignment", \
								useAbsolutePathAsPegasusFileName=False,\
								pegasusFileName=None, checkFileExistence=True)
				baiFile = self.registerOneInputFile(inputFname=os.path.join(self.data_dir, "%s.bai"%newAlignment.path), folderName="newAlignment", \
								useAbsolutePathAsPegasusFileName=False,\
								pegasusFileName=None, checkFileExistence=True)
				newAlignment.file = bamFile
				newAlignment.fileLs = [bamFile, baiFile]
				newAlignment.job = None
				newAlignment.jobLs = []
				alignmentData.newAlignment = newAlignment
				real_counter += 1
				if newAlignment.individual_sequence.individual.target_coverage>=self.minMediumCoverageThreshold:
					self.noOfMediumCoverageIndividuals +=1
		sys.stderr.write( '%s alignment(s) have new alignment found.\n'%(real_counter))
		
		return self

	def registerExecutables(self, workflow=None):
		"""
		"""
		if not workflow:
			workflow = self
		parentClass.registerExecutables(self, workflow=workflow)
		
		#2013.07.12
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'polymorphism/qc/mapper/MarkGenotypeMissingByAlignmentQuality.py'), \
									name='MarkGenotypeMissingByAlignmentQuality', \
									clusterSizeMultipler=0.5)

if __name__ == '__main__':
	main_class = MarkGenotypeMissingByAlignmentQualityWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()