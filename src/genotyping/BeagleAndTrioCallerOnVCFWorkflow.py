#!/usr/bin/env python
"""
Examples:
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "--noOfCallingJobsPerNode 3" controls clustering of calling programs.
	# "--clusters_size 30" controls clustering for other programs., "-S 447" dictates monkeys from VRC
	%s -u yh -a 524 -s 2 -z localhost -o dags/AlignmentToCall/AlignmentToTrioCallPipeline_VRC_top7559Contigs.xml -j hcondor -l hcondor 
		-N7559 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447 -e /u/home/eeskin/polyacti/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --needSSHDBTunnel
	
	# 2013.06.19
	%s -I /Network/Data/vervet/db/genotype_file/method_36 -u yh -a 524  -z localhost
		--pedigreeKinshipFilePath ...
		-o  dags/RefineCall/BeagleTrioCallerOnMethod36Contig96_100.xml -j hcondor -l hcondor
		--noOfCallingJobsPerNode 1 --clusters_size 1
		--data_dir /Network/Data/vervet/db/ --local_data_dir /Network/Data/vervet/db/
		--minContigID 96 --maxContigID 100
	
	# 2012.8.15 run TrioCaller on method 14 samtools calls, contig ID from 96 to 100 (--minContigID 96 --maxContigID 100)
	# sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2), onlyKeepBiAllelicSNP (--onlyKeepBiAllelicSNP) 
	# calling job clusters size=1, others =1 (--noOfCallingJobsPerNode 1 --clusters_size 1)
	# add -Y (not guess #loci from 1st number in filename) if the input VCF is not db affiliated
	# 3000 (-Z 3000) sites per unit, 500 overlapping between units (-U 500)
	# add "--treatEveryOneIndependent" if you want to treat everyone independent (no mendelian constraints from TrioCaller)
	%s --run_type 2 -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-u yh -a 524  -z localhost -o  dags/AlignmentToCall/TrioCallerOnMethod14Contig96_100.xml
		-j hcondor -l hcondor
		--noOfCallingJobsPerNode 1 --clusters_size 1 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--minContigID 96 --maxContigID 100 --sequence_filtered 1 --alignment_method_id 2  --onlyKeepBiAllelicSNP --needSSHDBTunnel
		# -Y -Z 3000 -U 500 --treatEveryOneIndependent
	
	#2013.6.25 test on ycondor
	%s --ref_ind_seq_id 3280 -I /Network/Data/vervet/db/genotype_file/method_101
		-u yh -z uclaOffice
		--pedigreeKinshipFilePath /Network/Data/vervet/Kinx2Aug2012.txt.gz 
		-o dags/RefineCall/BeagleTrioCallerOnMethod101Scaffold96_100.xml
		-j ycondor -l ycondor --noOfCallingJobsPerNode 1 --clusters_size 1
		--data_dir /Network/Data/vervet/db/ --local_data_dir /Network/Data/vervet/db/
		--minContigID 96 --maxContigID 100 --intervalOverlapSize 200 --intervalSize 2000
	
Description:
	2013.06.14
		a program which generates a 3-stage pedigree genotype-calling pegasus workflow dag (xml file).
			#. Beagle phasing/imputing on high-coverage members
			#. Beagle phasing/imputing on all members with step 1's output as reference panel
			#. TrioCaller on step 2's output.
		
		Sample IDs in --pedigreeKinshipFilePath are individual.ucla_id
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from pymodule import VCFFile
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow
from vervet.src.genotyping.AlignmentToTrioCallPipeline import AlignmentToTrioCallPipeline

parentClass=AlignmentToTrioCallPipeline


class BeagleAndTrioCallerOnVCFWorkflow(AbstractVervetWorkflow, parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetWorkflow.option_default_dict)
	option_default_dict.update(parentClass.commonCallPipelineOptionDict)
	#option_default_dict.pop(('ind_aln_id_ls', 0, ))
	#option_default_dict.pop(('ind_seq_id_ls', 0, ))
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						("minCoverageForRefPanel", 0, int): [8, '', 1, 'minimum coverage for a pedigree member to be in reference panel'],\
						("treatEveryOneIndependent", 0, int): [0, '', 0, 'toggle this to treat everyone in the pedigree independent (thus no replicates)'],\
						("run_type", 1, int): [1, '', 1, '1: not used right now'],\
						('pedigreeKinshipFilePath', 1, ): [None, '', 1, 'file that contains pairwise kinship between individuals (ID: ucla_id/code).\n\
	no header. coma-delimited 3-column file: individual1, individual2, kinsihp\n\
	The sampling will try to avoid sampling close pairs, kinship(i,j)<=maxPairwiseKinship'],\
						('maxPairwiseKinship', 0, float): [0.25, '', 1, 'maximum pairwise kinship allowed among the reference panel'],\
						("minProbForValidCall", 1, float): [0.6, '', 1, 'minimum probability for a call to be regarded as valid. DEPRECATED'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	#2012.8.26 run_type 2 lower the default. In the addTrioCallerJobsONVCFFiles(), these two numbers refer to the number of sites, not base pairs. 
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 1000
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
	def __init__(self,  **keywords):
		"""
		2013.05
		"""
		parentClass.__init__(self, **keywords)
		self.extractRefPanelSampleIDJob = None	#this job only needs to be run once during mapEachInterval(). 
		#can't put it inside preReduce because it requires a VCF with replicates among them
		
		#self.trioCallerPath =  self.insertHomePath(self.trioCallerPath, self.home_path)
		self.needSplitChrIntervalData = False	#2013.06.21 turn this off before setup_run() to not construct chr2IntervalDataLs
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.05.01
			1. a job that outputs the pedigree from db, with members from the VCF file. used by various filter programs and TrioCaller
			2. a job that extracts the high-coverage individuals from the VCF file
			3. figure out the existence of beagle unrelated cohort, trio cohort, pair/duo cohort for high-coverage group and all individuals
				need the pedigree graph, a VCF file => all sample IDs and only high-coverage individuals
				
		"""
		returnData = AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix, \
													passingData=passingData, transferOutput=transferOutput, **keywords)
		
		self.statDirJob = self.addMkDirJob(outputDir="%sStat"%(outputDirPrefix))
		self.highCoveragePanelDirJob = self.addMkDirJob(outputDir="%sHighCoveragePanel"%(outputDirPrefix))
		self.auxDirJob = self.addMkDirJob(outputDir="%sAuxilliary"%(outputDirPrefix))
		
		self.beagleReduceDirJob = self.addMkDirJob(outputDir="%sReduceBeagle"%(outputDirPrefix))
		# self.reduceOutputDirJob would contain non-replicate VCF files
		#this folder would store all the reduced VCF files with replicates among samles. 
		self.replicateVCFDirJob = self.addMkDirJob(outputDir="%sReplicateVCF"%(outputDirPrefix))
		
		self.pedigreeKinshipFile = self.registerOneInputFile(inputFname=self.pedigreeKinshipFilePath, \
										folderName='aux')
		
		inputFileBasenamePrefix = utils.getFileBasenamePrefixFromPath(self.firstVCFJobData.file.name)
		# output pedigree to get pedigree file (for TrioCaller etc. that requires pedigree to be split into trios/duos) and sampleID2FamilyCountF 
		#		(for ReplicateVCFGenotypeColumns job, setting TrioCaller up)
		pedigreeFileFormat = 2
		pedFile = File(os.path.join(self.auxDirJob.output, 'pedigree.replicates.%s.format%s.txt'%\
								(inputFileBasenamePrefix, pedigreeFileFormat)))
		sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'pedigree.replicates.sampleID2FamilyCount.%s.format%s.txt'%\
												(inputFileBasenamePrefix, pedigreeFileFormat)))
		self.outputReplicatePedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=self.firstVCFJobData.file, outputFile=pedFile, \
				sampleID2FamilyCountF=sampleID2FamilyCountF,\
				polymuttDatFile = None,\
				outputFileFormat=pedigreeFileFormat, \
				replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=self.firstVCFJobData.jobLs + [self.auxDirJob], \
				extraDependentInputLs=None, transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		#output pedigree, with no replicating certain individuals, no trio/duo splitting
		pedigreeFileFormat = 4
		pedFile = File(os.path.join(self.auxDirJob.output, 'pedigree.%s.format%s.txt'%\
								(inputFileBasenamePrefix, pedigreeFileFormat)))
		#sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'pedigree.sampleID2FamilyCount.%s.format%s.txt'%\
		#						(inputFileBasenamePrefix, pedigreeFileFormat)))
		self.outputPedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=self.firstVCFJobData.file, outputFile=pedFile, \
				sampleID2FamilyCountF=None,\
				polymuttDatFile = None,\
				outputFileFormat=pedigreeFileFormat, \
				replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=self.firstVCFJobData.jobLs + [self.auxDirJob], \
				extraDependentInputLs=None, transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		#ExtractSamplesFromVCF samples with coverage >=min_coverage
		# the input VCF does not contain replicates.
		outputFile = File(os.path.join(self.auxDirJob.output, '%s.minCoverage%s.sampleIDList.tsv'%\
									(inputFileBasenamePrefix, self.minCoverageForRefPanel)))
		extractRefPanelSampleIDJob = self.addExtractSampleIDJob(inputFile=self.firstVCFJobData.file, \
							outputFile=outputFile,\
							min_coverage=self.minCoverageForRefPanel, outputFormat=3,\
							returnData=returnData,\
							transferOutput=True, \
							parentJobLs=[self.firstVCFJobData.jobLs, self.auxDirJob])
		self.extractRefPanelSampleIDJob = extractRefPanelSampleIDJob
		
		
		# GATK SelectVariants: select High-coverage individuals out into a new VCF
		#	selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#	samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		#		?can it deal with Platypus output, which does not have AC/AF/DP?
		# selectHighCoverageSampleJob is needed here because a VCF file of high-coverage members is needed
		# 	for outputPedigreeOfHghCoverageSamplesJob
		#
		highCoverageSampleVCF = File(os.path.join(self.auxDirJob.output, '%s.minCoverage%s.vcf'%\
												(inputFileBasenamePrefix, self.minCoverageForRefPanel)))
		selectHighCoverageSampleJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=self.firstVCFJobData.file, \
				outputF=highCoverageSampleVCF, \
				refFastaFList=self.registerReferenceData.refFastaFList, \
				sampleIDKeepFile=self.extractRefPanelSampleIDJob.output,\
				parentJobLs=[self.auxDirJob, self.extractRefPanelSampleIDJob]+self.firstVCFJobData.jobLs, \
				extraDependentInputLs=[self.firstVCFJobData.tbi_F], transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000)
		
		# output a plink pedigree that contains these HC members only
		# output pedigree to get pedigree file (for GATK, TrioCaller, own programs) and sampleID2FamilyCountF (for ReplicateVCFGenotypeColumns job)
		# find a way to cache this job (used for same set of samples, but different chromosome intervals)
		pedigreeFileFormat = 4
		pedFile = File(os.path.join(self.auxDirJob.output, 'pedigree.minCoverage%s.%s.format%s.txt'%\
								(self.minCoverageForRefPanel, inputFileBasenamePrefix, pedigreeFileFormat)))
		#sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'pedigree.minCoverage%s.sampleID2FamilyCount.%s.format%s.txt'%\
		#									(self.minCoverageForRefPanel, inputFileBasenamePrefix, pedigreeFileFormat)))
		self.outputPedigreeOfHghCoverageSamplesJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=selectHighCoverageSampleJob.output, outputFile=pedFile, \
				sampleID2FamilyCountF=None,\
				polymuttDatFile = None,\
				outputFileFormat=pedigreeFileFormat, replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=[self.auxDirJob, selectHighCoverageSampleJob], \
				extraDependentInputLs=[], transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		#a job that outputs alignment coverage (alignment.read_group, median_depth)
		alignmentDepthFile = File(os.path.join(self.auxDirJob.folder, '%s.alignmentDepth.tsv'%(inputFileBasenamePrefix)))
		self.outputAlignmentDepthJob = self.addOutputVCFAlignmentDepthRangeJob(executable=self.OutputVCFAlignmentDepthRange, \
						inputFile=self.firstVCFJobData.file, \
						ref_ind_seq_id=self.ref_ind_seq_id, depthFoldChange=None, minGQ=None,\
						outputFile=alignmentDepthFile, outputFileFormat=1,\
						extraArgumentList=None,\
						parentJobLs=[self.auxDirJob]+self.firstVCFJobData.jobLs, \
						extraDependentInputLs=None, transferOutput=True, \
						job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		
		#a SelectDistantMembersFromGenotypeFile.py job to generate a ref panel for 2nd-round beagle
		# need the pedigree file
		# produces a list of samples
		phasedRefPanelSampleListFile = File(os.path.join(self.auxDirJob.folder, '%s.RefPanel.sampleList.maxPairwiseKinship%s.tsv'%\
														(inputFileBasenamePrefix, self.maxPairwiseKinship)))
		self.selectDistantMembersFromGenotypeFileJob = self.addGenericJob(executable=self.SelectDistantMembersFromGenotypeFile, \
						inputFile=selectHighCoverageSampleJob.output,
						outputFile=phasedRefPanelSampleListFile, outputArgumentOption="-o", \
						extraDependentInputLs=[self.pedigreeKinshipFile], \
						extraOutputLs=None, transferOutput=False, frontArgumentList=None, \
						extraArguments=None, \
						extraArgumentList=["--maxPairwiseKinship %s"%(self.maxPairwiseKinship), "--sampleSize 90", \
							"--pedigreeKinshipFile", self.pedigreeKinshipFile, \
							"--replicateIndividualTag", self.replicateIndividualTag,\
							"--individualAlignmentCoverageFname", self.outputAlignmentDepthJob.output, \
							"--pedigreeFname", self.outputPedigreeJob.output], \
						parentJobLs=[selectHighCoverageSampleJob, self.outputAlignmentDepthJob,  self.outputPedigreeJob,\
									self.auxDirJob],\
						no_of_cpus=None, job_max_memory = 4000, walltime= 120)
		
		"""
		
		#analyze the pedigree graph to figure out singletons, trios, duos
		self.alignmentLs = self.db.getAlignmentsFromVCFFile(inputFname=yh_pegasus.getAbsPathOutOfFile(self.firstVCFJobData.file))
		#2013.06.14 approach below does not work because pedigree of extracting-high-coverage + replication is different from that of replication + extracting-high-coverage (=reality).
		# some replicates might end up as singletons in the latter, while not so in the former.
		#
		self.highCoverageAlignmentLs = self.db.filterAlignments(alignmentLs=self.alignmentLs, min_coverage=self.minCoverageForRefPanel, \
			max_coverage=None, individual_site_id=None, \
			sequence_filtered=None, individual_site_id_set=None, \
			mask_genotype_method_id=None, parent_individual_alignment_id=None,\
			country_id_set=None, tax_id_set=None, excludeContaminant=False, excludeTissueIDSet=None,\
			local_realigned=None, reduce_reads=None, report=False)
		
		"""
		
		#a stat merge job (keeping track of how many mendel error sites were filtered)
		filterByRemoveMendelErrorSiteStatMergeFile = File(os.path.join(self.statDirJob.folder, 'filterByRemoveMendelErrorSiteStatMerge.tsv'))
		self.filterByRemoveMendelErrorSiteStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByRemoveMendelErrorSiteStatMergeFile, \
								transferOutput=False, parentJobLs=[self.statDirJob],\
								extraArguments="--keyColumnLs 1 --valueColumnLs 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
		returnData.jobDataLs.append(PassingData(jobLs=[self.filterByRemoveMendelErrorSiteStatMergeJob], \
											fileLs=[self.filterByRemoveMendelErrorSiteStatMergeJob.output]))
		#concordance stat reduce jobs
		#reduce the replicate concordance results from before TrioCaller (after beagle phasing)
		#
		"""
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.allSites.tsv'))
		reduceBeaglePhaseReplicateConcordanceJob_AllSites = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 2,3', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.homo.tsv'))
		reduceBeaglePhaseReplicateConcordanceJob_HomoOnly = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 5,6', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.tsv'))
		concatenateTwoBeaglePhaseConcordanceResultJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByMergeColumnsWithSameKey, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 2,3,4', transferOutput=False)
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoBeaglePhaseConcordanceResultJob, \
							parentJobLs=[reduceBeaglePhaseReplicateConcordanceJob_AllSites])
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoBeaglePhaseConcordanceResultJob, \
							parentJobLs=[reduceBeaglePhaseReplicateConcordanceJob_HomoOnly])
		returnData.jobDataLs.append(PassingData(jobLs=[concatenateTwoBeaglePhaseConcordanceResultJob], \
											fileLs=[concatenateTwoBeaglePhaseConcordanceResultJob.output]))
		#pass to self, as they will be used in reduceEachVCF()
		self.reduceBeaglePhaseReplicateConcordanceJob_AllSites = reduceBeaglePhaseReplicateConcordanceJob_AllSites
		self.reduceBeaglePhaseReplicateConcordanceJob_HomoOnly = reduceBeaglePhaseReplicateConcordanceJob_HomoOnly
		"""
		
		#reduce replicate concordance results from after-TrioCaller VCFs 
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.allSites.tsv'))
		reduceTrioCallerReplicateConcordanceJob_AllSites = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 2,3', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.homo.tsv'))
		reduceTrioCallerReplicateConcordanceJob_HomoOnly = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 5,6', transferOutput=False)
		
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.tsv'))
		concatenateTwoTrioCallerConcordanceResultJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByMergeColumnsWithSameKey, \
							outputF=outputFile, \
							extraArguments='--keyColumnLs 0,1 --valueColumnLs 2,3,4', transferOutput=False)
		
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoTrioCallerConcordanceResultJob, \
							parentJobLs=[reduceTrioCallerReplicateConcordanceJob_AllSites])
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoTrioCallerConcordanceResultJob, \
							parentJobLs=[reduceTrioCallerReplicateConcordanceJob_HomoOnly])
		returnData.jobDataLs.append(PassingData(jobLs=[concatenateTwoTrioCallerConcordanceResultJob], \
											fileLs=[concatenateTwoTrioCallerConcordanceResultJob.output]))
		#pass to self, as they will be used in reduceEachVCF()
		self.reduceTrioCallerReplicateConcordanceJob_AllSites = reduceTrioCallerReplicateConcordanceJob_AllSites
		self.reduceTrioCallerReplicateConcordanceJob_HomoOnly = reduceTrioCallerReplicateConcordanceJob_HomoOnly
		
		return returnData
	
	def mapEachInterval(self, workflow=None, intervalData=None, chromosome=None, \
					VCFJobData=None, passingData=None, 
					mapEachChromosomeData=None, transferOutput=False, \
					**keywords):
		"""
		2013.04.30
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		intervalFileBasenamePrefix = passingData.intervalFileBasenamePrefix
		span = passingData.span
		noOfIndividuals= passingData.noOfIndividuals
		SNPVCFFile = VCFJobData.file
		SNPVCFJobLs = VCFJobData.jobLs
		"""
		### 2013.06.19 intervalData does not exsit for input that is entirely composed of VCF files (SplitVCFFile job does not return intervals)
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFileBasenameSignature = intervalData.intervalFileBasenameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFileBasenameSignature = intervalData.overlapIntervalFnameSignature
		span = intervalData.span
		"""
		if chromosome is None:
			chromosome = getattr(passingData, 'chromosome', None)
		
		#noOfIndividuals
		realInputVolume = noOfIndividuals * span
		baseInputVolume = 600*2000	#600 individuals at 2000 sites
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=10000).value
		
		#splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		
		#### Part 1 generate high-quality reference panel through Beagle on high-coverage individuals
		# extractRefPanelSampleIDJob outputs sample IDs with replicate tags
		# select the high-coverage members
		outputVCF = File(os.path.join(self.highCoveragePanelDirJob.output, \
									'%s.minCoverageForRefPanel%s.vcf'%(intervalFileBasenamePrefix, self.minCoverageForRefPanel)))
		#selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		selectHighCoverageSampleJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=VCFJobData.file, outputF=outputVCF, \
				refFastaFList=self.registerReferenceData.refFastaFList, sampleIDKeepFile=self.extractRefPanelSampleIDJob.output,\
				parentJobLs=[self.highCoveragePanelDirJob, self.extractRefPanelSampleIDJob] + VCFJobData.jobLs, \
				extraDependentInputLs=[VCFJobData.tbi_F], transferOutput=False, \
				extraArguments=None, \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		# run Beagle 
		outputFnamePrefix = os.path.join(self.highCoveragePanelDirJob.folder, "%s.minCoverage%s.beagled"%\
								(intervalFileBasenamePrefix, self.minCoverageForRefPanel))
		beagleOnHighCoverageJob = self.addBeagle4Job(executable=self.BeagleOnHCMOnkeys, \
								inputFile=selectHighCoverageSampleJob.output, refPanelFile=None,\
								pedFile = self.outputPedigreeOfHghCoverageSamplesJob.output,\
								outputFnamePrefix=outputFnamePrefix, \
								burninIterations=7, phaseIterations=10, \
								noOfSamplingHaplotypesPerSample=4,\
								parentJobLs=[self.highCoveragePanelDirJob, selectHighCoverageSampleJob, self.outputPedigreeOfHghCoverageSamplesJob], \
								transferOutput=False, \
								extraArguments=None, extraArgumentList=None,\
								extraOutputLs=None, extraDependentInputLs=None, \
								no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=13000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		#index .vcf.gz, output of beagle, without index, GATK can't work on gzipped vcf
		tabixIndexFile = File('%s.tbi'%(beagleOnHighCoverageJob.output.name))
		tabixOnHighCoverageVCFJob = self.addGenericJob(executable=self.tabix, \
						inputFile=beagleOnHighCoverageJob.output, inputArgumentOption="",\
						outputFile=None, outputArgumentOption="-o", \
						extraDependentInputLs=None, \
						extraOutputLs=[tabixIndexFile], transferOutput=False, frontArgumentList=["-p vcf"], \
						extraArguments=None, \
						extraArgumentList=[], \
						parentJobLs=[beagleOnHighCoverageJob, self.highCoveragePanelDirJob],\
						no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=2000, maxJobPropertyValue=5000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=600).value)
		
		# select the high-coverage members
		outputVCF = File(os.path.join(self.highCoveragePanelDirJob.output, \
									'%s.minCoverage%s.maxPairwiseKinship%s.refPanel.beagled.vcf'%\
									(intervalFileBasenamePrefix, self.minCoverageForRefPanel, self.maxPairwiseKinship)))
		#selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		selectDistantMembersVariantsJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=beagleOnHighCoverageJob.output, outputF=outputVCF, \
				refFastaFList=self.registerReferenceData.refFastaFList, \
				sampleIDKeepFile=self.selectDistantMembersFromGenotypeFileJob.output,\
				parentJobLs=[self.highCoveragePanelDirJob, beagleOnHighCoverageJob, self.selectDistantMembersFromGenotypeFileJob,\
							tabixOnHighCoverageVCFJob], \
				extraDependentInputLs=[tabixOnHighCoverageVCFJob.output], transferOutput=False, \
				extraArguments=None, \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=7000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		##### Part 2 run Beagle on everyone with reference panel
		# run Beagle
		#refPanelFile=selectDistantMembersVariantsJob.output,\
		outputFnamePrefix = os.path.join(self.mapDirJob.folder, '%s.beagled'%(intervalFileBasenamePrefix))
		beagleJob = self.addBeagle4Job(executable=self.BeagleJava, \
						inputFile=VCFJobData.file, refPanelFile=None,\
						pedFile=self.outputPedigreeJob.output,\
						outputFnamePrefix=outputFnamePrefix, \
						burninIterations=7, phaseIterations=10, \
						noOfSamplingHaplotypesPerSample=4, duoscale=2, trioscale=2, \
						extraArguments=None, extraArgumentList=None,\
						parentJobLs=[self.mapDirJob, \
									self.outputPedigreeJob] + VCFJobData.jobLs, \
						transferOutput=False, no_of_cpus=None, \
						job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=13000).value,\
						walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value,\
						)
		returnData.beagleJob = beagleJob
		
		#index .vcf.gz, output of beagle, without index, GATK can't work on gzipped vcf
		tabixIndexFile = File('%s.tbi'%(beagleJob.output.name))
		tabixJob = self.addGenericJob(executable=self.tabix, \
						inputFile=beagleJob.output, inputArgumentOption="",\
						outputFile=None, outputArgumentOption="-o", \
						extraDependentInputLs=None, \
						extraOutputLs=[tabixIndexFile], transferOutput=False, frontArgumentList=["-p vcf"], \
						extraArguments=None, \
						extraArgumentList=[], \
						parentJobLs=[beagleJob, self.mapDirJob],\
						no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=2000, maxJobPropertyValue=4000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=600).value)
		
		#borrow PL to from pre-Beagle VCF to genotype 
		outputFile = File(os.path.join(self.mapDirJob.folder, '%s.beagled.withPL.vcf'%(intervalFileBasenamePrefix)))
		combineBeagleAndPreBeagleVariantsJob = self.addGATKJob(executable=self.CombineBeagleAndPreBeagleVariantsJava, \
					GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
					GATKAnalysisType="CombineBeagleAndPreBeagleVariants",\
					inputFile=None, inputArgumentOption=None, \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					inputFileList=None, argumentForEachFileInInputFileList="--variant",\
					interval=None, outputFile=outputFile, outputArgumentOption="--out", \
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=["--variant:first", beagleJob.output, "--variant:second", VCFJobData.file, \
								"-genotypeMergeOptions PRIORITIZE", "-priority first,second"], \
					extraOutputLs=None, \
					extraDependentInputLs=[beagleJob.output, VCFJobData.file] + tabixJob.outputLs, \
					parentJobLs=[beagleJob, tabixJob]+ VCFJobData.jobLs, transferOutput=False, \
					no_of_cpus=None, \
					key2ObjectForJob=None,\
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=2000, maxJobPropertyValue=4000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=600).value)
		#do not use "--variant:beagle" to name your vcf file as GATK would think it's in Beagle format
		
		#TrioCaller
		# 2013.06.11 replicate individuals who appear in more than 1 families
		round1_IndividualsReplicatedVCF = File( os.path.join(self.mapDirJob.folder, \
											'%s.replicate.vcf'%(intervalFileBasenamePrefix)))
		replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(\
					executable=self.ReplicateVCFGenotypeColumns, \
					inputF=combineBeagleAndPreBeagleVariantsJob.output, \
					sampleID2FamilyCountF=self.outputReplicatePedigreeJob.sampleID2FamilyCountF, \
					outputF=round1_IndividualsReplicatedVCF, \
					replicateIndividualTag=self.replicateIndividualTag,\
					parentJobLs=[self.outputReplicatePedigreeJob, self.mapDirJob, combineBeagleAndPreBeagleVariantsJob], \
					extraDependentInputLs=None, \
					transferOutput=False, \
					extraArguments=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value, \
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value,\
					)
		
		refineGenotypeOutputF = File(os.path.join(self.mapDirJob.folder, \
												'%s.trioCaller.vcf'%(intervalFileBasenamePrefix)))
		refineGenotypeJob = self.addTrioCallerJob(trioCallerWrapper=self.trioCallerWrapper, \
				trioCallerPath=self.trioCallerPath, \
				inputVCF=replicateVCFGenotypeColumnsJob.output,\
				pedFile=self.outputReplicatePedigreeJob.output, outputVCF=refineGenotypeOutputF, \
				inputPhased=True,\
				parentJobLs=[self.mapDirJob, replicateVCFGenotypeColumnsJob, self.outputReplicatePedigreeJob], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)	#1.2G memory for 12K loci
		
		returnData.refineGenotypeJob = refineGenotypeJob
		
		"""
		2013.07.10 the TrioCaller VCF has some info tags that are not described in VCF header
		"""
		outputFile = File(os.path.join(self.mapDirJob.folder, \
												'%s.extraInfoDesc.vcf'%(intervalFileBasenamePrefix)))
		addInfoDescJob = self.addGenericJob(executable=self.AddMissingInfoDescriptionToVCFHeader, \
					inputFile=refineGenotypeJob.output, \
					inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", \
					parentJobLs=[self.mapDirJob, refineGenotypeJob], \
					extraDependentInputLs=None, extraOutputLs=None, \
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, \
					transferOutput=False, sshDBTunnel=None, \
					key2ObjectForJob=None, objectWithDBArguments=None, \
					no_of_cpus=None, 
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
							minJobPropertyValue=1000, maxJobPropertyValue=3000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value,\
					max_walltime=None)
		
		# a CheckGenotypeConcordanceAmongReplicates.py job
		trioCallerReplicateConcordanceFile = File(os.path.join(self.statDirJob.folder, \
								'%s.trioCaller.concordance.tsv'%(intervalFileBasenamePrefix)))
		returnData.trioCallerReplicateConcordanceJob = self.addGATKJob(executable=self.CalculateConcordanceJava, \
					GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
					GATKAnalysisType="CalculateConcordanceAmongReplicates",\
					inputFile=refineGenotypeJob.output, inputArgumentOption="--variant", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					interval=None, \
					outputFile=trioCallerReplicateConcordanceFile, outputArgumentOption="--concordanceStatFname",\
					frontArgumentList=None, extraArguments="--replicateIndividualTag %s"%(self.replicateIndividualTag), \
					extraArgumentList=None, extraOutputLs=None, \
					parentJobLs=[self.statDirJob, refineGenotypeJob], \
					transferOutput=False, \
					no_of_cpus=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=6000, \
							minJobPropertyValue=9000, maxJobPropertyValue=16000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		#2013.06.14
		#merge replicates to generate consensus call
		# (not haplotype-based, as different recombination points across replicate haplotypes make it non-trivial )
		mergeReplicateOutputF = File(os.path.join(self.mapDirJob.folder, \
									'%s.replicatesMerged.vcf'%(intervalFileBasenamePrefix)))
		returnData.mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(\
							executable=self.MergeVCFReplicateHaplotypesJava,\
							GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
							inputF=addInfoDescJob.output, outputF=mergeReplicateOutputF, \
							replicateIndividualTag=self.replicateIndividualTag, \
							refFastaFList=self.registerReferenceData.refFastaFList, \
							parentJobLs=[self.mapDirJob, addInfoDescJob], \
							extraDependentInputLs=[], transferOutput=False, \
							extraArguments=None, \
							analysis_type='MergeVCFReplicateGenotypeColumns',\
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=5000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		return returnData
	
	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		2013.05.01
			#. concatenate all the sub-VCFs into one
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		
		refineGenotypeJobLs = [pdata.refineGenotypeJob for pdata in mapEachIntervalDataLs]
		mergeVCFReplicateColumnsJobLs = [pdata.mergeVCFReplicateColumnsJob for pdata in mapEachIntervalDataLs]
		
		
		realInputVolume = passingData.jobData.file.noOfIndividuals * passingData.jobData.file.noOfLoci
		baseInputVolume = 200*2000000
		#base is 4X coverage in 20Mb region => 120 minutes
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
							minJobPropertyValue=2000, maxJobPropertyValue=8000).value
		self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(passingData=passingData, \
						intervalJobLs=[pdata.beagleJob for pdata in mapEachIntervalDataLs],\
						outputDirJob=self.beagleReduceDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		
		self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(passingData=passingData, \
						intervalJobLs=refineGenotypeJobLs, outputDirJob=self.replicateVCFDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime, \
						**keywords)
		
		self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(passingData=passingData, \
						intervalJobLs=mergeVCFReplicateColumnsJobLs, outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		
		
		for pdata in mapEachIntervalDataLs:
			#add this output to the union job
			"""
			self.addInputToStatMergeJob(statMergeJob=self.reduceBeaglePhaseReplicateConcordanceJob_AllSites, \
							parentJobLs=[pdata.beaglePhasedReplicateConcordanceJob])
			self.addInputToStatMergeJob(statMergeJob=self.reduceBeaglePhaseReplicateConcordanceJob_HomoOnly, \
							parentJobLs=[pdata.beaglePhasedReplicateConcordanceJob])
			"""
			self.addInputToStatMergeJob(statMergeJob=self.reduceTrioCallerReplicateConcordanceJob_AllSites, \
							parentJobLs=[pdata.trioCallerReplicateConcordanceJob])
			self.addInputToStatMergeJob(statMergeJob=self.reduceTrioCallerReplicateConcordanceJob_HomoOnly, \
							parentJobLs=[pdata.trioCallerReplicateConcordanceJob])
		
		return returnData
		
		
	def registerCustomExecutables(self, workflow=None):
		"""
		2013.06
		"""
		parentClass.registerCustomExecutables(self, workflow)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='CalculateConcordanceJava', \
															clusterSizeMultipler=0.6)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='ProduceBeagleInputJava', \
															clusterSizeMultipler=0.6)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='PhaseByGATK',\
															clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='FindMendelErrorByGATK', \
															clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='RemoveMendelErrorSiteByGATK', \
															clusterSizeMultipler=0.5)
		#BeagleOnHCMOnkeys and BeagleJava are separate now because the former takes much less time than latter
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath,\
										name='BeagleOnHCMOnkeys', clusterSizeMultipler=0.1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/extractor/SelectDistantMembersFromGenotypeFile.py"), \
										name='SelectDistantMembersFromGenotypeFile', clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/splitter/SplitPedigreeVCFIntoBeagleTriosDuosFiles.py"), \
										name='SplitPedigreeVCFIntoBeagleTriosDuosFiles', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/converter/CombinePhasedBeagleOutputsIntoVCF.py"), \
										name='CombinePhasedBeagleOutputsIntoVCF', clusterSizeMultipler=1)		#2013.05.02 no clustering for OutputVRCPedigreeInTFAMGivenOrderFromFile
		self.setOrChangeExecutableClusterSize(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, clusterSizeMultipler=0)
		if hasattr(self, 'noOfCallingJobsPerNode') and self.noOfCallingJobsPerNode>=0:
			for executable in [self.BeagleJava]:
				#2013.2.26 use setOrChangeExecutableClusterSize to modify clusters size
				self.setOrChangeExecutableClusterSize(executable=executable, clusterSizeMultipler=self.noOfCallingJobsPerNode/float(self.clusters_size), \
													defaultClustersSize=self.clusters_size)
	
	
if __name__ == '__main__':
	main_class = BeagleAndTrioCallerOnVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()