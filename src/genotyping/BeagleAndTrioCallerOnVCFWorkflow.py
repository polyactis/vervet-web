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
				extraDependentInputLs=[], transferOutput=True, \
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
				extraDependentInputLs=[], transferOutput=True, \
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
				extraDependentInputLs=[], transferOutput=transferOutput, \
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
	
	def addBeagle3Job(self, workflow=None, executable=None, BeagleJar=None, \
					phasedBeagleInputFile=None,\
					likelihoodBeagleInputFile=None, triosBeagleInputFile=None, pairsBeagleInputFile=None,\
					unphasedBeagleInputFile=None,\
					markersBeagleInputFile=None,
					outputFnamePrefix=None, noOfIterations=None, noOfSamplingHaplotypesPerSample=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=120, **keywords):
		"""
		i.e.
			
		2013.06.13 a generic function to add Beagle version 3 jobs
		
		crocea@vervetNFS:~/bin/Beagle/example/imputation$ 
		prefix=method_36_Contig791_replicated_phasedByGATK_noMendelError_unphased_minProb0.65;
		java -Xmx15g -XX:PermSize=1000m -XX:MaxPermSize=6000m
			-jar ~/bin/Beagle/beagle.jar like=$prefix\_familySize1.bgl
			pairs=$prefix\_familySize2.bgl trios=$prefix\_familySize3.bgl
			markers=$prefix.markers out=$prefix\_beagled missing=?
			unphased=...bgl
			
		Note:
			"missing=? " can not be used if all input is likelihood data.
		
		Beagle help http://faculty.washington.edu/browning/beagle/beagle.html
		
		niterations=<number of iterations> where <number of iterations> is a positive even 
			integer giving the number of iterations of the phasing algorithm. If an odd integer is 
			specified, the next even integer is used. The niterations argument is optional. The default 
			value is niterations=10. The default value typically gives good accuracy
		
		nsamples=<number of samples> where <number of samples> is positive integer giving 
			the number of haplotype pairs to sample for each individual during each iteration of the 
			phasing algorithm. The nsamples argument is optional. The default value is nsamples=4.
			If you are phasing an extremely large sample (say > 4000 individuals), you may want to 
			use a smaller nsamples parameter (e.g. 1 or 2) to reduce computation time. If you are 
			phasing a small sample (say < 200 individuals), you may want to use a larger nsamples
			parameter (say 10 or 20) to increase accuracy.
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if BeagleJar is None:
			BeagleJar = self.BeagleJar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		#place holder for key output files, corresponding to these input like=, unphased=, trios=, pairs=
		key2ObjectForJob = {"likelihoodPhasedOutputFile":None,\
						"singletonPhasedOutputFile":None,\
						"pairsPhasedOutputFile":None,\
						"triosPhasedOutputFile":None,\
						}
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		frontArgumentList.extend([javaMemRequirement, '-jar', BeagleJar, "out=%s"%(outputFnamePrefix)])
		extraDependentInputLs.append(BeagleJar)
		
		if noOfIterations is not None:
			frontArgumentList.append("niterations=%s"%(noOfIterations))
		if noOfSamplingHaplotypesPerSample is not None:
			frontArgumentList.append("nsamples=%s"%(noOfSamplingHaplotypesPerSample))
		
		logFile=File('%s.log'%(outputFnamePrefix))
		extraOutputLs.append(logFile)
		key2ObjectForJob['logFile'] = logFile
		
		phasedBeagleFileSuffix = '.bgl.phased.gz'
		onlyLikelihoodInput = True	#determines whether missing=? should be added
		if likelihoodBeagleInputFile:
			frontArgumentList.extend(["like=%s"%(likelihoodBeagleInputFile.name)])
			extraDependentInputLs.append(likelihoodBeagleInputFile)
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(likelihoodBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['likelihoodPhasedOutputFile'] = phasedFile
			#additional output for the likelihood input
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if triosBeagleInputFile:
			frontArgumentList.extend(["trios=%s"%(triosBeagleInputFile.name)])
			extraDependentInputLs.append(triosBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(triosBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['triosPhasedOutputFile'] = phasedFile
		
		if pairsBeagleInputFile:
			frontArgumentList.extend(["pairs=%s"%(pairsBeagleInputFile.name)])
			extraDependentInputLs.append(pairsBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(pairsBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['pairsPhasedOutputFile'] = phasedFile
		
		if unphasedBeagleInputFile:
			frontArgumentList.extend(["unphased=%s"%(unphasedBeagleInputFile.name)])
			extraDependentInputLs.append(unphasedBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(unphasedBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['singletonPhasedOutputFile'] = phasedFile
			#additional likelihood-related output
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if phasedBeagleInputFile:	#reference panel
			frontArgumentList.extend(["phased=%s"%(phasedBeagleInputFile.name)])
			extraDependentInputLs.append(phasedBeagleInputFile)
			
		
		if markersBeagleInputFile:
			frontArgumentList.extend(["markers=%s"%(markersBeagleInputFile.name)])
			extraDependentInputLs.append(markersBeagleInputFile)
		if not onlyLikelihoodInput:
			frontArgumentList.append("missing=?")
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, \
					inputArgumentOption=None,  inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		
		return job
	
	def addBeagle4Job(self, workflow=None, executable=None, BeagleJar=None, \
					inputFile=None, refPanelFile=None, pedFile=None,\
					outputFnamePrefix=None, \
					burninIterations=None, phaseIterations=None, \
					noOfSamplingHaplotypesPerSample=None, singlescale=None, duoscale=None, trioscale=None,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					no_of_cpus=None, walltime=120, **keywords):
		"""
		i.e.
			
		2013.06.22 a generic function to add Beagle version 4 job
		
		java commandline example:
			...
		
		Beagle help http://faculty.washington.edu/browning/beagle/beagle.html
		
		gt=[file] specifies a VCF file containing a GT (genotype) format field for each marker.
		gl=[file] specifies a VCF file containing a GL or PL (genotype likelihood) format field for
			each marker. If both GL and PL format fields are present for a sample, the GL format will
			be used. See also the maxlr parameter.
		ref=[file] specifies a reference VCF file containing additional samples and phased
			genotypes for each marker. Use of an appropriate reference panel can increase accuracy.
		ped=[file] specifies a Linkage-format pedigree file for specifying family relationships.
			The pedigree file has one line per individual. The first 4 white-space delimited fields of
			each line are 1) pedigree ID, 2) individual ID, 3) father's ID, and 4) mother's ID. A "0" is
			used in column 3 or 4 if the father or mother is unknown. Beagle uses the data in columns
			2-4 to identify parent-offspring duos and trios. Any or all columns of the pedigree file
			after column 4 may be omitted. See also the duoscale and trioscale parameters.
		out=[prefix] specifies the output filename prefix. The prefix may be an absolute or
			relative filename, but it cannot be a directory name.
		excludesamples=[file] specifies a file containing non-reference samples (one sample per
			line) to be excluded from the analysis and output files.
		excluderefsamples=[file] specifies a file containing reference samples (one sample per
			line) to be excluded from the analysis.
		excludemarkers=[file] specifies a file containing markers (one marker per line) to be
			excluded from the analysis and the output files. An excluded marker identifier can either
			be an identifier from the VCF record's ID field or genomic coordinates in the format:
			CHROM:POS.
		chrom=[chrom:start-end] specifies a chromosome or chromosome interval using a
			chromosome identifier in the VCF file and the starting and ending positions of the
			interval. The entire chromosome, the beginning of the chromosome, and the end of a
			chromosome can be specified by chrom=[chrom], chrom=[chrom:-end], and chrom=[chrom:start-] respectively.
		
		window=[positive integer] (default: window=24000).
		overlap=[positive integer] (default: overlap=3000)
		gprobs=[true/false] (default: gprobs=true).
		usephase=[true/false] (default: usephase=false).
		singlescale=[positive number] (default: singlescale=1.0), change the scale to x, program samples x*x faster
		duoscale=[positive number]
		trioscale=[positive number] 
		burnin-its=[non-negative integer] (default: burnin=5).
		phase-its=[non-negative integer] (default: phase-its=5)
		
		
		Advanced options not recommended for general use:
		
		dump=[file] specifies a file containing sample identifiers (one identifier per line). For
			each marker window, all the sampled haplotypes for these individuals which are sampled
			after the burn-in iterations are printed to an output VCF files (dump.[window #].gz).
		nsamples=[positive integer] specifies the number of haplotype pairs to sample for each
			individual during each iteration of the algorithm (default: nsamples=4).
		buildwindow=[positive integer] specifies the number of markers used to build the
			haplotype frequency model at each locus (default: buildwindow=500).
		
		
		Three output files are created whose names begin with the output file prefix specified on
			the command line argument and whose names end with the suffixes: .log, .vcf.gz, and .ibd.

		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if BeagleJar is None:
			BeagleJar = self.Beagle4Jar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		key2ObjectForJob = {}
		
		extraArgumentList.extend(["out=%s"%(outputFnamePrefix)])
		if inputFile:
			extraArgumentList.append("gl=%s"%(inputFile.name))
			extraDependentInputLs.append(inputFile)
		if refPanelFile:
			extraArgumentList.append("ref=%s"%(refPanelFile.name))
			extraDependentInputLs.append(refPanelFile)
		if pedFile:
			extraArgumentList.append("ped=%s"%(pedFile.name))
			extraDependentInputLs.append(pedFile)
		if burninIterations is not None:
			extraArgumentList.append("burnin-its=%s"%(burninIterations))
		if phaseIterations is not None:
			extraArgumentList.append("phase-its=%s"%(phaseIterations))
		if noOfSamplingHaplotypesPerSample is not None:
			extraArgumentList.append("nsamples=%s"%(noOfSamplingHaplotypesPerSample))
		if singlescale is not None:
			extraArgumentList.append("singlescale=%s"%(singlescale))
		if duoscale is not None:
			extraArgumentList.append("duoscale=%s"%(duoscale))
		if trioscale is not None:
			extraArgumentList.append("trioscale=%s"%(trioscale))
		
		outputVCFFile = File("%s.vcf.gz"%(outputFnamePrefix))
		extraOutputLs.append(outputVCFFile)	#this would be accessible through job.output and job.vcfOutputFile 
		key2ObjectForJob['vcfOutputFile'] = outputVCFFile
		
		logFile=File('%s.log'%(outputFnamePrefix))
		extraOutputLs.append(logFile)
		key2ObjectForJob['logFile'] = logFile	#this would be accessible as job.logFile
		
		job = self.addGenericJavaJob(executable=executable, jarFile=BeagleJar, \
					inputFile=None, inputArgumentOption=None, \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=None, outputArgumentOption=None,\
					frontArgumentList=frontArgumentList, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraOutputLs=extraOutputLs, \
					extraDependentInputLs=extraDependentInputLs, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					key2ObjectForJob=key2ObjectForJob, \
					no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		return job
	
	def convertVCF2Beagle(self, workflow=None, VCFJobData=None, outputDirJob=None, \
						outputFileBasenamePrefix=None,\
						outputPedigreeJob=None, pedigreeSplitStructure=None,\
						transferOutput=False, \
						job_max_memory=None, walltime=None):
		"""
		2013.06.23 deprecated as Beagle v4 is used and it accepts VCF files
		2013.05.01
			convert VCF file into Beagle input format
		"""
		
		#replicate the individuals involved in more than one trio/duo
		#2012.4.2 replicate individuals who appear in more than 1 families
		#note: remove individuals in pedigree file but not in VCF file
		
		#GATK ProduceBeagleInput
		beagleLikelihoodFile = File(os.path.join(outputDirJob.folder, "%s.bgl"%(outputFileBasenamePrefix)))
		produceBeagleInputJob = self.addGATKJob(executable=self.ProduceBeagleInputJava, \
					GATKAnalysisType="ProduceBeagleInput", \
					inputFile=VCFJobData.file, inputArgumentOption="--variant:VCF", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					outputFile=beagleLikelihoodFile, \
					parentJobLs=VCFJobData.jobLs + [outputDirJob], \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=None, \
					extraOutputLs=None, extraDependentInputLs=None, \
					job_max_memory=job_max_memory, \
					no_of_cpus=None, walltime=walltime)
		
		#a SplitPedigreeVCFIntoBeagleTriosDuosFiles job
		#need pedigreeFile (from outputPedigreeJob)
		outputFnamePrefix = os.path.join(outputDirJob.folder, outputFileBasenamePrefix)
		
		key2File= {'size1File': None, 'size2File': None, 'size3File':None}	#to store output files
			#size1File: output for singletons (likelihood format)
			#size2File: output for duos (Beagle genotype format)
			#size3File: output for trios (Beagle genotype format)
		markerFile = File("%s.markers"%(outputFnamePrefix))
		extraOutputLs = []	#first is the markers file
		extraOutputLs.append(markerFile)
		key2File['markerFile'] = markerFile
		for familySize, familyLs in pedigreeSplitStructure.familySize2familyLs.iteritems():
			if familyLs:	#non-empty
				outputFile = File('%s_familySize%s.bgl'%(outputFnamePrefix, familySize))
				key2File['size%sFile'%(familySize)] = outputFile
				extraOutputLs.append(outputFile)
		
		splitPedigreeVCFIntoBeagleTriosDuosFilesJob = self.addGenericJob(executable=self.SplitPedigreeVCFIntoBeagleTriosDuosFiles, \
						inputFile=VCFJobData.file, \
						outputFile=None, \
						parentJobLs=VCFJobData.jobLs + [produceBeagleInputJob, outputDirJob, outputPedigreeJob], \
						extraDependentInputLs=[produceBeagleInputJob.output, outputPedigreeJob.output], \
						extraOutputLs=extraOutputLs,  transferOutput=transferOutput, \
						extraArguments=None, \
						extraArgumentList=["--gatkPrintBeagleFname", produceBeagleInputJob.output, \
										"--plinkPedigreeFname", outputPedigreeJob.output, \
										"--minProbForValidCall %s"%(self.minProbForValidCall), \
										"--dummyIndividualNamePrefix dummy", \
										"--outputFnamePrefix", outputFnamePrefix], \
						job_max_memory=job_max_memory, \
						sshDBTunnel=None, key2ObjectForJob=key2File, objectWithDBArguments=None, \
						no_of_cpus=None, walltime=walltime)
		
		#
		return splitPedigreeVCFIntoBeagleTriosDuosFilesJob
	
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
				extraDependentInputLs=[], transferOutput=False, \
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
		combineBeagleAndPreBeagleVariantsJob = self.addGATKJob(executable=self.CombineVariantsJava, \
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
						intervalJobLs=mergeVCFReplicateColumnsJobLs, outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(passingData=passingData, \
						intervalJobLs=refineGenotypeJobLs, outputDirJob=self.replicateVCFDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime, \
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