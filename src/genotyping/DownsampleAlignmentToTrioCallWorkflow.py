#!/usr/bin/env python
"""
Examples:
	#2012.1.10 downsample two trio to 1X,1X,4X (father,mother,child)
	%s  -u yh --ref_ind_seq_id 524 --site_type 2 --hostname localhost
		-o dags/DownsampleAlignment/DownsampleAln558To1_649To1_618To4_557To1_615To1_626To4_5sampling_top5Contigs.xml 
		-j hcondor -l hcondor --contigMaxRankBySize 5 --noOfCallingJobsPerNode 1 --clusters_size 30 --site_id_ls 447 -e ~/ 
		-t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/ 
		--alnId2targetDepth 558:1,649:1,618:4,557:1,615:1,626:4 --db_passwd secret --no_of_sampling 5
	
	#2013.10.03
	%s --ref_ind_seq_id 3488 --site_id_ls 447 --sequence_filtered 1 --local_realigned 0 --reduce_reads 0
		--completedAlignment 1 --excludeContaminant --alignment_method_id 6
		--ind_aln_id_ls 5372,5373,5374,5376,5382,5384,5386,5391,5412,5489,5490
		-o dags/DownsampleAlignment/DownsampleEverySingleOf103AlignmentsTo1X_Beagle_Chr28.xml
		-j hcondor -l hcondor --minContigID 28 --maxContigID 28 --noOfCallingJobsPerNode 1
		--clusters_size 15 -t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/
		--refineGenotypeCallerType 2 --GATKGenotypeCallerType HaplotypeCaller
		--heterozygosityForGATK 0.01 --alignmentDepthIntervalMethodShortName 10XVRCRef3488MinLength1000
		--no_of_sampling 1
		--originalVCFFolder ~/NetworkData/vervet/workflow/AlignmentToCall/AlignmentToTrioCall_VRC_103Alignments_Chr28AlnDepthMethod4_Ref3488.2013.Sep.18T174612/RefineSNP/
		--alnId2targetDepth 5382:1 --genotypeCallerType 0
		--originalTrioInconsistencyFname MendelInconsistency/TrioInconsistency_VRC_103Alignments_Chr28AlnDepthMethod4_Ref3488_Beagle.2013.Oct.3T084955/folderTrioInconsistency/trio_inconsistency.tsv
		--db_user yh --needSSHDBTunnel --hostname localhost 
	
	# 2013.09.16 run downsampling workflow on 103 alignments. "--noOfCallingJobsPerNode 3" controls clustering of calling programs.
	# "--clusters_size 30" controls clustering for other programs.
	%s  --ref_ind_seq_id 3488 --site_id_ls 447 --sequence_filtered 1 --local_realigned 0 --reduce_reads 0
		--completedAlignment 1 --excludeContaminant --alignment_method_id 6
		--ind_aln_id_ls 5372,5373,5374,5376,5382,5384,5386,5391,5412,5489,5490,5491,5492,...
		-o dags/DownsampleAlignment/DownsampleAln5382To1_3sampling_In103Alignments_Chr28.xml
		-j hcondor -l hcondor --minContigID 28 --maxContigID 28 --noOfCallingJobsPerNode 1
		--clusters_size 15 -t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/
		--GATKGenotypeCallerType HaplotypeCaller --heterozygosityForGATK 0.01
		--alignmentDepthIntervalMethodShortName 10XVRCCoverageOn3488 --no_of_sampling 3
		--alnId2targetDepth 5382:1 --genotypeCallerType 1 --db_user yh --needSSHDBTunnel --hostname localhost
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import random
from Pegasus.DAX3 import Executable, File, PFN, Job, Link
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from vervet.src.qc.CalculateTrioInconsistencyPipeline import CalculateTrioInconsistencyPipeline
from vervet.src.qc.CheckTwoVCFOverlapPipeline import CheckTwoVCFOverlapPipeline
from vervet.src.db import VervetDB
from AlignmentToTrioCallPipeline import AlignmentToTrioCallPipeline
from AlignmentToCallPipeline import AlignmentToCallPipeline

parentClass = AlignmentToTrioCallPipeline
class DownsampleAlignmentToTrioCallWorkflow(parentClass, CalculateTrioInconsistencyPipeline, CheckTwoVCFOverlapPipeline):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						("probToSample", 1, float): [0.1, '', 1, 'probability for a read in a bam to be included in down-sampled bam, overridden by alnId2targetDepth'],\
						("no_of_sampling", 1, int): [10, '', 1, 'how many samplings to run, deprecated'],\
						("alignmentMinCoverageToDownsample", 1, float): [10, '', 1, 'minimum coverage (median-depth) for an alignment to be chosen for downsampling'],\
						("alignmentMaxCoverageToDownsample", 1, float): [100, '', 1, 'maximum coverage (median-depth) for an alignment to be chosen for downsampling'],\
						("alnId2targetDepth", 0, ): [None, '', 1, '[deprecated] a coma-separated list, in the form of alignment_id:targetDepth. 620:1,632:1,648:1'],\
						("targetDepth", 0, int): [1, '', 1, 'the target depth for each downsampling'],\
						('minDepth', 0, float): [0, '', 1, 'minimum depth for a call to regarded as non-missing', ],\
						('originalVCFFolder', 1, ): [None, '', 1, 'the folder that contains the calls before downsampling', ],\
						('originalTrioInconsistencyFname', 1, ): [None, '', 1, 'file that contains trio inconsistency results from original calls', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	#no overlap, and use 2Mb as interval (deal with alignment, not 2M loci)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 2000000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		parentClass.__init__(self, **keywords)
		if self.alnId2targetDepth:
			alnId2targetDepthLs = getListOutOfStr(self.alnId2targetDepth, data_type=str)
			self.alnId2targetDepth = {}
			for alnIdTargetDepth in alnId2targetDepthLs:
				alnIdTargetDepth = alnIdTargetDepth.split(':')
				alnIdTargetDepth = map(int, alnIdTargetDepth)
				alnId, targetDepth = alnIdTargetDepth
				self.alnId2targetDepth[alnId] = targetDepth
		else:
			self.alnId2targetDepth = {}

	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		CalculateTrioInconsistencyPipeline.registerCustomExecutables(self, workflow=workflow)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
										name='DownsampleSamJava', clusterSizeMultipler=0.001)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'qc/AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey.py'), \
										name='AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey', clusterSizeMultipler=0.3)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'qc/AnalyzeTrioInconsistencyInRelationToPerturbedMonkey.py'), \
										name='AnalyzeTrioInconsistencyInRelationToPerturbedMonkey', clusterSizeMultipler=0.3)
		
	def registerJars(self, workflow=None):
		"""
		2012.1.6
			some custom jars
		"""
		parentClass.registerJars(self, workflow=workflow)
		
		self.registerOneJar(name="DownsampleSamJar", path=os.path.join(self.picard_path, 'DownsampleSam.jar'))
	
	def addDownsampleSamJob(self, workflow=None, executable=None, DownsampleSamJar=None, \
							inputFile=None, probToSample=0.1, outputFile=None, \
							extraDependentInputLs=None, \
							parentJobLs=None, job_max_memory=3000, walltime=600,\
							transferOutput=False, **keywords):
		"""
		2012.1.6
		"""
		if workflow is None:
			workflow = self
		#2012.1.10 Changing the RANDOM_SEED argument is important as otherwise it is identical series of random integers in java program
		#	given the same default RANDOM_SEED (=1). another option is to set it 'null'.
		extraArgumentList = ["P=%s"%(probToSample), "RANDOM_SEED=%s"%(int(random.random()*10000)), "VALIDATION_STRINGENCY=LENIENT"]
		job = self.addGenericJavaJob(executable=executable, jarFile=DownsampleSamJar, \
					inputFile=inputFile, inputArgumentOption="I=", \
					outputFile=outputFile, outputArgumentOption="O=",\
					frontArgumentList=None, \
					extraArguments=None, extraArgumentList=extraArgumentList, \
					extraOutputLs=None, \
					extraDependentInputLs=extraDependentInputLs, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					key2ObjectForJob=None, no_of_cpus=None, walltime=walltime, sshDBTunnel=None, **keywords)
		return job
	
	def addDownsampleJobToSelectedAlignment(self, workflow=None, alignmentDataLs=[], alnId2targetDepth={}, \
							downsampleDirJob=None, chr2IntervalDataLs=None, **keywords):
		"""
		2013.09.07 do chromosome selection to reduce memory usage
		2012.1.9
		"""
		if workflow is None:
			workflow = self
		
		sys.stderr.write("Adding downsampling jobs to %s selected alignments out of %s total alignments, job count=%s ..."%\
						(len(alnId2targetDepth), len(alignmentDataLs), self.no_of_jobs))
		returnData = []
		
		defaultDownsamplerMaxMemory = 6000	#in Mb
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			if alignment.id in alnId2targetDepth:
				targetDepth = alnId2targetDepth.get(alignment.id)
				if alignment.median_depth is None or alignment.median_depth==0:
					sys.stderr.write("Warning: alignment %s's median depth is %s and could not be down-sampled. Ignore.\n"%(alignment.id, alignment.median_depth))
					continue
				probToSample = float(targetDepth)/alignment.median_depth
				if probToSample<1:
					bamF = alignmentData.bamF
					baiF = alignmentData.baiF
					AlignmentJobAndOutputLs = []
					alignmentBasenamePrefix = os.path.splitext(os.path.basename(alignment.path))[0]
					
					for chromosome, intervalDataLs in chr2IntervalDataLs.iteritems():
						downsampleWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=alignment.median_depth, \
											baseInputVolume=6, baseJobPropertyValue=300, \
											minJobPropertyValue=120, maxJobPropertyValue=1320).value
						downsampleMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=alignment.median_depth, \
											baseInputVolume=6, baseJobPropertyValue=8000, \
											minJobPropertyValue=6000, maxJobPropertyValue=14000).value
					
						selectedBamFile = File(os.path.join(downsampleDirJob.output, "%s_%s.bam"%(alignmentBasenamePrefix, chromosome)))
						selectAlignmentJob, selectAlignmentIndexJob = self.addSelectAlignmentJob(inputFile=bamF, \
							outputFile=selectedBamFile, region=chromosome, \
							parentJobLs=[downsampleDirJob] + alignmentData.jobLs, \
							extraDependentInputLs=[baiF], transferOutput=False, \
							extraArguments=None, job_max_memory=3000, needBAMIndexJob=True)
					
						outputBamF = File(os.path.join(downsampleDirJob.output, "%s_%s_%sTo%s.bam"%(alignmentBasenamePrefix, chromosome, \
																								alignment.median_depth, targetDepth)))
						downsampleJob = self.addDownsampleSamJob(executable=workflow.DownsampleSamJava, DownsampleSamJar=self.DownsampleSamJar, \
										inputFile=selectAlignmentJob.output, \
										probToSample=probToSample, outputFile=outputBamF, \
										parentJobLs=[selectAlignmentJob, selectAlignmentIndexJob], \
										extraDependentInputLs=[selectAlignmentIndexJob.output], \
										job_max_memory=downsampleMaxMemory, walltime=downsampleWalltime,\
										transferOutput=False)
						index_sam_job = self.addBAMIndexJob(workflow, BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
							inputBamF=outputBamF, parentJobLs=[downsampleJob], transferOutput=False, javaMaxMemory=3500)
						AlignmentJobAndOutputLs.append(PassingData(jobLs=[downsampleJob, index_sam_job], file=downsampleJob.output,\
																fileLs=[downsampleJob.output, index_sam_job.output]))
					
					outputBamF = File(os.path.join(downsampleDirJob.output, "%s_%s.bam"%(alignmentBasenamePrefix, probToSample)))
					mergeAlignmentJob, mergeAlignmentIndexJob = self.addAlignmentMergeJob(AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
						outputBamFile=outputBamF,\
						MergeSamFilesJava=None, \
						BuildBamIndexFilesJava=None, \
						parentJobLs=None, transferOutput=False,\
						job_max_memory=8000, walltime=680)
					
					newAlignmentData = PassingData(alignment=alignment)
					#don't modify the old alignmentData as it will affect the original alignmentDataLs, which should remain same across different samplings
					newAlignmentData.jobLs = [mergeAlignmentJob, mergeAlignmentJob.bamIndexJob]	#downsampleJob has to be included otherwise its output (bamF) will be wiped out after index_sam_job is done
					newAlignmentData.bamF = mergeAlignmentJob.bamIndexJob.bamFile
					newAlignmentData.baiF = mergeAlignmentJob.bamIndexJob.baiFile
				else:
					newAlignmentData = alignmentData
					
			else:
				newAlignmentData = alignmentData
			
			returnData.append(newAlignmentData)
		sys.stderr.write(" job count=%s.\n"%(self.no_of_jobs))
		return returnData
	
	def setup(self, **keywords):
		"""
		2013.10.02
		"""
		parentClass.setup(self, **keywords)
		
		outputDirPrefix = "downSample"
		self.auxDirJob = self.addMkDirJob(outputDir="%sAuxilliary"%(outputDirPrefix))
		
		#output complete pedigree, use Individual.id as id
		# used by programs that analyze concordance and trio-inconsistency results
		pedigreeFileFormat = 5
		pedFile = File(os.path.join(self.auxDirJob.output, 'pedigree.format%s.txt'%(pedigreeFileFormat)))
		#sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'pedigree.sampleID2FamilyCount.%s.format%s.txt'%\
		#						(inputFileBasenamePrefix, pedigreeFileFormat)))
		self.wholePedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=None, outputFile=pedFile, \
				sampleID2FamilyCountF=None,\
				polymuttDatFile = None,\
				outputFileFormat=pedigreeFileFormat, \
				replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=[self.auxDirJob], \
				extraDependentInputLs=None, transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		
		#a job that outputs alignment coverage (alignment.read_group, median_depth)
		alignmentDepthFile = File(os.path.join(self.auxDirJob.folder, 'refIndSeq%s_AlignmentDepth.tsv'%(self.ref_ind_seq_id)))
		self.outputAlignmentDepthJob = self.addOutputVCFAlignmentDepthRangeJob(executable=self.OutputVCFAlignmentDepthRange, \
						inputFile=None, \
						ref_ind_seq_id=self.ref_ind_seq_id, depthFoldChange=None, minGQ=None,\
						outputFile=alignmentDepthFile, outputFileFormat=1,\
						data_dir=self.data_dir,\
						sequence_filtered=self.sequence_filtered, excludeContaminant=self.excludeContaminant, \
						local_realigned=self.local_realigned, reduce_reads=self.reduce_reads, \
						completedAlignment=self.completedAlignment, \
						alignment_method_id=self.alignment_method_id, \
						alignment_outdated_index=self.alignment_outdated_index,\
						extraArgumentList=None,\
						parentJobLs=[self.auxDirJob], \
						extraDependentInputLs=None, transferOutput=True, \
						job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
	def addConcordanceAnalysisJob(self, inputFile=None, outputFile=None,\
						perturbedMonkeyList=None, wholePedigreeJob=None, outputAlignmentDepthJob=None, \
						parentJobLs=None, transferOutput=True, job_max_memory=None, walltime=None):
		"""
		2013.10.03
		"""
		if parentJobLs is None:
			parentJobLs = []
		job = self.addGenericJob(executable=self.AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey, \
					inputFile=inputFile, outputFile=outputFile, \
					parentJobLs=parentJobLs + [wholePedigreeJob, outputAlignmentDepthJob], \
					extraArgumentList=["--perturbedMonkeyList %s"%(','.join(perturbedMonkeyList)), \
									'--individualAlignmentCoverageFname', outputAlignmentDepthJob.output, \
									"--pedigreeFname", wholePedigreeJob.output ], \
					extraDependentInputLs=[outputAlignmentDepthJob.output, wholePedigreeJob.output],\
					extraOutputLs=None,\
					transferOutput=transferOutput, \
					key2ObjectForJob=None, job_max_memory=job_max_memory, walltime=walltime)
		return job
	
	def addTrioInconsistencyAnalysisJob(self, inputFile=None, outputFile=None,\
						originalTrioInconsistencyFile=None,
						perturbedMonkeyList=None, wholePedigreeJob=None, outputAlignmentDepthJob=None, \
						parentJobLs=None, transferOutput=True, job_max_memory=None, walltime=None):
		"""
		2013.10.03
		"""
		if parentJobLs is None:
			parentJobLs = []
		job = self.addGenericJob(executable=self.AnalyzeTrioInconsistencyInRelationToPerturbedMonkey, \
					inputFile=inputFile, outputFile=outputFile, \
					parentJobLs=parentJobLs + [wholePedigreeJob, outputAlignmentDepthJob], \
					extraArgumentList=["--perturbedMonkeyList %s"%(','.join(perturbedMonkeyList)), \
									'--individualAlignmentCoverageFname', outputAlignmentDepthJob.output, \
									"--pedigreeFname", wholePedigreeJob.output,\
									"--originalTrioInconsistencyFname", originalTrioInconsistencyFile], \
					extraDependentInputLs=[outputAlignmentDepthJob.output, wholePedigreeJob.output, originalTrioInconsistencyFile],\
					extraOutputLs=None,\
					transferOutput=transferOutput, \
					key2ObjectForJob=None, job_max_memory=job_max_memory, walltime=walltime)
		return job
	
	def run(self):
		"""
		2011-7-11
		"""
		self.needSplitChrIntervalData = True	#2013.09.03 turned off in AlignmentToCallPipeline.py
		
		pdata = self.setup_run()
		workflow = pdata.workflow
		db_vervet = self.db
		
		self.setup()
		
		cumulativeMedianDepth = self.db.getCumulativeAlignmentMedianDepth(alignmentLs=pdata.alignmentLs, \
										defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		"""
		calculateTrioInconsistencyPipeline_ins = CalculateTrioInconsistencyPipeline(drivername=self.drivername, hostname=self.hostname, dbname=self.dbname, \
							schema=self.schema, db_user=self.db_user, db_passwd=self.db_passwd, ref_ind_seq_id=self.ref_ind_seq_id, \
							samtools_path=self.samtools_path, picard_path=self.picard_path, gatk_path=self.gatk_path,\
							vervetSrcPath=self.vervetSrcPath, home_path=self.home_path, tabixPath=self.tabixPath, javaPath=self.javaPath,\
							data_dir=self.data_dir, local_data_dir=self.local_data_dir, site_handler=self.site_handler,\
							input_site_handler=self.input_site_handler, clusters_size=self.clusters_size,
							outputFname=self.outputFname, checkEmptyVCFByReading=False,\
							debug=self.debug, report=self.report, inputDir="random", \
							maxContigID=self.maxContigID, minContigID=self.minContigID,\
							contigMaxRankBySize=self.contigMaxRankBySize, contigMinRankBySize=self.contigMinRankBySize, \
							ref_genome_tax_id=self.ref_genome_tax_id, ref_genome_version=self.ref_genome_version,\
							ref_genome_sequence_type_id=self.ref_genome_sequence_type_id,\
							ref_genome_outdated_index=self.ref_genome_outdated_index,\
							minDepth=self.minDepth, max_walltime=self.max_walltime)
							#checkEmptyVCFByReading is not used in addJobs()
		calculateTrioInconsistencyPipeline_ins.needSplitChrIntervalData = False	#no need for this as it's taking the jobs from addGenotypeCallJobs()
		calculateTrioInconsistencyPipeline_ins.setup_run()
		"""
		
		origAlignmentDataLs = self.alignmentDataLs
		
		"""
		#reduce the trio consistency by the same trio across different samplings 
		allTrioInconsistencyFile = File('trio_inconsistency_avg_all_samples.tsv')
		allTrioInconsistencyJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByAverageColumnsWithSameKey, \
						outputF=allTrioInconsistencyFile, extraArguments='-k 0 -v 1,2,3', parentJobLs=None, \
						extraDependentInputLs=None, transferOutput=True)
		"""
		trioLs = self.getDuoTrioFromAlignmentLs(db_vervet, self.alignmentLs)
		
		#add jobs to select only certain chromosome(s), (maxContigID, minContigID) out of alignment files
		originalVCFData = self.registerAllInputFiles(inputDir=self.originalVCFFolder, input_site_handler=None, \
					checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName='%sOriginalVCF'%(self.pegasusFolderName),\
					maxContigID=self.maxContigID, minContigID=self.minContigID, \
					db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLociInVCF=None, includeIndelVCF=True)
		originalTrioInconsistencyFile = self.registerOneInputFile(inputFname=self.originalTrioInconsistencyFname, folderName=self.pegasusFolderName, \
								checkFileExistence=True)
		
		#for i in xrange(self.no_of_sampling):
		noOfDownsampledMonkeys = 0
		for alignmentData in origAlignmentDataLs:	#downsample every single alignment
			if alignmentData.alignment.median_depth<self.alignmentMinCoverageToDownsample or alignmentData.alignment.median_depth>self.alignmentMaxCoverageToDownsample:
				#skip if original alignment's coverage is too low or too high
				continue
			noOfDownsampledMonkeys += 1
			perturbedMonkeyList = [alignmentData.alignment.read_group]
			downsampleDir="downsample_%s_%sTo%s"%(alignmentData.alignment.read_group, alignmentData.alignment.median_depth, self.targetDepth)
			downsampleDirJob = self.addMkDirJob(outputDir=downsampleDir)
			alnId2targetDepth = {alignmentData.alignment.id: self.targetDepth}
			alignmentDataLs = self.addDownsampleJobToSelectedAlignment(alignmentDataLs=origAlignmentDataLs, \
											alnId2targetDepth=alnId2targetDepth, downsampleDirJob=downsampleDirJob,\
											chr2IntervalDataLs=self.chr2IntervalDataLs)
			#--genotypeCallerType: '0: SAMtools, 1: GATK (--GATKGenotypeCallerType ...), 2: Platypus'
			genotypeCallJobData = self.addGenotypeCallJobs(alignmentDataLs=alignmentDataLs, chr2IntervalDataLs=self.chr2IntervalDataLs,\
						registerReferenceData=pdata.registerReferenceData, \
						site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						site_type=self.site_type, data_dir=self.data_dir, no_of_gatk_threads = 1,\
						outputDirPrefix="%s/"%(downsampleDir), genotypeCallerType=self.genotypeCallerType, \
						cumulativeMedianDepth=cumulativeMedianDepth, \
						sourceVCFFolder=None, extraArgumentPlatypus="")
			
			concordanceWithOriginalData = self.addCheckingVCFOverlapSubWorkflow(workflow=workflow, chr2size = self.chr2size, \
						inputVCFData1=genotypeCallJobData, \
						inputVCFData2=originalVCFData, \
						registerReferenceData=pdata.registerReferenceData,\
						outputDirPrefix="%s/CheckOverlapAndConcordance"%(downsampleDir))
			
			outputFile = File(os.path.join(downsampleDirJob.output, 'ConcordanceInRelationToPerturbedMonkey.tsv'))
			self.addConcordanceAnalysisJob(inputFile=concordanceWithOriginalData.perSampleMatchFractionReduceJob.output, \
								outputFile=outputFile, perturbedMonkeyList=perturbedMonkeyList, \
								wholePedigreeJob=self.wholePedigreeJob, outputAlignmentDepthJob=self.outputAlignmentDepthJob, \
								parentJobLs=[downsampleDirJob, concordanceWithOriginalData.perSampleMatchFractionReduceJob], \
								transferOutput=True,\
								job_max_memory=2000, walltime=120)
			
			#reduce small-interval jobs into bigger ones so that trioInconsistency jobs do not explode
				#one per trio,
				# bigIntervalSize is measured in bp
			intervalJobLs = [jobData.job for jobData in genotypeCallJobData.jobDataLs]
			bigIntervalSize = 40000000	#40Mb
			outputDirJob = self.addMkDirJob(outputDir="%s/bigInterval_%s_VCF"%(downsampleDir, bigIntervalSize))
			bigIntervalGenotypeCallJobData = self.reduceManySmallIntervalVCFIntoBigIntervalVCF(executable=self.CombineVariantsJava, \
							registerReferenceData=pdata.registerReferenceData, fileBasenamePrefix="", \
							intervalJobLs=intervalJobLs, outputDirJob=outputDirJob, bigIntervalSize=bigIntervalSize,
							transferOutput=False, job_max_memory=7000, walltime=300, needBGzipAndTabixJob=False)
			trioInconsistencyJobData = self.addCalculateTrioInconsistencyJobs(inputVCFData=bigIntervalGenotypeCallJobData, \
				trioLs=trioLs, \
				addTrioSpecificPlotJobs=None, addTrioContigSpecificPlotJobs=None,\
				outputDirPrefix="%s/"%downsampleDir)
			outputFile = File(os.path.join(downsampleDirJob.output, 'TrioInconsistencyInRelationToPerturbedMonkey.tsv'))
			self.addTrioInconsistencyAnalysisJob(inputFile=trioInconsistencyJobData.trioInconsistencySummaryJob.output, \
								outputFile=outputFile, \
								originalTrioInconsistencyFile=originalTrioInconsistencyFile, perturbedMonkeyList=perturbedMonkeyList, \
								wholePedigreeJob=self.wholePedigreeJob, outputAlignmentDepthJob=self.outputAlignmentDepthJob, \
								parentJobLs=[trioInconsistencyJobData.trioInconsistencySummaryJob], \
								transferOutput=True,\
								job_max_memory=2000, walltime=120)
			
			#add trio inconsistency summary output to reduction job
			#self.addInputToStatMergeJob(statMergeJob=allTrioInconsistencyJob, \
			#					inputF=trioInconsistencyJobData.trioInconsistencySummaryJob.output, \
			#					parentJobLs=[trioInconsistencyJobData.trioInconsistencySummaryJob])
		sys.stderr.write(" %s monkeys to be downsampled to %s.\n"%(noOfDownsampledMonkeys, self.targetDepth))
		self.end_run()

if __name__ == '__main__':
	main_class = DownsampleAlignmentToTrioCallWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
