#!/usr/bin/env python
"""
Examples:
	
	#2011-11-5 run it on hoffman2, need ssh tunnel for db (--needSSHDBTunnel)
	%s -a 524 -j hoffman2 -l hoffman2 -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_hoffman2.xml 
		--ind_aln_id_ls 552-661 -e /u/home/eeskin/polyacti/ --tmpDir /u/home/eeskin/polyacti/NetworkData/ 
		-J /u/local/apps/java/jre1.6.0_23/bin/java -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /Network/Data/vervet/db/
		--needSSHDBTunnel
	
	#2011-11-5 run on uschpc (input data is on uschpc), for each top contig as well
	%s -a 524 -j uschpc -l uschpc -u yh -z uclaOffice -o MarkDupAlnID552_661Pipeline_uschpc.xml
		--ind_aln_id_ls 552-661 --needPerContigJob -e /home/cmb-03/mn/yuhuang/ --tmpDir /home/cmb-03/mn/yuhuang/tmp/
		-J /usr/usc/jdk/default/bin/java -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -D /Network/Data/vervet/db/
	
	#2011-11-25 on hoffman2's condor pool, need ssh tunnel for db (--needSSHDBTunnel)
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost --contigMaxRankBySize 7559 -o InspectRefSeq524WholeAlignment.xml --clusters_size 30
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -J ~/bin/jdk/bin/java --needSSHDBTunnel
	
	#2012.4.3 change tmpDir (--tmpDir) for AddOrReplaceReadGroups, no job clustering (--clusters_size 1)
	%s -a 524 -j condorpool -l condorpool -u yh -z uclaOffice
		-o workflow/InspectAlignment/InspectAln1_To_661_RefSeq524Alignments.xml --ind_aln_id_ls 1-661
		--tmpDir /Network/Data/vervet/vervetPipeline/tmp/ --clusters_size 1
	
	#2012.5.8 do perContig depth estimation (--needPerContigJob) and skip alignments with stats in db already (--skipAlignmentWithStats)
	# need ssh tunnel for db (--needSSHDBTunnel)
	# add --individual_sequence_file_raw_id_type 2 (library-specific alignments, different libraries of one individual_sequence) 
	# add --individual_sequence_file_raw_id_type 3 (both all-library-fused and library-specific alignments)
	# add "--country_id_ls 135,136,144,148,151" to limit individuals from US,Barbados,StKitts,Nevis,Gambia (AND with -S, )
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost --contigMaxRankBySize 7559
		-o workflow/InspectAlignment/InspectAln1_To_1251_RefSeq524Alignments.xml
		--ind_aln_id_ls 1-1251 --clusters_size 1
		-e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java 
		--needPerContigJob --skipAlignmentWithStats --needSSHDBTunnel
		#--individual_sequence_file_raw_id_type 2 --country_id_ls 135,136,144,148,151 --tax_id_ls 60711 #sabaeus
		#--ind_seq_id_ls 632-3230 --site_id_ls 447 --sequence_filtered 1 --excludeContaminant	#VRC sequences
		#--sequence_filtered 1 --alignment_method_id  2
	
	#2013.10.03
	mLength=100;
	%s --country_id_ls 1,129,130,131,132,133,134,136,144,148,151,152 --tax_id_ls 460675
		--alignmentDepthIntervalMethodShortName 16CynosurusRef3488MinLength$mLength
		--sequence_filtered 1 --local_realigned 1 --reduce_reads 0 --completedAlignment 1
		--excludeContaminant --ind_seq_id_ls 632-5000 -a 3488 --ref_genome_tax_id 60711
		--ref_genome_sequence_type_id 1 --ref_genome_version 1 -j hcondor -l hcondor
		-u yh -z localhost --contigMaxRankBySize 3000
		-o dags/InspectAlignment/InspectCynosurusAlignment_RefSeq3488MinLength$mLength\_AlnMethod6.xml
		--clusters_size 1 --data_dir ~/NetworkData/vervet/db/ --local_data_dir ~/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java --skipAlignmentWithStats --needSSHDBTunnel --sequence_filtered 1
		--alignment_method_id 6 --completedAlignment 1 --min_segment_length $mLength
		
Description:
	2012.3.21
		use samtools flagstat
	2011-11-4
		a pegasus workflow that inspects no-of-reads-aligned, inferred insert size and etc.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])	#, sys.argv[0], sys.argv[0]

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import File, Link, Job
from pymodule import ProcessOptions, utils, PassingData
from pymodule.pegasus import yh_pegasus
from vervet.src import VervetDB
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow
from vervet.src.pegasus.AbstractVervetAlignmentWorkflow import AbstractVervetAlignmentWorkflow
parentClass = AbstractVervetAlignmentWorkflow

class InspectAlignmentPipeline(parentClass):
	__doc__ = __doc__
	commonOptionDict = parentClass.option_default_dict.copy()
	#commonOptionDict.pop(('inputDir', 0, ))
	commonOptionDict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	
	option_default_dict = commonOptionDict.copy()
	option_default_dict.update({
						('min_segment_length', 0, int): [100, '', 1, 'a parameter of segmentation algorithm used in segmenting the depth file', ],\
						("needPerContigJob", 0, int): [0, 'P', 0, 'toggle to add DepthOfCoverage and VariousReadCount jobs for each contig.'],\
						("skipAlignmentWithStats", 0, int): [0, 's', 0, 'If an alignment has depth stats filled, not DOC job will be run. similar for flagstat job.'],\
						("alignmentDepthIntervalMethodShortName", 0, ): [None, '', 1, 'AlignmentDepthIntervalMethod.short_name, used to store segmented depth intervals from all alignments into db. This portion of workflow will not run if this is not given.'],\
						})
	#	("fractionToSample", 0, float): [0.001, '', 1, 'fraction of loci to walk through for DepthOfCoverage walker.'],\
	option_default_dict[('completedAlignment', 0, int)][0]=1	#2013.05.03

	def __init__(self, **keywords):
		"""
		2011-11-4
		"""
		parentClass.__init__(self, **keywords)
		self.no_of_alns_with_depth_jobs = 0
		self.no_of_alns_with_flagstat_jobs = 0
		
		self.alignmentDepthJobDataList=[]	#2013.08.16 use to store 
		self.needSplitChrIntervalData = False	#2013.09.01 no need for this.
	
	def addDepthOfCoverageJob(self, workflow=None, DOCWalkerJava=None, GenomeAnalysisTKJar=None,\
							refFastaFList=None, bamF=None, baiF=None, DOCOutputFnamePrefix=None,\
							fractionToSample=None, minMappingQuality=20, minBaseQuality=20, \
							parentJobLs=None, extraArguments="", \
							transferOutput=False, \
							job_max_memory = 1000, walltime=None, **keywords):
		"""
		2013.06.12
			bugfix, instead of --minBaseQuality, it was --maxBaseQuality passed to GATK.
			set minMappingQuality (was 30) to 20.
		2013.06.09
			.sample_statistics is new GATK DOC output file (replacing the .sample_interval_summary file)
			ignore argument fractionToSample, not available
		2013.05.17
			re-activate this because "samtools depth" seems to have trouble working with local-realigned and BQSR-ed bam files
			use addGATKJob()
		2012.5.7
			no longer used, superceded by addSAMtoolsDepthJob()
		2012.4.17
			add --omitIntervalStatistics and --omitLocusTable to the walker
		2012.4.12
			add "--read_filter BadCigar" to GATK to avoid stopping because of malformed Cigar
				malformed: Read starting with deletion. Cigar: 1D65M299S 
		2012.4.3
			add argument fractionToSample
		2011-11-25
		"""
		sample_summary_file = File('%s.sample_summary'%(DOCOutputFnamePrefix))
		sample_statistics_file = File('%s.sample_statistics'%(DOCOutputFnamePrefix))
		extraOutputLs = [sample_summary_file, sample_statistics_file]
		extraArgumentList = ["-o", DOCOutputFnamePrefix, \
					"-pt sample", "--read_filter BadCigar", \
					"--omitDepthOutputAtEachBase", '--omitLocusTable', '--omitIntervalStatistics']
		if minMappingQuality is not None:
			extraArgumentList.append("--minMappingQuality %s"%(minMappingQuality))
		if minBaseQuality is not None:
			extraArgumentList.append("--minBaseQuality %s"%(minBaseQuality))
		
		#if fractionToSample and fractionToSample>0 and fractionToSample<=1:
		#	extraArgumentList.append("--fractionToSample %s"%(fractionToSample))
		extraDependentInputLs = [baiF]
		job = self.addGATKJob(workflow=workflow, executable=DOCWalkerJava, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
							GATKAnalysisType="DepthOfCoverage",\
					inputFile=bamF, inputArgumentOption="-I", refFastaFList=refFastaFList, inputFileList=None,\
					interval=None, outputFile=None, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraOutputLs=extraOutputLs, \
					extraDependentInputLs=extraDependentInputLs, no_of_cpus=1, walltime=walltime, **keywords)
		job.sample_summary_file = sample_summary_file
		job.sample_statistics_file = sample_statistics_file
		return job
	
	def addSAMtoolsDepthJob(self, workflow=None, samtoolsDepth=None, samtools_path=None,\
						bamF=None, outputFile=None, baiF=None, \
						parentJobLs=None, extraOutputLs=None, job_max_memory = 500, extraArguments=None, \
						transferOutput=False, minMappingQuality=None, minBaseQuality=None, walltime=120, **keywords):
		"""
		2013.08.27 default minMappingQuality and minBaseQuality are set to None
		2013.3.24 use addGenericJob()
		2012.5.7
			
		"""
		extraArgumentList = []
		if minMappingQuality is not None:
			extraArgumentList.append("%s"%minMappingQuality)
		if minBaseQuality is not None:
			extraArgumentList.append("%s"%minBaseQuality)
		job= self.addGenericJob(executable=samtoolsDepth, \
					frontArgumentList=[samtools_path],\
					inputFile=bamF, inputArgumentOption=None,\
					outputFile=outputFile, outputArgumentOption=None,\
				parentJobLs=parentJobLs, extraDependentInputLs=[baiF], \
				extraOutputLs=extraOutputLs, extraArguments=extraArguments, \
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, \
				key2ObjectForJob=None, job_max_memory=job_max_memory, \
				sshDBTunnel=None, walltime=walltime, **keywords)
		return job
	
	def addReadCountJob(self, workflow, VariousReadCountJava=None, GenomeAnalysisTKJar=None,\
					refFastaFList=None, bamF=None, baiF=None, readCountOutputF=None,\
					parentJobLs=[], job_max_memory = 1000, extraArguments="", \
					transferOutput=False):
		"""
		2011-11-25
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		refFastaF = refFastaFList[0]
		job = Job(namespace=workflow.namespace, name=VariousReadCountJava.name, version=workflow.version)
		job.addArguments(javaMemRequirement, '-jar', GenomeAnalysisTKJar, "-T", "VariousReadCount",\
			'-R', refFastaF, '-o', readCountOutputF, "-mmq 30")
		job.addArguments("-I", bamF)
		if extraArguments:
			job.addArguments(extraArguments)
		self.addJobUse(job, file=GenomeAnalysisTKJar, transfer=True, register=True, link=Link.INPUT)
		job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
		job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
		self.registerFilesAsInputToJob(job, refFastaFList)
		job.output = readCountOutputF
		job.uses(readCountOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		"""
		#2013.3.24 should use this
		job = self.addGATKJob(self, workflow=None, executable=None, GenomeAnalysisTKJar=None, GATKAnalysisType=None,\
					inputFile=None, inputArgumentOption=None, refFastaFList=None, inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					interval=None, outputFile=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=120, **keywords)
		"""
		return job
	
	def addReformatFlagstatOutputJob(self, workflow=None, executable=None, inputF=None, \
									alignmentID=None, outputF=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, walltime=20, **keywords):
		"""
		2013.3.24 use addGenericJob()
		2012.4.3
			input is output of "samtools flagstat"
		"""
		job= self.addGenericJob(executable=executable, \
					inputFile=inputF, inputArgumentOption='-i',\
					outputFile=outputF, outputArgumentOption='-o',\
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=None, extraArguments=extraArguments, \
				transferOutput=transferOutput, \
				extraArgumentList=['-a %s'%(alignmentID)], \
				key2ObjectForJob=None, job_max_memory=job_max_memory, \
				sshDBTunnel=None, walltime=walltime, **keywords)
		return job
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
			setup additional mkdir folder jobs, before mapEachAlignment, mapEachChromosome, mapReduceOneAlignment 
		"""
		if workflow is None:
			workflow = self
		returnData = parentClass.preReduce(self, workflow=workflow, passingData=passingData, \
							transferOutput=transferOutput)
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		self.logOutputDirJob = self.addMkDirJob(outputDir="%sLog"%(passingData.outputDirPrefix))

		depthOfCoverageOutputF = File(os.path.join(reduceOutputDirJob.output, 'DepthOfCoverage.tsv'))
		passingData.depthOfCoverageOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=depthOfCoverageOutputF, parentJobLs=[reduceOutputDirJob], transferOutput=True)
		
		if self.needPerContigJob:	#need for per-contig job
			depthOfCoveragePerChrOutputF = File(os.path.join(reduceOutputDirJob.output, 'DepthOfCoveragePerChr.tsv'))
			passingData.depthOfCoveragePerChrOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=depthOfCoveragePerChrOutputF,parentJobLs=[reduceOutputDirJob], transferOutput=True)
		else:
			passingData.depthOfCoveragePerChrOutputMergeJob = None
		
		flagStatOutputF = File(os.path.join(reduceOutputDirJob.output, 'FlagStat.tsv'))
		passingData.flagStatOutputMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=flagStatOutputF, parentJobLs=[reduceOutputDirJob], transferOutput=True)
		
		passingData.alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow=workflow, alignmentDataLs=passingData.alignmentDataLs, site_handler=self.site_handler, \
					input_site_handler=self.input_site_handler, \
					addOrReplaceReadGroupsJava=self.addOrReplaceReadGroupsJava, \
					AddOrReplaceReadGroupsJar=self.AddOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, BuildBamIndexJar=self.BuildBamIndexJar, \
					mv=self.mv, data_dir=self.data_dir, tmpDir=self.tmpDir)
		
		passingData.flagStatMapFolderJob = self.addMkDirJob(outputDir="%sFlagStatMap"%(passingData.outputDirPrefix))
		
		return returnData
	
	def mapEachAlignment(self, workflow=None, alignmentData=None,  passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		flagStatMapFolderJob = passingData.flagStatMapFolderJob
		
		refFastaF = passingData.refFastaFList[0]
		
		alignment = alignmentData.alignment
		parentJobLs = alignmentData.jobLs + [passingData.fastaDictJob, passingData.fastaIndexJob]
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		
		bamFileBasenamePrefix = alignment.getReadGroup()
		
		#4X coverage alignment => 120 minutes
		realInputVolume = alignment.individual_sequence.coverage
		if realInputVolume is None:	#default is 8X
			realInputVolume = 8
		jobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=4, baseJobPropertyValue=120, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value
		#base is 4X, => 5000M
		#2013.07.26 memory usage of DOC walker java job does not depend on coverage so much, more on ref genome size.
		jobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=4, baseJobPropertyValue=14000, \
							minJobPropertyValue=16000, maxJobPropertyValue=17000).value
		
		if self.skipAlignmentWithStats and alignment.median_depth is not None and alignment.mean_depth is not None and alignment.mode_depth is not None:
			pass
		else:
			"""
			depthOutputFile = File(os.path.join(topOutputDirJob.output, '%s_DOC.tsv.gz'%(alignment.id)))
			samtoolsDepthJob = self.addSAMtoolsDepthJob(workflow, samtoolsDepth=self.samtoolsDepth, \
						samtools_path=self.samtools_path,\
						bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
						parentJobLs=[topOutputDirJob] + alignmentData.jobLs, job_max_memory = 500, extraArguments="", \
						transferOutput=False)
			self.addRefFastaJobDependency(job=samtoolsDepthJob, refFastaF=passingData.refFastaF, \
						fastaDictJob=passingData.fastaDictJob, refFastaDictF=passingData.refFastaDictF,\
						fastaIndexJob = passingData.fastaIndexJob, refFastaIndexF=passingData.refFastaIndexF)
			meanMedianModeDepthFile = File(os.path.join(topOutputDirJob.output, "%s_meanMedianModeDepth.tsv"%(alignment.id)))
			meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(\
						executable=workflow.CalculateMedianMeanOfInputColumn, \
						inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, alignmentID=alignment.id, fractionToSample=self.fractionToSample, \
						noOfLinesInHeader=0, whichColumn=2, maxNumberOfSamplings=1E6,\
						parentJobLs=[topOutputDirJob, samtoolsDepthJob], job_max_memory = 500, extraArguments=None, \
						transferOutput=False)
			"""
			#2013.05.17 samtools depth + CalculateMedianMeanOfInputColumn is not working well for realigned and BQSRed alignments
			# use GATK DOC walker
			DOCOutputFnamePrefix = os.path.join(topOutputDirJob.output, '%s_DOC'%(alignment.id))
			DOCJob = self.addDepthOfCoverageJob(DOCWalkerJava=self.DOCWalkerJava, \
						refFastaFList=passingData.refFastaFList, bamF=bamF, baiF=baiF, \
						DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
						parentJobLs=parentJobLs + [topOutputDirJob], \
						transferOutput=False,\
						job_max_memory = jobMaxMemory, walltime=jobWalltime)	#1200 minutes is 20 hours
						#fractionToSample=self.fractionToSample, \
			depthOutputFile = DOCJob.sample_statistics_file
			meanMedianModeDepthFile = File(os.path.join(topOutputDirJob.output, "%s_meanMedianModeDepth.tsv"%(alignment.id)))
			meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(\
						executable=workflow.CalculateMedianMeanOfInputColumn, \
						inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, alignmentID=alignment.id, \
						parentJobLs=[topOutputDirJob, DOCJob], job_max_memory = 500, extraArguments="--inputFileFormat=2", \
						transferOutput=False)
			
			self.addInputToStatMergeJob(workflow, statMergeJob=passingData.depthOfCoverageOutputMergeJob, inputF=meanMedianModeDepthFile,\
						parentJobLs=[meanMedianModeDepthJob])
			self.no_of_alns_with_depth_jobs += 1
			
		if self.skipAlignmentWithStats and alignment.perc_reads_mapped is not None:
			pass
		else:
			#2013.05.17 GATK's flagstat, should be identical to samtools flagstat
			"""
		java -jar ~/script/gatk2/GenomeAnalysisTK.jar -T FlagStat
			-I ~/NetworkData/vervet/db/individual_alignment/3152_640_1985088_GA_vs_3280_by_method2_realigned1_reduced0_p2312_m87.bam
			-o ~/NetworkData/vervet/db/individual_alignment/3152_640_1985088_GA_vs_3280_by_method2_realigned1_reduced0_p2312_m87.flagstat.txt
			--reference_sequence ~/NetworkData/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta
			
		output (<4 hours) looks like:
			
			1119300506 in total
			0 QC failure
			186159065 duplicates
			1034122354 mapped (92.39%)
			1119300506 paired in sequencing
			559647234 read1
			559653272 read2
			859005395 properly paired (76.74%)
			949042688 with itself and mate mapped
			85079666 singletons (7.60%)
			80245327 with mate mapped to a different chr
			26716310 with mate mapped to a different chr (mapQ>=5)
			
			"""
			
			oneFlagStatOutputF = File(os.path.join(flagStatMapFolderJob.output, '%s_flagstat.txt.gz'%(alignment.id)))
			samtoolsFlagStatJob = self.addSamtoolsFlagstatJob(executable=self.samtoolsFlagStat, \
				samtoolsExecutableFile=self.samtoolsExecutableFile, inputFile=bamF, outputFile=oneFlagStatOutputF, \
				parentJobLs=parentJobLs + [flagStatMapFolderJob], extraDependentInputLs=[baiF], transferOutput=False, \
				extraArguments=None, job_max_memory=jobMaxMemory/3, walltime=jobWalltime/2)
			self.addRefFastaJobDependency(job=samtoolsFlagStatJob, refFastaF=passingData.refFastaF, \
						fastaDictJob=passingData.fastaDictJob, refFastaDictF=passingData.refFastaDictF,\
						fastaIndexJob = passingData.fastaIndexJob, refFastaIndexF=passingData.refFastaIndexF)
			reformatFlagStatOutputF = File(os.path.join(flagStatMapFolderJob.output, '%s_flagstat.tsv'%(alignment.id)))
			reformatFlagStatOutputJob = self.addReformatFlagstatOutputJob(executable=self.ReformatFlagstatOutput, \
								inputF=oneFlagStatOutputF, alignmentID=alignment.id, outputF=reformatFlagStatOutputF, \
								parentJobLs=[flagStatMapFolderJob, samtoolsFlagStatJob], extraDependentInputLs=[], \
								transferOutput=False, \
								extraArguments=None, job_max_memory=20, walltime=30)
			self.addInputToStatMergeJob(statMergeJob=passingData.flagStatOutputMergeJob, inputF=reformatFlagStatOutputJob.output, \
						parentJobLs=[reformatFlagStatOutputJob])
			self.no_of_alns_with_flagstat_jobs += 1
		
		if alignment.path_to_depth_file is None or not os.path.isfile(os.path.join(self.data_dir, alignment.path_to_depth_file)):
			depthOutputFile = File(os.path.join(topOutputDirJob.output, '%s_depth.tsv.gz'%(alignment.id)))
			samtoolsDepthJob = self.addSAMtoolsDepthJob(samtoolsDepth=self.samtoolsDepth, \
						samtools_path=self.samtools_path,\
						bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
						parentJobLs=[topOutputDirJob] + alignmentData.jobLs, job_max_memory = 2000, extraArguments=None, \
						transferOutput=False)
			self.addRefFastaJobDependency(job=samtoolsDepthJob, refFastaF=passingData.refFastaF, \
						fastaDictJob=passingData.fastaDictJob, refFastaDictF=passingData.refFastaDictF,\
						fastaIndexJob = passingData.fastaIndexJob, refFastaIndexF=passingData.refFastaIndexF)
			
			logFile = File(os.path.join(passingData.reduceOutputDirJob.output, "%s_depth_file_2DB.log"%(alignment.id)))
			outputFileRelativePath = "%s_depth.tsv.gz"%(os.path.splitext(alignment.path)[0])
			extraArgumentList = ["--db_entry_id %s"%(alignment.id), "--tableClassName IndividualAlignment", \
								"--filePathColumnName path_to_depth_file",\
								"--fileSizeColumnName depth_file_size", \
								"--outputFileRelativePath %s"%(outputFileRelativePath), "--data_dir %s"%(self.data_dir)]
			depthFile2DBJob = self.addPutStuffIntoDBJob(executable=self.AffiliateFile2DBEntry, \
					inputFile=samtoolsDepthJob.output, inputArgumentOption='-i',\
					logFile=logFile, commit=True, \
					parentJobLs=[samtoolsDepthJob, passingData.reduceOutputDirJob], \
					extraDependentInputLs=None, transferOutput=True, extraArguments=None, \
					extraArgumentList=extraArgumentList,\
					job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
			pdata = self.constructJobDataFromJob(samtoolsDepthJob)
		else:
			alignmentDepthFile = self.registerOneInputFile(inputFname=os.path.join(self.data_dir, alignment.path_to_depth_file), \
							input_site_handler=None, folderName=self.pegasusFolderName, useAbsolutePathAsPegasusFileName=False,\
							pegasusFileName=None, checkFileExistence=True)
			pdata = PassingData(job=None, jobLs=[], file=alignmentDepthFile, fileLs=[alignmentDepthFile])
		pdata.alignment = alignment
		self.alignmentDepthJobDataList.append(pdata)
		
		if self.needPerContigJob:	#need for per-contig job
			statOutputDir = 'perContigStatOfAlignment%s'%(alignment.id)
			passingData.statOutputDirJob = self.addMkDirJob(outputDir=statOutputDir)
		else:
			passingData.statOutputDirJob = None
		
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFJobData=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		if not self.needPerContigJob:	#no need for per-contig job
			return returnData
		
		alignment = alignmentData.alignment
		
		parentJobLs = alignmentData.jobLs
		bamF = alignmentData.bamF
		baiF = alignmentData.baiF
		
		bamFnamePrefix = alignment.getReadGroup()
		
		statOutputDirJob = passingData.statOutputDirJob
		
		depthOutputFile = File(os.path.join(statOutputDirJob.output, '%s_%s_DOC.tsv.gz'%(alignment.id, chromosome)))
		samtoolsDepthJob = self.addSAMtoolsDepthJob(workflow, samtoolsDepth=self.samtoolsDepth, \
												samtools_path=self.samtools_path,\
					bamF=bamF, outputFile=depthOutputFile, baiF=baiF, \
					parentJobLs=[statOutputDirJob]+alignmentData.jobLs, job_max_memory = 500, extraArguments=None, \
					transferOutput=False)
		self.addRefFastaJobDependency(job=samtoolsDepthJob, refFastaF=passingData.refFastaF, \
					fastaDictJob=passingData.fastaDictJob, refFastaDictF=passingData.refFastaDictF,\
					fastaIndexJob = passingData.fastaIndexJob, refFastaIndexF=passingData.refFastaIndexF)
		meanMedianModeDepthFile = File(os.path.join(statOutputDirJob.output, "%s_%s_meanMedianModeDepth.tsv"%(alignment.id, chromosome)))
		meanMedianModeDepthJob = self.addCalculateDepthMeanMedianModeJob(workflow, \
					executable=workflow.CalculateMedianMeanOfInputColumn, \
					inputFile=depthOutputFile, outputFile=meanMedianModeDepthFile, alignmentID="%s-%s"%(alignment.id, chromosome), \
					parentJobLs=[samtoolsDepthJob], job_max_memory = 500, extraArguments="-r %s"%(chromosome), \
					transferOutput=False)
		
		self.addInputToStatMergeJob(workflow, statMergeJob=passingData.depthOfCoveragePerChrOutputMergeJob, \
					inputF=meanMedianModeDepthFile,\
					parentJobLs=[meanMedianModeDepthJob])
		"""
		DOCOutputFnamePrefix = os.path.join(statOutputDir, '%s_%s_DOC'%(alignment.id, chromosome))
		DOCJob = self.addDepthOfCoverageJob(workflow, DOCWalkerJava=ContigDOCWalkerJava, \
				GenomeAnalysisTKJar=GenomeAnalysisTKJar,\
				refFastaFList=refFastaFList, bamF=bamF, baiF=baiF, \
				DOCOutputFnamePrefix=DOCOutputFnamePrefix,\
				parentJobLs=[statOutputDirJob]+alignmentData.jobLs, job_max_memory = perContigJobMaxMemory*3, extraArguments="-L %s"%(chromosome),\
				transferOutput=False,\
				fractionToSample=self.fractionToSample)
		
		reduceDepthOfCoverageJob.addArguments(DOCJob.sample_statistics_file)
		reduceDepthOfCoverageJob.uses(DOCJob.sample_statistics_file, transfer=True, register=True, link=Link.INPUT)
		workflow.depends(parent=DOCJob, child=reduceDepthOfCoverageJob)
		"""
		
		"""
		#2012.4.3 no more VariousReadCountJava job
		readCountOutputF = File(os.path.join(statOutputDir, '%s_%s_variousReadCount.tsv'%(alignment.id, chromosome)))
		readCountJob = self.addReadCountJob(workflow, VariousReadCountJava=ContigVariousReadCountJava, \
					GenomeAnalysisTKJar=GenomeAnalysisTKJar, refFastaFList=refFastaFList, \
					bamF=bamF, baiF=baiF, readCountOutputF=readCountOutputF,\
					parentJobLs=statOutputDirJob, job_max_memory = perContigJobMaxMemory, extraArguments="-L %s"%(chromosome), \
					transferOutput=False)
		
		reduceVariousReadCountJob.addArguments(readCountOutputF)
		reduceVariousReadCountJob.uses(readCountOutputF, transfer=True, register=True, link=Link.INPUT)
		workflow.depends(parent=readCountJob, child=reduceVariousReadCountJob)
		"""
		
		return returnData
	
	def reduceAfterEachAlignment(self, workflow=None, passingData=None, mapEachChromosomeDataLs=None,\
								reduceAfterEachChromosomeDataLs=None,\
								transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceAfterEachChromosomeDataLs = reduceAfterEachChromosomeDataLs
		return returnData

	def reduce(self, workflow=None, passingData=None, reduceAfterEachAlignmentDataLs=None,
			transferOutput=True, **keywords):
		"""
		2013.08.14 add 2DB jobs only when their input is not empty
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.reduceAfterEachAlignmentDataLs = reduceAfterEachAlignmentDataLs
		
		reduceOutputDirJob = passingData.reduceOutputDirJob
		
		if passingData.flagStatOutputMergeJob.inputLs:	#the merge job's input is not empty or None
			flagStat2DBLogFile = File(os.path.join(reduceOutputDirJob.output, "flagStat2DB.log"))
			flagStat2DBJob = self.addPutStuffIntoDBJob(workflow, executable=self.PutFlagstatOutput2DB, \
						inputFileList=[passingData.flagStatOutputMergeJob.output], \
						logFile=flagStat2DBLogFile, commit=True, \
						parentJobLs=[reduceOutputDirJob, passingData.flagStatOutputMergeJob], \
						extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
						job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
		if passingData.depthOfCoverageOutputMergeJob.inputLs:
			DOC2DBLogFile = File(os.path.join(reduceOutputDirJob.output, "DOC2DB.log"))
			DOC2DBJob = self.addPutStuffIntoDBJob(workflow, executable=self.PutDOCOutput2DB, \
					inputFileList=[passingData.depthOfCoverageOutputMergeJob.output], \
					logFile=DOC2DBLogFile, commit=True, \
					parentJobLs=[reduceOutputDirJob, passingData.depthOfCoverageOutputMergeJob], \
					extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=10, sshDBTunnel=self.needSSHDBTunnel)
		if self.alignmentDepthJobDataList and self.alignmentDepthIntervalMethodShortName:
			if not self.min_segment_length:
				sys.stderr.write("alignmentDepthIntervalMethodShortName=%s is given but min_segment_length (%s) is not.\n"%\
								(self.alignmentDepthIntervalMethodShortName, self.min_segment_length))
				sys.exit(4)
			#2013.08.16
			alignmentIDList = [pdata.alignment.id for pdata in self.alignmentDepthJobDataList]
			alignmentIDListInStr = utils.getSuccinctStrOutOfList(alignmentIDList)
			#job to add an AlignmentDepthIntervalMethod
			logFile = File(os.path.join(self.logOutputDirJob.output, 'AddAlignmentDepthIntervalMethod2DB.log'))
			addMethod2DBJob = self.addGenericFile2DBJob(executable=self.AddAlignmentDepthIntervalMethod2DB, \
					inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", \
					data_dir=self.data_dir, logFile=logFile, commit=True,\
					parentJobLs=[self.logOutputDirJob], extraDependentInputLs=None, extraOutputLs=None, \
					transferOutput=True, extraArguments=None, \
					extraArgumentList=["--methodShortName %s"%(self.alignmentDepthIntervalMethodShortName), \
									"--alignmentIDList %s"%(alignmentIDListInStr),\
									"--min_segment_length %s"%(self.min_segment_length)], \
					job_max_memory=2000, walltime=30,  sshDBTunnel=self.needSSHDBTunnel)
			
			logFile = File(os.path.join(self.logOutputDirJob.output, 'updateMethodNoOfIntervals.log'))
			updateMethodNoOfIntervalsJob = self.addGenericFile2DBJob(executable=self.UpdateAlignmentDepthIntervalMethodNoOfIntervals, \
					data_dir=self.data_dir, logFile=logFile, commit=True,\
					parentJobLs=[self.logOutputDirJob], extraDependentInputLs=None, extraOutputLs=None, \
					transferOutput=True, extraArguments=None, \
					extraArgumentList=["--methodShortName %s"%(self.alignmentDepthIntervalMethodShortName) ], \
					job_max_memory=2000, walltime=30, sshDBTunnel=self.needSSHDBTunnel)
			
			for chromosome, chromosomeSize in self.chr2size.iteritems():
				#add a ReduceSameChromosomeAlignmentDepthFiles job
				outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_alignments_chr_%s_depth.tsv.gz'%(len(self.alignmentDepthJobDataList), chromosome)))
				reduceSameChromosomeAlignmentDepthFilesJob = self.addGenericJob(executable=self.ReduceSameChromosomeAlignmentDepthFiles, \
									inputFile=None, outputFile=outputFile, \
									parentJobLs=[reduceOutputDirJob], extraDependentInputLs=None, \
									extraArgumentList=["-w 2 --chromosomePositionColumnIndex 1 --chromosomeSize %s"%(chromosomeSize)], extraOutputLs=None,\
									transferOutput=False, \
									key2ObjectForJob=None, job_max_memory=2000, walltime=60)
				for alignmentDepthJobData in self.alignmentDepthJobDataList:
					#add a chromosome selection job
					outputFile = File(os.path.join(passingData.topOutputDirJob.output, \
												'%s_chr_%s.tsv.gz'%(utils.getFileBasenamePrefixFromPath(alignmentDepthJobData.file.name), chromosome)))
					selectRowsFromMatrixCCJob = self.addGenericJob(executable=self.SelectRowsFromMatrixCC, \
									inputFile=alignmentDepthJobData.file, outputFile=outputFile, \
									parentJobLs=alignmentDepthJobData.jobLs + [passingData.topOutputDirJob], extraDependentInputLs=None, \
									extraArgumentList=["--inputFileSortMode 1 -w 0 --whichColumnValue %s"%(chromosome)], extraOutputLs=None,\
									transferOutput=False, \
									key2ObjectForJob=None, job_max_memory=1000, walltime=60)
					self.addInputToStatMergeJob(statMergeJob=reduceSameChromosomeAlignmentDepthFilesJob, inputF=selectRowsFromMatrixCCJob.output, \
											inputArgumentOption="-i", parentJobLs=[selectRowsFromMatrixCCJob], \
											extraDependentInputLs=None)
				#add GADA job
				# add segmentation jobs to figure out intervals at similar
				outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_alignments_%s_depth_GADAOut_minSegLength%s.tsv.gz'%\
											(len(self.alignmentDepthJobDataList), chromosome, self.min_segment_length)))
				#adjust memory based on chromosome size, 135Mb => 21.4g memory
				realInputVolume = chromosomeSize
				jobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
									baseInputVolume=60000000, baseJobPropertyValue=600, \
									minJobPropertyValue=60, maxJobPropertyValue=2400).value
				#base is 135M, => 21G
				jobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
									baseInputVolume=135000000, baseJobPropertyValue=25000, \
									minJobPropertyValue=11000, maxJobPropertyValue=29000).value
				GADAJob = self.addGenericJob(executable=self.GADA, \
									inputFile=reduceSameChromosomeAlignmentDepthFilesJob.output, outputFile=outputFile, \
									parentJobLs=[reduceOutputDirJob, reduceSameChromosomeAlignmentDepthFilesJob], extraDependentInputLs=None, \
									extraArgumentList=["--MinSegLen %s"%(self.min_segment_length), '--debug -T 10 -a 0.5'], extraOutputLs=None,\
									transferOutput=False, \
									key2ObjectForJob=None, job_max_memory=jobMaxMemory, walltime=jobWalltime)
				"""
				GADAJob = self.addGenericJob(executable=self.GADA, \
									inputFile=reduceSameChromosomeAlignmentDepthFilesJob.output, outputFile=outputFile, \
									parentJobLs=[reduceOutputDirJob, reduceSameChromosomeAlignmentDepthFilesJob], extraDependentInputLs=None, \
									extraArgumentList=["-M %s"%(self.min_segment_length)], extraOutputLs=None,\
									transferOutput=False, \
									key2ObjectForJob=None, job_max_memory=10000, walltime=200)
				"""
				#job that adds AlignmentDepthIntervalFile
				logFile = File(os.path.join(self.logOutputDirJob.output, 'AddAlignmentDepthIntervalFile2DB_chr_%s.log'%(chromosome)))
				addFile2DBJob = self.addGenericFile2DBJob(executable=self.AddAlignmentDepthIntervalFile2DB, \
					inputFile=GADAJob.output, \
					inputArgumentOption="-i", \
					inputFileList=None, argumentForEachFileInInputFileList=None,\
					outputFile=None, outputArgumentOption="-o", \
					data_dir=self.data_dir, logFile=logFile, commit=True,\
					parentJobLs=[GADAJob, addMethod2DBJob, self.logOutputDirJob], \
					extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					extraArguments=None, \
					extraArgumentList=["--methodShortName %s"%(self.alignmentDepthIntervalMethodShortName), \
									"--alignmentIDList %s"%(alignmentIDListInStr), '--chromosome %s'%(chromosome),\
									"--format tsv"], \
					job_max_memory=2000, walltime=30, sshDBTunnel=self.needSSHDBTunnel)
				workflow.depends(parent=addFile2DBJob, child=updateMethodNoOfIntervalsJob)
		sys.stderr.write(" %s jobs, %s alignments with depth jobs, %s alignments with flagstat jobs.\n"%(self.no_of_jobs, \
							self.no_of_alns_with_depth_jobs, self.no_of_alns_with_flagstat_jobs))
		return returnData
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-25
			split out of run()
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		#2013.08.23
		#self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.variationSrcPath, '../bin/GADA/src/testGADA.py'), \
		#								name='GADA', clusterSizeMultipler=0.1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.variationSrcPath, '../bin/GADA/src/GADA'), \
										name='GADA', clusterSizeMultipler=0.1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/AddAlignmentDepthIntervalMethod2DB.py'), \
										name='AddAlignmentDepthIntervalMethod2DB', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/AddAlignmentDepthIntervalFile2DB.py'), \
										name='AddAlignmentDepthIntervalFile2DB', clusterSizeMultipler=0.3)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/UpdateAlignmentDepthIntervalMethodNoOfIntervals.py'), \
										name='UpdateAlignmentDepthIntervalMethodNoOfIntervals', clusterSizeMultipler=0)
		#2013.08.08
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/AffiliateFile2DBEntry.py'), \
										name='AffiliateFile2DBEntry', clusterSizeMultipler=0.1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'reducer/ReduceDepthOfCoverage.py'), \
										name='ReduceDepthOfCoverage', clusterSizeMultipler=0)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'reducer/ReduceVariousReadCount.py'), \
										name='ReduceVariousReadCount', clusterSizeMultipler=0)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
										name='ContigDOCWalkerJava', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
										name='ContigVariousReadCountJava', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "mapper/ReformatFlagstatOutput.py"), \
										name='ReformatFlagstatOutput', clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "shell/samtoolsDepth.sh"), \
										name='samtoolsDepth', clusterSizeMultipler=0.1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, "mapper/CalculateMedianModeFromSAMtoolsDepthOutput.py"), \
										name='CalculateMedianModeFromSAMtoolsDepthOutput', clusterSizeMultipler=1)
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/PutFlagstatOutput2DB.py'), \
										name='PutFlagstatOutput2DB', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.vervetSrcPath, 'db/input/PutDOCOutput2DB.py'), \
										name='PutDOCOutput2DB', clusterSizeMultipler=0)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, 'shell/pipeCommandOutput2File.sh'), \
										name='samtoolsFlagStat', clusterSizeMultipler=1)

	
if __name__ == '__main__':
	main_class = InspectAlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
