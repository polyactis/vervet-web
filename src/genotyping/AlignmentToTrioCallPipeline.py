#!/usr/bin/env python
"""
Examples:
	#2012.04.02 test-run on 40 VRC sequences (at least >1 trios)
	%s  -i 15-50 -u yh -a 524 -s 2 -z 10.8.0.10 -o AlignmentToTrioCallPipeline_VRC_ISQ_15_50_ReplicateParentsTop2Contigs.xml
		-j condorpool -l condorpool -N2 --clusters_size 1 --noOfCallingJobsPerNode1
		-S 447
	
	# 2011.12.14 on all VRC (-S 447) monkeys), top 25 contigs
	%s -u yh -a 524 -s 2 -z 10.8.0.10 -o AlignmentToTrioCallPipeline_VRC_top25Contigs.xml 
		-j condorpool -l condorpool -N25 --noOfCallingJobsPerNode 5 -S 447
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "--noOfCallingJobsPerNode 3" controls clustering of calling programs.
	# "--clusters_size 30" controls clustering for other programs., "-S 447" dictates monkeys from VRC
	%s -u yh -a 524 -s 2 -z localhost -o dags/AlignmentToCall/AlignmentToTrioCallPipeline_VRC_top7559Contigs.xml -j hcondor -l hcondor 
		-N7559 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447 --data_dir ~/NetworkData/vervet/db/
		--local_data_dir ~/NetworkData/vervet/db/ --needSSHDBTunnel
		
	# 2012.4.13 run VRC with 5 SK + 5 Nevis, sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2)
	%s -u yh -a 524 -s 2  -z localhost -o dags/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_VRC_SK_Nevis_FilteredSeq_top1000Contigs.xml 
		-j hcondor -l hcondor -N1000 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447,417,420,427,431,432,435,437,439,440,442 
		--data_dir ~/NetworkData/vervet/db/
		--local_data_dir ~/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 --needSSHDBTunnel
	
	# 2012.4.13 run 5 SK + 5 Nevis, sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2), 
	# TrioCaller fails because 10 is too small sample size.
	# 2012.6.13 supply the alignment depth stat file (-q), maxSNPMissingRate (-L 0.30), onlyKeepBiAllelicSNP (--onlyKeepBiAllelicSNP) 
	%s -u yh -a 524 -s 2 -z localhost -o dags/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_SK_Nevis_FilteredSeq_top1000Contigs.xml
		-j hcondor -l hcondor -N1000 --noOfCallingJobsPerNode 2 --clusters_size 10 -S 417,420,427,431,432,435,437,439,440,442 
		--data_dir ~/NetworkData/vervet/db/
		--local_data_dir ~/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 -q aux/alnStatForFilter.2012.6.13.tsv --onlyKeepBiAllelicSNP -L 0.30 --needSSHDBTunnel
		
	# 2012.8.15 run TrioCaller on method 14 samtools calls, contig ID from 96 to 100 (--minContigID 96 --maxContigID 100)
	# sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2), onlyKeepBiAllelicSNP (--onlyKeepBiAllelicSNP) 
	# calling job clusters size=1, others =1 (--noOfCallingJobsPerNode 1 --clusters_size 1)
	# add --notToUseDBToInferVCFNoOfLoci (not guess #loci from the 1st number in filename) if input VCF files are not db affiliated
	# 3000 (--intervalSize 3000) sites per unit, 500 overlapping between units (--intervalOverlapSize 500)
	# add "--treatEveryOneIndependent" if you want to treat everyone independent (no mendelian constraints from TrioCaller)
	%s --run_type 2 -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-u yh -a 524  -z localhost -o  dags/AlignmentToCall/TrioCallerOnMethod14Contig96_100.xml
		-j hcondor -l hcondor
		--noOfCallingJobsPerNode 1 --clusters_size 1
		--data_dir ~/NetworkData/vervet/db/ --local_data_dir ~/NetworkData/vervet/db/
		--minContigID 96 --maxContigID 100 --sequence_filtered 1 --alignment_method_id 2  --onlyKeepBiAllelicSNP --needSSHDBTunnel
		# --notToUseDBToInferVCFNoOfLoci --intervalSize 3000 --intervalOverlapSize 500 --treatEveryOneIndependent
	
	#2013.1.4 run polymutt on method 36, no overlapping between intervals, --intervalOverlapSize 0
	#(because polymutt runs on loci one by one, no dependency)
	%s --run_type 3 -I ~/NetworkData/vervet/db/genotype_file/method_36/ -u yh -a 524
		-z localhost -o dags/AlignmentToCall/PolymuttOnMethod36Contig96_100.xml
		-j hcondor -l hcondor --noOfCallingJobsPerNode 1 --clusters_size 1 
		--data_dir ~/NetworkData/vervet/db/ --local_data_dir ~/NetworkData/vervet/db/
		--minContigID 96 --maxContigID 100 --sequence_filtered 1 --alignment_method_id 2  --onlyKeepBiAllelicSNP --needSSHDBTunnel
		--intervalOverlapSize 0 --intervalSize 2000
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
		run_type 1: TrioCaller on alignments; 2: TrioCaller on VCF files; 3: polymutt on VCF files
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from pymodule import VCFFile
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow
from vervet.src.genotyping.AlignmentToCallPipeline import AlignmentToCallPipeline

parentClass=AlignmentToCallPipeline
class AlignmentToTrioCallPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetWorkflow.option_default_dict)
	option_default_dict.update(parentClass.commonCallPipelineOptionDict)
	#option_default_dict.pop(('ind_aln_id_ls', 0, ))
	#option_default_dict.pop(('ind_seq_id_ls', 0, ))
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						("onlyKeepBiAllelicSNP", 0, int): [0, 'c', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 0, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('depthFoldChange', 1, float): [2.0, 'y', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						("treatEveryOneIndependent", 0, int): [0, '', 0, 'toggle this to treat everyone in the pedigree independent and also no replicates'],\
						("run_type", 1, int): [1, 'n', 1, '1: TrioCaller on alignments; 2: TrioCaller on VCF files; 3: Polymutt on VCF files'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	#2012.8.26 run_type 2 lower the default. In the addTrioCallerJobsONVCFFiles(), these two numbers refer to the number of sites, not base pairs. 
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 1000
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		parentClass.__init__(self, **keywords)
		#self.trioCallerPath = self.trioCallerPath%self.home_path
		#self.trioCallerPath =  self.insertHomePath(self.trioCallerPath, self.home_path)
	
	
	def addTrioCallerJob(self, workflow=None, trioCallerWrapper=None, trioCallerPath=None, \
						inputVCF=None, pedFile=None, outputVCF=None, \
						inputPhased=False,\
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, walltime=None, **keywords):
		"""
		2013.06.27 added argument inputPhased
		2012.12.5 call addGenericJob()
		2012.1.20
			increase the rounds from 30 to 40
			add --burnin 20
		2011-12-4
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = [outputVCF]
		key2ObjectForJob = {}
		extraArgumentList.extend([trioCallerPath, "--vcf", inputVCF, "--pedfile", pedFile, \
						"--states 50 --rounds 40 --burnin 20",\
						"--prefix %s"%(os.path.splitext(outputVCF.name)[0])])
		if inputPhased:	#2013.06.27
			extraArgumentList.append("--inputPhased")
		else:
			extraArgumentList.append("--randomPhase")
		extraDependentInputLs.extend([inputVCF, pedFile])
		
		job = self.addGenericJob(executable=trioCallerWrapper, inputFile=None, inputArgumentOption="", \
					outputFile=None, outputArgumentOption="", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, walltime=walltime, **keywords)
		return job
	
	def addPolyMuttJob(self, workflow=None, executable=None, inputVCF=None, pedFile=None, datFile=None, \
					outputVCF=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, walltime=None, **keywords):
		"""
		2012.12.5 datFile could be empty file.
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = [outputVCF]
		key2ObjectForJob = {}
		extraArgumentList.extend(["-p", pedFile, "-d", datFile])
		extraDependentInputLs.extend([pedFile, datFile])
		
		job = self.addGenericJob(executable=executable, inputFile=inputVCF, inputArgumentOption="--in_vcf", \
					outputFile=outputVCF, outputArgumentOption="--out_vcf", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, walltime=walltime, **keywords)
		return job
	
	def addRefineSNPGenotypeJob(self, inputFile=None, vcfBaseFname=None, outputDirJob=None, statDirJob=None, \
					refFastaFList=None, intervalData=None,\
					baseInputVolume=450*2000000, realInputVolume=None,\
					parentJobLs=None, \
					transferOutput=False, \
					no_of_cpus=None, job_max_memory=2000, walltime=180, \
					max_walltime=None, **keywords):
		"""
		2013.09.02
			#memory/walltime setting
			#base for Platypus is 450X coverage in 2Mb region => 80 minutes
		"""

		
		returnData = PassingData()
		if getattr(self, "outputPedigreeJob", None) is None:
			outputFileFormat = 2	#trioCaller
			pedFile = File(os.path.join(outputDirJob.output, 'pedigree_outputFileFormat%s.txt'%(outputFileFormat)))
			sampleID2FamilyCountF = File(os.path.join(outputDirJob.output, 'sampleID2FamilyCount_outputFileFormat%s.txt'%(outputFileFormat)))
			if outputFileFormat==3:	#for polymutt
				polymuttDatFile = File(os.path.join(outputDirJob.output, 'datFile_outputFileFormat%s.txt'%(outputFileFormat)))
			else:
				polymuttDatFile = None
			self.outputPedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
					inputFile=inputFile, \
					outputFile=pedFile, sampleID2FamilyCountF=sampleID2FamilyCountF,\
					polymuttDatFile = polymuttDatFile,\
					outputFileFormat=outputFileFormat, replicateIndividualTag=self.replicateIndividualTag,\
					treatEveryOneIndependent=False,\
					parentJobLs=parentJobLs + [outputDirJob], extraDependentInputLs=None, transferOutput=True, \
					extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		SNPOnlyOutputF = File(os.path.join(outputDirJob.folder, '%s.SNP.vcf'%vcfBaseFname))
		selectSNPJob = self.addGATKJob(executable=self.SelectVariantsJava, GATKAnalysisType="SelectVariants",\
			inputFile=inputFile, inputArgumentOption="--variant", refFastaFList=refFastaFList, inputFileList=None,\
			argumentForEachFileInInputFileList=None,\
			interval=None, outputFile=SNPOnlyOutputF, \
			parentJobLs=[outputDirJob] + parentJobLs, transferOutput=False, job_max_memory=job_max_memory,\
			frontArgumentList=None, extraArguments="-selectType SNP", \
			extraArgumentList=None, extraOutputLs=None, \
			extraDependentInputLs=None, no_of_cpus=None, walltime=80)
		
		# 2013.06.11 replicate individuals who appear in more than 1 families
		round1_IndividualsReplicatedVCF = File( os.path.join(outputDirJob.folder, \
											'%s.replicate.vcf'%(vcfBaseFname)))
		replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(\
					executable=self.ReplicateVCFGenotypeColumns, \
					inputF=selectSNPJob.output, \
					sampleID2FamilyCountF=self.outputPedigreeJob.sampleID2FamilyCountF, \
					outputF=round1_IndividualsReplicatedVCF, \
					replicateIndividualTag=self.replicateIndividualTag,\
					parentJobLs=[self.outputPedigreeJob, outputDirJob, selectSNPJob], \
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
		
		trioCallerOutputFile = File(os.path.join(outputDirJob.folder, \
												'%s.trioCaller.vcf'%(vcfBaseFname)))
		trioCallerJob = self.addTrioCallerJob(trioCallerWrapper=self.trioCallerWrapper, \
				trioCallerPath=self.trioCallerPath, \
				inputVCF=replicateVCFGenotypeColumnsJob.output,\
				pedFile=self.outputPedigreeJob.output, outputVCF=trioCallerOutputFile, \
				inputPhased=False,\
				parentJobLs=[outputDirJob, replicateVCFGenotypeColumnsJob, self.outputPedigreeJob], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)	#1.2G memory for 12K loci
		
		returnData.trioCallerJob = trioCallerJob
		
		"""
		2013.07.10 the TrioCaller VCF has some info tags that are not described in VCF header
		"""
		outputFile = File(os.path.join(outputDirJob.folder, \
												'%s.extraInfoDesc.vcf'%(vcfBaseFname)))
		addInfoDescJob = self.addGenericJob(executable=self.AddMissingInfoDescriptionToVCFHeader, \
					inputFile=trioCallerJob.output, \
					inputArgumentOption="-i", \
					outputFile=outputFile, outputArgumentOption="-o", \
					parentJobLs=[outputDirJob, trioCallerJob], \
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
		
		if statDirJob:
			# a CheckGenotypeConcordanceAmongReplicates.py job
			trioCallerReplicateConcordanceFile = File(os.path.join(self.statDirJob.folder, \
									'%s.trioCaller.concordance.tsv'%(vcfBaseFname)))
			returnData.trioCallerReplicateConcordanceJob = self.addGATKJob(executable=self.CalculateConcordanceJava, \
						GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
						GATKAnalysisType="CalculateConcordanceAmongReplicates",\
						inputFile=trioCallerJob.output, inputArgumentOption="--variant", \
						refFastaFList=self.registerReferenceData.refFastaFList, \
						interval=None, \
						outputFile=trioCallerReplicateConcordanceFile, outputArgumentOption="--concordanceStatFname",\
						frontArgumentList=None, extraArguments="--replicateIndividualTag %s"%(self.replicateIndividualTag), \
						extraArgumentList=None, extraOutputLs=None, \
						parentJobLs=[self.statDirJob, trioCallerJob], \
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
		mergeReplicateOutputF = File(os.path.join(outputDirJob.folder, \
									'%s.replicatesMerged.vcf'%(vcfBaseFname)))
		mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(\
							executable=self.MergeVCFReplicateHaplotypesJava,\
							GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
							inputF=addInfoDescJob.output, outputF=mergeReplicateOutputF, \
							replicateIndividualTag=self.replicateIndividualTag, \
							refFastaFList=self.registerReferenceData.refFastaFList, \
							parentJobLs=[outputDirJob, addInfoDescJob], \
							extraDependentInputLs=[], transferOutput=False, \
							extraArguments=None, \
							analysis_type='MergeVCFReplicateGenotypeColumns',\
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=5000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		returnData.refineGenotypeJob = mergeVCFReplicateColumnsJob	#the final gentoype job
		returnData.refineGenotypeJob.intervalData = intervalData	#attached so that it could be used by downstream jobs
		return returnData
	
	def addTrioCallerJobsONVCFFiles(self, workflow=None, alignmentLs=None, inputData=None, samtools=None, \
				genotyperJava=None, SelectVariantsJava=None, GenomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				MergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexJar=None,\
				mv=None, CallVariantBySamtools=None,\
				trioCallerPath=None, trioCallerWrapper=None, \
				replicateIndividualTag="copy", treatEveryOneIndependent=False,\
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None, ligateVcf=None, ligateVcfExecutableFile=None,\
				registerReferenceData=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, intervalOverlapSize=100000, site_type=1, data_dir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", \
				maxSNPMissingRate=None, alnStatForFilterF=None, onlyKeepBiAllelicSNP=True, \
				cumulativeMedianDepth=5000, job_max_memory = 2000, vcf_job_max_memory = 1000,\
				run_type=2, transferOutput=True, **keywords):
		"""
		2012.12.5 added argument run_type (same as self.run_type) 2: TrioCaller; 3: polymutt
		2012.8.15
		"""
		sys.stderr.write("Adding trioCaller jobs for  %s vcf files ..."%(len(inputData.jobDataLs)))
		if workflow is None :
			workflow = self
		refFastaFList = registerReferenceData.refFastaFList
		refFastaF = refFastaFList[0]
		
		if needFastaDictJob or registerReferenceData.needPicardFastaDictJob:
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
												refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = registerReferenceData.refPicardFastaDictF
		
		if needFastaIndexJob or registerReferenceData.needSAMtoolsFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = registerReferenceData.refSAMtoolsFastaIndexF
		
		trioCallerOutputDir = "%sRefinedCalls"%(outputDirPrefix)
		trioCallerOutputDirJob = self.addMkDirJob(outputDir=trioCallerOutputDir)
		round1CallDir = "%sPreRefinedCalls"%(outputDirPrefix)
		round1CallDirJob = self.addMkDirJob(outputDir=round1CallDir)
		
		outputPedigreeJob = None
		
		# add merge jobs for every reference
		returnData = PassingData()
		returnData.jobDataLs = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputF = jobData.vcfFile
			inputFBaseName = os.path.basename(inputF.name)
			chr_id = self.getChrFromFname(inputFBaseName)
			commonPrefix = inputFBaseName.split('.')[0]
			
			overlapInterval = chr_id
			#split VCF job
			outputFnamePrefix = os.path.join(round1CallDirJob.folder, '%s_splitVCF'%commonPrefix)
			splitVCFJob = self.addSplitVCFFileJob(executable=self.SplitVCFFile, inputFile=inputF, outputFnamePrefix=outputFnamePrefix, \
					noOfOverlappingSites=intervalOverlapSize, noOfSitesPerUnit=intervalSize, noOfTotalSites=inputF.noOfLoci, \
					parentJobLs=jobData.jobLs+[round1CallDirJob], \
					extraDependentInputLs=[jobData.tbi_F], \
					extraArguments=None, transferOutput=False, job_max_memory=job_max_memory)
			
			#ligate vcf job (different segments of a chromosome into one chromosome)
			concatTrioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.vcf'%chr_id)
			concatTrioCallerOutputF = File(concatTrioCallerOutputFname)
			trioCallerWholeContigConcatJob = self.addLigateVcfJob(executable=ligateVcf, ligateVcfExecutableFile=ligateVcfExecutableFile, \
										outputFile=concatTrioCallerOutputF, \
										parentJobLs=[trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
										extraArguments=None, job_max_memory=vcf_job_max_memory)
			
			#bgzip and tabix the trio caller output
			bgzip_concatTrioCallerOutputF = File("%s.gz"%concatTrioCallerOutputFname)
			bgzip_tabix_concatTrioCallerOutput_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
					parentJob=trioCallerWholeContigConcatJob, inputF=concatTrioCallerOutputF, outputF=bgzip_concatTrioCallerOutputF, \
					transferOutput=transferOutput)
			
			returnData.jobDataLs.append(PassingData(vcfFile=bgzip_concatTrioCallerOutputF, jobLs=[bgzip_tabix_concatTrioCallerOutput_job]))
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			
			noOfUnits = max(1, utils.getNoOfUnitsNeededToCoverN(N=inputF.noOfLoci, s=intervalSize, o=intervalOverlapSize)-1)
			for unitNumber in xrange(1, noOfUnits+1):
				splitVCFFile = getattr(splitVCFJob, 'unit%sFile'%(unitNumber))
				
				#2012.4.2
				tranferIntermediateFilesForDebug=False
				overlapIntervalFnameSignature = '%s_%s'%(commonPrefix, unitNumber)
				
				#selectVariants would generate AC, AF so that TrioCaller could read it. (samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
				round1_VCF4OutputFname = os.path.join(round1CallDirJob.folder, '%s.niceformat.vcf'%overlapIntervalFnameSignature)
				round1_VCF4OutputF = File(round1_VCF4OutputFname)
				round1_vcf_convert_job = self.addSelectVariantsJob(SelectVariantsJava=SelectVariantsJava, \
						inputF=splitVCFFile, outputF=round1_VCF4OutputF, \
						refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob, splitVCFJob], \
						extraDependentInputLs=None, transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=job_max_memory, interval=overlapInterval)
				
				if outputPedigreeJob is None:	#2013.1.2
					outputFileFormat = run_type
					pedFile = File(os.path.join(trioCallerOutputDirJob.output, 'pedigree_outputFileFormat%s.txt'%(outputFileFormat)))
					sampleID2FamilyCountF = File(os.path.join(trioCallerOutputDirJob.output, 'sampleID2FamilyCount_outputFileFormat%s.txt'%(outputFileFormat)))
					if outputFileFormat==3:	#for polymutt
						polymuttDatFile = File(os.path.join(trioCallerOutputDirJob.output, 'datFile_outputFileFormat%s.txt'%(outputFileFormat)))
					else:
						polymuttDatFile = None
					outputPedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
							inputFile=round1_VCF4OutputF, outputFile=pedFile, sampleID2FamilyCountF=sampleID2FamilyCountF,\
							polymuttDatFile = polymuttDatFile,\
							outputFileFormat=outputFileFormat, replicateIndividualTag=replicateIndividualTag,\
							treatEveryOneIndependent=treatEveryOneIndependent,\
							parentJobLs=[round1_vcf_convert_job, trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=True, \
							extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
				
				#2012.4.2 replicate individuals who appear in more than 1 families
				round1_IndividualsReplicatedVCF = File( os.path.join(round1CallDirJob.folder, \
													'%s.replicate.vcf'%overlapIntervalFnameSignature))
				replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(workflow, \
							executable=workflow.ReplicateVCFGenotypeColumns, inputF=round1_VCF4OutputF, \
							sampleID2FamilyCountF=outputPedigreeJob.sampleID2FamilyCountF, outputF=round1_IndividualsReplicatedVCF, \
							replicateIndividualTag=replicateIndividualTag,\
							parentJobLs=[round1_vcf_convert_job, outputPedigreeJob], extraDependentInputLs=[], \
							transferOutput=tranferIntermediateFilesForDebug, \
							extraArguments=None, job_max_memory=vcf_job_max_memory)
				
				if run_type==2:
					#TrioCaller job
					refineGenotypeOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature))
					refineGenotypeJob = self.addTrioCallerJob(workflow, trioCallerWrapper=trioCallerWrapper, trioCallerPath=trioCallerPath, \
							inputVCF=round1_IndividualsReplicatedVCF,\
							pedFile=outputPedigreeJob.output, outputVCF=refineGenotypeOutputF, \
							parentJobLs=[trioCallerOutputDirJob, replicateVCFGenotypeColumnsJob, outputPedigreeJob], \
							extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
							extraArguments=None, job_max_memory=vcf_job_max_memory)	#1.2G memory for 12K loci
				elif run_type==3:
					#add polymutt jobs
					#TrioCaller job
					refineGenotypeOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.polymutt.vcf'%overlapIntervalFnameSignature))
					refineGenotypeJob = self.addPolyMuttJob(executable=self.polymutt, inputVCF=round1_IndividualsReplicatedVCF, \
										pedFile=outputPedigreeJob.output, datFile=outputPedigreeJob.polymuttDatFile, \
										outputVCF=refineGenotypeOutputF, parentJobLs=[trioCallerOutputDirJob, replicateVCFGenotypeColumnsJob, outputPedigreeJob], \
										extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
										extraArguments=None, job_max_memory=vcf_job_max_memory*2)
				
				"""
				2013.07.10 the TrioCaller VCF has some info tags that are not described in VCF header
				"""
				outputFile = File(os.path.join(trioCallerOutputDirJob.folder, \
														'%s.extraInfoDesc.vcf'%(overlapIntervalFnameSignature)))
				addInfoDescJob = self.addGenericJob(executable=self.AddMissingInfoDescriptionToVCFHeader, \
							inputFile=refineGenotypeJob.output, \
							inputArgumentOption="-i", \
							outputFile=outputFile, outputArgumentOption="-o", \
							parentJobLs=[trioCallerOutputDirJob]+ refineGenotypeJob, \
							extraDependentInputLs=None, extraOutputLs=None, \
							frontArgumentList=None, extraArguments=None, extraArgumentList=None, \
							transferOutput=False, sshDBTunnel=None, \
							key2ObjectForJob=None, objectWithDBArguments=None, \
							no_of_cpus=None, 
							job_max_memory=vcf_job_max_memory/4, \
							walltime=None,\
							max_walltime=None)
				#2012.4.2
				mergeReplicateOutputF = File(os.path.join(trioCallerOutputDirJob.folder, \
												'%s.noReplicate.vcf'%overlapIntervalFnameSignature))
				noOfAlignments= len(alignmentLs)
				
				no_of_loci =  getattr(inputF, 'no_of_loci', None)
				if no_of_loci is None:
					no_of_loci = 26000.	#assume 26K
				memoryRequest = min(10000, max(4000, int(6000*(noOfAlignments/323.0)*(no_of_loci/26000.0))) )
					#extrapolates from 4000Mb memory for a 323-sample with 26K loci, (2.6 Megabase)
					#upper bound is 10g. lower bound is 4g.
				mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(workflow, \
									executable=self.MergeVCFReplicateHaplotypesJava,\
									GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
									inputF=addInfoDescJob.output, outputF=mergeReplicateOutputF, \
									replicateIndividualTag=replicateIndividualTag, \
									refFastaFList=refFastaFList, parentJobLs=[addInfoDescJob], \
									extraDependentInputLs=None, transferOutput=tranferIntermediateFilesForDebug, \
									extraArguments=None, job_max_memory=memoryRequest,\
									analysis_type='MergeVCFReplicateGenotypeColumns')
				
				#convert to vcf4 so that other vcftools software could be used.
				vcf4_trioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.noreplicte.v4.vcf'%overlapIntervalFnameSignature)
				vcf4_trioCallerOutputF = File(vcf4_trioCallerOutputFname)
				vcf_convert_TrioCallerOutputJob = self.addVCFFormatConvertJob(vcf_convert=vcf_convert, \
							parentJob=mergeVCFReplicateColumnsJob, inputF=mergeVCFReplicateColumnsJob.output, \
							outputF=vcf4_trioCallerOutputF, transferOutput=False)
				
				
				#bgzip and tabix the trio caller output
				trioGzipOutputF = File("%s.gz"%vcf4_trioCallerOutputFname)
				trioGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_trioCallerOutputFname)
				bgzip_tabix_trioOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=vcf_convert_TrioCallerOutputJob, inputF=vcf4_trioCallerOutputF, outputF=trioGzipOutputF, \
						transferOutput=False)
				
				lisOfJobsThatNeedRefIndexFiles = [round1_vcf_convert_job, replicateVCFGenotypeColumnsJob, \
												mergeVCFReplicateColumnsJob]
				for job in lisOfJobsThatNeedRefIndexFiles:
					self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
				#add this output to the union job
				# 2012.6.1 done it through addInputToStatMergeJob()
				self.addInputToStatMergeJob(statMergeJob=trioCallerWholeContigConcatJob, inputF=trioGzipOutputF, \
							parentJobLs=[bgzip_tabix_trioOutputF_job], \
							extraDependentInputLs=[trioGzipOutput_tbi_F])
				
		sys.stderr.write(" %s jobs. \n"%(self.no_of_jobs))
		return returnData
	
	def addReplicateVCFGenotypeColumnsJob(self, workflow=None, executable=None, inputF=None, \
					sampleID2FamilyCountF=None, outputF=None, \
					replicateIndividualTag=None,\
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, walltime=None, **keywords):
		"""
		2012.12.5 call addGenericJob()
		2012.4.2
			a job that replicates the genotype columns of individuals who appear in >1 families.
			This is to be run before TrioCaller is applied.
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		extraArgumentList.extend(['--sampleID2FamilyCountFname', sampleID2FamilyCountF, \
								'--replicateIndividualTag', replicateIndividualTag])
		extraDependentInputLs.append(sampleID2FamilyCountF)
		
		job = self.addGenericJob(executable=executable, inputFile=inputF, inputArgumentOption="-i", \
					outputFile=outputF, outputArgumentOption="-o", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, walltime=walltime, **keywords)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow is None:
			workflow = self
		parentClass.registerCustomExecutables(self, workflow)
		
	
	def run(self):
		"""
		2011-7-11
		"""
		if self.run_type!=1:
			self.needSplitChrIntervalData = False	#2013.06.21 turn this off before setup_run() to not construct chr2IntervalDataLs
		pdata = self.setup_run()
		workflow = pdata.workflow
		db_vervet = self.db
		
		if self.run_type==1:
			alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, data_dir=self.local_data_dir)
		elif self.run_type in [2,3]:
			inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
									checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
									pegasusFolderName=self.pegasusFolderName,\
									maxContigID=self.maxContigID, \
									minContigID=self.minContigID,  db_vervet=db_vervet, \
									needToKnowNoOfLoci=abs(1-self.notToKnowNoOfLoci),\
									minNoOfLociInVCF=self.minNoOfLociInVCF)	#ignore files with too few loci
			inputF = inputData.jobDataLs[0].vcfFile
			vcfFile = VCFFile(inputFname=inputF.abspath)
			alignmentLs = db_vervet.getAlignmentsFromVCFSampleIDList(vcfFile.getSampleIDList())
			del vcfFile
		
		alignmentLs = db_vervet.filterAlignments(alignmentLs=alignmentLs, sequence_filtered=self.sequence_filtered, \
												individual_site_id_set=set(self.site_id_ls),\
												mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls),\
									excludeContaminant=self.excludeContaminant)
		cumulativeMedianDepth = db_vervet.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
															defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		
		registerReferenceData = pdata.registerReferenceData
		
		
		if self.run_type==1:
			alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow=workflow, alignmentLs=alignmentLs, data_dir=self.data_dir)
			chr2size = self.chr2size
			#chr2size = set(['Contig149'])	#temporary when testing Contig149
			#chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
			chrLs = chr2size.keys()
			chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=chr2size, \
															intervalSize=self.intervalSize, \
															intervalOverlapSize=self.intervalOverlapSize)
			# 2012.8.2 if maxContigID/minContigID is not well defined. restrictContigDictionry won't do anything.
			chr2IntervalDataLs = self.restrictContigDictionry(dc=chr2IntervalDataLs, \
													maxContigID=self.maxContigID, minContigID=self.minContigID)
			#2012.6.12
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=30)	#minGQ doesn't matter anymore.
			alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
			
			self.addGenotypeCallJobs(workflow=workflow, alignmentDataLs=alignmentDataLs, chr2IntervalDataLs=chr2IntervalDataLs, \
						registerReferenceData=registerReferenceData, \
						site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
						site_type=self.site_type, data_dir=self.data_dir,\
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, maxSNPMissingRate=self.maxSNPMissingRate,\
						alnStatForFilterF=alnStatForFilterF, cumulativeMedianDepth=cumulativeMedianDepth,\
						transferOutput=True)
		elif self.run_type in [2, 3]:
			self.addTrioCallerJobsONVCFFiles(workflow=workflow, alignmentLs=alignmentLs, inputData=inputData, \
						samtools=workflow.samtools, \
						genotyperJava=workflow.genotyperJava,  SelectVariantsJava=workflow.SelectVariantsJava, \
						GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
						addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, AddOrReplaceReadGroupsJar=workflow.AddOrReplaceReadGroupsJar, \
						CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
						MergeSamFilesJar=workflow.MergeSamFilesJar, \
						BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
						mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools, \
						trioCallerPath=self.trioCallerPath, trioCallerWrapper=workflow.trioCallerWrapper, \
						replicateIndividualTag=self.replicateIndividualTag, treatEveryOneIndependent=self.treatEveryOneIndependent,\
						bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, \
						vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
						concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
						ligateVcf=self.ligateVcf, ligateVcfExecutableFile=self.ligateVcfExecutableFile,\
						registerReferenceData=registerReferenceData, \
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
						site_type=self.site_type, data_dir=self.data_dir,\
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, maxSNPMissingRate=self.maxSNPMissingRate,\
						alnStatForFilterF=None, cumulativeMedianDepth=cumulativeMedianDepth,\
						run_type=self.run_type, transferOutput=True)
		
		self.end_run()

if __name__ == '__main__':
	main_class = AlignmentToTrioCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
