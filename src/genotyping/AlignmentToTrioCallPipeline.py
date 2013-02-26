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
		-N7559 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447 -e /u/home/eeskin/polyacti/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --needSSHDBTunnel
		
	# 2012.4.13 run VRC with 5 SK + 5 Nevis, sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2)
	%s -u yh -a 524 -s 2  -z localhost -o dags/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_VRC_SK_Nevis_FilteredSeq_top1000Contigs.xml 
		-j hcondor -l hcondor -N1000 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447,417,420,427,431,432,435,437,439,440,442 
		-e /u/home/eeskin/polyacti/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 --needSSHDBTunnel
	
	# 2012.4.13 run 5 SK + 5 Nevis, sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2), 
	# TrioCaller fails because 10 is too small sample size.
	# 2012.6.13 supply the alignment depth stat file (-q), maxSNPMissingRate (-L 0.30), onlyKeepBiAllelicSNP (--onlyKeepBiAllelicSNP) 
	%s -u yh -a 524 -s 2 -z localhost -o dags/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_SK_Nevis_FilteredSeq_top1000Contigs.xml
		-j hcondor -l hcondor -N1000 --noOfCallingJobsPerNode 2 --clusters_size 10 -S 417,420,427,431,432,435,437,439,440,442 
		-e /u/home/eeskin/polyacti/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--sequence_filtered 1 --alignment_method_id 2 -q aux/alnStatForFilter.2012.6.13.tsv --onlyKeepBiAllelicSNP -L 0.30 --needSSHDBTunnel
		
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
	
	#2013.1.4 run polymutt on method 36, no overlapping between intervals, --intervalOverlapSize 0
	#(because polymutt runs on loci one by one, no dependency)
	%s --run_type 3 -I ~/NetworkData/vervet/db/genotype_file/method_36/ -u yh -a 524
		-z localhost -o dags/AlignmentToCall/PolymuttOnMethod36Contig96_100.xml
		-j hcondor -l hcondor --noOfCallingJobsPerNode 1 --clusters_size 1 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
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

class AlignmentToTrioCallPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetWorkflow.option_default_dict)
	option_default_dict.update(AlignmentToCallPipeline.commonCallPipelineOptionDict)
	#option_default_dict.pop(('ind_aln_id_ls', 0, ))
	#option_default_dict.pop(('ind_seq_id_ls', 0, ))
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						("onlyKeepBiAllelicSNP", 0, int): [0, 'c', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 0, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('depthFoldChange', 1, float): [2.0, 'y', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						("notToKnowNoOfLoci", 0, int): [0, 'Y', 0, 'toggle this to not guess the number of loci in one VCF file (for TrioCaller) from the 1st number in its filename'],\
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
		AlignmentToCallPipeline.__init__(self, **keywords)
		#self.trioCallerPath = self.trioCallerPath%self.home_path
		#self.trioCallerPath =  self.insertHomePath(self.trioCallerPath, self.home_path)
	
	
	def addTrioCallerJob(self, workflow=None, trioCallerWrapper=None, trioCallerPath=None, inputVCF=None, pedFile=None, outputVCF=None, \
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
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
						"--states 50 --randomPhase --rounds 40 --burnin 20",\
						"--prefix %s"%(os.path.splitext(outputVCF.name)[0])])
		extraDependentInputLs.extend([inputVCF, pedFile])
		
		job = self.addGenericJob(executable=trioCallerWrapper, inputFile=None, inputArgumentOption="", \
					outputFile=None, outputArgumentOption="", \
					parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def addPolyMuttJob(self, workflow=None, executable=None, inputVCF=None, pedFile=None, datFile=None, \
					outputVCF=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, **keywords):
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
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def addGenotypeCallJobs(self, workflow=None, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				genotyperJava=None, SelectVariantsJava=None, GenomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				MergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexJar=None,\
				mv=None, CallVariantBySamtools=None,\
				trioCallerPath=None, trioCallerWrapper=None, \
				replicateIndividualTag="copy", treatEveryOneIndependent=False,\
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None, ligateVcf=None, ligateVcfPerlPath=None,\
				refFastaFList=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, intervalOverlapSize=100000, site_type=1, data_dir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", \
				maxSNPMissingRate=None, alnStatForFilterF=None, onlyKeepBiAllelicSNP=True, \
				cumulativeMedianDepth=5000, job_max_memory = 2000, vcf_job_max_memory = 1000,\
				transferOutput=True, **keywords):
		"""
		2012.8.15
			use chr2IntervalDataLs as guide to partition jobs
		2012.6.12
			add argument maxSNPMissingRate, alnStatForFilterF, onlyKeepBiAllelicSNP to filter SNPs after first round.
		2012.1.9
			add outputDirPrefix to differentiate one run from another if multiple trio call workflows are run simultaneously
			use alignmentDataLs instead of alignmentLs
		2011-9-22
			add argument concatGATK, concatSamtools.
		2011-9-15
			bamListF is now useless. samtools_job could accept variable-length list of bam input files
		2011-9-14
			argument intervalSize determines how many sites gatk/samtools works on at a time
		"""
		sys.stderr.write("Adding genotype call jobs for %s chromosomes/contigs ..."%(len(chr2IntervalDataLs)))
		refFastaF = refFastaFList[0]
		no_of_jobs = 0
		
		if needFastaDictJob:	# the .dict file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
												CreateSequenceDictionaryJar=CreateSequenceDictionaryJar, refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
			no_of_jobs += 1
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
			no_of_jobs += 1
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
		trioCallerOutputDir = "%strioCaller"%(outputDirPrefix)
		trioCallerOutputDirJob = self.addMkDirJob(outputDir=trioCallerOutputDir)
		round1CallDir = "%spreTrioCaller"%(outputDirPrefix)
		round1CallDirJob = self.addMkDirJob(outputDir=round1CallDir)
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow, alignmentDataLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexJar=BuildBamIndexJar, \
					mv=mv, namespace=namespace, version=version, data_dir=data_dir)
		
		# add merge jobs for every reference
		returnData = PassingData()
		returnData.jobDataLs = []
		outputPedigreeJob = None	#2013.1.4
		for chr, intervalDataLs in chr2IntervalDataLs.iteritems():
			#reduce the number of chunks 1 below needed. last trunk to reach the end of contig
			#however set it to 1 for contigs smaller than intervalSize 	
			concatTrioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.vcf'%chr)
			concatTrioCallerOutputF = File(concatTrioCallerOutputFname)
			trioCallerWholeContigConcatJob = self.addLigateVcfJob(executable=ligateVcf, ligateVcfPerlPath=ligateVcfPerlPath, \
										outputFile=concatTrioCallerOutputF, \
										parentJobLs=[trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
										extraArguments=None, job_max_memory=vcf_job_max_memory)
			"""
			wholeContigConcatJob = self.addVCFConcatJob(workflow, concatExecutable=vcf_concat, \
							parentDirJob=trioCallerOutputDirJob, \
							outputF=concatTrioCallerOutputF, namespace=namespace, version=version, transferOutput=transferOutput, \
							vcf_job_max_memory=vcf_job_max_memory)
			"""
			#bgzip and tabix the trio caller output
			bgzip_concatTrioCallerOutputF = File("%s.gz"%concatTrioCallerOutputFname)
			bgzip_concatTrioCallerOutput_tbi_F = File("%s.gz.tbi"%concatTrioCallerOutputFname)
			bgzip_tabix_concatTrioCallerOutput_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
					parentJob=trioCallerWholeContigConcatJob, inputF=concatTrioCallerOutputF, outputF=bgzip_concatTrioCallerOutputF, \
					transferOutput=transferOutput)
			
			returnData.jobDataLs.append(PassingData(vcfFile=bgzip_concatTrioCallerOutputF, jobLs=[bgzip_tabix_concatTrioCallerOutput_job]))
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			
			#2011-9-22 union of all 1st-rounding calling intervals for one contig
			round1_VCFConcatOutputFname = os.path.join(round1CallDirJob.folder, '%s.vcf.gz'%chr)
			round1_VCFConcatOutputF = File(round1_VCFConcatOutputFname)
			round1_VCFConcatJob = self.addVCFConcatJob(workflow, concatExecutable=concatGATK, parentDirJob=round1CallDirJob, \
							outputF=round1_VCFConcatOutputF, namespace=namespace, version=version, transferOutput=transferOutput, \
							vcf_job_max_memory=vcf_job_max_memory)
			
			no_of_jobs += 3
			
			for intervalData in intervalDataLs:
				if intervalData.file:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.file
				else:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.interval
				intervalFnameSignature = intervalData.intervalFnameSignature
				overlapInterval = intervalData.overlapInterval
				overlapIntervalFnameSignature = intervalData.overlapIntervalFnameSignature
				overlapStart = intervalData.overlapStart
				overlapStop = intervalData.overlapStop
				
				
				#1st-round genotype calling
				round1CallOutputFname = os.path.join(round1CallDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature)
				round1CallOutputF = File(round1CallOutputFname)
				indelVCFOutputFname = "%s.indel.vcf"%(round1CallOutputFname)
				indelVCFOutputF = File(indelVCFOutputFname)
				preTrioCallerCallJob = self.addSAMtoolsCallJob(workflow, CallVariantBySamtools=CallVariantBySamtools, \
					samtoolsOutputF=round1CallOutputF, indelVCFOutputF=indelVCFOutputF, \
					refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=vcf_job_max_memory, site_type=site_type, mpileupInterval=mpileupInterval,\
					bcftoolsInterval=bcftoolsInterval, maxDP=cumulativeMedianDepth*5)
				
				#2012.6.12 filter via depth
				vcf1AfterDepthFilter = File(os.path.join(round1CallDirJob.folder, '%s.depthFiltered.vcf'%(overlapIntervalFnameSignature)))
				vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
						GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
						refFastaFList=refFastaFList, inputVCFF=round1CallOutputF, outputVCFF=vcf1AfterDepthFilter, \
						parentJobLs=[preTrioCallerCallJob], \
						alnStatForFilterF=alnStatForFilterF, \
						extraDependentInputLs=[], \
						onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
				
				#convert to vcf4 so that vcftools (it only recognizes VCF4) could be used.
				round1_VCF4OverlapOutputFname = os.path.join(round1CallDirJob.folder, '%s.depthFiltered.v4.vcf'%overlapIntervalFnameSignature)
				round1_VCF4OverlapOutputF = File(round1_VCF4OverlapOutputFname)
				round1OverlapVCFconvert_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=vcf1FilterByDepthJob, inputF=vcf1FilterByDepthJob.output, outputF=round1_VCF4OverlapOutputF, \
							transferOutput=False)
				
				#2012.6.12 filter 1st-round calls via missing percentage
				outputFnamePrefix = os.path.join(round1CallDirJob.folder, '%s.filter_by_vcftools'%(overlapIntervalFnameSignature))
				vcf1FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=round1_VCF4OverlapOutputF, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[round1OverlapVCFconvert_job], \
						snpMisMatchStatFile=None, \
						minMAC=None, minMAF=None, \
						maxSNPMissingRate=maxSNPMissingRate,\
						extraDependentInputLs=[], extraArguments="--recode-INFO-all")
				
				
				"""
				#GATK as 1st round caller
				#GATK produces "./." for missing genotype and TrioCaller has trouble passing that.
				round1CallOutputFname = os.path.join(round1CallDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature)
				round1CallOutputF = File(round1CallOutputFname)
				gatkIDXOutputFname = os.path.join(round1CallDirJob.folder, '%s.vcf.idx'%(overlapIntervalFnameSignature))
				gatkIDXOutput = File(gatkIDXOutputFname)
				
				preTrioCallerCallJob= self.addGATKCallJob(workflow, genotyperJava=genotyperJava, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
						round1CallOutputF=round1CallOutputF, gatkIDXOutput=gatkIDXOutput, refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob], \
						extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
						job_max_memory=job_max_memory, no_of_gatk_threads=no_of_gatk_threads, site_type=site_type, \
						interval=overlapInterval)
				"""
				# following few steps is prepare 1st-round calls to be concatenated (not for 2nd-round caller)
				#select the variants to get rid of overlap, so that union of whole contig doesn't have overlap
				round1_NonOverlapOutputF = File(os.path.join(round1CallDirJob.folder, '%s.nonoverlap.vcf'%intervalFnameSignature))
				round1SelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						GenomeAnalysisTKJar=GenomeAnalysisTKJar, inputF=vcf1FilterByvcftoolsJob.output, outputF=round1_NonOverlapOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[vcf1FilterByvcftoolsJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=mpileupInterval)
				
				#convert to vcf4 so that other vcftools software could be used.
				round1_VCF4NonOverlapOutputFname = os.path.join(round1CallDirJob.folder, '%s.v4.vcf'%intervalFnameSignature)
				round1_VCF4NonOverlapOutputF = File(round1_VCF4NonOverlapOutputFname)
				round1NonOverlapVCFconvert_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=round1SelectVariantJob, inputF=round1SelectVariantJob.output, outputF=round1_VCF4NonOverlapOutputF, \
							namespace=namespace, version=version, transferOutput=False)
				
				round1VCFGzipOutputF = File("%s.gz"%round1_VCF4NonOverlapOutputFname)
				round1VCFGzipOutput_tbi_F = File("%s.gz.tbi"%round1_VCF4NonOverlapOutputFname)
				round1_bgzip_tabix_VCFOutputJOb = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=round1NonOverlapVCFconvert_job, inputF=round1_VCF4NonOverlapOutputF, outputF=round1VCFGzipOutputF, \
						namespace=namespace, version=version, transferOutput=False)
				
				#add this output to a GATK union job
				# 2012.6.1 done it through addInputToStatMergeJob()
				self.addInputToStatMergeJob(statMergeJob=round1_VCFConcatJob, inputF=round1VCFGzipOutputF, \
							parentJobLs=[round1_bgzip_tabix_VCFOutputJOb], \
							extraDependentInputLs=[round1VCFGzipOutput_tbi_F])
				
				"""
				#convert to vcf4 so that TrioCaller could read it. (samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
				round1_VCF4OutputFname = os.path.join(round1CallDirJob.folder, '%s.v4.vcf'%overlapIntervalFnameSignature)
				round1_VCF4OutputF = File(round1_VCF4OutputFname)
				round1_vcf_convert_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=preTrioCallerCallJob, inputF=preTrioCallerCallJob.output, outputF=round1_VCF4OutputF, \
							namespace=namespace, version=version, transferOutput=False)
				"""
				
				#2012.4.2
				tranferIntermediateFilesForDebug=False
				
				#selectVariants would generate AC, AF so that TrioCaller could read it. (samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
				round1_VCF4OutputFname = os.path.join(round1CallDirJob.folder, '%s.niceformat.vcf'%overlapIntervalFnameSignature)
				round1_VCF4OutputF = File(round1_VCF4OutputFname)
				round1_vcf_convert_job = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						GenomeAnalysisTKJar=GenomeAnalysisTKJar, inputF=vcf1FilterByvcftoolsJob.output, outputF=round1_VCF4OutputF, \
						refFastaFList=refFastaFList, parentJobLs=[vcf1FilterByvcftoolsJob], \
						extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=job_max_memory, interval=overlapInterval)
				if outputPedigreeJob is None:
					outputFileFormat=2	#trioCaller
					pedFile = File(os.path.join(trioCallerOutputDirJob.output, 'pedigree_outputFileFormat%s.txt'%(outputFileFormat)))
					sampleID2FamilyCountF = File(os.path.join(trioCallerOutputDirJob.output, 'sampleID2FamilyCount_outputFileFormat%s.txt'%(outputFileFormat)))
					outputPedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
							inputFile=round1_VCF4OutputF, outputFile=pedFile, sampleID2FamilyCountF=sampleID2FamilyCountF,\
							polymuttDatFile = None,\
							outputFileFormat=outputFileFormat, replicateIndividualTag=replicateIndividualTag,\
							treatEveryOneIndependent=treatEveryOneIndependent,\
							parentJobLs=[round1_vcf_convert_job, trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=True, \
							extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
				
				
				#2012.4.2 replicate individuals who appear in more than 1 families
				round1_IndividualsReplicatedVCF = File( os.path.join(round1CallDirJob.folder, '%s.replicate.vcf'%overlapIntervalFnameSignature))
				replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(workflow, \
								executable=workflow.ReplicateVCFGenotypeColumns, inputF=round1_VCF4OutputF, \
								sampleID2FamilyCountF=outputPedigreeJob.sampleID2FamilyCountF, outputF=round1_IndividualsReplicatedVCF, \
								replicateIndividualTag=replicateIndividualTag,\
								parentJobLs=[round1_vcf_convert_job, outputPedigreeJob], extraDependentInputLs=[], \
								transferOutput=tranferIntermediateFilesForDebug, \
								extraArguments=None, job_max_memory=500)
				
				#TrioCaller job
				refineGenotypeOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature))
				refineGenotypeJob = self.addTrioCallerJob(workflow, trioCallerWrapper=trioCallerWrapper, trioCallerPath=trioCallerPath, \
						inputVCF=round1_IndividualsReplicatedVCF,\
						pedFile=outputPedigreeJob.output, outputVCF=refineGenotypeOutputF, \
						parentJobLs=[trioCallerOutputDirJob, replicateVCFGenotypeColumnsJob, outputPedigreeJob], \
						extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=vcf_job_max_memory)
				
				#2012.4.2
				mergeReplicateOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.noReplicate.vcf'%overlapIntervalFnameSignature))
				noOfAlignments= len(alignmentDataLs)
				entireLength = overlapStop - overlapStart + 1	#could be very small for shorter reference contigs
				memoryRequest = min(42000, max(4000, int(20000*(noOfAlignments/323.0)*(entireLength/2600000.0))) )
					#extrapolates (20000Mb memory for a 323-sample + 2.6Mbase reference length/26K loci)
					#upper bound is 42g. lower bound is 4g.
				mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(workflow, \
									executable=workflow.MergeVCFReplicateHaplotypesJava,\
									GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
									inputF=refineGenotypeOutputF, outputF=mergeReplicateOutputF, \
									replicateIndividualTag=replicateIndividualTag, \
									refFastaFList=refFastaFList, parentJobLs=[refineGenotypeJob], \
									extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
									extraArguments=None, job_max_memory=memoryRequest)
				
				"""
				#2012.6.1 commented out, the overlap is key for ligateVcf.pl to ligate them
				#select the variants to get rid of overlap
				nonOverlapTrioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.nonoverlap.vcf'%intervalFnameSignature))
				trioCallerSelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						GenomeAnalysisTKJar=GenomeAnalysisTKJar, inputF=mergeVCFReplicateColumnsJob.output, \
						outputF=nonOverlapTrioCallerOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[mergeVCFReplicateColumnsJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=interval)
				"""
				
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
				
				#add this output to the union job
				# 2012.6.1 done it through addInputToStatMergeJob()
				self.addInputToStatMergeJob(statMergeJob=trioCallerWholeContigConcatJob, inputF=trioGzipOutputF, \
							parentJobLs=[bgzip_tabix_trioOutputF_job], \
							extraDependentInputLs=[trioGzipOutput_tbi_F])
				
				
				lisOfJobs = [preTrioCallerCallJob]
				for job in lisOfJobs:
					self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
				
				self.addAlignmentAsInputToJobLs(workflow, alignmentDataLs, jobLs=[preTrioCallerCallJob], jobInputOption="")
				no_of_jobs +=10
		
		sys.stderr.write(" %s jobs. \n"%(no_of_jobs))
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
				concatGATK=None, concatSamtools=None, ligateVcf=None, ligateVcfPerlPath=None,\
				refFastaFList=None, \
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
		refFastaF = refFastaFList[0]
		
		if needFastaDictJob:	# the .dict file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, \
												refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		
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
			trioCallerWholeContigConcatJob = self.addLigateVcfJob(executable=ligateVcf, ligateVcfPerlPath=ligateVcfPerlPath, \
										outputFile=concatTrioCallerOutputF, \
										parentJobLs=[trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
										extraArguments=None, job_max_memory=vcf_job_max_memory)
			
			#bgzip and tabix the trio caller output
			bgzip_concatTrioCallerOutputF = File("%s.gz"%concatTrioCallerOutputFname)
			bgzip_concatTrioCallerOutput_tbi_F = File("%s.gz.tbi"%concatTrioCallerOutputFname)
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
				round1_vcf_convert_job = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						GenomeAnalysisTKJar=GenomeAnalysisTKJar, inputF=splitVCFFile, outputF=round1_VCF4OutputF, \
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
									executable=workflow.MergeVCFReplicateHaplotypesJava,\
									GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
									inputF=refineGenotypeJob.output, outputF=mergeReplicateOutputF, \
									replicateIndividualTag=replicateIndividualTag, \
									refFastaFList=refFastaFList, parentJobLs=[refineGenotypeJob], \
									extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
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
					extraArguments=None, job_max_memory=2000, **keywords):
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
					key2ObjectForJob=key2ObjectForJob, **keywords)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		AlignmentToCallPipeline.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = self.db_vervet
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
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
									minNoOfLoci=10)	#ignore files with too few loci
			inputF = inputData.jobDataLs[0].vcfFile
			vcfFile = VCFFile(inputFname=inputF.abspath)
			alignmentLs = db_vervet.getAlignmentsFromVCFSampleIDList(vcfFile.getSampleIDList())
			del vcfFile
		
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
												individual_site_id_set=set(self.site_id_ls),\
												mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls),\
									excludeContaminant=self.excludeContaminant)
		cumulativeMedianDepth = db_vervet.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
															defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.data_dir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		
		if self.run_type==1:
			alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, data_dir=self.data_dir)
			chr2size = self.getTopNumberOfContigs(self.topNumberOfContigs, contigMinRankBySize=self.contigMinRankBySize)
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
						ligateVcf=self.ligateVcf, ligateVcfPerlPath=self.ligateVcfPerlPath,\
						refFastaFList=refFastaFList, \
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
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
						ligateVcf=self.ligateVcf, ligateVcfPerlPath=self.ligateVcfPerlPath,\
						refFastaFList=refFastaFList, \
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
						site_type=self.site_type, data_dir=self.data_dir,\
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, maxSNPMissingRate=self.maxSNPMissingRate,\
						alnStatForFilterF=None, cumulativeMedianDepth=cumulativeMedianDepth,\
						run_type=self.run_type, transferOutput=True)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)

if __name__ == '__main__':
	main_class = AlignmentToTrioCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
