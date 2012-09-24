#!/usr/bin/env python
"""
Examples:
	#2012.04.02 test-run on 40 VRC sequences (at least >1 trios)
	%s  -i 15-50 -u yh -a 524 -s 2 -z 10.8.0.10 -o AlignmentToTrioCallPipeline_VRC_ISQ_15_50_ReplicateParentsTop2Contigs.xml
		-j condorpool -l condorpool -N2  -w pedigreeVRC_ISQ_15_50_Merlin.txt -C 1 -O1
		-n pedigreeVRC_VRC_ISQ_15_50_sampleID2FamilyCount.tsv   -S 447
	
	# 2011.12.14 on all VRC (-S 447) monkeys), top 25 contigs
	%s -u yh -a 524 -s 2 -z 10.8.0.10 -o AlignmentToTrioCallPipeline_VRC_top25Contigs.xml 
		-j condorpool -l condorpool -N25 -w pedigreeVRCMerlin.txt -O 5 -S 447
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "-O 3" controls clustering of calling programs.
	# "-C 30" controls clustering for other programs., "-S 447" dictates monkeys from VRC
	%s -u yh -a 524 -s 2 -z localhost -o workflow/AlignmentToCall/AlignmentToTrioCallPipeline_VRC_top7559Contigs.xml -j hcondor -l hcondor 
		-N7559 -w pedigreeVRCMerlin.txt -O 3 -C 30 -S 447 -e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		
	# 2012.4.13 run VRC with 5 SK + 5 Nevis, sequenced filtered (-Q1), alignment by method 2 (-G2)
	%s -u yh -a 524 -s 2  -z localhost -o workflow/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_VRC_SK_Nevis_FilteredSeq_top1000Contigs.xml 
		-j hcondor -l hcondor -N1000 -w aux/pedigreeVRC_SK_Nevis_Merlin.2012.4.13.txt -O 3 -C 30 -S 447,417,420,427,431,432,435,437,439,440,442 
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-Q1 -G2  -n aux/VRC_SK_Nevis_sampleID2FamilyCount.2012.5.4.tsv
	
	# 2012.4.13 run 5 SK + 5 Nevis, sequenced filtered (-Q1), alignment by method 2 (-G2), TrioCaller fails because 10 is too small sample size.
	# 2012.6.13 supply the alignment depth stat file (-q), maxSNPMissingRate (-L 0.30), onlyKeepBiAllelicSNP (-c) 
	%s -u yh -a 524 -s 2 -z localhost -o workflow/AlignmentToCall/AlignmentToTrioCall_ReplicateIndividual_SK_Nevis_FilteredSeq_top1000Contigs.xml
		-j hcondor -l hcondor -N1000 -w aux/pedigreeSK_NV_Merlin.txt -O 2 -C 10 -S 417,420,427,431,432,435,437,439,440,442 
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-Q1 -G2 -n aux/10SKNV_sampleID2FamilyCount.tsv -q aux/alnStatForFilter.2012.6.13.tsv -c -L 0.30
		
	# 2012.8.15 run TrioCaller on method 14 samtools calls, contig ID from 96 to 100 (-V 96 -x 100)
	# sequenced filtered (-Q1), alignment by method 2 (-G2), onlyKeepBiAllelicSNP (-c) 
	# calling job clusters size=1, others =1 (-O 1 -C 1)
	# add -Y (not guess #loci from 1st number in filename) if the input VCF is not db affiliated
	# 3000 (-Z 3000) sites per unit, 500 overlapping between units (-U 500)
	# add "--treatEveryOneIndependent" if you want to treat everyone independent (no mendelian constraints from TrioCaller)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-u yh -a 524  -z localhost -o  workflow/AlignmentToCall/TrioCallerOnMethod14Contig96_100.xml
		-j hcondor -l hcondor  -w aux/pedigreeVRC_2012.8.15T1246.txt
		-O 1 -C 1 -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-V 96 -x 100 -Q1 -G2 -n aux/sampleID2FamilyCount_VRC_2012.8.15T1248.tsv  -c
		# -Y -Z 3000 -U 500 --treatEveryOneIndependent
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB, csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from Pegasus.DAX3 import *
from AlignmentToCallPipeline import AlignmentToCallPipeline
from AbstractVervetWorkflow import AbstractVervetWorkflow
from pymodule.VCFFile import VCFFile

class AlignmentToTrioCallPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update(AlignmentToCallPipeline.commonCallPipelineOptionDict)
	option_default_dict.pop(('ind_aln_id_ls', 0, ))
	option_default_dict.pop(('ind_seq_id_ls', 0, ))
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('pedigreeOutputFname', 1, ): ['', 'w', 1, 'the file which would contain the pedigree file in Merlin format'],\
						("trioCallerPath", 1, ): ["%s/script/vervet/bin/trioCaller/TrioCaller", 'P', 1, 'path to TrioCaller binary'],\
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						('sampleID2FamilyCountFname', 1, ): ['', 'n', 1, 'a tab-delimited file that records how many families in which each individual occurs'],\
						("onlyKeepBiAllelicSNP", 0, int): [0, 'c', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 0, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('depthFoldChange', 1, float): [2.0, 'y', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						("notToKnowNoOfLoci", 0, int): [0, 'Y', 0, 'toggle this to not guess the number of loci in one VCF file (for TrioCaller) from the 1st number in its filename'],\
						("treatEveryOneIndependent", 0, int): [0, '', 0, 'toggle this to treat everyone in the pedigree independent and also no replicates'],\
						
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
		self.trioCallerPath =  self.insertHomePath(self.trioCallerPath, self.home_path)
		
		#2012.8.15 change it to run after SAMtools
		self.run_type = 2
	
	def writeOneLineToPedigreeFile(self, writer, family_id, individual_alignment=None, father_id=0, mother_id=0,
								sex=None, replicateIndividualTag='copy', sampleID2FamilyCount=None):
		"""
		2012.3.29
		2011.12.5
			used by outputPedgreeOfAlignmentsInMerlinFormat()
		"""
		if type(individual_alignment)==str:
			sample_id = individual_alignment
			sexByGuess = 'x'
		else:
			#sample_id = individual_alignment.getReadGroup()
			sample_id = self.getSampleIDWithReplicateCount(individual_alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
			sexByGuess = individual_alignment.individual_sequence.individual.codeSexInNumber()
		if sex is None:
			sex = sexByGuess
		
		data_row = [family_id, sample_id, father_id, mother_id, sex]
		writer.writerow(data_row)
	
	def getSampleIDWithReplicateCount(self, alignment, replicateIndividualTag=None, sampleID2FamilyCount=None):
		"""
		2012.3.29
			use alignment.getReadGroup() + someReplicateIndicator to identify each VCF column
		"""
		hashedID = '%s%s%s'%(alignment.getReadGroup(), replicateIndividualTag, \
							sampleID2FamilyCount.get(alignment.getReadGroup()))
		return hashedID

	def adjustSampleID2FamilyCountFromAlignmentDBEntry(self, alignment, sampleID2FamilyCount=None):
		"""
		2012.3.29
		"""
		sample_id = alignment.getReadGroup()
		if sample_id not in sampleID2FamilyCount:
			sampleID2FamilyCount[sample_id] = 0
		sampleID2FamilyCount[sample_id] += 1
		
	def outputPedgreeOfAlignmentsInMerlinFormat(self, db_vervet, alignmentLs, pedigreeOutputFname=None, replicateIndividualTag='copy',\
											treatEveryOneIndependent=False):
		"""
		2012.8.28
			add argument treatEveryOneIndependent, which removes the family structure and also no replicates
		2012.3.29
			every individual gets replicated by the number of families it appears in.
		2011-11-27
			
		2011-11-3
			find all trios, a list of "father,mother,child", all of which are identified by IndividualSequence.id.
			
			for each individual_alignment, -> ind_seq -> individual
				get its parents from ind2ind
				check if both parents have also been sequenced and aligned.
		"""
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		DG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		
		#2012.8.29
		if treatEveryOneIndependent:
			removeFamilyFromGraph = True
		else:
			removeFamilyFromGraph = False
		#find trios first
		trioLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=3, removeFamilyFromGraph=removeFamilyFromGraph)
		#find duos
		duoLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=2, removeFamilyFromGraph=removeFamilyFromGraph)
		#find singletons (familySize=1 => noOfIncomingEdges=0, noOfIncomingEdges=0 => will not be parents of others)
		singletonLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=1, removeFamilyFromGraph=removeFamilyFromGraph, noOfOutgoingEdges=0)
		#[[node] for node in DG.nodes()]. this will be good only when removeFamilyFromGraph=True for all sentences above.
		
		familyLs = trioLs + duoLs + singletonLs
		sys.stderr.write("Outputting %s families in units of trios (Merlin format)..."%(len(familyLs)))
		#time to output
		sampleID2FamilyCount = {}	#2012.3.29
		writer = csv.writer(open(pedigreeOutputFname, 'w'), delimiter=' ')
		for family in familyLs:
			familySize = len(family)
			for memberID in family:	#2012.3.29 make sure to record the occurrence of each member before writing them out
				# as the occurrence is used in output 
				alignment = individual_id2alignmentLs.get(memberID)[0]
				self.adjustSampleID2FamilyCountFromAlignmentDBEntry(alignment, sampleID2FamilyCount)
				
			if familySize==1:
				individual_id = family[0]
				alignment = individual_id2alignmentLs.get(individual_id)[0]
				
				sampleID = self.getSampleIDWithReplicateCount(alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				family_id = "F%s_%s"%(individual_id, sampleID2FamilyCount[alignment.getReadGroup()])
				self.writeOneLineToPedigreeFile(writer, family_id, alignment, father_id=0, mother_id=0,\
										replicateIndividualTag=replicateIndividualTag, \
										sampleID2FamilyCount=sampleID2FamilyCount)
			elif familySize==2:
				parent1ID, offspring_id = family[:2]
				#output the single parent first
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				parent1sampleID = self.getSampleIDWithReplicateCount(parent1Alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				family_id = "F%s_%s"%(parent1ID, sampleID2FamilyCount[parent1Alignment.getReadGroup()])	#2012.3.29 add count to family ID
				self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0,\
										replicateIndividualTag=replicateIndividualTag, \
										sampleID2FamilyCount=sampleID2FamilyCount)
				if treatEveryOneIndependent:	#2012.8.29
					#output the offspring
					childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
					father_id = 0
					mother_id = 0
					self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				else:
					#output a fake 2nd parent, with opposite sex
					sndParentID = 'dummy'
					sndParentSex = 1-(parent1Alignment.individual_sequence.individual.codeSexInNumber()-1)+1
					self.writeOneLineToPedigreeFile(writer, family_id, sndParentID, father_id=0, mother_id=0, sex=sndParentSex,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
					#output the offspring
					childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
					if sndParentSex==1:
						father_id = sndParentID
						mother_id = parent1sampleID
					else:
						father_id = parent1sampleID
						mother_id = sndParentID
					self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
			elif familySize==3:
				parent1ID, parent2ID, offspring_id = family[:3]
				#output one parent
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				parent1SampleID = self.getSampleIDWithReplicateCount(parent1Alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				family_id = "F%s_%s"%(parent1ID, sampleID2FamilyCount[parent1Alignment.getReadGroup()])
				self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0,\
										replicateIndividualTag=replicateIndividualTag, \
										sampleID2FamilyCount=sampleID2FamilyCount)
				#output 2nd parent
				parent2Alignment = individual_id2alignmentLs.get(parent2ID)[0]
				parent2SampleID = self.getSampleIDWithReplicateCount(parent2Alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				self.writeOneLineToPedigreeFile(writer, family_id, parent2Alignment, father_id=0, mother_id=0,\
										replicateIndividualTag=replicateIndividualTag, \
										sampleID2FamilyCount=sampleID2FamilyCount)
				#output offspring
				if treatEveryOneIndependent:	#2012.8.29
					#output the offspring
					childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
					father_id = 0
					mother_id = 0
					self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
				else:
					childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
					parent2Sex = parent2Alignment.individual_sequence.individual.codeSexInNumber()
					if parent2Sex==1:
						father_id = parent2SampleID
						mother_id = parent1SampleID
					else:
						father_id = parent1SampleID
						mother_id = parent2SampleID
					self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
		del writer
		sys.stderr.write("Done.\n")
		return sampleID2FamilyCount
		
	def outputSampleID2FamilyCount(self, sampleID2FamilyCount, outputFname=None):
		"""
		2012.3.29
		"""
		sys.stderr.write("Outputting sampleID2FamilyCount (%s entries) to %s ..."%(len(sampleID2FamilyCount), outputFname))
		header= ['individualID', 'familyCount']
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(header)
		count = 0
		for individualID, familyCount in sampleID2FamilyCount.iteritems():
			data_row = [individualID, familyCount]
			writer.writerow(data_row)
			count +=1 
		del writer
		sys.stderr.write("%s entries outputted.\n"%(len(sampleID2FamilyCount)))
		
	
	def addTrioCallerJob(self, workflow, trioCallerWrapper=None, trioCallerPath=None, inputVCF=None, pedFile=None, outputVCF=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.1.20
			increase the rounds from 30 to 40
			add --burnin 20
		2011-12-4
		"""
		job = Job(namespace=workflow.namespace, name=trioCallerWrapper.name, version=workflow.version)
		job.addArguments(trioCallerPath, "--vcf", inputVCF, "--pedfile", pedFile, \
						"--states 50 --randomPhase --rounds 40 --burnin 20",\
						"--prefix %s"%(os.path.splitext(outputVCF.name)[0]))
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputVCF, transfer=True, register=True, link=Link.INPUT)
		job.uses(pedFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputVCF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputVCF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addGenotypeCallJobs(self, workflow=None, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				genotyperJava=None, SelectVariantsJava=None, genomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				mergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, CallVariantBySamtools=None,\
				trioCallerPath=None, trioCallerWrapper=None, pedFile=None, \
				sampleID2FamilyCountF=None,\
				replicateIndividualTag="copy", 
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None, ligateVcf=None, ligateVcfPerlPath=None,\
				refFastaFList=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, intervalOverlapSize=100000, site_type=1, dataDir=None, no_of_gatk_threads = 1, \
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
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
												createSequenceDictionaryJar=createSequenceDictionaryJar, refFastaF=refFastaF)
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
		trioCallerOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=trioCallerOutputDir)
		round1CallDir = "%spreTrioCaller"%(outputDirPrefix)
		round1CallDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=round1CallDir)
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow, alignmentDataLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		
		# add merge jobs for every reference
		returnData = PassingData()
		returnData.jobDataLs = []
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
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
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
				
				preTrioCallerCallJob= self.addGATKCallJob(workflow, genotyperJava=genotyperJava, genomeAnalysisTKJar=genomeAnalysisTKJar, \
						round1CallOutputF=round1CallOutputF, gatkIDXOutput=gatkIDXOutput, refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob], \
						extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
						job_max_memory=job_max_memory, no_of_gatk_threads=no_of_gatk_threads, site_type=site_type, \
						interval=overlapInterval)
				"""
				# following few steps is prepare 1st-round calls to be concatenated (not for 2nd-round caller)
				#select the variants to get rid of overlap, so that union of whole contig doesn't have overlap
				round1_NonOverlapOutputF = File(os.path.join(round1CallDirJob.folder, '%s.nonoverlap.vcf'%intervalFnameSignature))
				round1SelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=vcf1FilterByvcftoolsJob.output, outputF=round1_NonOverlapOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[vcf1FilterByvcftoolsJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=interval)
				
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
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=vcf1FilterByvcftoolsJob.output, outputF=round1_VCF4OutputF, \
						refFastaFList=refFastaFList, parentJobLs=[vcf1FilterByvcftoolsJob], \
						extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=job_max_memory, interval=overlapInterval)
				
				#2012.4.2 replicate individuals who appear in more than 1 families
				round1_IndividualsReplicatedVCF = File( os.path.join(round1CallDirJob.folder, '%s.replicate.vcf'%overlapIntervalFnameSignature))
				replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(workflow, executable=workflow.ReplicateVCFGenotypeColumns, inputF=round1_VCF4OutputF, \
										sampleID2FamilyCountF=sampleID2FamilyCountF, outputF=round1_IndividualsReplicatedVCF, \
										replicateIndividualTag=replicateIndividualTag,\
										parentJobLs=[round1_vcf_convert_job], extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
										extraArguments=None, job_max_memory=500)
				
				#TrioCaller job
				trioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature))
				trioCallerJob = self.addTrioCallerJob(workflow, trioCallerWrapper=trioCallerWrapper, trioCallerPath=trioCallerPath, \
						inputVCF=round1_IndividualsReplicatedVCF,\
						pedFile=pedFile, outputVCF=trioCallerOutputF, \
						parentJobLs=[trioCallerOutputDirJob, replicateVCFGenotypeColumnsJob], \
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
									genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
									inputF=trioCallerOutputF, outputF=mergeReplicateOutputF, \
									replicateIndividualTag=replicateIndividualTag, \
									refFastaFList=refFastaFList, parentJobLs=[trioCallerJob], \
									extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
									extraArguments=None, job_max_memory=memoryRequest)
				
				"""
				#2012.6.1 commented out, the overlap is key for ligateVcf.pl to ligate them
				#select the variants to get rid of overlap
				nonOverlapTrioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.nonoverlap.vcf'%intervalFnameSignature))
				trioCallerSelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=mergeVCFReplicateColumnsJob.output, \
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
				genotyperJava=None, SelectVariantsJava=None, genomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				mergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, CallVariantBySamtools=None,\
				trioCallerPath=None, trioCallerWrapper=None, pedFile=None, \
				sampleID2FamilyCountF=None,\
				replicateIndividualTag="copy", 
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None, ligateVcf=None, ligateVcfPerlPath=None,\
				refFastaFList=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, intervalOverlapSize=100000, site_type=1, dataDir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", \
				maxSNPMissingRate=None, alnStatForFilterF=None, onlyKeepBiAllelicSNP=True, \
				cumulativeMedianDepth=5000, job_max_memory = 2000, vcf_job_max_memory = 1000,\
				transferOutput=True, **keywords):
		"""
		2012.8.15
		"""
		sys.stderr.write("Adding trioCaller jobs for  %s vcf files ..."%(len(inputData.jobDataLs)))
		if workflow is None :
			workflow = self
		refFastaF = refFastaFList[0]
		no_of_jobs = 0
		
		if needFastaDictJob:	# the .dict file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
												refFastaF=refFastaF)
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
		trioCallerOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=trioCallerOutputDir)
		round1CallDir = "%spreTrioCaller"%(outputDirPrefix)
		round1CallDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=round1CallDir)
		
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
			
			outputFnamePrefix = os.path.join(round1CallDirJob.folder, '%s_splitVCF'%commonPrefix)
			splitVCFJob = self.addSplitVCFFileJob(executable=self.SplitVCFFile, inputFile=inputF, outputFnamePrefix=outputFnamePrefix, \
					noOfOverlappingSites=intervalOverlapSize, noOfSitesPerUnit=intervalSize, noOfTotalSites=inputF.noOfLoci, \
					parentJobLs=jobData.jobLs+[round1CallDirJob], \
					extraDependentInputLs=[jobData.tbi_F], \
					extraArguments=None, transferOutput=False, job_max_memory=job_max_memory)
			no_of_jobs +=1
			
			concatTrioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.vcf'%chr_id)
			concatTrioCallerOutputF = File(concatTrioCallerOutputFname)
			trioCallerWholeContigConcatJob = self.addLigateVcfJob(executable=ligateVcf, ligateVcfPerlPath=ligateVcfPerlPath, \
										outputFile=concatTrioCallerOutputF, \
										parentJobLs=[trioCallerOutputDirJob], extraDependentInputLs=[], transferOutput=False, \
										extraArguments=None, job_max_memory=vcf_job_max_memory)
			no_of_jobs +=1
			
			#bgzip and tabix the trio caller output
			bgzip_concatTrioCallerOutputF = File("%s.gz"%concatTrioCallerOutputFname)
			bgzip_concatTrioCallerOutput_tbi_F = File("%s.gz.tbi"%concatTrioCallerOutputFname)
			bgzip_tabix_concatTrioCallerOutput_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
					parentJob=trioCallerWholeContigConcatJob, inputF=concatTrioCallerOutputF, outputF=bgzip_concatTrioCallerOutputF, \
					transferOutput=transferOutput)
			
			returnData.jobDataLs.append(PassingData(vcfFile=bgzip_concatTrioCallerOutputF, jobLs=[bgzip_tabix_concatTrioCallerOutput_job]))
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			no_of_jobs +=1
			
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
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=splitVCFFile, outputF=round1_VCF4OutputF, \
						refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob, splitVCFJob], \
						extraDependentInputLs=None, transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=job_max_memory, interval=overlapInterval)
				no_of_jobs +=1
				
				#2012.4.2 replicate individuals who appear in more than 1 families
				round1_IndividualsReplicatedVCF = File( os.path.join(round1CallDirJob.folder, \
													'%s.replicate.vcf'%overlapIntervalFnameSignature))
				replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(workflow, \
							executable=workflow.ReplicateVCFGenotypeColumns, inputF=round1_VCF4OutputF, \
							sampleID2FamilyCountF=sampleID2FamilyCountF, outputF=round1_IndividualsReplicatedVCF, \
							replicateIndividualTag=replicateIndividualTag,\
							parentJobLs=[round1_vcf_convert_job], extraDependentInputLs=[], \
							transferOutput=tranferIntermediateFilesForDebug, \
							extraArguments=None, job_max_memory=vcf_job_max_memory)
				no_of_jobs +=1
				
				#TrioCaller job
				trioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.orig.vcf'%overlapIntervalFnameSignature))
				trioCallerJob = self.addTrioCallerJob(workflow, trioCallerWrapper=trioCallerWrapper, trioCallerPath=trioCallerPath, \
						inputVCF=round1_IndividualsReplicatedVCF,\
						pedFile=pedFile, outputVCF=trioCallerOutputF, \
						parentJobLs=[trioCallerOutputDirJob, replicateVCFGenotypeColumnsJob], \
						extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
						extraArguments=None, job_max_memory=vcf_job_max_memory)	#1.2G memory for 12K loci
				no_of_jobs +=1
				
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
									genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
									inputF=trioCallerOutputF, outputF=mergeReplicateOutputF, \
									replicateIndividualTag=replicateIndividualTag, \
									refFastaFList=refFastaFList, parentJobLs=[trioCallerJob], \
									extraDependentInputLs=[], transferOutput=tranferIntermediateFilesForDebug, \
									extraArguments=None, job_max_memory=memoryRequest,\
									analysis_type='MergeVCFReplicateGenotypeColumns')
				no_of_jobs +=1
				
				#ligate job
				
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
				
				no_of_jobs +=2
				
				#add this output to the union job
				# 2012.6.1 done it through addInputToStatMergeJob()
				self.addInputToStatMergeJob(statMergeJob=trioCallerWholeContigConcatJob, inputF=trioGzipOutputF, \
							parentJobLs=[bgzip_tabix_trioOutputF_job], \
							extraDependentInputLs=[trioGzipOutput_tbi_F])
				
		sys.stderr.write(" %s jobs. \n"%(no_of_jobs))
		return returnData
	
	def addReplicateVCFGenotypeColumnsJob(self, workflow, executable=None, inputF=None, sampleID2FamilyCountF=None, outputF=None, \
					replicateIndividualTag=None,\
					parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.4.2
			a job that replicates the genotype columns of individuals who appear in >1 families.
			This is to be run before TrioCaller is applied.
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-n', sampleID2FamilyCountF, '-T', replicateIndividualTag, '-i', inputF, '-o', outputF)
		job.uses(sampleID2FamilyCountF, transfer=True, register=True, link=Link.INPUT)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
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
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
		if self.run_type==1:
			alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, dataDir=self.localDataDir)
		elif self.run_type==2:
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
												individual_site_id_set=set(self.site_id_ls))
		cumulativeMedianDepth = db_vervet.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
															defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		sampleID2FamilyCount = self.outputPedgreeOfAlignmentsInMerlinFormat(db_vervet, alignmentLs, self.pedigreeOutputFname,\
																treatEveryOneIndependent=self.treatEveryOneIndependent)
		
		self.outputSampleID2FamilyCount(sampleID2FamilyCount, outputFname=self.sampleID2FamilyCountFname)
		
		
		pedFile = yh_pegasus.registerFile(workflow, self.pedigreeOutputFname)
		sampleID2FamilyCountF = self.registerOneInputFile(workflow, self.sampleID2FamilyCountFname, folderName="")
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		
		if self.run_type==1:
			alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
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
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
						createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
						mergeSamFilesJar=workflow.mergeSamFilesJar, \
						BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
						mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools, \
						trioCallerPath=self.trioCallerPath, trioCallerWrapper=workflow.trioCallerWrapper, pedFile=pedFile, \
						sampleID2FamilyCountF=sampleID2FamilyCountF,\
						replicateIndividualTag=self.replicateIndividualTag, 
						bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
						concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
						ligateVcf=self.ligateVcf, ligateVcfPerlPath=self.ligateVcfPerlPath,\
						refFastaFList=refFastaFList, \
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
						site_type=self.site_type, dataDir=self.dataDir,\
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, maxSNPMissingRate=self.maxSNPMissingRate,\
						alnStatForFilterF=alnStatForFilterF, cumulativeMedianDepth=cumulativeMedianDepth,\
						transferOutput=True)
		elif self.run_type==2:
			self.addTrioCallerJobsONVCFFiles(workflow=workflow, alignmentLs=alignmentLs, inputData=inputData, \
						samtools=workflow.samtools, \
						genotyperJava=workflow.genotyperJava,  SelectVariantsJava=workflow.SelectVariantsJava, \
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
						createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
						mergeSamFilesJar=workflow.mergeSamFilesJar, \
						BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
						mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools, \
						trioCallerPath=self.trioCallerPath, trioCallerWrapper=workflow.trioCallerWrapper, pedFile=pedFile, \
						sampleID2FamilyCountF=sampleID2FamilyCountF,\
						replicateIndividualTag=self.replicateIndividualTag, 
						bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
						concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
						ligateVcf=self.ligateVcf, ligateVcfPerlPath=self.ligateVcfPerlPath,\
						refFastaFList=refFastaFList, \
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
						site_type=self.site_type, dataDir=self.dataDir,\
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, maxSNPMissingRate=self.maxSNPMissingRate,\
						alnStatForFilterF=None, cumulativeMedianDepth=cumulativeMedianDepth,\
						transferOutput=True)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)

if __name__ == '__main__':
	main_class = AlignmentToTrioCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
