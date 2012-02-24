#!/usr/bin/env python
"""
Examples:
	#2011.12.14 test on one trio
	%s -I 620,632,648   -u yh -a 524 -s 2 -z 10.8.0.10 
		-o AlignmentToTrioCallPipeline_620_632_648_top2Contigs.xml -j condorpool -l condorpool -N2 
		-F pedigreeTrio620_632_648.txt -C 1 -O 1
	
	# 2011.12.14 on all VRC (-S 447) monkeys), top 25 contigs
	%s -u yh -a 524 -s 2 -z 10.8.0.10 -o AlignmentToTrioCallPipeline_VRC_top25Contigs.xml 
		-j condorpool -l condorpool -N25 -F pedigreeVRCMerlin.txt -O 5 -S 447
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "-O 3" controls clustering of calling programs.
	# "-C 30" controls clustering for other programs., "-S 447" dictates monkeys from VRC
	%s -u yh -a 524 -s 2 -z localhost -o AlignmentToTrioCallPipeline_VRC_top7559Contigs.xml -j hcondor -l hcondor 
		-N7559 -F pedigreeVRCMerlin.txt -O 3 -C 30 -S 447 -e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
		If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB, csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from AlignmentToCallPipeline import AlignmentToCallPipeline
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class AlignmentToTrioCallPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						('pedigreeOutputFname', 1, ): ['', 'F', 1, 'the file which would contain the pedigree file in Merlin format'],\
						("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\
						("topNumberOfContigs", 1, int): [156, 'N', 1, 'number of contigs'],\
						("needFastaIndexJob", 0, int): [0, '', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, '', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						("site_type", 1, int): [2, 's', 1, '1: all genome sites, 2: variants only'],\
						("noOfCallingJobsPerNode", 1, int): [1, 'O', 1, 'this parameter controls how many genotype calling jobs should be clustered together to run on one node. \
									Increase it to above 1 only when your average genotyping job is short and the number of input bam files are short.'],\
						("site_id", 1, int): [447, 'S', 1, 'individuals must come from this site.'],\
						("trioCallerPath", 1, ): ["%s/script/vervet/bin/trioCaller/TrioCaller", 'P', 1, 'path to TrioCaller binary'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AlignmentToCallPipeline.__init__(self, **keywords)
		self.trioCallerPath = self.trioCallerPath%self.home_path
	
	def writeOneLineToPedigreeFile(self, writer, family_id, individual_alignment, father_id=0, mother_id=0,
								sex=None):
		"""
		2011.12.5
			used by outputPedgreeOfAlignmentsInMerlinFormat()
		"""
		if type(individual_alignment)==str:
			individual_id = individual_alignment
			sexByGuess = 'x'
		else:
			individual_id = individual_alignment.getReadGroup()
			sexByGuess = individual_alignment.ind_sequence.individual.codeSexInNumber()
		if sex is None:
			sex = sexByGuess
		data_row = [family_id, individual_id, father_id, mother_id, sex]
		writer.writerow(data_row)
	
	def outputPedgreeOfAlignmentsInMerlinFormat(self, db_vervet, alignmentLs, pedigreeOutputFname):
		"""
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
		
		#find trios first
		trioLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=3, removeFamilyFromGraph=True)
		#find duos
		duoLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=2, removeFamilyFromGraph=True)
		#remaining are singletons
		singletonLs = [[node] for node in DG.nodes()]
		
		familyLs = trioLs + duoLs + singletonLs
		sys.stderr.write("Outputting %s families in units of trios (Merlin format)..."%(len(familyLs)))
		#time to output
		writer = csv.writer(open(pedigreeOutputFname, 'w'), delimiter=' ')
		for family in familyLs:
			familySize = len(family)
			if familySize==1:
				individual_id = family[0]
				family_id = "F%s"%(individual_id)
				alignment = individual_id2alignmentLs.get(individual_id)[0]
				self.writeOneLineToPedigreeFile(writer, family_id, alignment, father_id=0, mother_id=0)
			elif familySize==2:
				parent1ID, offspring_id = family[:2]
				family_id = "F%s"%(parent1ID)
				#output the single parent first
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0)
				#output a fake 2nd parent, with opposite sex
				sndParentID = 'dummy'
				sndParentSex = 1-(parent1Alignment.ind_sequence.individual.codeSexInNumber()-1)+1
				self.writeOneLineToPedigreeFile(writer, family_id, sndParentID, father_id=0, mother_id=0, sex=sndParentSex)
				#output the offspring
				childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
				if sndParentSex==1:
					father_id = sndParentID
					mother_id = parent1Alignment.getReadGroup()
				else:
					father_id = parent1Alignment.getReadGroup()
					mother_id = sndParentID
				self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id)
			elif familySize==3:
				parent1ID, parent2ID, offspring_id = family[:3]
				family_id = "F%s"%(parent1ID)
				#output one parent
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0)
				#output 2nd parent
				parent2Alignment = individual_id2alignmentLs.get(parent2ID)[0]
				self.writeOneLineToPedigreeFile(writer, family_id, parent2Alignment, father_id=0, mother_id=0)
				#output offspring
				childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
				parent2Sex = parent2Alignment.ind_sequence.individual.codeSexInNumber()
				if parent2Sex==1:
					father_id = parent2Alignment.getReadGroup()
					mother_id = parent1Alignment.getReadGroup()
				else:
					father_id = parent1Alignment.getReadGroup()
					mother_id = parent2Alignment.getReadGroup()
				self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id)
		del writer
		sys.stderr.write("Done.\n")
	
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
		job.addArguments(trioCallerPath, "--shotgun", inputVCF, "--pedfile", pedFile, \
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
	
	def addGenotypeCallJobs(self, workflow, alignmentDataLs, refName2size, samtools=None, \
				genotyperJava=None, SelectVariantsJava=None, genomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				mergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, CallVariantBySamtools=None,\
				trioCallerPath=None, trioCallerWrapper=None, pedFile=None, \
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None,\
				refFastaFList=None, \
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, intervalOverlapSize=100000, site_type=1, dataDir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", **keywords):
		"""
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
		sys.stderr.write("Adding genotype call jobs for %s chromosomes/contigs ..."%(len(refName2size)))
		job_max_memory = 2000	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		vcf_job_max_memory = 1000
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
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow, alignmentDataLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		
		# add merge jobs for every reference
		refName2mergedBamCallJob = {}
		returnData = PassingData()
		returnData.jobDataLs = []
		for refName, refSize in refName2size.iteritems():
			no_of_intervals = max(1, int(math.ceil(refSize/float(intervalSize)))-1)
			#reduce the number of chunks 1 below needed. last trunk to reach the end of contig
			#however set it to 1 for contigs smaller than intervalSize 	
			concatTrioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.vcf.gz'%refName)
			concatTrioCallerOutputF = File(concatTrioCallerOutputFname)
			wholeContigConcatJob = self.addVCFConcatJob(workflow, concatExecutable=vcf_concat, \
							parentDirJob=trioCallerOutputDirJob, \
							outputF=concatTrioCallerOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			
			returnData.jobDataLs.append(PassingData(vcfFile=concatTrioCallerOutputF, jobLs=[wholeContigConcatJob]))
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			
			#2011-9-22 union of all 1st-rounding calling intervals for one contig
			round1_VCFConcatOutputFname = os.path.join(round1CallDirJob.folder, '%s.vcf.gz'%refName)
			round1_VCFConcatOutputF = File(round1_VCFConcatOutputFname)
			round1_VCFConcatJob = self.addVCFConcatJob(workflow, concatExecutable=concatGATK, parentDirJob=round1CallDirJob, \
							outputF=round1_VCFConcatOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			
			no_of_jobs += 2
			
			
			for i in range(no_of_intervals):
				originalStartPos = i*intervalSize + 1
				#to render adjacent intervals overlapping because trioCaller uses LD
				startPos = max(1, originalStartPos-intervalOverlapSize)
				if i<no_of_intervals-1:
					originalStopPos = min((i+1)*intervalSize, refSize)
				else:	#last chunk, include bp till the end
					originalStopPos = refSize
				#to render adjacent intervals overlapping because trioCaller uses LD
				stopPos = min(refSize, originalStopPos+intervalOverlapSize)
				
				nonOverlapInterval = "%s:%s-%s"%(refName, originalStartPos, originalStopPos)
				interval = "%s:%s-%s"%(refName, startPos, stopPos)
				nonOverlapVCFBaseFname = '%s_%s_%s'%(refName, originalStartPos, originalStopPos)
				vcfBaseFname = '%s_%s_%s'%(refName, startPos, stopPos)
				
				#samtools part
				round1CallOutputFname = os.path.join(round1CallDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
				round1CallOutputF = File(round1CallOutputFname)
				indelVCFOutputFname = "%s.indel.vcf"%(round1CallOutputFname)
				indelVCFOutputF = File(indelVCFOutputFname)
				preTrioCallerCallJob = self.addSAMtoolsCallJob(workflow, CallVariantBySamtools=CallVariantBySamtools, \
					samtoolsOutputF=round1CallOutputF, indelVCFOutputF=indelVCFOutputF, \
					refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=vcf_job_max_memory, site_type=site_type, interval=interval)
				"""
				#GATK job
				#GATK produces "./." for missing genotype and TrioCaller has trouble passing that.
				round1CallOutputFname = os.path.join(round1CallDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
				round1CallOutputF = File(round1CallOutputFname)
				gatkIDXOutputFname = os.path.join(round1CallDirJob.folder, '%s.vcf.idx'%(vcfBaseFname))
				gatkIDXOutput = File(gatkIDXOutputFname)
				
				preTrioCallerCallJob= self.addGATKCallJob(workflow, genotyperJava=genotyperJava, genomeAnalysisTKJar=genomeAnalysisTKJar, \
						round1CallOutputF=round1CallOutputF, gatkIDXOutput=gatkIDXOutput, refFastaFList=refFastaFList, parentJobLs=[round1CallDirJob], \
						extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
						job_max_memory=job_max_memory, no_of_gatk_threads=no_of_gatk_threads, site_type=site_type, \
						interval=interval)
				"""
				#select the variants to get rid of overlap, so that union of whole contig doesn't have overlap
				nonOverlapGATKOutputF = File(os.path.join(round1CallDirJob.folder, '%s.nonoverlap.vcf'%nonOverlapVCFBaseFname))
				round1SelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=round1CallOutputF, outputF=nonOverlapGATKOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[preTrioCallerCallJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=nonOverlapInterval)
				
				#convert to vcf4 so that other vcftools software could be used.
				round1_VCF4NonOverlapOutputFname = os.path.join(round1CallDirJob.folder, '%s.v4.vcf'%nonOverlapVCFBaseFname)
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
				round1_VCFConcatJob.addArguments(round1VCFGzipOutputF)
				round1_VCFConcatJob.uses(round1VCFGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				round1_VCFConcatJob.uses(round1VCFGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=round1_bgzip_tabix_VCFOutputJOb, child=round1_VCFConcatJob)
				
				"""
				#convert to vcf4 so that TrioCaller could read it. (samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
				round1_VCF4OutputFname = os.path.join(round1CallDirJob.folder, '%s.v4.vcf'%vcfBaseFname)
				round1_VCF4OutputF = File(round1_VCF4OutputFname)
				round1_vcf_convert_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=preTrioCallerCallJob, inputF=preTrioCallerCallJob.output, outputF=round1_VCF4OutputF, \
							namespace=namespace, version=version, transferOutput=False)
				"""
				#selectVariants would generate AC, AF so that TrioCaller could read it. (samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
				round1_VCF4OutputFname = os.path.join(round1CallDirJob.folder, '%s.niceformat.vcf'%vcfBaseFname)
				round1_VCF4OutputF = File(round1_VCF4OutputFname)
				round1_vcf_convert_job = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=preTrioCallerCallJob.output, outputF=round1_VCF4OutputF, \
						refFastaFList=refFastaFList, parentJobLs=[preTrioCallerCallJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=interval)
				
				
				#TrioCaller job
				trioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.orig.vcf'%vcfBaseFname))
				trioCallerJob = self.addTrioCallerJob(workflow, trioCallerWrapper=trioCallerWrapper, trioCallerPath=trioCallerPath, \
						inputVCF=round1_VCF4OutputF,\
						pedFile=pedFile, outputVCF=trioCallerOutputF, \
						parentJobLs=[trioCallerOutputDirJob, round1_vcf_convert_job], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=vcf_job_max_memory)
				
				#select the variants to get rid of overlap
				nonOverlapTrioCallerOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.nonoverlap.vcf'%nonOverlapVCFBaseFname))
				trioCallerSelectVariantJob = self.addSelectVariantsJob(workflow, SelectVariantsJava=SelectVariantsJava, \
						genomeAnalysisTKJar=genomeAnalysisTKJar, inputF=trioCallerOutputF, outputF=nonOverlapTrioCallerOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[trioCallerJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=job_max_memory, interval=nonOverlapInterval)
				
				#convert to vcf4 so that other vcftools software could be used.
				vcf4_trioCallerOutputFname = os.path.join(trioCallerOutputDirJob.folder, '%s.v4.vcf'%nonOverlapVCFBaseFname)
				vcf4_trioCallerOutputF = File(vcf4_trioCallerOutputFname)
				vcf_convert_TrioCallerOutputJob = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=trioCallerSelectVariantJob, inputF=trioCallerSelectVariantJob.output, \
							outputF=vcf4_trioCallerOutputF, \
							namespace=namespace, version=version, transferOutput=False)
				
				
				#bgzip and tabix the trio caller output
				trioGzipOutputF = File("%s.gz"%vcf4_trioCallerOutputFname)
				trioGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_trioCallerOutputFname)
				bgzip_tabix_trioOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=vcf_convert_TrioCallerOutputJob, inputF=vcf4_trioCallerOutputF, outputF=trioGzipOutputF, \
						transferOutput=False)
				
				
				#add this output to the union job
				wholeContigConcatJob.addArguments(trioGzipOutputF)
				wholeContigConcatJob.uses(trioGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				wholeContigConcatJob.uses(trioGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=bgzip_tabix_trioOutputF_job, child=wholeContigConcatJob)
				
				
				lisOfJobs = [preTrioCallerCallJob]
				for job in lisOfJobs:
					self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
				
				self.addAlignmentAsInputToJobLs(workflow, alignmentDataLs, jobLs=[preTrioCallerCallJob], jobInputOption="")
				no_of_jobs +=10
		
		sys.stderr.write(" %s jobs. \n"%(no_of_jobs))
		return returnData
	
	
	def registerCustomExecutables(self, workflow):
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
		
		refName2size = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refName2size = set(['Contig149'])	#temporary when testing Contig149
		#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = refName2size.keys()
		
		alignmentLs = self.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2, dataDir=self.localDataDir)
		alignmentLs = self.filterAlignments(alignmentLs, individual_site_id=self.site_id)
		self.outputPedgreeOfAlignmentsInMerlinFormat(db_vervet, alignmentLs, self.pedigreeOutputFname)
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		pedFile = yh_pegasus.registerFile(workflow, self.pedigreeOutputFname)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
		
		self.addGenotypeCallJobs(workflow, alignmentDataLs, refName2size, samtools=workflow.samtools, \
					genotyperJava=workflow.genotyperJava,  SelectVariantsJava=workflow.SelectVariantsJava, \
					genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
					createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
					mergeSamFilesJar=workflow.mergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools, \
					trioCallerPath=self.trioCallerPath, trioCallerWrapper=workflow.trioCallerWrapper, pedFile=pedFile, \
					bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
					concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
					refFastaFList=refFastaFList, \
					namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
					needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
					intervalSize=2000000, intervalOverlapSize=100000, site_type=self.site_type, dataDir=self.dataDir)
		
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = AlignmentToTrioCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
