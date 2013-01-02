#!/usr/bin/env python
"""
Examples:
	# scaffold (120) is used as reference in alignment. get genome sequence id from 1 to 8.
	%s -o workflow.xml -a 120 -i 1-8 -y2 -s2
	
	# 1Mb-BAC (9) is used as reference.
	%s -o workflow.xml -a 9 -i 1-4 -y2 -s2
	
	# 8 genomes versus top 156 contigs
	%s -o workflow_8GenomeVsTop156Contigs.xml -u yh  -a 128 -i 1-8 -N 156 -y2 -s2
	
	# 2011-7-21 use GATK + coverage filter, top 2 contigs
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8 -N 2 -y1  -s2
	
	# 2011-7-21 use GATK + coverage filter on hoffman2 and site_handler, top 5 contigs
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8
		-N 5 -y1 -l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db  -s2
		-J /u/home/eeskin/polyacti/bin/jdk/bin/java
	
	#2011-8-31 work 10 VWP and one VRC ref monkeys, variants only
	%s -a 9 -I 495,498-507 -u yh  
		-l condorpool -y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC.xml -s2 -q /tmp/all_isq_coverage.tsv
	
	%s -a 120 -I 34,38 -u yh -l hoffman2
	-y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC_hoffman2.xml  -s2 -e /u/home/eeskin/polyacti
	-t /u/home/eeskin/polyacti/NetworkData/vervet/db -N 4 -q /tmp/all_isq_coverage.tsv
	
	#2011-9-14 top 25 contigs, variants only, run on uschpc cluster
	%s -I 559-656 -j uschpc -l uschpc -u yh -a 524 -s 2 -e /home/cmb-03/mn/yuhuang -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/
		-z 10.8.0.10 -o ./AlignmentToCallPipeline_559_656_vs_524_top_25Contigs_uschpc.xml
		-D ~/mnt/hpc-cmb_home/NetworkData/vervet/db/ -N25
	
	# 2011-11-4 run GATK/samtools on single-sample at a time, for 4 high-coverage VRC monkeys, top 804 contigs
	%s -a 524 -i 15-18 -u yh -l condorpool -j condorpool -s 2 -N 804
		-o AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_single_sample_condor.xml -z uclaOffice -n 2
		-O 5
	
	# 2012.6.26 run on hoffman2 condor pool (hcondor), filtered sequences (-Q 1), alignment method 2 (-G 2)
	# no site ID filtering (-S ""), clustering size =5 for calling jobs (-O 5).
	# with 2million bp interval (-Z 2000000).
	%s -a 524 -i 633,1495,...,1524,1459,1505,1478,1486,1442,1472,1516,1453
		-u yh -z localhost -N 7559 -S "" -Q1 -G 2 -l hcondor -j hcondor -e /u/home/eeskin/polyacti 
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db 
		-J /u/home/eeskin/polyacti/bin/jdk/bin/java -O 5 -Z 2000000
		-o workflow/AlignmentToCall_130PoplationVervets_vs_524_top7559Contigs.xml
	
	# 2012.7.30 genotype-call 723 alignments on method 7 sites (-R ...). "-N ..." (top number of contigs) doesn't matter here.
	# 2000 method 7 sites for each calling job (-K 2000)
	%s -a 524 -S 447 -u yh -z localhost -Q1 -G2
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db
		-J /u/home/eeskin/polyacti/bin/jdk/bin/java
		-O 5  -o workflow/AlignmentToCall_AllVRC_vs_524_top1000Contigs.xml -R method7_BED.tsv -K 2000
	
	#2012.8.2 testing calling at known sites (-R method7_BED.n10k.tsv, Contig731 and partial Contig645) 
	#			with 500 regions for each job (-K 500),
	# run GATK along with samtools (-T), no clustering for any job (-O 1 -C1)
	# only on four alignments of isq-id (-i 643-646)
	%s -a 524 -i 643-646 -S 447 -u yh -z localhost  -Q1 -G2 
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db
		-J /u/home/eeskin/polyacti/bin/jdk/bin/java
		-O 1 -C1
		-o workflow/AlignmentToCall_ISQ643_646_vs_524_method7n10kSites.xml -R method7_BED.n10k.tsv -T -K 500
	# 2012.8.2 part of the test above, now run multi-sample on Contig731 on all sites
	# use -x 731 (maxContigID) -V 731 (minContigID) to restrict the top 1000 contigs to only Contig731.
	# run on intervals of 200kb (-Z), with both GATK (-T) and SAMtools
	# add --individual_sequence_file_raw_id_type 2 (library-specific alignments, different libraries of one individual_sequence) 
	# add --individual_sequence_file_raw_id_type 3 (both all-library-fused and library-specific alignments)
	# add "--country_id_ls 135,136,144,148,151" to limit individuals from US,Barbados,StKitts,Nevis,Gambia (AND with -S, )
	%s -a 524 -i 643-646
		#-S 447
		-u yh -z localhost -Q1 -G2
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db -D /u/home/eeskin/polyacti/NetworkData/vervet/db
		-J /u/home/eeskin/polyacti/bin/jdk/bin/java -O 1 -C1
		-o workflow/AlignmentToCall/AlignmentToCall_ISQ643_646_vs_524_Contig731.xml -T -N 1000 -x 731 -V 731 -Z 200000
		#--individual_sequence_file_raw_id_type 2 --country_id_ls 135,136,144,148,151 --tax_id_ls 60711 #sabaeus
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
			It will also stage out every output file.
		If needFastaIndexJob is off, the reference fasta file and its affiliated files will not be staged in.
			If on, the reference fasta file will be staged in and affiliated index/dict files will be created by a job.
		If selectedRegionFname is not set, the calls are made throughout all chosen contigs/chromosomes interval by interval.
			start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
			However, if the former is not null, workflow will make calls in those selected regions in blocks (=10000 regions).
	2012.8.2 --intervalOverlapSize has no effect. intervals are not overlapped for SAMtools/GATK. 
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0],\
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from vervet.src import AbstractAlignmentAndVCFWorkflow, AbstractVervetWorkflow, VervetDB

parentClass = AbstractAlignmentAndVCFWorkflow
class AlignmentToCallPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('inputDir', 0, ))
	
	
	commonCallPipelineOptionDict= {
						#('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("site_type", 1, int): [2, 's', 1, '1: all genome sites, 2: variants only'],\
						("noOfCallingJobsPerNode", 1, int): [1, 'O', 1, 'this parameter controls how many genotype calling jobs should be clustered together to run on one node. \
									Increase it to above 1 only when your average genotyping job is short and the number of input bam files are short.'],\
						}
	
	commonCallPipelineOptionDict.update(parentClass.partitionWorkflowOptionDict.copy())
	commonCallPipelineOptionDict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(commonCallPipelineOptionDict.copy())
	option_default_dict.update({
						("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\
						("run_type", 1, int): [1, 'n', 1, '1: multi-sample calling, 2: single-sample one by one'],\
						("addGATKJobs", 0, int): [0, 'T', 0, 'toggle if you want to have GATK genotyping jobs as well.'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.9.17 call parentClass.__init__() directly
		2011-7-11
		"""
		parentClass.__init__(self, **keywords)
	
	
	def addRefFastaFileSplitJobs(self, workflow, refFastaF, selectAndSplitFasta, chrLs, mkdir=None, samtools=None,
								java=None, createSequenceDictionaryJar=None,\
								site_handler=None, namespace='workflow', version='1.0',):
		"""
		2011-7-25
			split the whole fasta file into files, each containing one fasta record (from chrLs)
			return the data
		"""
		
		# Add a mkdir job for the call directory.
		#letting numerou genotype call jobs detect&create this directory runs into race condition.
		mkDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
		fastaOutputDir = "fasta"
		mkDirJob.addArguments(fastaOutputDir)
		#mkDirJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
		workflow.addJob(mkDirJob)
		
		selectAndSplitFastaJob = Job(namespace=namespace, name=selectAndSplitFasta.name, version=version)
		selectAndSplitFastaJob.addArguments('-i', refFastaF, "-o", fastaOutputDir)
		#selectAndSplitFastaJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		workflow.addJob(selectAndSplitFastaJob)
		workflow.depends(parent=mkDirJob, child=selectAndSplitFastaJob)
		
		chr2jobDataLs = {}
		for chr in chrLs:
			if chr not in chr2jobDataLs:
				chr2jobDataLs[chr] = []
			selectAndSplitFastaJob.addArguments(chr)
			
			fastaFname = os.path.join(fastaOutputDir, '%s.fasta'%(chr))
			fastaFile = File(fastaFname)
			
			# add the index job
			fai_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			fai_index_job.addArguments("faidx", fastaFname)
			fastaFAIIndexFname = os.path.join(fastaOutputDir, '%s.fasta.fai'%(chr))
			fastaFAIIndexFile = File(fastaFAIIndexFname)
			fai_index_job.uses(fastaFAIIndexFile, transfer=True, register=False, link=Link.OUTPUT)	#this file is input & output
			#fai_index_job.uses(output, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
			#fai_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(fai_index_job)
			workflow.depends(parent=selectAndSplitFastaJob, child=fai_index_job)
			
			# the add-read-group job
			createSeqDictJob = Job(namespace=namespace, name=java.name, version=version)
			dictFname = os.path.join(fastaOutputDir, '%s.dict'%(chr))
			#outputRGFname = '%s_%s.RG.bam'%(inputFileBaseNamePrefix, chr)
			dictFile = File(dictFname)
			createSeqDictJob.addArguments('-jar', createSequenceDictionaryJar, \
								"R=", fastaFile, 'O=', dictFile)
			createSeqDictJob.uses(dictFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			createSeqDictJob.uses(fastaFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			workflow.addJob(createSeqDictJob)
			workflow.depends(parent=selectAndSplitFastaJob, child=createSeqDictJob)
			
			
			chr2jobDataLs[chr] = [fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile]
		return PassingData(chr2jobDataLs=chr2jobDataLs, workflow=workflow)
	
	def addSelectAndSplitBamJobs(self, db_vervet, workflow, alignmentLs, site_handler, maxContigID, chrLs, samtools=None, \
							java=None, addOrReplaceReadGroupsAndCleanSQHeaderJar=None, BuildBamIndexFilesJar=None, mkdir=None, namespace="workflow",\
							version="1.0", mkCallDirJob=None,\
							addRGExecutable=None, dataDir=None):
		"""
		2011-7-14
			1. select reference out of whole-alignment
			2. index
			3. add read groups
			4. index
		"""
		sys.stderr.write("Adding Bam select and split jobs for %s alignments ..."%(len(alignmentLs)))
		chr2jobDataLs = {}
		for alignment in alignmentLs:
			# Add input file to the DAX-level replica catalog
			inputFname = os.path.join(dataDir, alignment.path)
			input = File(inputFname)
			input.addPFN(PFN("file://" + inputFname, site_handler))
			workflow.addFile(input)
			
			inputFileBaseNamePrefix = os.path.splitext(os.path.basename(alignment.path))[0]
			outputDir = inputFileBaseNamePrefix
			
			
			# add RG to this bam
			sequencer = alignment.individual_sequence.sequencer
			read_group = alignment.getReadGroup()
			if sequencer=='454':
				platform_id = 'LS454'
			elif sequencer=='GA':
				platform_id = 'ILLUMINA'
			else:
				platform_id = 'ILLUMINA'
			# Add a mkdir job
			mkdirJob = Job(namespace=namespace, name=mkdir.name, version=version)
			
			mkdirJob.addArguments(outputDir)
			#mkdir.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
			workflow.addJob(mkdirJob)
			workflow.depends(parent=mkCallDirJob, child=mkdirJob)
			
			for chr in chrLs:
				if chr not in chr2jobDataLs:
					chr2jobDataLs[chr] = []
				#select reads that are aligned to one reference name
				selectRefJob = Job(namespace=namespace, name=samtools.name, version=version)
				
				outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix, chr))
				#outputFname = '%s_%s.bam'%(inputFileBaseNamePrefix, chr)
				output = File(outputFname)
				selectRefJob.addArguments('view', '-h', input, chr, "-o", output, "-b", "-u")	# -b -u forces uncompressed bam output
				#selectRefJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
				workflow.addJob(selectRefJob)
				workflow.depends(parent=mkdirJob, child=selectRefJob)
				
				# add the index job
				index_1_job = Job(namespace=namespace, name=samtools.name, version=version)
				index_1_job.addArguments("index", output)
				#index_1_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				index_1_job.uses(output, transfer=False, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputFname)
				index_1_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(index_1_job)
				workflow.depends(parent=selectRefJob, child=index_1_job)
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=java.name, version=version)
				#2011-7-27 somehow AddOrReplaceReadGroupsAndCleanSQHeader.jar couldn't output a un-corrupted bam file. so sam first.
				outputRGSAMFname = os.path.join(outputDir, '%s_%s.RG.sam'%(inputFileBaseNamePrefix, chr))
				outputRGSAM = File(outputRGSAMFname)
				"""
				tmpRGFname = os.path.join(outputDir, '%s_%s.RG.txt'%(inputFileBaseNamePrefix, chr))
				addRGJob.addArguments(read_group, platform_id, output, tmpRGFname, outputRGSAM)
				"""
				addRGJob.addArguments('-jar', addOrReplaceReadGroupsAndCleanSQHeaderJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SQName=%s'%(chr))	#'SORT_ORDER=coordinate', (adding this is useless)
				"""
				
				#addRGJob.addArguments('-jar', addOrReplaceReadGroupsJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM)	#'SORT_ORDER=coordinate', (adding this is useless)
				"""
				addRGJob.uses(output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(outputRGSAM, transfer=False, register=True, link=Link.OUTPUT)	#time to discard them
				workflow.addJob(addRGJob)
				workflow.depends(parent=index_1_job, child=addRGJob)
				#output.addPFN(PFN("file://" + outputFname, site_handler))
				#selectAndSplitJob.uses(output, transfer=True, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
				
				samToBamJob = Job(namespace=namespace, name=samtools.name, version=version)
				outputRGFname = os.path.join(outputDir, '%s_%s.RG.bam'%(inputFileBaseNamePrefix, chr))
				outputRG = File(outputRGFname)
				samToBamJob.addArguments('view', '-F4', '-Sbh', "-o", outputRG, "-u", outputRGSAM)	# -b -u forces uncompressed bam output
				samToBamJob.uses(outputRGSAM, transfer=False, register=True, link=Link.INPUT)
				#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
				workflow.addJob(samToBamJob)
				workflow.depends(parent=addRGJob, child=samToBamJob)
				
				# add the index job
				index_sam_job = Job(namespace=namespace, name=java.name, version=version)
				bai_output = File('%s.bai'%outputRGFname)
				index_sam_job.addArguments("-Xms128m", "-Xmx2500m", "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
								"INPUT=", outputRG, \
								"OUTPUT=", bai_output)
				index_sam_job.uses(outputRG, transfer=False, register=False, link=Link.INPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
				index_sam_job.uses(bai_output, transfer=False, register=False, link=Link.OUTPUT)
				
				"""
				samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
				samtools_index_job.addArguments("index", outputRG)
				#samtools_index_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				samtools_index_job.uses(outputRG, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputRGFname)
				samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				"""
				workflow.addJob(index_sam_job)
				workflow.depends(parent=samToBamJob, child=index_sam_job)
				chr2jobDataLs[chr].append((outputRG, index_sam_job, bai_output))
				
				"""
				#2011-9-1 temporary addition to make sure it's sorted for Vasily
				# input bam doesn't need to be indexed.
				outputRGSortSAMFnamePrefix = os.path.join(outputDir, '%s_%s.RG.sorted'%(inputFileBaseNamePrefix, chr))
				outputRGSortSAMFname = File('%s.bam'%(outputRGSortSAMFnamePrefix))
				sam_sort_job = Job(namespace=namespace, name=samtools.name, version=version)
				sam_sort_job.addArguments('sort', '-m', '2000000000', outputRG, outputRGSortSAMFnamePrefix)	#maxMemory is down to 2G
				sam_sort_job.uses(outputRG, transfer=False, register=False, link=Link.INPUT)
				sam_sort_job.uses(outputRGSortSAMFname, transfer=True, register=False, link=Link.OUTPUT)
				job_max_memory = 2000	#in MB
				sam_sort_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				sam_sort_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				workflow.addJob(sam_sort_job)
				workflow.depends(parent=samToBamJob, child=sam_sort_job)
				"""
				
		sys.stderr.write(".Done\n")
		return PassingData(chr2jobDataLs=chr2jobDataLs, workflow=workflow)
	
	def addMergeAlignmentAndGenotypeCallJobs(self, workflow, chr2jobDataLs, chrLs, samtools, \
				java, createSequenceDictionaryJar=None, GenotypeCallByCoverage=None, refFastaF=None, \
				namespace='workflow', version="1.0", callOutputDir = "call", genotypeCallerType=1,\
				mergeSamFilesJar=None, genomeAnalysisTKJar=None, calcula=None, chr2splitFastaJobDataLs=None, \
				seqCoverageF=None, needFastaIndexJob=False, \
				needFastaDictJob=False, site_type=1):
		"""
		2011-9-2
			add argument seqCoverageF
		2011-8-26
			add argument needFastaDictJob, needFastaIndexJob, site_type
		2011-7-26
			add chr2splitFastaJobDataLs, for GATK
		2011-7-14
			1. merge alignments based on same reference from different genomes into one
			2. run GenotypeCallByCoverage on each merged alignment
		"""
		sys.stderr.write("Adding alignment merge and genotype call jobs for %s references ..."%(len(chr2jobDataLs)))
		if needFastaDictJob:	# the .dict file is required for GATK
			refFastaDictFname = '%s.dict'%(os.path.splitext(refFastaF.name)[0])
			refFastaDictF = File(refFastaDictFname)
			#not os.path.isfile(refFastaDictFname) or 
			fastaDictJob = Job(namespace=namespace, name=java.name, version=version)
			fastaDictJob.addArguments('-jar', createSequenceDictionaryJar, \
					'REFERENCE=', refFastaF, 'OUTPUT=', refFastaDictF)
			fastaDictJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaDictJob.uses(refFastaDictF, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaDictJob)
		else:
			fastaDictJob = None
			refFastaDictF = None
		
		if needFastaIndexJob:
			refFastaIndexFname = '%s.fai'%(refFastaF.name)	# the .fai file is required for GATK
			refFastaIndexF = File(refFastaIndexFname)
			#not os.path.isfile(refFastaIndexFname)
			fastaIndexJob = Job(namespace=namespace, name=samtools.name, version=version)
			fastaIndexJob.addArguments("faidx", refFastaF)
			fastaIndexJob.uses(refFastaF, register=False, link=Link.INPUT)
			fastaIndexJob.uses(refFastaIndexFname, transfer=True, register=False, link=Link.OUTPUT)
			workflow.addJob(fastaIndexJob)
		else:
			fastaIndexJob = None
			refFastaIndexF = None
		#workflow.depends(parent=fastaDictJob, child=fastaIndexJob)	#no dependency between these two jobs
		
		# add merge jobs for every reference
		chr2mergedBamCallJob = {}
		for chr, jobDataLs in chr2jobDataLs.iteritems():
			# add the index job
			picard_job = Job(namespace=namespace, name=java.name, version=version)
			outputFname = '%s.bam'%(chr)
			picard_output = File(outputFname)
			picard_job.addArguments('-jar', mergeSamFilesJar, \
				'USE_THREADING=true', 'SORT_ORDER=coordinate', 'ASSUME_SORTED=false', 'OUTPUT=', picard_output)
			#picard_job.uses(picard_output, transfer=False, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
			workflow.addJob(picard_job)
			
			for jobData in jobDataLs:
				bamFile, samtools_index_job, bamFileBai = jobData[:3]
				picard_job.addArguments('INPUT=', bamFile)
				picard_job.uses(bamFileBai, transfer=True, register=True, link=Link.INPUT)	#register them here to be deleted 
				picard_job.uses(bamFile, transfer=True, register=True, link=Link.INPUT)
				#this picard merge job depends on a bunch of prior samtools index jobs
				workflow.depends(parent=samtools_index_job, child=picard_job)
			
			# add the index job on the merged bam file
			samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			samtools_index_job.addArguments("index", picard_output)
			samtools_index_job.uses(picard_output, transfer=True, register=True, link=Link.OUTPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
			bai_output = File('%s.bai'%outputFname)
			samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(samtools_index_job)
			workflow.depends(parent=picard_job, child=samtools_index_job)
			
			
			if genotypeCallerType==2:
				coverageFilterParentJob = samtools_index_job
				coverageFilterParentOutput = picard_output
				coverageFilterParentOutput_bai = bai_output
			else:
				splitFastaJobDataLs = chr2splitFastaJobDataLs.get(chr)
				fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile = splitFastaJobDataLs[:5]
				
				gatk_job = Job(namespace=namespace, name=java.name, version=version)
				gatkOutputFname = os.path.join(callOutputDir, '%s.vcf'%(chr))
				gatk_output = File(gatkOutputFname)
				gatkIDXOutputFname = os.path.join(callOutputDir, '%s.vcf.idx'%(chr))
				gatkIDXOutput = File(gatkIDXOutputFname)
				
				gatk_job.addArguments('-jar', genomeAnalysisTKJar, \
					"-I", picard_output, "-R", fastaFile, "-T", "UnifiedGenotyper","--out", gatk_output,\
					'-U', '-S SILENT', "-nt 4")
				if site_type==1:
					gatk_job.addArguments('--output_mode EMIT_ALL_SITES')	#2011-8-24 new GATK no longers ues "-all_bases"
				gatk_job.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
				gatk_job.uses(picard_output, transfer=False, register=True, link=Link.INPUT)
				
				gatk_job.uses(fastaFAIIndexFile, transfer=True, register=True, link=Link.INPUT)
				gatk_job.uses(fastaFile, transfer=True, register=True, link=Link.INPUT)
				gatk_job.uses(dictFile, transfer=False, register=True, link=Link.INPUT)
				
				gatk_job.uses(gatk_output, transfer=True, register=True, link=Link.OUTPUT)
				gatk_job.uses(gatkIDXOutput, transfer=True, register=True, link=Link.OUTPUT)
				
				
				gatk_job_max_memory = 1500	#in MB
				yh_pegasus.setJobProperRequirement(gatk_job, job_max_memory=gatk_job_max_memory)
				
				workflow.addJob(gatk_job)
				workflow.depends(parent=samtools_index_job, child=gatk_job)
				workflow.depends(parent=fai_index_job, child=gatk_job)
				workflow.depends(parent=createSeqDictJob, child=gatk_job)
				if fastaIndexJob:	#2011-7-22 if job doesn't exist, don't add it. means this job isn't necessary to run.
					workflow.depends(parent=fastaIndexJob, child=gatk_job)
				if fastaDictJob:
					workflow.depends(parent=fastaDictJob, child=gatk_job)
				
				coverageFilterParentJob = gatk_job
				coverageFilterParentOutput = gatk_output
				coverageFilterParentOutput_bai = None
				
			
			#add the cover filter (or filter+call) job after index is done
			GenotypeCallByCoverage_job = Job(namespace=namespace, name=GenotypeCallByCoverage.name, version=version)
			genotypeCallOutputFname = os.path.join(callOutputDir, '%s.call'%(chr))	#GenotypeCallByCoverage_job would create directory "call".
			genotypeCallOutput = File(genotypeCallOutputFname)
			GenotypeCallByCoverage_job.addArguments("-i", coverageFilterParentOutput, "-n", str(len(jobDataLs)), \
					"-o", genotypeCallOutput, '-e', refFastaF, '-y', str(genotypeCallerType), \
					'-s', repr(self.site_type))
			if seqCoverageF:
				GenotypeCallByCoverage_job.addArguments("-q", seqCoverageF)
				GenotypeCallByCoverage_job.uses(seqCoverageF, transfer=True, register=True, link=Link.INPUT)
			if coverageFilterParentOutput_bai:
				GenotypeCallByCoverage_job.uses(coverageFilterParentOutput_bai, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
			GenotypeCallByCoverage_job.uses(coverageFilterParentOutput, transfer=True, register=True, link=Link.INPUT)
			GenotypeCallByCoverage_job.uses(genotypeCallOutput, transfer=True, register=True, link=Link.OUTPUT)
			if self.site_type==1:	#all sites require additional memory
				job_max_memory = 5000	#in MB
			else:	#variants only, less memory
				job_max_memory=1000	#in MB. 
			yh_pegasus.setJobProperRequirement(GenotypeCallByCoverage_job, job_max_memory=job_max_memory)
			workflow.addJob(GenotypeCallByCoverage_job)
			workflow.depends(parent=coverageFilterParentJob, child=GenotypeCallByCoverage_job)
			
			
			#add the pairwise distance matrix job after filter is done
			calcula_job = Job(namespace=namespace, name=calcula.name, version=version)
			
			convertHetero2NA = True
			max_NA_rate = 0.4
			min_MAF = 0
			calculaOutputFname ='%s.pairwiseDist.convertHetero2NA%s.minMAF%s.maxNA%s.tsv'%(os.path.splitext(genotypeCallOutputFname)[0], \
								convertHetero2NA, min_MAF, max_NA_rate)
			calculaOutput = File(calculaOutputFname)
			calcula_job.addArguments("-i", genotypeCallOutput, "-n", str(min_MAF), \
								"-o", calculaOutput, '-m', repr(max_NA_rate), '-c', str(1))
			calcula_job.uses(genotypeCallOutput, transfer=True, register=False, link=Link.INPUT)
			calcula_job.uses(calculaOutput, transfer=True, register=False, link=Link.OUTPUT)
			
			workflow.addJob(calcula_job)
			workflow.depends(parent=GenotypeCallByCoverage_job, child=calcula_job)
			
			chr2mergedBamCallJob[chr] = [genotypeCallOutput, GenotypeCallByCoverage_job]
		sys.stderr.write(".Done\n")
		return PassingData(chr2mergedBamCallJob=chr2mergedBamCallJob)
	
	def addAddRG2BamJobsAsNeeded(self, workflow, alignmentDataLs, site_handler, input_site_handler=None, \
							addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
							BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
							mv=None, namespace="workflow", version="1.0", \
							dataDir=None, tmpDir="/tmp"):
		"""
		2012.4.5
			fix some bugs here
		2011-9-15
			add a read group only when the alignment doesn't have it according to db record
			DBVervet.pokeBamReadGroupPresence() from misc.py helps to fill in db records if it's unclear.
		2011-9-14
			The read-group adding jobs will have a "move" part that overwrites the original bam&bai if site_handler and input_site_handler is same.
			For those alignment files that don't need to. It doesn't matter. pegasus will transfer/symlink them.
		"""
		sys.stderr.write("Adding add-read-group2BAM jobs for %s alignments if read group is not detected ..."%(len(alignmentDataLs)))
		job_max_memory = 3500	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		indexJobMaxMem=2500
		
		addRG2BamDir = None
		addRG2BamDirJob = None
		
		no_of_rg_jobs = 0
		returnData = []
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			if alignment.read_group_added!=1:
				if addRG2BamDir is None:
					addRG2BamDir = "addRG2Bam"
					addRG2BamDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=addRG2BamDir)
				
				# add RG to this bam
				sequencer = alignment.individual_sequence.sequencer
				#read_group = '%s_%s_%s_%s_vs_%s'%(alignment.id, alignment.ind_seq_id, alignment.individual_sequence.individual.code, \
				#						sequencer, alignment.ref_ind_seq_id)
				read_group = alignment.getReadGroup()	##2011-11-02
				if sequencer=='454':
					platform_id = 'LS454'
				elif sequencer=='GA':
					platform_id = 'ILLUMINA'
				else:
					platform_id = 'ILLUMINA'
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=addOrReplaceReadGroupsJava.name, version=version)
				outputRGSAM = File(os.path.join(addRG2BamDir, os.path.basename(alignment.path)))
				
				addRGJob.addArguments(javaMemRequirement, '-jar', addOrReplaceReadGroupsJar, \
									"INPUT=", bamF,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SORT_ORDER=coordinate', "VALIDATION_STRINGENCY=LENIENT")
									#(adding the SORT_ORDER doesn't do sorting but it marks the header as sorted so that BuildBamIndexFilesJar won't fail.)
				if tmpDir:
					addRGJob.addArguments("TMP_DIR=%s"%tmpDir)
				addRGJob.uses(bamF, transfer=True, register=True, link=Link.INPUT)
				addRGJob.uses(baiF, transfer=True, register=True, link=Link.INPUT)
				addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
				yh_pegasus.setJobProperRequirement(addRGJob, job_max_memory=job_max_memory)
				for parentJob in parentJobLs:
					if parentJob:
						workflow.depends(parent=parentJob, child=addRGJob)
				workflow.addJob(addRGJob)
				
				
				index_sam_job = self.addBAMIndexJob(workflow, BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					inputBamF=outputRGSAM, parentJobLs=[addRGJob], transferOutput=True, javaMaxMemory=2000)
				newAlignmentData = PassingData(alignment=alignment)
				newAlignmentData.jobLs = [index_sam_job, addRGJob]
				newAlignmentData.bamF = index_sam_job.bamFile
				newAlignmentData.baiF = index_sam_job.baiFile
				"""
				# add the index job to the bamF (needs to be re-indexed)
				index_sam_job = Job(namespace=namespace, name=BuildBamIndexFilesJava.name, version=version)
				
				if input_site_handler==site_handler:	#on the same site. overwrite the original file without RG
					mvJob = Job(namespace=namespace, name=mv.name, version=version)
					mvJob.addArguments(outputRGSAM, inputFname)	#watch, it's inputFname, not input. input is in relative path.
					#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
					workflow.addJob(mvJob)
					workflow.depends(parent=addRGJob, child=mvJob)
					bai_output = File('%s.bai'%inputFname)	#in absolute path, don't register it to the job
				else:
					##on different site, input for index should be outputRGSAM and register it as well
					mvJob = addRGJob
					bamF = outputRGSAM	
					addRGJob.uses(outputRGSAM, transfer=True, register=True, link=Link.OUTPUT)
					bai_output = File('%s.bai'%outputRGSAMFname)
					index_sam_job.uses(bai_output, transfer=True, register=False, link=Link.OUTPUT)
				
				index_sam_job.addArguments("-Xms128m -Xmx%sm"%indexJobMaxMem, "-jar", BuildBamIndexFilesJar, "VALIDATION_STRINGENCY=LENIENT", \
								"INPUT=", bamF, \
								"OUTPUT=", bai_output)
				# bamF is in relative path, either symlink to inputFname (same site) or outputRGSAM (different site)
				index_sam_job.uses(bamF, transfer=True, register=False, link=Link.INPUT)
				
				yh_pegasus.setJobProperRequirement(index_sam_job, job_max_memory=indexJobMaxMem)
				
				workflow.addJob(index_sam_job)
				workflow.depends(parent=mvJob, child=index_sam_job)
				alignmentId2RGJobDataLs[alignment.id]= [index_sam_job, input, bai_output]
				"""
				no_of_rg_jobs += 1
			else:
				newAlignmentData = alignmentData
			returnData.append(newAlignmentData)
		sys.stderr.write(" %s alignments need read-group addition. Done\n"%(no_of_rg_jobs))
		return returnData
	
	def addGenotypeCallJobs(self, workflow, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				genotyperJava=None,  genomeAnalysisTKJar=None, \
				addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				mergeSamFilesJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, CallVariantBySamtools=None,\
				bgzip_tabix=None, vcf_convert=None, vcf_isec=None, vcf_concat=None, \
				concatGATK=None, concatSamtools=None,\
				GenotypeCallByCoverage=None, refFastaFList=None, bamListF=None, \
				callOutputDirJob =None, gatkDirJob=None, samtoolsDirJob=None, unionDirJob=None, intersectionDirJob=None,\
				namespace='workflow', version="1.0", site_handler=None, input_site_handler=None,\
				seqCoverageF=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				intervalSize=2000000, site_type=1, dataDir=None, no_of_gatk_threads = 1,\
				outputDirPrefix="", addGATKJobs=False, cumulativeMedianDepth=5000, **keywords):
		"""
		2012.7.31
			add argument addGATKJobs (default=False)
		2011-9-22
			add argument concatGATK, concatSamtools.
		2011-9-15
			bamListF is now useless. samtools_job could accept variable-length list of bam input files
		2011-9-14
			argument intervalSize determines how many sites gatk/samtools works on at a time
		"""
		sys.stderr.write("Adding genotype call jobs for %s chromosomes/contigs ..."%(len(chr2IntervalDataLs)))
		job_max_memory = 2000	#in MB
		javaMemRequirement = "-Xms128m -Xmx%sm"%job_max_memory
		vcf_job_max_memory = 1000
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
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(workflow, alignmentDataLs, site_handler, input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
					mv=mv, namespace=namespace, version=version, dataDir=dataDir)
		#alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		
		# add merge jobs for every reference
		chr2mergedBamCallJob = {}
		returnData = PassingData()
		for chr, intervalDataLs in chr2IntervalDataLs.iteritems():
			#reduce the number of chunks 1 below needed. last trunk to reach the end of contig
			#however set it to 1 for contigs smaller than intervalSize 	
			callOutputFname = os.path.join(callOutputDirJob.folder, '%s.vcf.gz'%chr)
			callOutputF = File(callOutputFname)
			"""
			#2012.7.21 no more gatk-samtools intersection job
			wholeRefUnionOfIntersectionJob = self.addVCFConcatJob(workflow, concatExecutable=vcf_concat, \
							parentDirJob=callOutputDirJob, \
							outputF=callOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			"""
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			if addGATKJobs:
				#2011-9-22 union of all GATK intervals for one contig
				gatkUnionOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf.gz'%chr)
				gatkUnionOutputF = File(gatkUnionOutputFname)
				gatkUnionJob = self.addVCFConcatJob(workflow, concatExecutable=concatGATK, parentDirJob=gatkDirJob, \
								outputF=gatkUnionOutputF, namespace=namespace, version=version, transferOutput=True, \
								vcf_job_max_memory=vcf_job_max_memory)
			
			
			#2011-9-22 union of all samtools intervals for one contig
			samtoolsUnionOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf.gz'%chr)
			samtoolsUnionOutputF = File(samtoolsUnionOutputFname)
			samtoolsUnionJob = self.addVCFConcatJob(workflow, concatExecutable=concatSamtools, parentDirJob=samtoolsDirJob, \
							outputF=samtoolsUnionOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			
			
			#2011-9-22 union of all samtools intervals for one contig
			samtoolsIndelUnionOutputFname = os.path.join(samtoolsDirJob.folder, '%s.indel.vcf.gz'%chr)
			samtoolsIndelUnionOutputF = File(samtoolsIndelUnionOutputFname)
			samtoolsIndelUnionJob = self.addVCFConcatJob(workflow, concatExecutable=concatSamtools, parentDirJob=samtoolsDirJob, \
							outputF=samtoolsIndelUnionOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			no_of_jobs += 4
			for intervalData in intervalDataLs:
				if intervalData.file:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.file
				else:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.interval
				vcfBaseFname = intervalData.intervalFnameSignature
				
				if addGATKJobs:
					#GATK job
					gatkOutputFname = os.path.join(gatkDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
					gatkOutputF = File(gatkOutputFname)
					gatkIDXOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf.idx'%(vcfBaseFname))
					gatkIDXOutput = File(gatkIDXOutputFname)
					
					gatk_job= self.addGATKCallJob(workflow, genotyperJava=genotyperJava, genomeAnalysisTKJar=genomeAnalysisTKJar, \
							gatkOutputF=gatkOutputF, gatkIDXOutput=gatkIDXOutput, refFastaFList=refFastaFList, \
							parentJobLs=[gatkDirJob]+intervalData.jobLs, \
							extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
							job_max_memory=job_max_memory*max(1, len(alignmentDataLs)/100), no_of_gatk_threads=no_of_gatk_threads, site_type=site_type, \
							interval=bcftoolsInterval)
					
					vcf4_gatkOutputFname = os.path.join(gatkDirJob.folder, '%s.vcf'%vcfBaseFname)
					vcf4_gatkOutputF = File(vcf4_gatkOutputFname)
					vcf_convert_gatkOutputF_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
								parentJob=gatk_job, inputF=gatkOutputF, outputF=vcf4_gatkOutputF, \
								namespace=namespace, version=version, transferOutput=False)
					
					gatkGzipOutputF = File("%s.gz"%vcf4_gatkOutputFname)
					gatkGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_gatkOutputFname)
					bgzip_tabix_gatkOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
							parentJob=vcf_convert_gatkOutputF_job, inputF=vcf4_gatkOutputF, outputF=gatkGzipOutputF, \
							namespace=namespace, version=version, transferOutput=False)
					
					#add this output to a GATK union job
					gatkUnionJob.addArguments(gatkGzipOutputF)
					gatkUnionJob.uses(gatkGzipOutputF, transfer=False, register=True, link=Link.INPUT)
					gatkUnionJob.uses(gatkGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
					workflow.depends(parent=bgzip_tabix_gatkOutputF_job, child=gatkUnionJob)
				
				#samtools part
				samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
				samtoolsOutputF = File(samtoolsOutputFname)
				indelVCFOutputFname = "%s.indel.vcf"%(samtoolsOutputFname)
				indelVCFOutputF = File(indelVCFOutputFname)
				samtools_job= self.addSAMtoolsCallJob(workflow, CallVariantBySamtools=CallVariantBySamtools, \
					samtoolsOutputF=samtoolsOutputF, indelVCFOutputF=indelVCFOutputF, \
					refFastaFList=refFastaFList, parentJobLs=[samtoolsDirJob]+intervalData.jobLs, \
					extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=job_max_memory*max(1, len(alignmentDataLs)/200), \
					site_type=site_type, maxDP=cumulativeMedianDepth*5, mpileupInterval=mpileupInterval, \
					bcftoolsInterval=bcftoolsInterval)
				
				#deal with samtools's indel VCF
				"""
				# convert to VCF4. vcf-convert requires "-r refFastaF" but still fails in sanity check.
				# complaining ref sequence from this vcf mismatches from the reference sequence. so comment this out. 
				samtoolsIndelVCF4Fname = os.path.join(samtoolsDirJob.folder, '%s.indel.vcf'%vcfBaseFname)
				samtoolsIndelVCF4F = File(samtoolsIndelVCF4Fname)
				samtoolsIndelVCF4F_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=samtools_job, inputF=indelVCFOutputF, outputF=samtoolsIndelVCF4F, \
							namespace=namespace, version=version, transferOutput=False)
				"""
				samtoolsIndelGzipOutputF = File("%s.gz"%indelVCFOutputFname)
				samtoolsIndelGzipOutput_tbi_F = File("%s.gz.tbi"%indelVCFOutputFname)
				samtoolsIndelGzipOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=samtools_job, inputF=indelVCFOutputF, outputF=samtoolsIndelGzipOutputF, \
						namespace=namespace, version=version, transferOutput=False)
				
				
				samtoolsIndelUnionJob.addArguments(samtoolsIndelGzipOutputF)
				samtoolsIndelUnionJob.uses(samtoolsIndelGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				samtoolsIndelUnionJob.uses(samtoolsIndelGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=samtoolsIndelGzipOutputF_job, child=samtoolsIndelUnionJob)
				
				#deal with samtools's snp VCF
				vcf4_samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf'%vcfBaseFname)
				vcf4_samtoolsOutputF = File(vcf4_samtoolsOutputFname)
				vcf_convert_samtoolsOutputF_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
							parentJob=samtools_job, inputF=samtoolsOutputF, outputF=vcf4_samtoolsOutputF, \
							namespace=namespace, version=version, transferOutput=False)
				
				samtoolsGzipOutputF = File("%s.gz"%vcf4_samtoolsOutputFname)
				samtoolsGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_samtoolsOutputFname)
				bgzip_tabix_samtoolsOutputF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
						parentJob=vcf_convert_samtoolsOutputF_job, inputF=vcf4_samtoolsOutputF, outputF=samtoolsGzipOutputF, \
						namespace=namespace, version=version, transferOutput=False)
				
				#add this output to a samtools union job
				samtoolsUnionJob.addArguments(samtoolsGzipOutputF)
				samtoolsUnionJob.uses(samtoolsGzipOutputF, transfer=False, register=True, link=Link.INPUT)
				samtoolsUnionJob.uses(samtoolsGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=bgzip_tabix_samtoolsOutputF_job, child=samtoolsUnionJob)
				
				"""
				#2012.7.21 no more intersection
				# intersection of GATK and samtools
				isectJob = Job(namespace=namespace, name=vcf_isec.name, version=version)
				intersectionOutputFname = os.path.join(intersectionDirJob.folder, '%s.isec.vcf.gz'%vcfBaseFname)
				intersectionOutputF = File(intersectionOutputFname)
				intersectionOutput_tbi_F = File("%s.tbi"%intersectionOutputFname)
				
				#put samtools in front of GATK cuz the former outputs better-quality data based on inconsistent rates
				#	estimated from trios
				isectJob.addArguments(intersectionOutputF, samtoolsGzipOutputF, gatkGzipOutputF)
				isectJob.uses(samtoolsGzipOutputF, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(samtoolsGzipOutput_tbi_F, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(gatkGzipOutputF, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(gatkGzipOutput_tbi_F, transfer=True, register=True, link=Link.INPUT)
				isectJob.uses(intersectionOutputF, transfer=False, register=True, link=Link.OUTPUT)
				isectJob.uses(intersectionOutput_tbi_F, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(isectJob)
				workflow.depends(parent=intersectionDirJob, child=isectJob)
				workflow.depends(parent=bgzip_tabix_gatkOutputF_job, child=isectJob)
				workflow.depends(parent=bgzip_tabix_samtoolsOutputF_job, child=isectJob)
				
				yh_pegasus.setJobProperRequirement(isectJob, job_max_memory=vcf_job_max_memory)
				#union of intersected-intervals from one contig
				wholeRefUnionOfIntersectionJob.addArguments(intersectionOutputF)
				wholeRefUnionOfIntersectionJob.uses(intersectionOutputF, transfer=False, register=True, link=Link.INPUT)
				wholeRefUnionOfIntersectionJob.uses(intersectionOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
				workflow.depends(parent=isectJob, child=wholeRefUnionOfIntersectionJob)
				"""
				
				"""
				java -jar /path/to/dist/GenomeAnalysisTK.jar \
				  -T CombineVariants \
				  -R /path/to/reference.fasta \
				  -o /path/to/output.vcf \
				  -V:foo /path/to/variants1.vcf \
				  -V:bar /path/to/variants2.vcf \
				  -genotypeMergeOptions PRIORITIZE \
				
				java -Xmx2g -jar dist/GenomeAnalysisTK.jar -T SelectVariants 
				-R ~/Desktop/broadLocal/localData/human_g1k_v37.fasta 
				-L 1:1-1,000,000 --variant union.vcf 
				-select 'set == "Intersection"' -o intersect.vcf
				
				#CombineVariants has problem working with samtools vcf output in which GT is not placed in the first position.
				
				genotypeUnionJob = Job(namespace=namespace, name=java.name, version=version)
				unionOutputFname = os.path.join(unionDirJob.folder, '%s.vcf'%vcfBaseFname)
				unionOutputF = File(unionOutputFname)
				genotypeUnionJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "CombineVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:gatk,VCF", gatkOutputF, "-B:samtools,VCF", samtoolsOutputF, \
					"--out", unionOutputF, \
					'-U', '-S LENIENT', "-genotypeMergeOptions PRIORITIZE")	#no support for multi-threading
				genotypeUnionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeUnionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				genotypeUnionJob.uses(gatkOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(samtoolsOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(unionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeUnionJob)
				workflow.depends(parent=unionDirJob, child=genotypeUnionJob)
				workflow.depends(parent=samtools_job, child=genotypeUnionJob)
				workflow.depends(parent=gatk_job, child=genotypeUnionJob)
				
				genotypeIntersectionJob = Job(namespace=namespace, name=java.name, version=version)
				intersectionOutputFname = os.path.join(intersectionDirJob.folder, '%s.vcf'%vcfBaseFname)
				intersectionOutputF = File(intersectionOutputFname)
				genotypeIntersectionJob.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "SelectVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:%s,VCF"%vcfBaseFname, unionOutputF, \
					"-select 'set == "Intersection"'", "--out", intersectionOutputF, \
					'-U', '-S LENIENT')	#no support for multi-threading
				genotypeIntersectionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeIntersectionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				genotypeIntersectionJob.uses(unionOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeIntersectionJob.uses(intersectionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeIntersectionJob)
				workflow.depends(parent=intersectionDirJob, child=genotypeIntersectionJob)
				workflow.depends(parent=genotypeUnionJob, child=genotypeIntersectionJob)
				
				#add this intersected output to the union job for the whole chr
				wholeRefUnionOfIntersectionJob.addArguments("-B:%s,VCF"%vcfBaseFname, intersectionOutputF)
				workflow.depends(parent=genotypeIntersectionJob, child=wholeRefUnionOfIntersectionJob)
				"""
				no_of_jobs += 8
				jobAndInputOptionList = [(samtools_job, "")]
				if addGATKJobs:
					jobAndInputOptionList.append((gatk_job, "-I"))
				for job, jobInputOption in jobAndInputOptionList:
					self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
					self.addAlignmentAsInputToJobLs(workflow, alignmentDataLs, jobLs=[job], jobInputOption=jobInputOption)
		
		
		sys.stderr.write(" %s jobs.\n"%(no_of_jobs))
		return returnData
	
	def addAlignmentAsInputToJobLs(self, workflow, alignmentDataLs, jobLs=[], jobInputOption=""):
		"""
		2012.1.9
			used in addGenotypeCallJobs() to add alignment files as input to calling jobs
		"""
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			for job in jobLs:
				if jobInputOption:
					job.addArguments(jobInputOption)
				job.addArguments(bamF)
				#it's either symlink or stage-in
				job.uses(bamF, transfer=True, register=True, link=Link.INPUT)
				job.uses(baiF, transfer=True, register=True, link=Link.INPUT)
				for parentJob in parentJobLs:
					if parentJob:
						workflow.depends(parent=parentJob, child=job)
		
	
	def outputSeqCoverage(self, outputFname, ind_seq_id_ls=[]):
		"""
		2011-9-2
		"""
		sys.stderr.write("Outputting sequence coverage to %s ..."%outputFname)
		import csv
		counter = 0
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['isq.id', 'coverage'])
		TableClass = VervetDB.IndividualSequence
		query = TableClass.query.filter(TableClass.coverage!=None)
		if ind_seq_id_ls:
			query = query.filter(TableClass.id.in_(ind_seq_id_ls))
		query = query.order_by(TableClass.id)
		for row in query:
			writer.writerow([row.id, row.coverage])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
	def outputListOfBam(self, alignmentLs, dataDir=None, outputFname=None, workflow=None, site_handler=None):
		"""
		2011-9-15
			deprecated. CallVariantBySamtools.sh doesn't need this file anymore.
			using relative path, alignments are symlinked or transferred
		2011-9-14
			samtools mpileup uses this file.
		"""
		sys.stderr.write("Outputting path to %s alignments to %s ..."%(len(alignmentLs), outputFname))
		outf = open(outputFname, 'w')
		for alignment in alignmentLs:
			#abs_path = os.path.join(dataDir, alignment.path)
			outf.write("%s\n"%alignment.path)	#using relative path, alignments are symlinked or transferred
		del outf
		bamListF = File(os.path.basename(outputFname))
		bamListF.addPFN(PFN("file://" + os.path.abspath(outputFname), site_handler))
		workflow.addFile(bamListF)
		sys.stderr.write("Done.\n")
		return bamListF
	
	def addMkDirJob(self, workflow, mkdir=None, outputDir=None, namespace=None, version=None):
		"""
		2011-9-14
		"""
		return yh_pegasus.addMkDirJob(workflow, mkdir=mkdir, outputDir=outputDir, namespace=namespace, version=version)
		"""
		# Add a mkdir job for any directory.
		mkDirJob = Job(namespace=namespace, name=mkdir.name, version=version)
		mkDirJob.addArguments(outputDir)
		mkDirJob.folder = outputDir	#custom attribute
		workflow.addJob(mkDirJob)
		return mkDirJob
		"""
	
	def addGATKCallJob(self, workflow, genotyperJava=None, genomeAnalysisTKJar=None, gatkOutputF=None, gatkIDXOutput=None, \
					refFastaFList=[], parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
					extraArguments=None, job_max_memory=2000, no_of_gatk_threads=1, site_type=2, interval=None, **keywords):
		"""
		2012.7.30
			interval could be a BED file, rather than just a string (chr:start-stop).
			start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
			
		2011-12-4
		"""
		#GATK job
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		#= "-Xms%sm -Xmx%sm"%(job_max_memory/2, job_max_memory)
		refFastaF = refFastaFList[0]
		job = Job(namespace=workflow.namespace, name=genotyperJava.name, version=workflow.version)
		job.addArguments(javaMemRequirement, '-jar', genomeAnalysisTKJar, "-T", "UnifiedGenotyper",\
			"-mbq 20", "-R", refFastaF, "--out", gatkOutputF,\
			self.defaultGATKArguments, "-nt %s"%no_of_gatk_threads, "--baq CALCULATE_AS_NECESSARY")
		# 2011-11-22 "-mmq 30",  is no longer included
		if site_type==1:
			job.addArguments('--output_mode EMIT_ALL_SITES')	#2011-8-24 new GATK no longers ues "-all_bases"
		if extraArguments:
			job.addArguments(extraArguments)
		if interval is not None:
			if isinstance(interval, File):	#2012.7.30
				job.uses(interval, transfer=True, register=True, link=Link.INPUT)
				job.addArguments("-L:bed", interval)
			else:
				job.addArguments("-L", interval)
		for refFastaFile in refFastaFList:
			job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		job.uses(gatkOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.uses(gatkIDXOutput, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = gatkOutputF
		job.gatkIDXOutput = gatkIDXOutput
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, no_of_cpus=no_of_gatk_threads)
		workflow.addJob(job)
		
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=parentJob, child=job)
		if extraDependentInputLs:
			for input in extraDependentInputLs:
				if input:
					job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def addSAMtoolsCallJob(self, workflow, CallVariantBySamtools=None, samtoolsOutputF=None, indelVCFOutputF=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, site_type=2, maxDP=15000, mpileupInterval=None, \
					bcftoolsInterval=None, **keywords):
		"""
		2012.8.7
			split argument interval into mpileupInterval and bcftoolsInterval
		2012.8.6 add argument maxDP
		2012.7.30
			mpileupInterval could be a BED file, rather than just a string (chr:start-stop). CallVariantBySamtools would adjust itself.
			start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
		2011-12-4
		"""
		refFastaF = refFastaFList[0]
		job = Job(namespace=workflow.namespace, name=CallVariantBySamtools.name, version=workflow.version)
		job.addArguments(refFastaF, mpileupInterval, bcftoolsInterval, samtoolsOutputF, repr(site_type), "%s"%(maxDP))
		if isinstance(mpileupInterval, File):
			job.uses(mpileupInterval, transfer=True, register=True, link=Link.INPUT)
		if isinstance(bcftoolsInterval, File) and mpileupInterval!=bcftoolsInterval:	#do not add the same file as INPUT twice
			job.uses(bcftoolsInterval, transfer=True, register=True, link=Link.INPUT)
		for refFastaFile in refFastaFList:
			job.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		#job.uses(bamListF, transfer=True, register=True, link=Link.INPUT)
		job.uses(samtoolsOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.uses(indelVCFOutputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = samtoolsOutputF
		job.indelVCFOutputF = indelVCFOutputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		workflow.addJob(job)
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=parentJob, child=job)
		if extraDependentInputLs:
			for input in extraDependentInputLs:
				if input:
					job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2012.6.1
			self has become a workflow. added MergeVCFReplicateHaplotypesJava, ligateVcf
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableList = []
		genotypingExecutableSet = set([self.genotyperJava, self.samtools, self.CallVariantBySamtools, self.MergeVCFReplicateHaplotypesJava])
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		trioCallerWrapper = Executable(namespace=namespace, name="trioCallerWrapper", version=version, os=operatingSystem,\
									arch=architecture, installed=True)
		trioCallerWrapper.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/trioCallerWrapper.sh"), site_handler))
		executableClusterSizeMultiplierList.append((trioCallerWrapper, 0))
		genotypingExecutableSet.add(trioCallerWrapper)
		
		MergeVCFReplicateGenotypeColumnsJava = Executable(namespace=namespace, name="MergeVCFReplicateGenotypeColumnsJava", \
											version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		MergeVCFReplicateGenotypeColumnsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		executableClusterSizeMultiplierList.append((MergeVCFReplicateGenotypeColumnsJava, 0))
		genotypingExecutableSet.add(MergeVCFReplicateGenotypeColumnsJava)
		
		ligateVcf = Executable(namespace=namespace, name="ligateVcf", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		ligateVcf.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/ligateVcf.sh"), site_handler))
		executableClusterSizeMultiplierList.append((ligateVcf, 1))
		
		ReplicateVCFGenotypeColumns = Executable(namespace=namespace, name="ReplicateVCFGenotypeColumns", version=version, os=operatingSystem,\
									arch=architecture, installed=True)
		ReplicateVCFGenotypeColumns.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ReplicateVCFGenotypeColumns.py"), site_handler))
		executableClusterSizeMultiplierList.append((ReplicateVCFGenotypeColumns, 1))
	
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		
		if self.noOfCallingJobsPerNode>1:
			for executable in genotypingExecutableSet:
				clusterinProf = Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.noOfCallingJobsPerNode)
				if not executable.hasProfile(clusterinProf):
					executable.addProfile(clusterinProf)
	
	def addLigateVcfJob(self, executable=None, ligateVcfPerlPath=None, outputFile=None, \
						parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.6.1
			ligateVcf ligates overlapping VCF files.
		"""
		extraArgumentList = [ligateVcfPerlPath, outputFile]
		if extraArguments:
			extraArgumentList.append(extraArguments)
		#do not pass outputFile as argument to addGenericJob() because addGenericJob will add "-o" in front of it.
		return self.addGenericJob(executable=executable, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, extraOutputLs=[outputFile], \
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
	
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
		
		alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, dataDir=self.localDataDir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
									individual_site_id_set=set(self.site_id_ls),\
									mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls),\
									excludeContaminant=self.excludeContaminant)
		
		cumulativeMedianDepth = db_vervet.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
										defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		vervetSrcPath = self.vervetSrcPath
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
		
		if self.selectedRegionFname:
			chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitBEDFile(intervalFname=self.selectedRegionFname, \
														noOfLinesPerUnit=self.maxNoOfRegionsPerJob, \
														folderName=self.pegasusFolderName, parentJobLs= None)
		else:
			chr2size = self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, contigMinRankBySize=self.contigMinRankBySize)
			#chr2size = set(['Contig149'])	#temporary when testing Contig149
			#chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
			chrLs = chr2size.keys()
			chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=chr2size, \
														intervalSize=self.intervalSize, \
														intervalOverlapSize=self.intervalOverlapSize)
		# 2012.8.2 if maxContigID/minContigID is not well defined. restrictContigDictionry won't do anything.
		chr2IntervalDataLs = self.restrictContigDictionry(dc=chr2IntervalDataLs, \
												maxContigID=self.maxContigID, minContigID=self.minContigID)
		"""
		#2011-9-2
		self.outputSeqCoverage(self.seqCoverageFname)
		seqCoverageF = File(os.path.basename(self.seqCoverageFname))
		seqCoverageF.addPFN(PFN("file://" + os.path.abspath(self.seqCoverageFname), \
											self.input_site_handler))
		workflow.addFile(seqCoverageF)
		"""
		
		if self.run_type==1:	#multi-sample calling
			dirPrefix = ""
			# Add a mkdir job for the call directory.
			callOutputDir = "%scall"%(dirPrefix)
			callOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=callOutputDir)
			gatkDir = "%sgatk"%(dirPrefix)
			gatkDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=gatkDir)
			samtoolsDir = "%ssamtools"%(dirPrefix)
			samtoolsDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=samtoolsDir)
			unionDir = "%sgatk_samtools_union"%(dirPrefix)
			unionDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=unionDir)
			intersectionDir = "%sgatk_samtools_intersection"%(dirPrefix)
			intersectionDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=intersectionDir)
			
			
			self.addGenotypeCallJobs(workflow, alignmentDataLs, chr2IntervalDataLs=chr2IntervalDataLs, samtools=workflow.samtools, \
					genotyperJava=workflow.genotyperJava, genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
					createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
					mergeSamFilesJar=workflow.mergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools,\
					bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
					concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
					GenotypeCallByCoverage=workflow.GenotypeCallByCoverage, refFastaFList=refFastaFList, bamListF=None, \
					callOutputDirJob =callOutputDirJob, gatkDirJob=gatkDirJob, samtoolsDirJob=samtoolsDirJob, unionDirJob=unionDirJob, intersectionDirJob=intersectionDirJob,\
					namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
					seqCoverageF=None, \
					needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
					intervalSize=self.intervalSize, site_type=self.site_type, dataDir=self.dataDir, addGATKJobs=self.addGATKJobs,\
					cumulativeMedianDepth=cumulativeMedianDepth)
		elif self.run_type==2:
			#2011-11-4 for single-alignment-calling pipeline, adjust the folder name so that they are unique from each other
			for alignmentData in alignmentDataLs:
				alignment = alignmentData.alignment
				dirPrefix = "%s"%(alignment.id)
				# Add a mkdir job for the call directory.
				callOutputDir = "%s_call"%(dirPrefix)
				callOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=callOutputDir)
				gatkDir = "%s_gatk"%(dirPrefix)
				gatkDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=gatkDir)
				samtoolsDir = "%s_samtools"%(dirPrefix)
				samtoolsDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=samtoolsDir)
				unionDir = "%s_gatk_samtools_union"%(dirPrefix)
				unionDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=unionDir)
				intersectionDir = "%s_gatk_samtools_intersection"%(dirPrefix)
				intersectionDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=intersectionDir)
				
				self.addGenotypeCallJobs(workflow, [alignment], chr2IntervalDataLs, samtools=workflow.samtools, \
						genotyperJava=workflow.genotyperJava, genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
						createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
						mergeSamFilesJar=workflow.mergeSamFilesJar, \
						BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
						mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools,\
						bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
						concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
						GenotypeCallByCoverage=workflow.GenotypeCallByCoverage, refFastaFList=refFastaFList, bamListF=None, \
						callOutputDirJob =callOutputDirJob, gatkDirJob=gatkDirJob, samtoolsDirJob=samtoolsDirJob, unionDirJob=unionDirJob, intersectionDirJob=intersectionDirJob,\
						namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						seqCoverageF=None, \
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						intervalSize=self.intervalSize, site_type=self.site_type, dataDir=self.dataDir)
				
		
		"""
		addSelectAndSplitBamReturnData = self.addSelectAndSplitBamJobs(db_vervet, workflow, alignmentLs, site_handler, \
							self.maxContigID, chrLs, samtools=samtools, \
							java=java, addOrReplaceReadGroupsAndCleanSQHeaderJar=workflow.addOrReplaceReadGroupsAndCleanSQHeaderJar, \
							BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar,\
							mkdir=workflow.mkdirWrap, namespace=workflow.namespace, \
							version=workflow.version, mkCallDirJob=callOutputDirJob,\
							addRGExecutable=workflow.addRGToBAM, dataDir=self.dataDir)
		chr2jobDataLs = addSelectAndSplitBamReturnData.chr2jobDataLs
		
		returnData3 = self.addRefFastaFileSplitJobs(workflow, refFastaF, selectAndSplitFasta, chrLs, mkdir=workflow.mkdirWrap, samtools=workflow.samtools,
								java=workflow.java, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar,\
								site_handler=workflow.site_handler, namespace=workflow.namespace, version=workflow.version)
		chr2splitFastaJobDataLs = returnData3.chr2jobDataLs
		
		returnData2 = self.addMergeAlignmentAndGenotypeCallJobs(workflow, chr2jobDataLs, chrLs, samtools, \
							java, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, GenotypeCallByCoverage=workflow.GenotypeCallByCoverage, \
							refFastaF=refFastaF, \
							namespace=workflow.namespace, version=workflow.version, callOutputDir = callOutputDir, \
							genotypeCallerType=self.genotypeCallerType, \
							mergeSamFilesJar=workflow.mergeSamFilesJar, genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, calcula=workflow.calcula, \
							chr2splitFastaJobDataLs=chr2splitFastaJobDataLs, seqCoverageF=seqCoverageF, \
							needFastaIndexJob=self.needFastaIndexJob, \
							needFastaDictJob=self.needFastaDictJob, site_type=self.site_type)
		
		chr2mergedBamCallJob =returnData2.chr2mergedBamCallJob
		
		#merge all genotype call files
		mergeGenotypeMatrix_job = Job(namespace=workflow.namespace, name=workflow.mergeGenotypeMatrix.name, version=workflow.version)
		finalCallOutputFname = '%s_genomes_vs_top%sReferences_call.tsv'%(len(alignmentLs), len(chrLs))
		finalCallOutput = File(finalCallOutputFname)
		mergeGenotypeMatrix_job.addArguments("-o", finalCallOutput)
		mergeGenotypeMatrix_job.uses(finalCallOutput, transfer=True, link=Link.OUTPUT, register=True)
		workflow.addJob(mergeGenotypeMatrix_job)
		inputFnameLs = []
		for chr, callJobData in chr2mergedBamCallJob.iteritems():
			callJobOutput, callJob = callJobData[:2]
			#callJobOutput = callJob.used_files[1]
			
			mergeGenotypeMatrix_job.uses(callJobOutput, transfer=False, register=True, link=Link.INPUT)
			inputFnameLs.append(callJobOutput.name)
			workflow.depends(parent=callJob, child=mergeGenotypeMatrix_job)
		
			mergeGenotypeMatrix_job.addArguments(callJobOutput)
		"""
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = AlignmentToCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
