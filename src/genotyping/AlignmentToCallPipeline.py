#!/usr/bin/env python
"""
Examples:
	#2013.10.07 call population sabaeus, only at certain intervals
	tax_id=60711 
	%s -a 3488 --country_id_ls 1,129,130,131,132,133,134,136,144,148,151,152
		--tax_id_ls $tax_id --sequence_filtered 1 --local_realigned 1 --reduce_reads 0
		--completedAlignment 1 --excludeContaminant
		-o dags/AlignmentToCall/AlignmentToCall_VervetTax$tax_id\_LocalRealigned1_vs_3488_ChrRank0_2000_AlnDepthIntervalMethod8_Min1Kb_FreeBayes.xml
		-u yh --hostname localhost --alignment_method_id 6 --data_dir ~/NetworkData/vervet/db/
		--local_data_dir ~/NetworkData/vervet/db -J ~/bin/jdk/bin/java
		--noOfCallingJobsPerNode 3 --clusters_size 40 --contigMinRankBySize 0
		--contigMaxRankBySize 2000 --intervalSize 1000000 --intervalOverlapSize 0
		--ref_genome_tax_id 60711 --ref_genome_sequence_type_id 1
		--ref_genome_version 1 --genotypeCallerType 3 --heterozygosityPrior 0.01
		--GATKGenotypeCallerType HaplotypeCaller --mask_genotype_method_id 227
		--minAlignmentDepthIntervalLength 1000 --alignmentDepthIntervalMethodShortName PopSabaeusRef3488MinLength100
		--alignmentDepthMinFold 0.01 --alignmentDepthMaxFold 2 -l hcondor -j hcondor
	
	
	#2011-9-14 top 25 contigs, variants only, run on uschpc cluster
	%s --ind_aln_id_ls 559-656 --input_site_handler uschpc --site_handler uschpc -u yh -a 524 -s 2 -e /home/cmb-03/mn/yuhuang
		--data_dir /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/
		--hostname 10.8.0.10 -o ./AlignmentToCallPipeline_559_656_vs_524_top_25Contigs_uschpc.xml
		--local_data_dir ~/mnt/hpc-cmb_home/NetworkData/vervet/db/ -N25
	
	# 2011-11-4 run GATK/samtools on single-sample at a time, for 4 high-coverage VRC monkeys, top 804 contigs
	%s -a 524 --ind_seq_id_ls 15-18 -u yh --site_handler condorpool --input_site_handler condorpool -s 2 -N 804
		-o AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_single_sample_condor.xml --hostname uclaOffice -n 2
		--noOfCallingJobsPerNode 5
	
	# 2012.6.26 run on hoffman2 condor pool (hcondor), filtered sequences (-Q 1), alignment method 2 (-G 2)
	# no site ID filtering (-S ""), clustering size =5 for calling jobs (--noOfCallingJobsPerNode 5).
	# with 2million bp interval (--intervalSize 2000000).
	%s -a 524 --ind_seq_id_ls 633,1495,...,1524,1459,1505,1478,1486,1442,1472,1516,1453
		-u yh --hostname localhost -N 7559 -S "" --sequence_filtered 1 -G 2 --site_handler hcondor --input_site_handler hcondor -e ~ 
		--data_dir ~/NetworkData/vervet/db --local_data_dir ~/NetworkData/vervet/db 
		-J ~/bin/jdk/bin/java --noOfCallingJobsPerNode 5 --intervalSize 2000000 --intervalOverlapSize 0
		-o dags/AlignmentToCall_130PoplationVervets_vs_524_top7559Contigs.xml
	
	# 2012.7.30 genotype-call 723 alignments on method 7 sites (--selectedRegionFname ...). "-N ..." (top number of contigs) doesn't matter here.
	# 2000 method 7 sites for each calling job (--maxNoOfRegionsPerJob 2000)
	%s -a 524 -S 447 -u yh --hostname localhost --sequence_filtered 1 --alignment_method_id  2
		--site_handler hcondor --input_site_handler hcondor
		--data_dir ~/NetworkData/vervet/db --local_data_dir ~/NetworkData/vervet/db
		-J ~/bin/jdk/bin/java
		--noOfCallingJobsPerNode 5  -o dags/AlignmentToCall_AllVRC_vs_524_top1000Contigs.xml
		--sourceVCFFolder ~/NetworkData/vervet/db/genotype/method_7/ --selectedRegionFname aux/method7_BED.tsv --maxNoOfRegionsPerJob 2000
	
	#2012.8.2 testing calling at known sites (--selectedRegionFname method7_BED.n10k.tsv, Contig731 and partial Contig645) 
	#			with 500 regions for each job (--maxNoOfRegionsPerJob 500),
	# run GATK (--genotypeCallerType 1), no clustering for any job (--noOfCallingJobsPerNode 1 --clusters_size1)
	# only on four alignments of isq-id (--ind_seq_id_ls 643-646)
	%s -a 524 --ind_seq_id_ls 643-646 -S 447 -u yh --hostname localhost  --sequence_filtered 1 --alignment_method_id  2 
		--site_handler hcondor --input_site_handler hcondor
		--data_dir ~/NetworkData/vervet/db --local_data_dir ~/NetworkData/vervet/db
		-J ~/bin/jdk/bin/java
		--noOfCallingJobsPerNode 1 --clusters_size 1
		-o dags/AlignmentToCall_ISQ643_646_vs_524_method7n10kSites.xml
		--genotypeCallerType 1
		--sourceVCFFolder ~/NetworkData/vervet/db/genotype/method_7/ --selectedRegionFname aux/method7_BED.n10k.tsv --maxNoOfRegionsPerJob 5000
	
	# 2012.8.2 part of the test above, now run multi-sample on Contig731 on all sites
	# use --maxContigID 731 (maxContigID) -V 731 (minContigID) to restrict the top 1000 contigs to only Contig731.
	# run on intervals of 200kb (--intervalSize), with GATK (--genotypeCallerType 1)
	# add --individual_sequence_file_raw_id_type 2 (library-specific alignments, different libraries of one individual_sequence) 
	# add --individual_sequence_file_raw_id_type 3 (both all-library-fused and library-specific alignments)
	# add "--country_id_ls 135,136,144,148,151" to limit individuals from US,Barbados,StKitts,Nevis,Gambia (AND with -S, )
	%s -a 524 --ind_seq_id_ls 643-646
		#-S 447
		-u yh --hostname localhost --sequence_filtered 1 --alignment_method_id  2
		--site_handler hcondor --input_site_handler hcondor
		-e ~ --data_dir ~/NetworkData/vervet/db
		--local_data_dir ~/NetworkData/vervet/db
		-J ~/bin/jdk/bin/java --noOfCallingJobsPerNode 1 --clusters_size 1
		-o dags/AlignmentToCall/AlignmentToCall_ISQ643_646_vs_524_Contig731.xml
		--genotypeCallerType 1
		--contigMaxRankBySize 1000 --maxContigID 731 --minContigID 731 --intervalSize 200000 --intervalOverlapSize 0
		#--individual_sequence_file_raw_id_type 2 --country_id_ls 135,136,144,148,151 --tax_id_ls 60711 #sabaeus
		#--heterozygosityPrior 0.01
		# 2013.2.25 add this to try HaplotypeCaller
		#--minPruningForGATKHaplotypCaller 2 --GATKGenotypeCallerType HaplotypeCaller
	
Description:
	2013.06.03 "--selectedRegionFname" is required for GATK & SAMtools genotyping mode, not for Platypus.
		"--sourceVCFFolder" is used to fill "--selectedRegionFname". If the latter file already exists,
			and you do not want to overwrite it, do not supply "--sourceVCFFolder".
	
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
				sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import Executable, File, PFN, Job, Link
from pymodule import ProcessOptions, PassingData
from pymodule.pegasus import yh_pegasus
from vervet.src import AbstractVervetWorkflow, VervetDB, AbstractVervetAlignmentAndVCFWorkflow

parentClass = AbstractVervetAlignmentAndVCFWorkflow
class AlignmentToCallPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('inputDir', 0, ))
	
	
	commonCallPipelineOptionDict= {
						#('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("site_type", 1, int): [2, 's', 1, '1: all genome sites, 2: variants only'],\
						("noOfCallingJobsPerNode", 1, int): [1, 'O', 1, 'this parameter controls how many genotype calling jobs should be clustered together to run on one node. \
									Increase it to above 1 only when your average genotyping job is short and the number of input bam files are short.'],\
						("polymuttPath", 1, ): ["%s/bin/polymutt", '', 1, 'path to the polymutt binary'],\
						("trioCallerPath", 1, ): ["%s/script/vervet/bin/trioCaller/TrioCaller", '', 1, 'path to TrioCaller binary'],\
						("GATKGenotypeCallerType", 1, ): ['UnifiedGenotyper', '', 1, 'passed to -T of GATK. alternative is HaplotypeCaller'],\
						("heterozygosityPrior", 1, float): [0.005, '', 1, 'heterozygosity prior for GATK. GATK default is 0.001 '],\
						("minPruningForGATKHaplotypCaller", 1, int): [2, '', 1, 'The minimum allowed pruning factor in assembly graph.\n\
	Paths with <= X supporting kmers are pruned from the graph.\n\
	default in GATK2 is 1.'],\
						("genotypeCallerType", 0, int): [0, '', 1, '0: SAMtools, 1: GATK (--GATKGenotypeCallerType ...), 2: Platypus, 3: FreeBayes'],\
						}
	
	commonCallPipelineOptionDict.update(parentClass.partitionWorkflowOptionDict.copy())
	commonCallPipelineOptionDict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(commonCallPipelineOptionDict.copy())
	option_default_dict.update({
						("run_type", 1, int): [1, 'n', 1, '1: multi-sample calling, 2: single-sample one by one'],\
						("sourceVCFFolder", 0, ): [None, '', 1, 'a VCF folder which provides source vcf files for Platypus to run in genotyping mode (no snp discovery)'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
						#("genotypeCallerType", 1, int): [1, 'y', 1, '1: GATK + coverage filter; 2: ad-hoc coverage based caller; 3: samtools + coverage filter'],\

	def __init__(self,  **keywords):
		"""
		2012.9.17 call parentClass.__init__() directly
		2011-7-11
		"""
		self.pathToInsertHomePathList.extend(['polymuttPath', 'trioCallerPath'])
		
		parentClass.__init__(self, **keywords)
	
	
	def addRefFastaFileSplitJobs(self, workflow, refFastaF, selectAndSplitFasta, chrLs, mkdir=None, samtools=None,
								java=None, CreateSequenceDictionaryJar=None,\
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
		for chromosome in chrLs:
			if chromosome not in chr2jobDataLs:
				chr2jobDataLs[chromosome] = []
			selectAndSplitFastaJob.addArguments(chromosome)
			
			fastaFname = os.path.join(fastaOutputDir, '%s.fasta'%(chromosome))
			fastaFile = File(fastaFname)
			
			# add the index job
			fai_index_job = Job(namespace=namespace, name=samtools.name, version=version)
			fai_index_job.addArguments("faidx", fastaFname)
			fastaFAIIndexFname = os.path.join(fastaOutputDir, '%s.fasta.fai'%(chromosome))
			fastaFAIIndexFile = File(fastaFAIIndexFname)
			fai_index_job.uses(fastaFAIIndexFile, transfer=True, register=False, link=Link.OUTPUT)	#this file is input & output
			#fai_index_job.uses(output, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
			#fai_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
			workflow.addJob(fai_index_job)
			workflow.depends(parent=selectAndSplitFastaJob, child=fai_index_job)
			
			# the add-read-group job
			createSeqDictJob = Job(namespace=namespace, name=java.name, version=version)
			dictFname = os.path.join(fastaOutputDir, '%s.dict'%(chromosome))
			#outputRGFname = '%s_%s.RG.bam'%(inputFileBaseNamePrefix, chromosome)
			dictFile = File(dictFname)
			createSeqDictJob.addArguments('-jar', CreateSequenceDictionaryJar, \
								"R=", fastaFile, 'O=', dictFile)
			self.addJobUse(createSeqDictJob, file=CreateSequenceDictionaryJar, transfer=True, register=True, link=Link.INPUT)
			createSeqDictJob.uses(dictFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			createSeqDictJob.uses(fastaFile, transfer=True, register=True, link=Link.OUTPUT)	#time to discard them
			workflow.addJob(createSeqDictJob)
			workflow.depends(parent=selectAndSplitFastaJob, child=createSeqDictJob)
			
			chr2jobDataLs[chromosome] = [fai_index_job, fastaFAIIndexFile, createSeqDictJob, fastaFile, dictFile]
			self.no_of_jobs += 1
		return PassingData(chr2jobDataLs=chr2jobDataLs, workflow=workflow)
	
	def addSelectAndSplitBamJobs(self, db_vervet, workflow, alignmentLs, site_handler, maxContigID, chrLs, samtools=None, \
							java=None, AddOrReplaceReadGroupsAndCleanSQHeaderJar=None, BuildBamIndexJar=None, mkdir=None, namespace="workflow",\
							version="1.0", mkCallDirJob=None,\
							addRGExecutable=None, data_dir=None):
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
			inputFname = os.path.join(data_dir, alignment.path)
			inputFile = File(inputFname)
			inputFile.addPFN(PFN("file://" + inputFname, site_handler))
			workflow.addFile(inputFile)
			
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
			#mkdir.uses(inputFile, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
			workflow.addJob(mkdirJob)
			workflow.depends(parent=mkCallDirJob, child=mkdirJob)
			
			for chromosome in chrLs:
				if chromosome not in chr2jobDataLs:
					chr2jobDataLs[chromosome] = []
				#select reads that are aligned to one reference name
				selectRefJob = Job(namespace=namespace, name=samtools.name, version=version)
				
				outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix, chromosome))
				#outputFname = '%s_%s.bam'%(inputFileBaseNamePrefix, chromosome)
				output = File(outputFname)
				selectRefJob.addArguments('view', '-h', inputFile, chromosome, "-o", output, "-b", "-u")	# -b -u forces uncompressed bam output
				#selectRefJob.uses(inputFile, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
				workflow.addJob(selectRefJob)
				workflow.depends(parent=mkdirJob, child=selectRefJob)
				
				# add the index job
				index_1_job = Job(namespace=namespace, name=samtools.name, version=version)
				index_1_job.addArguments("index", output)
				#index_1_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is inputFile & output
				index_1_job.uses(output, transfer=False, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputFname)
				index_1_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(index_1_job)
				workflow.depends(parent=selectRefJob, child=index_1_job)
				
				# the add-read-group job
				#addRGJob = Job(namespace=namespace, name=addRGExecutable.name, version=version)
				addRGJob = Job(namespace=namespace, name=java.name, version=version)
				#2011-7-27 somehow AddOrReplaceReadGroupsAndCleanSQHeader.jar couldn't output a un-corrupted bam file. so sam first.
				outputRGSAMFname = os.path.join(outputDir, '%s_%s.RG.sam'%(inputFileBaseNamePrefix, chromosome))
				outputRGSAM = File(outputRGSAMFname)
				"""
				tmpRGFname = os.path.join(outputDir, '%s_%s.RG.txt'%(inputFileBaseNamePrefix, chromosome))
				addRGJob.addArguments(read_group, platform_id, output, tmpRGFname, outputRGSAM)
				"""
				addRGJob.addArguments('-jar', AddOrReplaceReadGroupsAndCleanSQHeaderJar, \
									"INPUT=", output,\
									'RGID=%s'%(read_group), 'RGLB=%s'%(platform_id), 'RGPL=%s'%(platform_id), \
									'RGPU=%s'%(read_group), 'RGSM=%s'%(read_group),\
									'OUTPUT=', outputRGSAM, 'SQName=%s'%(chromosome))	#'SORT_ORDER=coordinate', (adding this is useless)
				self.addJobUse(addRGJob, file=AddOrReplaceReadGroupsAndCleanSQHeaderJar, transfer=True, register=True, link=Link.INPUT)
				
				addRGJob.uses(output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#time to discard them
				addRGJob.uses(outputRGSAM, transfer=False, register=True, link=Link.OUTPUT)	#time to discard them
				workflow.addJob(addRGJob)
				workflow.depends(parent=index_1_job, child=addRGJob)
				#output.addPFN(PFN("file://" + outputFname, site_handler))
				#selectAndSplitJob.uses(output, transfer=True, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
				
				samToBamJob = Job(namespace=namespace, name=samtools.name, version=version)
				outputRGFname = os.path.join(outputDir, '%s_%s.RG.bam'%(inputFileBaseNamePrefix, chromosome))
				outputRG = File(outputRGFname)
				samToBamJob.addArguments('view', '-F4', '-Sbh', "-o", outputRG, "-u", outputRGSAM)	# -b -u forces uncompressed bam output
				samToBamJob.uses(outputRGSAM, transfer=False, register=True, link=Link.INPUT)
				#samToBamJob.uses(outputRG, transfer=False, register=True, link=Link.OUTPUT)	#don't register it here
				workflow.addJob(samToBamJob)
				workflow.depends(parent=addRGJob, child=samToBamJob)
				
				# add the index job
				index_sam_job = Job(namespace=namespace, name=java.name, version=version)
				bai_output = File('%s.bai'%outputRGFname)
				index_sam_job.addArguments("-Xms128m", "-Xmx2500m", "-jar", BuildBamIndexJar, "VALIDATION_STRINGENCY=LENIENT", \
								"INPUT=", outputRG, \
								"OUTPUT=", bai_output)
				index_sam_job.uses(outputRG, transfer=False, register=False, link=Link.INPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
				index_sam_job.uses(bai_output, transfer=False, register=False, link=Link.OUTPUT)
				
				"""
				samtools_index_job = Job(namespace=namespace, name=samtools.name, version=version)
				samtools_index_job.addArguments("index", outputRG)
				#samtools_index_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is inputFile & output
				samtools_index_job.uses(outputRG, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputRGFname)
				samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				"""
				workflow.addJob(index_sam_job)
				workflow.depends(parent=samToBamJob, child=index_sam_job)
				chr2jobDataLs[chromosome].append((outputRG, index_sam_job, bai_output))
				
				"""
				#2011-9-1 temporary addition to make sure it's sorted for Vasily
				# input bam doesn't need to be indexed.
				outputRGSortSAMFnamePrefix = os.path.join(outputDir, '%s_%s.RG.sorted'%(inputFileBaseNamePrefix, chromosome))
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
	
	
	def addCombineChromosomeVCFJob(self, outputDirPrefix="", callerName="", chromosome=None, refFastaFList=None,\
								job_max_memory=5000, transferOutput=True):
		"""
		2013.09.04
		"""
		callOutputDirJobAttrName = "%s%sCallOutputDirJob"%(outputDirPrefix, callerName)
		snpOutputDirJobAttrName = "%s%sSNPOutputDirJob"%(outputDirPrefix, callerName)
		indelOutputDirJobAttrName = "%s%sIndelOutputDirJob"%(outputDirPrefix, callerName)
		
		callOutputDirJob = getattr(self, callOutputDirJobAttrName, None)
		if callOutputDirJob is None:
			callOutputDirJob = self.addMkDirJob(outputDir="%s%s"%(outputDirPrefix, callerName))
			setattr(self, callOutputDirJobAttrName, callOutputDirJob)
		
		snpOutputDirJob = getattr(self, snpOutputDirJobAttrName, None)
		if snpOutputDirJob is None:
			snpOutputDirJob = self.addMkDirJob(outputDir="%s%sSNP"%(outputDirPrefix, callerName))
			setattr(self, snpOutputDirJobAttrName, snpOutputDirJob)
		
		indelOutputDirJob = getattr(self, indelOutputDirJobAttrName, None)
		if indelOutputDirJob is None:
			indelOutputDirJob = self.addMkDirJob(outputDir="%s%sIndel"%(outputDirPrefix, callerName))
			setattr(self, indelOutputDirJobAttrName, indelOutputDirJob)
		
		#2011-9-22 union of all GATK intervals for one contig
		combineCallOutputFname = os.path.join(callOutputDirJob.folder, '%s.vcf'%chromosome)
		combineCallOutputF = File(combineCallOutputFname)
		combineCallJob = self.addGATKCombineVariantsJob(refFastaFList=refFastaFList, outputFile=combineCallOutputF, \
													genotypeMergeOptions='UNSORTED', \
			parentJobLs=[callOutputDirJob], \
			extraArguments=None, extraArgumentList=None, extraDependentInputLs=None,\
			transferOutput=False, job_max_memory=job_max_memory,walltime=180)
		
		gatkGzipUnionOutputF = File("%s.gz"%combineCallOutputFname)
		bgzip_tabix_combineCallOutputF_job = self.addBGZIP_tabix_Job(\
				parentJob=combineCallJob, inputF=combineCallJob.output, outputF=gatkGzipUnionOutputF, \
				transferOutput=transferOutput)
		
		indelUnionOutputF = File(os.path.join(indelOutputDirJob.folder, '%s.indel.vcf'%chromosome))
		selectIndelJob = self.addGATKJob(executable=self.SelectVariantsJava, GATKAnalysisType="SelectVariants",\
			inputFile=combineCallJob.output, inputArgumentOption="--variant", refFastaFList=refFastaFList, inputFileList=None,\
			argumentForEachFileInInputFileList=None,\
			interval=None, outputFile=indelUnionOutputF, \
			parentJobLs=[indelOutputDirJob, combineCallJob], transferOutput=False, job_max_memory=job_max_memory,\
			frontArgumentList=None, extraArguments="-selectType INDEL -selectType MNP", \
			extraArgumentList=None, extraOutputLs=None, \
			extraDependentInputLs=None, no_of_cpus=None, walltime=80)
		
		bgzip_tabix_job = self.addBGZIP_tabix_Job(parentJob=selectIndelJob, \
				inputF=selectIndelJob.output, outputF=File("%s.gz"%indelUnionOutputF.name), \
				transferOutput=transferOutput)
		
		SNPOnlyOutputF = File(os.path.join(snpOutputDirJob.folder, '%s.SNP.vcf'%chromosome))
		selectSNPJob = self.addGATKJob(executable=self.SelectVariantsJava, GATKAnalysisType="SelectVariants",\
			inputFile=combineCallJob.output, inputArgumentOption="--variant", refFastaFList=refFastaFList, inputFileList=None,\
			argumentForEachFileInInputFileList=None,\
			interval=None, outputFile=SNPOnlyOutputF, \
			parentJobLs=[snpOutputDirJob, combineCallJob], transferOutput=False, job_max_memory=job_max_memory,\
			frontArgumentList=None, extraArguments="-selectType SNP", \
			extraArgumentList=None, extraOutputLs=None, \
			extraDependentInputLs=None, no_of_cpus=None, walltime=80)
		bgzip_tabix_job = self.addBGZIP_tabix_Job(parentJob=selectSNPJob, \
				inputF=selectSNPJob.output, outputF=File("%s.gz"%SNPOnlyOutputF.name), \
				transferOutput=transferOutput)
		return PassingData(callUnionJob=combineCallJob, callOutputDirJob=callOutputDirJob, selectIndelJob=selectIndelJob, selectSNPJob=selectSNPJob)
	
	def addGenotypeCallJobs(self, workflow=None, alignmentDataLs=None, chr2IntervalDataLs=None, \
				registerReferenceData=None, \
				site_handler=None, input_site_handler=None,\
				needFastaIndexJob=False, needFastaDictJob=False, \
				site_type=1, data_dir=None, no_of_gatk_threads = 1,\
				outputDirPrefix="", genotypeCallerType=0, cumulativeMedianDepth=5000, \
				sourceVCFFolder=None, extraArgumentPlatypus="", **keywords):
		"""
		2013.06.03 added argument extraArgumentPlatypus
		2013.05.20 added argument sourceVCFFolder
		2013.05.16
			argument addGATKJobs renamed to genotypeCallerType
				0: SAMtools, 1: GATK (--GATKGenotypeCallerType ...), 2: Platypus
		2012.7.31
			add argument addGATKJobs (default=False)
		2011-9-22
			add argument concatGATK, concatSamtools.
		2011-9-15
			bamListF is now useless. samtools_job could accept variable-length list of bam input files
		2011-9-14
			argument intervalSize determines how many sites gatk/samtools works on at a time
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding genotype call jobs for %s chromosomes/contigs ..."%(len(chr2IntervalDataLs)))
		
		job_max_memory = 5000	#in MB
		vcf_job_max_memory = 1000
		refFastaFList = registerReferenceData.refFastaFList
		refFastaF = refFastaFList[0]
		
		if needFastaDictJob or registerReferenceData.needPicardFastaDictJob:	# the .dict, .fai file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, CreateSequenceDictionaryJava=self.CreateSequenceDictionaryJava, \
												CreateSequenceDictionaryJar=self.CreateSequenceDictionaryJar, refFastaF=refFastaF)
			refFastaDictF = fastaDictJob.refFastaDictF
		else:
			fastaDictJob = None
			refFastaDictF = registerReferenceData.refPicardFastaDictF
		
		if needFastaIndexJob or registerReferenceData.needSAMtoolsFastaIndexJob:
			fastaIndexJob = self.addRefFastaFaiIndexJob(workflow, samtools=self.samtools, refFastaF=refFastaF)
			refFastaIndexF = fastaIndexJob.refFastaIndexF
		else:
			fastaIndexJob = None
			refFastaIndexF = registerReferenceData.refSAMtoolsFastaIndexF
		
		# Add a mkdir job for the call directory.
		#callOutputDirJob = self.addMkDirJob(outputDir="%sCall"%(outputDirPrefix))
		samtoolsDirJob = self.addMkDirJob(outputDir="%sSAMtools"%(outputDirPrefix))
		samtoolsIndelDirJob = self.addMkDirJob(outputDir="%sSAMtoolsIndel"%(outputDirPrefix))
		#unionDirJob = self.addMkDirJob(outputDir="%s%s_SAMtools_union"%(outputDirPrefix, callerName))
		#intersectionDirJob = self.addMkDirJob(outputDir="%s%s_SAMtools_intersection"%(outputDirPrefix, callerName))
		
		alignmentDataLs = self.addAddRG2BamJobsAsNeeded(alignmentDataLs=alignmentDataLs, site_handler=site_handler, \
					input_site_handler=input_site_handler, \
					addOrReplaceReadGroupsJava=self.addOrReplaceReadGroupsJava, AddOrReplaceReadGroupsJar=self.AddOrReplaceReadGroupsJar, \
					BuildBamIndexFilesJava=self.BuildBamIndexFilesJava, BuildBamIndexJar=self.BuildBamIndexJar, \
					mv=self.mv, data_dir=data_dir)
		#alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		
		if sourceVCFFolder and genotypeCallerType==2:	#2013.05.29 only for Platypus
			sourceVCFDirJob = self.addMkDirJob(outputDir="%sSourceVCF"%(outputDirPrefix))
			sourceVCFFileID2object = self.getVCFFileID2object(self.sourceVCFFolder)
		else:
			sourceVCFDirJob = None
			sourceVCFFileID2object = None
		
		# add merge jobs for every reference
		returnData = PassingData(jobDataLs=[])
		chromosome2jobData = {}	#2013.10.22 
		no_of_chromosomes = 0
		no_of_intervals = 0
		for chromosome, intervalDataLs in chr2IntervalDataLs.iteritems():
			if sourceVCFFileID2object and chromosome not in sourceVCFFileID2object:	#2013.05.20 genotyping mode and this chromosome is not in the source folder
				continue
			#reduce the number of chunks 1 below needed. last trunk to reach the end of contig
			#however set it to 1 for contigs smaller than intervalSize 	
			"""
			callOutputFname = os.path.join(callOutputDirJob.folder, '%s.vcf.gz'%chromosome)
			callOutputF = File(callOutputFname)
			#2012.7.21 no more gatk-samtools intersection job
			wholeRefUnionOfIntersectionJob = self.addVCFConcatJob(workflow, concatExecutable=vcf_concat, \
							parentDirJob=callOutputDirJob, \
							outputF=callOutputF, namespace=namespace, version=version, transferOutput=True, \
							vcf_job_max_memory=vcf_job_max_memory)
			"""
			#self.addRefFastaJobDependency(workflow, wholeRefUnionOfIntersectionJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
			#							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
			if chromosome not in chromosome2jobData:
				chromosome2jobData[chromosome] = PassingData(genotypingJobLs=[], refineGenotypeJobLs=[])
			if genotypeCallerType>=1:
				#2011-9-22 union of all GATK intervals for one contig
				combineChromosomeJobData = self.addCombineChromosomeVCFJob(outputDirPrefix=outputDirPrefix, callerName="CallerType%s"%(genotypeCallerType), \
								chromosome=chromosome, refFastaFList=refFastaFList,\
								job_max_memory=5000, transferOutput=True)
			
			elif genotypeCallerType==0:
				#2011-9-22 union job for SNPs
				samtoolsUnionOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf.gz'%chromosome)
				samtoolsUnionOutputF = File(samtoolsUnionOutputFname)
				samtoolsUnionJob = self.addVCFConcatJob(concatExecutable=self.concatSamtools, parentDirJob=samtoolsDirJob, \
								outputF=samtoolsUnionOutputF, transferOutput=True, \
								vcf_job_max_memory=vcf_job_max_memory)
				
				#2011-9-22 union job for indels
				samtoolsIndelUnionOutputFname = os.path.join(samtoolsIndelDirJob.folder, '%s.indel.vcf.gz'%chromosome)
				samtoolsIndelUnionOutputF = File(samtoolsIndelUnionOutputFname)
				samtoolsIndelUnionJob = self.addVCFConcatJob(concatExecutable=self.concatSamtools, parentDirJob=samtoolsIndelDirJob, \
								outputF=samtoolsIndelUnionOutputF, transferOutput=True, \
								vcf_job_max_memory=vcf_job_max_memory)
			if hasattr(self, "addRefineSNPGenotypeJob"):
				combineChromosomeRefineSNPJobData = self.addCombineChromosomeVCFJob(outputDirPrefix=outputDirPrefix, callerName="Refine", \
								chromosome=chromosome, refFastaFList=refFastaFList,\
								job_max_memory=5000, transferOutput=True)
			no_of_chromosomes += 1
			for intervalData in intervalDataLs:
				no_of_intervals += 1
				if intervalData.file:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.file
				else:
					mpileupInterval = intervalData.interval
					bcftoolsInterval = intervalData.interval
				vcfBaseFname = intervalData.intervalFileBasenameSignature
				
				span = intervalData.span
				readSpace = cumulativeMedianDepth * span
				#base for Platypus is 450X coverage in 4Mb region => 80 minutes
				genotypingJobWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=readSpace, \
									baseInputVolume=450*2000000, baseJobPropertyValue=160, \
									minJobPropertyValue=120, maxJobPropertyValue=500).value
				#base for Platypus is 450X, => 5.2g
				if genotypeCallerType==1:	#GATK
					maxGenotypingJobMemory = 10000
				else:
					maxGenotypingJobMemory = 8000
				genotypingJobMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=cumulativeMedianDepth, \
									baseInputVolume=450, baseJobPropertyValue=5200, \
									minJobPropertyValue=4000, maxJobPropertyValue=maxGenotypingJobMemory).value
									
				if genotypeCallerType>=1:	#non-SAMtools
					#extra caller
					callOutputFname = os.path.join(combineChromosomeJobData.callOutputDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
					callOutputFile = File(callOutputFname)
					
					if genotypeCallerType==1:	#GATK
						#2013.05.09 add downsample setting for gatk job 
						genotypingJob= self.addGATKCallJob(genotyperJava=self.genotyperJava, GenomeAnalysisTKJar=self.GenomeAnalysisTK2Jar, \
							outputFile=callOutputFile, refFastaFList=refFastaFList, \
							parentJobLs=[combineChromosomeJobData.callOutputDirJob]+intervalData.jobLs, \
							extraArguments=" --downsample_to_coverage 250 --downsampling_type BY_SAMPLE", \
							extraDependentInputLs=None, transferOutput=False, \
							job_max_memory=genotypingJobMaxMemory, \
							no_of_gatk_threads=no_of_gatk_threads, \
							site_type=site_type, \
							interval=bcftoolsInterval,\
							heterozygosityPrior=self.heterozygosityPrior, \
							minPruningForGATKHaplotypCaller=self.minPruningForGATKHaplotypCaller, \
							GATKGenotypeCallerType = self.GATKGenotypeCallerType,\
							walltime=genotypingJobWalltime,\
							)
					elif genotypeCallerType==2:
						#2013.05.16 Platypus job
						sourceVCFFile = None
						sourceVCFIndexFile = None
						if sourceVCFFileID2object:
							sourceVCFFileObject = sourceVCFFileID2object.get(chromosome)
							sourceVCFFilePath = sourceVCFFileObject.vcfFilePath
							vcfIndexFilePath = sourceVCFFileObject.vcfIndexFilePath
							
							if sourceVCFFilePath:
								if sourceVCFFilePath[-3:]!='.gz':
									sys.stderr.write("Error: source VCF file for Platypus genotyping mode must be bgzipped. This one %s is not.\n"%\
													(sourceVCFFilePath))
									sys.exit(3)
								sourceVCFFile = self.registerOneInputFile(inputFname=sourceVCFFilePath, \
															folderName=sourceVCFDirJob.output)
								if vcfIndexFilePath:
									sourceVCFIndexFile = self.registerOneInputFile(inputFname=vcfIndexFilePath, \
															folderName=sourceVCFDirJob.output)
								else:
									sys.stderr.write("Error: index file (%s) for VCF file %s does not exist.\n"%\
													(vcfIndexFilePath, sourceVCFFilePath))
									sys.exit(3)
						#--bufferSize=BUFFERSIZE, Data will be buffered in regions of this size
						#--maxReads=MAXREADS   Maximium coverage in window
						bufferSize=1000	#used to be 15kb, too much for 725 monkeys
						maxReads=bufferSize/100*cumulativeMedianDepth*4	#100 is read length. maxReads is 4X median number of reads in this buffer window
						genotypingJob = self.addPlatypusCallJob(executable=self.Platypus, outputFile=callOutputFile, \
							refFastaFList=refFastaFList, \
							sourceVCFFile=sourceVCFFile, interval=mpileupInterval, skipRegionsFile=None,\
							site_type=site_type, \
							extraArguments="--bufferSize=%s --maxReads=%s --maxReadLength=4268 %s"%\
								(bufferSize, maxReads, extraArgumentPlatypus), \
							job_max_memory=genotypingJobMaxMemory, no_of_cpus=1, \
							walltime=genotypingJobWalltime, \
							parentJobLs=[combineChromosomeJobData.callOutputDirJob], extraDependentInputLs=[sourceVCFIndexFile], transferOutput=False)
					elif genotypeCallerType==3:
						#2013.10.07 freebayes
						genotypingJob = self.addFreeBayesJob(executable=self.freebayes, outputFile=callOutputFile, \
							refFastaFList=refFastaFList, \
							interval=bcftoolsInterval, heterozygosity=self.heterozygosityPrior, minCoverage=2,\
							site_type=site_type, \
							extraArguments=None, job_max_memory=genotypingJobMaxMemory, no_of_cpus=1, \
							walltime=genotypingJobWalltime, \
							parentJobLs=[combineChromosomeJobData.callOutputDirJob]+intervalData.jobLs, \
							extraDependentInputLs=None, transferOutput=False)
					chromosome2jobData[chromosome].genotypingJobLs.append(genotypingJob)
					#add this output to a GATK union job
					#self.addInputToStatMergeJob(statMergeJob=combineChromosomeJobData.callUnionJob, parentJobLs=[genotypingJob], \
					#						inputArgumentOption="--variant")
					
				else:	#SAMtools
					#samtools part
					samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.orig.vcf'%vcfBaseFname)
					samtoolsOutputF = File(samtoolsOutputFname)
					indelVCFOutputFname = "%s.indel.vcf"%(samtoolsOutputFname)
					indelVCFOutputF = File(indelVCFOutputFname)
					genotypingJob= self.addSAMtoolsCallJob(workflow, CallVariantBySamtools=self.CallVariantBySamtools, \
						samtoolsOutputF=samtoolsOutputF, indelVCFOutputF=indelVCFOutputF, \
						refFastaFList=refFastaFList, parentJobLs=[samtoolsDirJob]+intervalData.jobLs, \
						extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=genotypingJobMaxMemory, \
						site_type=site_type, maxDP=cumulativeMedianDepth*5, mpileupInterval=mpileupInterval, \
						bcftoolsInterval=bcftoolsInterval, walltime=genotypingJobWalltime)
					
					#deal with samtools's indel VCF
					"""
					# convert to VCF4. vcf-convert requires "-r refFastaF" but still fails in sanity check.
					# complaining ref sequence from this vcf mismatches from the reference sequence. so comment this out. 
					samtoolsIndelVCF4Fname = os.path.join(samtoolsDirJob.folder, '%s.indel.vcf'%vcfBaseFname)
					samtoolsIndelVCF4F = File(samtoolsIndelVCF4Fname)
					samtoolsIndelVCF4F_job = self.addVCFFormatConvertJob(workflow, vcf_convert=vcf_convert, \
								parentJob=genotypingJob, inputF=indelVCFOutputF, outputF=samtoolsIndelVCF4F, \
								namespace=namespace, version=version, transferOutput=False)
					"""
					#indels
					samtoolsIndelGzipOutputF = File("%s.gz"%indelVCFOutputFname)
					samtoolsIndelGzipOutput_tbi_F = File("%s.gz.tbi"%indelVCFOutputFname)
					samtoolsIndelGzipOutputF_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix, \
							parentJob=genotypingJob, inputF=indelVCFOutputF, outputF=samtoolsIndelGzipOutputF, \
							transferOutput=False)
					
					
					samtoolsIndelUnionJob.addArguments(samtoolsIndelGzipOutputF)
					samtoolsIndelUnionJob.uses(samtoolsIndelGzipOutputF, transfer=False, register=True, link=Link.INPUT)
					samtoolsIndelUnionJob.uses(samtoolsIndelGzipOutput_tbi_F, transfer=False, register=True, link=Link.INPUT)
					workflow.depends(parent=samtoolsIndelGzipOutputF_job, child=samtoolsIndelUnionJob)
					
					#SNPs 
					vcf4_samtoolsOutputFname = os.path.join(samtoolsDirJob.folder, '%s.vcf'%vcfBaseFname)
					vcf4_samtoolsOutputF = File(vcf4_samtoolsOutputFname)
					vcf_convert_samtoolsOutputF_job = self.addVCFFormatConvertJob(vcf_convert=self.vcf_convert, \
								parentJob=genotypingJob, inputF=samtoolsOutputF, outputF=vcf4_samtoolsOutputF, \
								transferOutput=False)
					
					samtoolsGzipOutputF = File("%s.gz"%vcf4_samtoolsOutputFname)
					samtoolsGzipOutput_tbi_F = File("%s.gz.tbi"%vcf4_samtoolsOutputFname)
					bgzip_tabix_samtoolsOutputF_job = self.addBGZIP_tabix_Job(bgzip_tabix=self.bgzip_tabix, \
							parentJob=vcf_convert_samtoolsOutputF_job, inputF=vcf4_samtoolsOutputF, outputF=samtoolsGzipOutputF, \
							transferOutput=False)
					
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
				genotypeUnionJob.addArguments(javaMemRequirement, '-jar', GenomeAnalysisTKJar, "-T", "CombineVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:gatk,VCF", callOutputFile, "-B:samtools,VCF", samtoolsOutputF, \
					"--out", unionOutputF, \
					'-U', '-S LENIENT', "-genotypeMergeOptions PRIORITIZE")	#no support for multi-threading
				genotypeUnionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeUnionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				self.addJobUse(DOCJob, file=GenomeAnalysisTKJar, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(callOutputFile, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(samtoolsOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeUnionJob.uses(unionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeUnionJob)
				workflow.depends(parent=unionDirJob, child=genotypeUnionJob)
				workflow.depends(parent=genotypingJob, child=genotypeUnionJob)
				workflow.depends(parent=genotypingJob, child=genotypeUnionJob)
				
				genotypeIntersectionJob = Job(namespace=namespace, name=java.name, version=version)
				intersectionOutputFname = os.path.join(intersectionDirJob.folder, '%s.vcf'%vcfBaseFname)
				intersectionOutputF = File(intersectionOutputFname)
				genotypeIntersectionJob.addArguments(javaMemRequirement, '-jar', GenomeAnalysisTKJar, "-T", "SelectVariants",\
					"-L", interval, "-R", refFastaF, \
					"-B:%s,VCF"%vcfBaseFname, unionOutputF, \
					"-select 'set == "Intersection"'", "--out", intersectionOutputF, \
					'-U', '-S LENIENT')	#no support for multi-threading
				genotypeIntersectionJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				genotypeIntersectionJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%))
				self.addJobUse(DOCJob, file=GenomeAnalysisTKJar, transfer=True, register=True, link=Link.INPUT)
				genotypeIntersectionJob.uses(unionOutputF, transfer=True, register=True, link=Link.INPUT)
				genotypeIntersectionJob.uses(intersectionOutputF, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(genotypeIntersectionJob)
				workflow.depends(parent=intersectionDirJob, child=genotypeIntersectionJob)
				workflow.depends(parent=genotypeUnionJob, child=genotypeIntersectionJob)
				
				#add this intersected output to the union job for the whole chromosome
				wholeRefUnionOfIntersectionJob.addArguments("-B:%s,VCF"%vcfBaseFname, intersectionOutputF)
				workflow.depends(parent=genotypeIntersectionJob, child=wholeRefUnionOfIntersectionJob)
				"""
				
				if genotypeCallerType in [0,1,3]:
					if genotypeCallerType==1:	#GATK
						jobAndInputOptionList=[(genotypingJob, "-I")]
					elif genotypeCallerType==3:	#freebayes
						jobAndInputOptionList=[(genotypingJob, "--bam")]
					else:	#SAMtools
						jobAndInputOptionList = [(genotypingJob, "")]
					for job, jobInputOption in jobAndInputOptionList:
						self.addRefFastaJobDependency(workflow, job, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
								refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
						self.addAlignmentAsInputToJobLs(workflow, alignmentDataLs, jobLs=[job], jobInputOption=jobInputOption)
				elif genotypeCallerType==2:	#Platypus
					self.addRefFastaJobDependency(workflow, genotypingJob, refFastaF=refFastaF, fastaDictJob=fastaDictJob, \
							refFastaDictF=refFastaDictF, fastaIndexJob = fastaIndexJob, refFastaIndexF = refFastaIndexF)
					self.addAlignmentAsInputToPlatypusJobLs(alignmentDataLs=alignmentDataLs, jobLs=[genotypingJob], \
														jobInputOption="--bamFiles")
				genotypingJob.intervalData = intervalData
				if hasattr(self, "addRefineSNPGenotypeJob"):
					refineGenotypeJobData = self.addRefineSNPGenotypeJob(inputFile=genotypingJob.output, vcfBaseFname=vcfBaseFname, \
								outputDirJob=combineChromosomeRefineSNPJobData.callOutputDirJob, statDirJob=None, \
								refFastaFList=refFastaFList, intervalData=intervalData,\
								baseInputVolume=450*2000000, realInputVolume=intervalData.span*cumulativeMedianDepth,\
								parentJobLs=[genotypingJob], \
								transferOutput=False, \
								no_of_cpus=None, job_max_memory=2000, walltime=180, \
								max_walltime=None)
					chromosome2jobData[chromosome].refineGenotypeJobLs.append(refineGenotypeJobData.refineGenotypeJob)
					returnData.jobDataLs.append(self.constructJobDataFromJob(refineGenotypeJobData.refineGenotypeJob))
				else:
					returnData.jobDataLs.append(self.constructJobDataFromJob(genotypingJob))
			#combine intervals hierarchically, to avoid open file number limit on cluster
			self.hierarchicalCombineIntervalJobIntoOneChromosome(unionJob=combineChromosomeJobData.callUnionJob, \
					refFastaFList=refFastaFList, fileBasenamePrefix="subUnion",\
					intervalJobLs=chromosome2jobData[chromosome].genotypingJobLs, \
					maxNoOfIntervalJobsInOneUnion=250,\
					outputDirJob=combineChromosomeJobData.callOutputDirJob,
					transferOutput=False, job_max_memory=None, walltime=120, needBGzipAndTabixJob=False)
					#transferOutput = False and needBGzipAndTabixJob=False because both are for intermediate files.
					# unionJob will transfer and bgzip final chromosomal level output
			if hasattr(self, "addRefineSNPGenotypeJob"):
				self.hierarchicalCombineIntervalJobIntoOneChromosome(unionJob=combineChromosomeRefineSNPJobData.callUnionJob, \
					refFastaFList=refFastaFList, fileBasenamePrefix="subUnionRefineSNP",\
					intervalJobLs=chromosome2jobData[chromosome].refineGenotypeJobLs, \
					maxNoOfIntervalJobsInOneUnion=250,\
					outputDirJob=combineChromosomeRefineSNPJobData.callOutputDirJob,
					transferOutput=False, job_max_memory=None, walltime=120, needBGzipAndTabixJob=False)
		sys.stderr.write(" %s jobs covering %s chromosomes, %s intervals.\n"%\
						(self.no_of_jobs, no_of_chromosomes, no_of_intervals))
		return returnData
	
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
	
	def outputListOfBam(self, alignmentLs, data_dir=None, outputFname=None, workflow=None, site_handler=None):
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
			#abs_path = os.path.join(data_dir, alignment.path)
			outf.write("%s\n"%alignment.path)	#using relative path, alignments are symlinked or transferred
		del outf
		bamListF = File(os.path.basename(outputFname))
		bamListF.addPFN(PFN("file://" + os.path.abspath(outputFname), site_handler))
		workflow.addFile(bamListF)
		sys.stderr.write("Done.\n")
		return bamListF
	
	def addPlatypusCallJob(self, workflow=None, executable=None, outputFile=None, \
					refFastaFList=None, \
					sourceVCFFile=None, interval=None, skipRegionsFile=None,\
					site_type=2, \
					extraArguments=None, job_max_memory=4000, no_of_cpus=1, \
					walltime=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					**keywords):
		"""
		2013.05.23 added "--assemble=1" to Platypus
		2013.05.16
			
			"--useCortex=..." is optional. the two results are of no difference.
			~/bin/Platypus/Platypus.py callVariants
				--bamFiles=individual_alignment/2311_639_1987079_GA_vs_3280_by_method2_realigned0.bam,..bam,...bam
				--regions=Scaffold56:8000001-11975416 --output=GATK/Scaffold56_8000001_11975416.orig.2.vcf
				--refFile=reference/3280_vervet_ref_6.0.3.fasta
				--useCortex=1
		
		use this argument to focus on certain regions
			--skipRegionsFile=SKIPREGIONSFILE
				region as comma-separated list of chromosome:start-end, or
				just list of chromosome, or nothing
			--outputRefCalls=OUTPUTREFCALLS, If 1, output block reference calls.

		"""
		
		if workflow is None:
			workflow = self
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		refFastaF = refFastaFList[0]
		frontArgumentList = ["callVariants", "--regions", interval, "--refFile", refFastaF]	#, "--assemble=1"
		if sourceVCFFile:	#only calling variants that exist in this vcf file
			#sourceVCFFile must be compressed by bgzip and tabix
			frontArgumentList.extend(["--source", sourceVCFFile, "--minPosterior=0", "--getVariantsFromBAMs=0"])
			extraDependentInputLs.append(sourceVCFFile)
		if skipRegionsFile:
			frontArgumentList.extend(["--skipRegionsFile", skipRegionsFile])
			extraDependentInputLs.append(skipRegionsFile)
		if site_type==1:	#all sites, not just variant sites
			frontArgumentList.append("--outputRefCalls=1")
		
		for refFastaFile in refFastaFList:
			extraDependentInputLs.append(refFastaFile)
		#job.uses(bamListF, transfer=True, register=True, link=Link.INPUT)
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, \
					outputFile=outputFile, outputArgumentOption='--output', inputFileList=None, \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=None, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		return job
		
	def addFreeBayesJob(self, workflow=None, executable=None, outputFile=None, \
					refFastaFList=None, \
					interval=None, heterozygosity=0.001, minCoverage=2,\
					site_type=2, \
					extraArguments=None, job_max_memory=4000, no_of_cpus=1, \
					walltime=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					**keywords):
		"""
		2013.10.07 https://github.com/ekg/freebayes
		
		"""
		
		if workflow is None:
			workflow = self
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		refFastaF = refFastaFList[0]
		
		frontArgumentList = ["--fasta-reference", refFastaF, "--min-coverage %s"%minCoverage,\
							"--ploidy 2", "--theta %s"%heterozygosity, ]
		if isinstance(interval, File):
			#Limit analysis to targets listed in the BED-format FILE.
			frontArgumentList.extend(["--targets", interval])
			extraDependentInputLs.append(interval)
		else:
			# Limit analysis to the specified region, 0-base coordinates,
			#end_position included.  Either '-' or '..' maybe used as a separator.
			frontArgumentList.extend(["--region", interval])
		if site_type==1:	#all sites, not just variant sites
			frontArgumentList.append("--report-monomorphic")
		
		for refFastaFile in refFastaFList:
			extraDependentInputLs.append(refFastaFile)
		#job.uses(bamListF, transfer=True, register=True, link=Link.INPUT)
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, \
					outputFile=outputFile, outputArgumentOption='--vcf', inputFileList=None, \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=None, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=None, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		return job
	
	def addGATKCallJob(self, workflow=None, genotyperJava=None, GenomeAnalysisTKJar=None, outputFile=None, \
					refFastaFList=None, \
					interval=None, \
					heterozygosityPrior=0.001, minPruningForGATKHaplotypCaller=None, GATKGenotypeCallerType = 'UnifiedGenotyper', \
					site_type=2, \
					extraArguments=None, job_max_memory=4000, no_of_gatk_threads=1, \
					walltime=None, \
					parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					**keywords):
		"""
		2013.2.25 added argument GATKGenotypeCallerType, (default is UnifiedGenotyper, could try HaplotypeCaller)
			added argument heterozygosityPrior, minPruningForGATKHaplotypCaller
			use addGenericJob
		2012.7.30
			interval could be a BED file, rather than just a string (chromosome:start-stop).
			start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
			
		2011-12-4
		"""
		if workflow is None:
			workflow = self
		if GenomeAnalysisTKJar is None:
			GenomeAnalysisTKJar = workflow.GenomeAnalysisTK2Jar
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		extraArgumentList = []
		
		# 2011-11-22 "-mmq 30",  is no longer included
		if site_type==1:
			extraArgumentList.append('--output_mode EMIT_ALL_SITES')	#2011-8-24 new GATK no longers ues "-all_bases"
		if heterozygosityPrior:
			extraArgumentList.append("--heterozygosity %s"%(heterozygosityPrior))
		if minPruningForGATKHaplotypCaller and GATKGenotypeCallerType=='HaplotypeCaller':
			extraArgumentList.append("--minPruning %s"%(minPruningForGATKHaplotypCaller))
		if GATKGenotypeCallerType=='UnifiedGenotyper':	#2013.2.25 only for UnifiedGenotyper, not for HaplotypeCaller
			extraArgumentList.append("-mbq 20 --baq CALCULATE_AS_NECESSARY ")
		
		job = self.addGATKJob(workflow=workflow, executable=genotyperJava, GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
							GATKAnalysisType=GATKGenotypeCallerType,\
					inputFile=None, inputArgumentOption=None, refFastaFList=refFastaFList, inputFileList=None,\
					interval=interval, outputFile=outputFile, \
					parentJobLs=parentJobLs, transferOutput=transferOutput, job_max_memory=job_max_memory,\
					frontArgumentList=None, extraArguments=extraArguments, extraArgumentList=extraArgumentList, \
					extraOutputLs=None, \
					extraDependentInputLs=extraDependentInputLs, no_of_cpus=no_of_gatk_threads, walltime=walltime, **keywords)
		return job
	
	def addSAMtoolsCallJob(self, workflow=None, CallVariantBySamtools=None, samtoolsOutputF=None, indelVCFOutputF=None, \
					refFastaFList=[], parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
					extraArguments=None, job_max_memory=2000, site_type=2, maxDP=15000, mpileupInterval=None, \
					bcftoolsInterval=None, walltime=None, no_of_cpus=None, **keywords):
		"""
		2013.2.25 use addGenericJob
		2012.8.7
			split argument interval into mpileupInterval and bcftoolsInterval
		2012.8.6 add argument maxDP
		2012.7.30
			mpileupInterval could be a BED file, rather than just a string (chromosome:start-stop). CallVariantBySamtools would adjust itself.
			start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.
		2011-12-4
		"""
		if workflow is None:
			workflow = self
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		refFastaF = refFastaFList[0]
		frontArgumentList = [refFastaF, mpileupInterval, bcftoolsInterval, samtoolsOutputF, repr(site_type), "%s"%(maxDP)]
		if isinstance(mpileupInterval, File):
			extraDependentInputLs.append(mpileupInterval)
		if isinstance(bcftoolsInterval, File) and mpileupInterval!=bcftoolsInterval:	#do not add the same file as INPUT twice
			extraDependentInputLs.append(bcftoolsInterval)
		for refFastaFile in refFastaFList:
			extraDependentInputLs.append(refFastaFile)
		#job.uses(bamListF, transfer=True, register=True, link=Link.INPUT)

		
		job = self.addGenericJob(workflow=workflow, executable=CallVariantBySamtools, inputFile=None, \
					outputFile=None, inputFileList=None, \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=[samtoolsOutputF, indelVCFOutputF], \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=None, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=None, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		job.output = samtoolsOutputF
		job.indelVCFOutputF = indelVCFOutputF
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
		
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "shell/trioCallerWrapper.sh"), \
												name="trioCallerWrapper", clusterSizeMultipler=0.001)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, \
												name="MergeVCFReplicateGenotypeColumnsJava", clusterSizeMultipler=0.001)
		
		ReplicateVCFGenotypeColumns = Executable(namespace=namespace, name="ReplicateVCFGenotypeColumns", version=version, os=operatingSystem,\
									arch=architecture, installed=True)
		ReplicateVCFGenotypeColumns.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ReplicateVCFGenotypeColumns.py"), site_handler))
		executableClusterSizeMultiplierList.append((ReplicateVCFGenotypeColumns, 1))
	
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		genotypingExecutableSet = set([self.genotyperJava, self.samtools, self.CallVariantBySamtools, \
									self.MergeVCFReplicateHaplotypesJava, self.Platypus, self.trioCallerWrapper, \
									self.MergeVCFReplicateGenotypeColumnsJava, self.freebayes])
		#2012.1.3. added polymutt
		if hasattr(self, 'polymuttPath'):	#some inherited classes do not have polymuttPath
			polymutt = self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.expanduser(self.polymuttPath), \
												name="polymutt", clusterSizeMultipler=0)
			genotypingExecutableSet.add(polymutt)
		
		if hasattr(self, 'noOfCallingJobsPerNode') and self.noOfCallingJobsPerNode>=0:
			for executable in genotypingExecutableSet:
				#2013.2.26 use setOrChangeExecutableClusterSize to modify clusters size
				self.setOrChangeExecutableClusterSize(executable=executable, clusterSizeMultipler=self.noOfCallingJobsPerNode/float(self.clusters_size), \
													defaultClustersSize=self.clusters_size)
	
	def addJobsToDeriveGenotypeCallingIntervals(self, workflow=None, alignmentLs=None):
		"""
		2013.08.11
		#. check/get depth.tsv.gz for each alignment
		#. add all depth.tsv.gz up and output it into a GADA input format
		#. segment it (GADA)
		#. output segmented intervals into one .bed file  (staged out in the end) and also small temporary .bed files (evenly by region sizes), under 5-fold threshold
			#.  only one chromosome intervals in one .bed file
			#.  keep track of chromosomes for each .bed file
		"""
		
		returnData = []
		for alignment in alignmentLs:
			alignmentDepthFile = self.registerOneInputFile(inputFname=alignment.path_to_depth_file, \
							input_site_handler=None, folderName="", useAbsolutePathAsPegasusFileName=False,\
							pegasusFileName=None, checkFileExistence=True)
		
		
		#. a job that adds all depth.tsv.gz up and output it into a GADA input format for each chromosome
		#. for each chromosome
		#.	run GADA job
		
		#.	output segmented intervals into one .bed file
		#.	filter out segments that beyond 2-fold depth range, or too small (bad alignments),
		#.		then output them in a fixed number of BED files for each chromosome, 
		
		sys.stderr.write(" %s jobs.\n"%(self.no_of_jobs))
		
	
	def run(self):
		"""
		2011-7-11
		"""
		
		self.needSplitChrIntervalData = False	#2013.09.01 as it'll be handled in run()
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		if self.selectedRegionFname and self.genotypeCallerType!=2:	#2013.05.29 not for Platypus
			if self.sourceVCFFolder:
				if os.path.isfile(self.selectedRegionFname):
					sys.stderr.write("Warning: selectedRegionFname %s already exists, will be overwritten.\n"%(self.selectedRegionFname))
					
				self.extractVCFSitesIntoBEDFile(inputVCFFolder=self.sourceVCFFolder, outputFname=self.selectedRegionFname)
			if os.path.isfile(self.selectedRegionFname):
				chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitBEDFile(intervalFname=self.selectedRegionFname, \
														noOfLinesPerUnit=self.maxNoOfRegionsPerJob, \
														folderName=self.pegasusFolderName, parentJobLs= None)
			
		else:
			chr2size = self.chr2size
			#self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, contigMinRankBySize=self.contigMinRankBySize)
			#chr2size = set(['Contig149'])	#temporary when testing Contig149
			#chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
			if self.alignmentDepthIntervalMethodShortName:	#2013.09.01 fetch intervals from db
				chr2IntervalDataLs = self.getChr2IntervalDataLsFromDBAlignmentDepthInterval(db=self.db, \
									intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize,\
									alignmentDepthIntervalMethodShortName=self.alignmentDepthIntervalMethodShortName, \
									alignmentDepthMinFold=self.alignmentDepthMinFold, alignmentDepthMaxFold=self.alignmentDepthMaxFold, \
									minAlignmentDepthIntervalLength=self.minAlignmentDepthIntervalLength,\
									maxContigID=self.maxContigID, minContigID=self.minContigID)
			else: #split evenly using chromosome size
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
		cumulativeMedianDepth = self.db.getCumulativeAlignmentMedianDepth(alignmentLs=pdata.alignmentLs, \
										defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		
		if self.run_type==1:	#multi-sample calling
			"""
			for extraArgumentPlatypus, outputDirPrefix in [("--maxHaplotypes=1500", "maxHaplotypes1500"), ("--maxGOF=100", "maxGOF100"),\
														("--sbThreshold=0", "sbThreshold0") ,("--abThreshold=0", "abThreshold0"),\
														("--rmsmqThreshold=0", "rmsmqThreshold0"), ("--qdThreshold=0", "qdThreshold0"),\
														("--hapScoreThreshold=0", "hapScoreThreshold1000")]:
			"""
			outputDirPrefix = ""
			extraArgumentPlatypus = ""
			self.addGenotypeCallJobs(workflow=workflow, alignmentDataLs=pdata.alignmentDataLs, \
						chr2IntervalDataLs=chr2IntervalDataLs,
						registerReferenceData=pdata.registerReferenceData, \
						site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						site_type=self.site_type, data_dir=self.data_dir, \
						outputDirPrefix=outputDirPrefix,\
						genotypeCallerType=self.genotypeCallerType,\
						cumulativeMedianDepth=cumulativeMedianDepth, sourceVCFFolder=self.sourceVCFFolder,\
						extraArgumentPlatypus=extraArgumentPlatypus)
		elif self.run_type==2:
			#2011-11-4 for single-alignment-calling pipeline, adjust the folder name so that they are unique from each other
			for alignmentData in pdata.alignmentDataLs:
				alignment = alignmentData.alignment
				
				self.addGenotypeCallJobs(workflow=workflow, alignmentDataLs=[alignmentData], \
						chr2IntervalDataLs=chr2IntervalDataLs, \
						registerReferenceData=pdata.registerReferenceData, \
						namespace=workflow.namespace, version=workflow.version, \
						site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
						needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
						site_type=self.site_type, data_dir=self.data_dir, \
						outputDirPrefix="alignmentID%d"%(alignment.id),
						genotypeCallerType=self.genotypeCallerType)
				
		
		self.end_run()
		


	
if __name__ == '__main__':
	main_class = AlignmentToCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
