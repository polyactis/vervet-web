#!/usr/bin/env python
"""
Examples:
	# 2011-8-30 workflow on condor, always commit (--commit)
	%s --ind_seq_id_ls 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool
		-n1 -z dl324b-1.cmb.usc.edu --commit --needSSHDBTunnel
	
	# 2011-8-30 a workflow with 454 long-read and short-read PE. need a ref index job (-n1). 
	%s --ind_seq_id_ls 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9
		-e /u/home/eeskin/polyacti -l hoffman2 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db -n1
		-z dl324b-1.cmb.usc.edu --commit
		--tmpDir /work/ --needSSHDBTunnel
	
	# 2011-8-30 output a workflow to run alignments on hoffman2's condor pool (--local_data_dir changes local_data_dir. --data_dir changes data_dir.)
	# 2012.3.20 use /work/ or /u/scratch/p/polyacti/tmp as TMP_DIR for MarkDuplicates.jar (/tmp is too small for 30X genome)
	# 2012.5.4 cluster 4 alignment jobs (before merging) as a unit (--alignmentJobClustersSizeFraction 0.2), skip done alignment (--skipDoneAlignment)
	# 2012.9.21 add "--needSSHDBTunnel" because AddAlignmentFile2DB need db conneciton
	# 2012.9.21 add "--alignmentPerLibrary" to also get alignment for each library within one individual_sequence
	# 2013.3.15 add "--coreAlignmentJobWallTimeMultiplier 0.5" to reduce wall time for core-alignment (bwa/stampy) jobs by half 
	ref=3280; %s --ind_seq_id_ls 632-3230 --sequence_min_coverage 15 --sequence_max_coverage 80 --site_id_ls 447 --sequence_filtered 1
		--excludeContaminant -a $ref -o dags/ShortRead2Alignment/ShortRead2AlignmentPipeline_VRCPart1_vs_$ref\_AlnMethod2.xml
		-u yh -l hcondor -j hcondor -z localhost -u yh --commit --tmpDir /work/
		--home_path /u/home/eeskin/polyacti --no_of_aln_threads 1 --skipDoneAlignment
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--clusters_size 20 --alignment_method_name bwaShortRead
		--coreAlignmentJobWallTimeMultiplier 0.5
		--alignmentJobClustersSizeFraction 0.2 
		--needSSHDBTunnel --ref_genome_version 2 --needRefIndexJob --db_passwd secret
		#--alignmentPerLibrary
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files. need ref index job (--needRefIndexJob).
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s --ind_seq_id_ls 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  --commit --needSSHDBTunnel
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files.
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s --ind_seq_id_ls 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  --commit --needSSHDBTunnel
		
	# 2011-8-30 a workflow to run on uschpc, with ref index job. Note the site_handler and input_site_handler.
	# to enable replica-transfer.
	%s --ind_seq_id_ls 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j local -l uschpc -n1 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10 -p secret  --commit --needSSHDBTunnel

	# 2011-8-30 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler and "--data_dir ..." to enable symlink
	%s --ind_seq_id_ls 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j uschpc -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4 --needSSHDBTunnel
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		--data_dir /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ --javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
	
	# 2011-11-16 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler. this will stage in all input and output (--notStageOutFinalOutput).
	%s --ind_seq_id_ls 391-397,456,473,493
		-o dags/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_local2usc.xml -u yh -a 9
		-j local -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		--javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
		--needSSHDBTunnel
	
	
	#2011-9-13 no ref index job, staging input files from localhost to uschpc, stage output files back to localhost
	# modify the refFastaFile's path in xml manually
	%s --ind_seq_id_ls 1-3 -o ShortRead2Alignment_1_3_vs_524_local2uschpc.xml -u yh -a 524
		-j local -l uschpc --needRefIndexJob -p secret --commit --no_of_aln_threads 4 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
		--data_dir /Network/Data/vervet/db/
		--needSSHDBTunnel
	
	# 2011-8-31 output the same workflow above but for condorpool
	%s --ind_seq_id_ls 391-397,456,473,493, -o dags/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z 10.8.0.10  -p secret  --commit --alignmentPerLibrary
	
	# 2012-4-5 new alignment method, stampy (--alignment_method_name)
	%s --ind_seq_id_ls 167,176,178,182,183,207-211,391-397,456,473,493
		-o dags/ShortRead2Alignment/ShortRead2Alignment_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
		-u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  --commit --alignment_method_name stampy
	
	# 2013.2.28 use the new alignment-method: bwaShortReadHighMismatches
	#double the core alignment (bwa aln) job walltime (=23 hrs) (--coreAlignmentJobWallTimeMultiplier) because it takes much longer
	# set max walltime for any job to be 1 day (--max_walltime 1440)
	ref=3231; %s --ind_seq_id_ls 638 -a $ref
		-o dags/ShortRead2Alignment/ShortRead2AlignmentPipeline_Aethiops_vs_$ref\_AlnMethod5.xml
		-u yh -l hcondor -j hcondor -z localhost -u yh --commit --tmpDir /work/ --home_path /u/home/eeskin/polyacti
		--no_of_aln_threads 1 --skipDoneAlignment -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ --clusters_size 1 --alignment_method_name bwaShortReadHighMismatches
		--coreAlignmentJobWallTimeMultiplier 2  --needSSHDBTunnel
		--max_walltime 1440
	
	# 2013.04.04 use the new alignment-method: bwa-mem, get rid of "-q 20" by --additionalArguments " " as mem doesn't support -q.
	# 23*0.1 hrs walltime for the core alignment (bwa mem) jobs (--coreAlignmentJobWallTimeMultiplier 0.1) because it's much faster
	# set max walltime for any job to be 1 day (--max_walltime 1440)
	ref=3280; %s --ind_seq_id_ls 638 -a $ref --additionalArguments " "
		-o dags/ShortRead2Alignment/ShortRead2AlignmentPipeline_Aethiops_vs_$ref\_AlnMethod6.xml
		-u yh -l hcondor -j hcondor -z localhost -u yh --commit --tmpDir /work/ --home_path /u/home/eeskin/polyacti
		--no_of_aln_threads 1 --skipDoneAlignment -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ --clusters_size 1 --alignment_method_name mem
		--coreAlignmentJobWallTimeMultiplier 0.1  --needSSHDBTunnel
		--max_walltime 1440
Description:
	2013.04.07
		A program which generates a pegasus workflow dag (xml file) which does the alignment for all available sequences.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
		It will also stage out every output file.
		Be careful about -R, only toggle it if you know every input individual_sequence_file is not empty.
			Empty read files would fail alignment jobs and thus no final alignment for a few indivdiuals.
		Use "--alignmentJobClustersSizeFraction ..." to cluster the alignment jobs if the input read file is small enough (~1Million reads for bwa, ~300K for stampy).
		The arguments related to how many chromosomes/contigs do not matter unless local_realigned=1. 
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus import yh_pegasus
from pymodule.pegasus.ShortRead2AlignmentWorkflow import ShortRead2AlignmentWorkflow
from vervet.src import VervetDB

class ShortRead2AlignmentPipeline(ShortRead2AlignmentWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(ShortRead2AlignmentWorkflow.option_default_dict)
	option_default_dict.pop(('refSequenceFname', 1, ))
	option_default_dict.update({
						('skipDoneAlignment', 0, int):[0, '', 0, 'skip alignment whose path is a valid file'],\
						("alignmentPerLibrary", 0, int): [0, '', 0, 'toggle to run alignment for each library of one individual_sequence'],\
						})
	option_default_dict[('local_realigned', 0, int)][0] = 0
	
	def __init__(self,  **keywords):
		"""
		2012.3.29
			default to stage out final output.
			Argument stageOutFinalOutput morphs into notStageOutFinalOutput.
		2011-7-11
		"""
		ShortRead2AlignmentWorkflow.__init__(self, **keywords)
		
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
	
	def addAllAlignmentJobs(self, db_vervet=None, individualSequenceID2FilePairLs=None, \
					data_dir=None, \
					isqLs=None,\
					refSequence=None, registerReferenceData=None, \
					chr2IntervalDataLs=None,\
					workflow=None, bwa=None, additionalArguments=None, samtools=None, mkdirWrap=None, mv=None,\
					java=None, MergeSamFilesJava=None, MergeSamFilesJar=None, \
					MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir='/tmp',\
					BuildBamIndexFilesJava=None, BuildBamIndexJar=None, \
					SortSamFilesJava=None, SortSamJar=None, \
					addOrReplaceReadGroupsJava=None, AddOrReplaceReadGroupsJar=None,\
					alignment_method_name='bwaShortRead', alignment_format='bam',\
					namespace='workflow', version='1.0', transferOutput=False,\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					no_of_aln_threads=3, stampy=None, skipDoneAlignment=False, \
					alignmentPerLibrary=False, outputDirPrefix="", **keywords):
		"""
		2013.03.29 scale MarkDuplicates and MergeSam jobs memory & walltime to the sequence coverage
		2013.3.13 individual_sequence.sequencer and individual_sequence.sequence_type has changed its meaning
		2012.9.19
			isq_id2LibrarySplitOrder2FileLs is replaced by isqLs.
			add argument alignmentPerLibrary
		2012.4.20
			bugfix, pass alignment_method.short_name instead of alignment_method_name to db_vervet.getAlignment()
				because alignment_method might be changed according to sequencer regardless of alignment_method_name.
		2012.4.12
			add argument skipDoneAlignment
		2012.4.5
			fetch the alignment_method directly based on the alignment_method_name, except for 454 sequences.
				officially merges "bwa-short-read-SR" (single-read) into "bwa-short-read"
		2012.2.24
			pass data_dir to db_vervet.getAlignment()
			add stampy part
		2011-9-15
			adjust alignment_method_name according to individual_sequence.sequencer and individual_sequence.sequence_type
			only when this is not possible, value of argument alignment_method_name is used.
		2011-9-14
			give different names to different java jobs according to jars
		2011-8-28
		"""
		sys.stderr.write("Adding alignment jobs for %s individual sequences ..."%(len(isqLs)))
		
		no_of_alignment_jobs = 0
		no_of_merging_jobs = 0
		
		alignmentFolder = "%sAlignment"%(outputDirPrefix)
		alignmentFolderJob = self.addMkDirJob(outputDir=alignmentFolder)
		
		oneLibraryAlignmentFolder = "%sOneLibAlignment"%(outputDirPrefix)
		oneLibraryAlignmentFolderJob = self.addMkDirJob(outputDir=oneLibraryAlignmentFolder)
		
		refFastaFList = registerReferenceData.refFastaFList
		refIndexJob = None
		if self.needRefIndexJob or registerReferenceData.needBWARefIndexJob or registerReferenceData.needStampyRefIndexJob:
			if self.alignment_method_name.find('bwa')>=0 and registerReferenceData.needBWARefIndexJob:
				refIndexJob = self.addBWAReferenceIndexJob(workflow, refFastaFList=refFastaFList, \
												refSequenceBaseCount=refSequence.base_count, bwa=workflow.bwa)
			elif self.alignment_method_name.find('stampy')>=0 and registerReferenceData.needStampyRefIndexJob:
				refIndexJob = self.addStampyGenomeIndexHashJob(workflow, executable=workflow.stampy, refFastaFList=refFastaFList, \
						parentJobLs=None, job_max_memory=3000, walltime = 200, \
						extraDependentInputLs=None, \
						transferOutput=True)
		
		for individual_sequence  in isqLs:
			if individual_sequence is not None and individual_sequence.format=='fastq':
				library2Data = individual_sequence.library2Data
				AlignmentJobAndOutputLs = []
				# get or add alignment method
				if individual_sequence.sequencer.short_name=='454' or \
						(individual_sequence.sequence_type and individual_sequence.sequence_type.read_length_mean is not None \
						and individual_sequence.sequence_type.read_length_mean>150):
					alignment_method = db_vervet.getAlignmentMethod("bwaLongRead")
				else:
					alignment_method = db_vervet.getAlignmentMethod(alignment_method_name)
				"""
				#2012.4.5 have all this commented out
				elif individual_sequence.sequencer.short_name=='GA':
					if individual_sequence.sequence_type=='SR':	#single-end
						alignment_method = db_vervet.getAlignmentMethod("bwa-short-read-SR")
					else:	#default is PE
						alignment_method = db_vervet.getAlignmentMethod("bwa-short-read")
				"""
				#alignment for the whole individual_sequence
				individual_alignment = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
												individual_sequence_id=individual_sequence.id,\
									path_to_original_alignment=None, sequencer_id=individual_sequence.sequencer_id, \
									sequence_type_id=individual_sequence.sequence_type_id, sequence_format=individual_sequence.format, \
									ref_individual_sequence_id=refSequence.id, \
									alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
									individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
									data_dir=data_dir, local_realigned=self.local_realigned)
				skipIndividualAlignment = False
				if individual_alignment.path:
					if skipDoneAlignment and db_vervet.isThisAlignmentComplete(individual_alignment=individual_alignment, data_dir=data_dir):
						skipIndividualAlignment = True
				#2012.3.29 this folder will store the alignment output by the alignment jbos
				tmpOutputDir = os.path.basename(individual_sequence.path)
				# add a mkdir job
				mkdirJob = None
				for library, pdata in library2Data.iteritems():
					minIsqFileRawID = min(pdata.isqFileRawID2Index.keys())
					splitOrder2Index = pdata.splitOrder2Index
					fileObjectPairLs = pdata.fileObjectPairLs
					
					oneLibraryAlignmentJobAndOutputLs = []
					splitOrderLs = splitOrder2Index.keys()
					splitOrderLs.sort()
					oneLibraryCumulativeBaseCount = 0
					if alignmentPerLibrary:
						#alignment for this library of the individual_sequence
						oneLibraryAlignmentEntry = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
												individual_sequence_id=individual_sequence.id,\
									path_to_original_alignment=None, sequencer_id=individual_sequence.sequencer_id, \
									sequence_type_id=individual_sequence.sequence_type_id, sequence_format=individual_sequence.format, \
									ref_individual_sequence_id=refSequence.id, \
									alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
									individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
									data_dir=data_dir, individual_sequence_file_raw_id=minIsqFileRawID,\
									local_realigned=self.local_realigned)
						skipLibraryAlignment = False
						if oneLibraryAlignmentEntry.path:
							if skipDoneAlignment and db_vervet.isThisAlignmentComplete(individual_alignment=individual_alignment, data_dir=data_dir):
								# 2013.3.22
								# file_size is updated in the last of AddAlignmentFile2DB.py.
								# if it fails in the middle of copying, file_size would be None.
								skipLibraryAlignment = True
					else:
						skipLibraryAlignment = True
					if skipIndividualAlignment and skipLibraryAlignment:	#2012.9.19 if both skip flags are true, then yes
						continue
					for splitOrder in splitOrderLs:
						splitOrderIndex = splitOrder2Index[splitOrder]
						fileObjectLs = fileObjectPairLs[splitOrderIndex]
						oneLibraryCumulativeBaseCount += sum([fileObject.db_entry.base_count for fileObject in fileObjectLs])
						if mkdirJob is None:	#now it's time to add the mkdirJob
							# add a mkdir job
							mkdirJob = self.addMKDIRJob(workflow, mkdirWrap=mkdirWrap, dirName=tmpOutputDir)
						#newFilePair = self.registerFileToWorkflow(filePair, workflow, data_dir=data_dir)
						newFileObjLs = self.registerISQFileObjLsToWorkflow(fileObjectLs=fileObjectLs, workflow=workflow)
						#2012.9.19 individual_alignment is passed as None so that ReadGroup addition job is not added in addAlignmentJob()
						alignmentJob, alignmentOutput = self.addAlignmentJob(workflow=workflow, fileObjectLs=newFileObjLs, \
																			individual_alignment=None, \
							data_dir=data_dir, refFastaFList=refFastaFList, bwa=bwa, \
							additionalArguments=additionalArguments, samtools=samtools, \
							refIndexJob=refIndexJob, parentJobLs=[refIndexJob, mkdirJob], \
							alignment_method=alignment_method, \
							outputDir=tmpOutputDir, namespace=namespace, version=version,\
							PEAlignmentByBWA=PEAlignmentByBWA, ShortSEAlignmentByBWA=ShortSEAlignmentByBWA, \
							LongSEAlignmentByBWA=LongSEAlignmentByBWA,\
							java=java, SortSamFilesJava=SortSamFilesJava, SortSamJar=SortSamJar,\
							addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
							no_of_aln_threads=no_of_aln_threads, stampy=stampy)
						no_of_alignment_jobs += 1
						
						fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
														os.path.basename(alignmentOutput.name))[0]
						if not skipIndividualAlignment:
							#2012.9.19 add a AddReadGroup job
							outputRGBAM = File(os.path.join(tmpOutputDir, "%s.isq_RG.bam"%(fname_prefix)))
							addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=individual_alignment, \
												inputBamFile=alignmentJob.output, \
												outputBamFile=outputRGBAM,\
												addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, \
												AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
												parentJobLs=[alignmentJob, mkdirJob], extraDependentInputLs=None, \
												extraArguments=None, job_max_memory = 2500, transferOutput=False)
							AlignmentJobAndOutputLs.append(PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
						if not skipLibraryAlignment:
							#2012.9.19 add a AddReadGroup job for the library bam file
							outputRGBAM = File(os.path.join(tmpOutputDir, "%s.isq_library_%s_RG.bam"%(fname_prefix, library)))
							addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=oneLibraryAlignmentEntry, \
												inputBamFile=alignmentJob.output, \
												outputBamFile=outputRGBAM,\
												addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, \
												AddOrReplaceReadGroupsJar=AddOrReplaceReadGroupsJar,\
												parentJobLs=[alignmentJob, mkdirJob], extraDependentInputLs=None, \
												extraArguments=None, job_max_memory = 2500, transferOutput=False)
							oneLibraryAlignmentJobAndOutputLs.append(PassingData(parentJobLs=[addRGJob], file=addRGJob.output))
					if alignmentPerLibrary and not skipLibraryAlignment and oneLibraryAlignmentJobAndOutputLs:	#2012.9.19
						baseCoverage = 4*3000000000	#baseline
						minMergeAlignmentWalltime = 240	#in minutes, 4 hours, when coverage is defaultCoverage
						maxMergeAlignmentWalltime = 2880	#in minutes, 2 days
						minMergeAlignmentMaxMemory = 5000	#in MB, when coverage is defaultCoverage
						maxMergeAlignmentMaxMemory = 12000	#in MB
						
						mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=oneLibraryCumulativeBaseCount, \
												baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentWalltime, \
												minJobPropertyValue=minMergeAlignmentWalltime, maxJobPropertyValue=maxMergeAlignmentWalltime).value
						mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=oneLibraryCumulativeBaseCount, \
												baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentMaxMemory, \
												minJobPropertyValue=minMergeAlignmentMaxMemory, maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
						markDuplicateWalltime= mergeAlignmentWalltime
						markDuplicateMaxMemory = mergeAlignmentMaxMemory
						fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
														os.path.basename(oneLibraryAlignmentEntry.constructRelativePath()))[0]
						mergedBamFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_merged.bam'%(fname_prefix, library)))
						alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, \
											AlignmentJobAndOutputLs=oneLibraryAlignmentJobAndOutputLs, \
											outputBamFile=mergedBamFile, \
											samtools=samtools, java=java, \
											MergeSamFilesJava=MergeSamFilesJava, MergeSamFilesJar=MergeSamFilesJar, \
											BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, \
											BuildBamIndexJar=BuildBamIndexJar, \
											mv=mv, namespace=namespace, version=version, \
											transferOutput=False, \
											job_max_memory=mergeAlignmentMaxMemory, walltime=mergeAlignmentWalltime, \
											parentJobLs=[oneLibraryAlignmentFolderJob])
						
						finalBamFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_dupMarked.bam'%(fname_prefix, library)))
						markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[alignmentMergeJob, bamIndexJob], \
								inputBamF=alignmentMergeJob.output, \
								inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
								MarkDuplicatesJava=MarkDuplicatesJava, MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
								BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=BuildBamIndexJar, \
								namespace=namespace, version=version, job_max_memory=markDuplicateMaxMemory, \
								walltime=markDuplicateWalltime, transferOutput=False)
						no_of_merging_jobs += 1
						
						if self.local_realigned:
							alignmentData = PassingData(jobLs=[markDupJob, markDupBamIndexJob], bamF=markDupJob.output, \
													baiF=markDupBamIndexJob.output, alignment=oneLibraryAlignmentEntry)
							#2013.03.31
							preDBAlignmentJob, preDBAlignmentIndexJob = self.addLocalRealignmentSubWorkflow(workflow=workflow, \
								chr2IntervalDataLs=chr2IntervalDataLs, \
								registerReferenceData=registerReferenceData, \
								alignmentData=alignmentData,\
								inputBamF=markDupJob.output, \
								outputBamF=None, \
								parentJobLs=[markDupJob, markDupBamIndexJob], \
								outputDirPrefix='%s_%s_localRealignment'%(fname_prefix, library), transferOutput=False)
						else:
							preDBAlignmentJob = markDupJob
							preDBAlignmentIndexJob = markDupBamIndexJob
						#2012.9.19 add/copy the alignment file to db-affliated storage
						#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
						logFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_2db.log'%(fname_prefix, library)))
						alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
											inputFile=preDBAlignmentJob.output, \
											otherInputFileList=[], \
											individual_alignment_id=oneLibraryAlignmentEntry.id, \
											individual_sequence_file_raw_id=minIsqFileRawID,\
											format=None, local_realigned=self.local_realigned,\
											logFile=logFile, data_dir=data_dir, \
											parentJobLs=[preDBAlignmentJob, preDBAlignmentIndexJob], \
											extraDependentInputLs=[preDBAlignmentIndexJob.output,], \
											extraArguments=None, transferOutput=transferOutput, \
											job_max_memory=2000, walltime=180, \
											sshDBTunnel=self.needSSHDBTunnel, commit=True)
						
				if AlignmentJobAndOutputLs and not skipIndividualAlignment:
					baseCoverage = 4	#baseline
					actualCoverage = getattr(individual_sequence, 'coverage', baseCoverage)
					minMergeAlignmentWalltime = 240	#in minutes, 4 hours, when coverage is defaultCoverage
					maxMergeAlignmentWalltime = 2880	#in minutes, 2 days
					minMergeAlignmentMaxMemory = 5000	#in MB, when coverage is defaultCoverage
					maxMergeAlignmentMaxMemory = 12000	#in MB
					
					mergeAlignmentWalltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=actualCoverage, \
											baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentWalltime, \
											minJobPropertyValue=minMergeAlignmentWalltime, maxJobPropertyValue=maxMergeAlignmentWalltime).value
					mergeAlignmentMaxMemory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=actualCoverage, \
											baseInputVolume=baseCoverage, baseJobPropertyValue=minMergeAlignmentMaxMemory, \
											minJobPropertyValue=minMergeAlignmentMaxMemory, maxJobPropertyValue=maxMergeAlignmentMaxMemory).value
					markDuplicateWalltime= mergeAlignmentWalltime
					markDuplicateMaxMemory = mergeAlignmentMaxMemory
					
					#2012.3.29	merge alignment output only when there is something to merge!
					fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
													os.path.basename(individual_alignment.constructRelativePath()))[0]
					mergedBamFile = File(os.path.join(alignmentFolder, '%s_merged.bam'%(fname_prefix)))
					alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
										outputBamFile=mergedBamFile, \
										samtools=samtools, java=java, \
										MergeSamFilesJava=MergeSamFilesJava, MergeSamFilesJar=MergeSamFilesJar, \
										BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=BuildBamIndexJar, \
										mv=mv, namespace=namespace, version=version, \
										parentJobLs=[alignmentFolderJob],\
										transferOutput=False, job_max_memory=mergeAlignmentMaxMemory, \
										walltime=mergeAlignmentWalltime)
					#relative path in the scratch
					finalBamFile = File(os.path.join(alignmentFolder, '%s_dupMarked.bam'%(fname_prefix)))
					
					markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[alignmentMergeJob, bamIndexJob],
							inputBamF=alignmentMergeJob.output, \
							inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
							MarkDuplicatesJava=MarkDuplicatesJava, MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
							BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexJar=BuildBamIndexJar, \
							namespace=namespace, version=version, job_max_memory=markDuplicateMaxMemory, \
							walltime=markDuplicateWalltime, \
							transferOutput=False)
					no_of_merging_jobs += 1
					
					
					if self.local_realigned:
						#2013.03.31
						alignmentData = PassingData(jobLs=[markDupJob, markDupBamIndexJob], bamF=markDupJob.output, \
													baiF=markDupBamIndexJob.output, alignment=individual_alignment)
						preDBAlignmentJob, preDBAlignmentIndexJob = self.addLocalRealignmentSubWorkflow(workflow=workflow, \
							chr2IntervalDataLs=chr2IntervalDataLs, registerReferenceData=registerReferenceData, \
							alignmentData=alignmentData,\
							inputBamF=markDupJob.output, \
							outputBamF=None, \
							parentJobLs=[markDupJob, markDupBamIndexJob], \
							outputDirPrefix='%s_localRealignment'%(fname_prefix), transferOutput=False)
					else:
						preDBAlignmentJob = markDupJob
						preDBAlignmentIndexJob = markDupBamIndexJob
					#2012.9.19 add/copy the alignment file to db-affliated storage
					#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
					logFile = File(os.path.join(alignmentFolder, '%s_2db.log'%(fname_prefix)))
					alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
										inputFile=preDBAlignmentJob.output, \
										otherInputFileList=[],\
										individual_alignment_id=individual_alignment.id, \
										format=None, local_realigned=self.local_realigned,\
										logFile=logFile, data_dir=data_dir, \
										parentJobLs=[preDBAlignmentJob, preDBAlignmentIndexJob], \
										extraDependentInputLs=[preDBAlignmentIndexJob.output], \
										extraArguments=None, transferOutput=transferOutput, \
										job_max_memory=2000, walltime=max(180, markDuplicateWalltime/4), \
										sshDBTunnel=self.needSSHDBTunnel, commit=True)
					
		sys.stderr.write("%s alignment jobs; %s merge alignment jobs; %s jobs in total.\n"%(no_of_alignment_jobs, no_of_merging_jobs,\
																				self.no_of_jobs))
	
	def registerFileToWorkflow(self, filePair, workflow=None, data_dir=None):
		'''
		2011-8-30
		'''
		
		newFilePair = []
		for fileRecord in filePair:
			relativePath = fileRecord[0]
			fastqF = File(relativePath)
			fastqF.addPFN(PFN("file://" + os.path.join(data_dir, relativePath), self.input_site_handler))
			workflow.addFile(fastqF)
			newFileRecord = [fastqF] + fileRecord[1:]
			newFilePair.append(newFileRecord)
		return newFilePair
	
	def addMKDIRJob(self, workflow=None, mkdirWrap=None, dirName=None, namespace='workflow', version='1.0'):
		"""
		2012.10.10 call yh_pegasus.addMkDirJob
		"""
		if workflow is None:
			workflow = self
		return self.addMkDirJob(outputDir=dirName)
		#mkdirJob = Job(namespace=namespace, name=mkdirWrap.name, version=version)
		#mkdirJob.addArguments(dirName)
		#workflow.addJob(mkdirJob)
		#return mkdirJob
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, db_user=self.db_user,
					db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		session = db_vervet.session
		session.begin()
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		# Create a abstract dag
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerCustomJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=self.chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		#individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, data_dir=self.local_data_dir)
		isqLs = db_vervet.getISQDBEntryLsForAlignment(self.ind_seq_id_ls, data_dir=self.data_dir, \
												filtered=self.sequence_filtered, ignoreEmptyReadFile=self.ignoreEmptyReadFile)
		isqLs = db_vervet.filterIndividualSequenceList(individual_sequence_list=isqLs, min_coverage=self.sequence_min_coverage,\
						max_coverage=self.sequence_max_coverage, \
						individual_site_id_set=set(self.site_id_ls), individual_id_set=None, \
						sequence_type_id_set=set(self.sequence_type_id_ls),\
						sequencer_id_set=set(self.sequencer_id_ls), sequence_filtered=self.sequence_filtered,\
						sequence_batch_id_set=set(self.sequence_batch_id_ls), parent_individual_sequence_id_set=None, \
						version_set=set(self.version_ls),\
						country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls), \
						excludeContaminant=self.excludeContaminant, \
						report=True)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		
		#2011-11-16 new way of registering reference fasta file. but still dont' want to trasnfer 7Gb of data
		refFastaFname = os.path.join(self.data_dir, refSequence.path)
		registerReferenceData = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		
		self.addAllAlignmentJobs(db_vervet=db_vervet, individualSequenceID2FilePairLs=None, \
					isqLs = isqLs,\
					data_dir=self.data_dir,\
					refSequence=refSequence, registerReferenceData=registerReferenceData, \
					chr2IntervalDataLs=chr2IntervalDataLs,\
					workflow=workflow, bwa=workflow.bwa, additionalArguments=self.additionalArguments, \
					samtools=workflow.samtools, \
					mkdirWrap=workflow.mkdirWrap, mv=workflow.cp, \
					java=workflow.java, \
					MergeSamFilesJava=workflow.MergeSamFilesJava, MergeSamFilesJar=workflow.MergeSamFilesJar, \
					MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=self.tmpDir,\
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar, \
					SortSamFilesJava=workflow.SortSamFilesJava, SortSamJar=workflow.SortSamJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, AddOrReplaceReadGroupsJar=workflow.AddOrReplaceReadGroupsJar,\
					alignment_method_name=self.alignment_method_name, alignment_format='bam',\
					
					namespace=workflow.namespace, version=workflow.version, transferOutput=self.stageOutFinalOutput,\
					PEAlignmentByBWA=workflow.PEAlignmentByBWA, ShortSEAlignmentByBWA=workflow.ShortSEAlignmentByBWA, \
					LongSEAlignmentByBWA=workflow.LongSEAlignmentByBWA, no_of_aln_threads=self.no_of_aln_threads,\
					stampy=workflow.stampy, skipDoneAlignment=self.skipDoneAlignment,\
					alignmentPerLibrary=self.alignmentPerLibrary, outputDirPrefix="")
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)
		
		if self.commit:
			session.commit()
		else:
			session.rollback()
	
if __name__ == '__main__':
	main_class = ShortRead2AlignmentPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
