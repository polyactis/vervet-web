#!/usr/bin/env python
"""
Examples:
	# 2011-8-30 workflow on condor, always commit (-c)
	%s -i 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9 -l condorpool
		-n1 -z dl324b-1.cmb.usc.edu -c -H
	
	# 2011-8-30 a workflow with 454 long-read and short-read PE. need a ref index job (-n1). 
	%s -i 165-167 -o ShortRead2Alignment_isq_id_165_167_vs_9.xml -u yh -a 9
		-e /u/home/eeskin/polyacti -l hoffman2 --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db -n1
		-z dl324b-1.cmb.usc.edu -c
		--tmpDir /work/ -H
	
	# 2011-8-30 output a workflow to run alignments on hoffman2's condor pool (--local_data_dir changes local_data_dir. --data_dir changes data_dir.)
	# 2012.3.20 use /work/ or /u/scratch/p/polyacti/tmp as TMP_DIR for MarkDuplicates.jar (/tmp is too small for 30X genome)
	# 2012.5.4 cluster 10 alignment jobs (before merging) as a unit (--cluster_size_for_aln_jobs 10), skip done alignment (-K)
	# 2012.9.21 add "-H" because AddAlignmentFile2DB need db conneciton
	# 2012.9.21 add "--alignmentPerLibrary" to also get alignment for each library within one individual_sequence
	%s  --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-l hcondor -j hcondor 
		-z localhost -u yh -c
		-i 631-700 -o workflow/ShortRead2Alignment_Isq_631-700_vs_524_hcondor.xml  -a 524 
		--tmpDir /work/ -e /u/home/eeskin/polyacti  --cluster_size_for_aln_jobs 10 -K  -H --alignmentPerLibrary
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files. need ref index job (--needRefIndexJob).
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  -c -H
	
	# 2011-8-30 a workflow to run on condorpool, no ref index job. Note the site_handler and input_site_handler are both condorpool
	# to enable symlink of input files.
	# If input_site_handler is "local", pegasus will report error saying it doesn't know how to replica-transfer input files.
	%s -i 176,178-183,207-211
		-o ShortRead2Alignment_8VWP_vs_9_condor_no_refIndex.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z dl324b-1.cmb.usc.edu -p secret  -c -H
		
	# 2011-8-30 a workflow to run on uschpc, with ref index job. Note the site_handler and input_site_handler.
	# to enable replica-transfer.
	%s -i 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j local -l uschpc -n1 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10 -p secret  -c -H

	# 2011-8-30 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler and "--data_dir ..." to enable symlink
	%s -i 391-397,456,473,493
		-o ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_uschpc.xml -u yh -a 9
		-j uschpc -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4 -H
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		--data_dir /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ --javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java
	
	# 2011-11-16 a workflow to run on uschpc, Need ref index job (--needRefIndexJob), and 4 threads for each alignment job 
	# Note the site_handler, input_site_handler. this will stage in all input and output (--notStageOutFinalOutput).
	%s -i 391-397,456,473,493
		-o workflow/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_local2usc.xml -u yh -a 9
		-j local -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4
		-e /home/cmb-03/mn/yuhuang -z 10.8.0.10 
		--javaPath /home/cmb-03/mn/yuhuang/bin/jdk/bin/java  -H
	
	
	#2011-9-13 no ref index job, staging input files from localhost to uschpc, stage output files back to localhost
	# modify the refFastaFile's path in xml manually
	%s -i 1-3 -o ShortRead2Alignment_1_3_vs_524_local2uschpc.xml -u yh -a 524
		-j local -l uschpc --needRefIndexJob -p secret -c --no_of_aln_threads 4 -e /home/cmb-03/mn/yuhuang -z 10.8.0.10
		--data_dir /Network/Data/vervet/db/  -H
	
	# 2011-8-31 output the same workflow above but for condorpool
	%s -i 391-397,456,473,493, -o workflow/ShortRead2Alignment/ShortRead2Alignment_4DeepVRC_6LowCovVRC_392_397_vs_9_condorpool.xml
		-u yh -a 9 -j condorpool -l condorpool --needRefIndexJob -z 10.8.0.10  -p secret  -c --alignmentPerLibrary  -H
	
	# 2012-4-5 new alignment method, stampy (--alignment_method_name)
	%s -i 167,176,178,182,183,207-211,391-397,456,473,493
		-o workflow/ShortRead2Alignment/ShortRead2Alignment_10VWP_4DeepVRC_6LowCovVRC_392_397_vs_508_condorpool.xml
		-u yh -a 508 -j condorpool -l condorpool -n1 -z 10.8.0.10  -p secret  -c --alignment_method_name stampy  -H
	
Description:
	2012.5.3
		a program which generates a pegasus workflow dag (xml file) which does the alignment for all available sequences.
		The workflow will stage in (or symlink if site_handler and input_site_handler are same.) every input file.
		It will also stage out every output file.
		Be careful about -R, only toggle it if you know every input individual_sequence_file is not empty.
			Empty read files would fail alignment jobs and thus no final alignment for a few indivdiuals.
		Use "--cluster_size_for_aln_jobs 10" to cluster the alignment jobs if the input read file is small enough (~1Million reads for bwa, ~300K for stampy). 
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, yh_pegasus
from pymodule.pegasus.ShortRead2AlignmentWorkflow import ShortRead2AlignmentWorkflow
from vervet.src import VervetDB


class ShortRead2AlignmentPipeline(ShortRead2AlignmentWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(ShortRead2AlignmentWorkflow.option_default_dict)
	option_default_dict.pop(('refSequenceFname', 1, ))
	option_default_dict.update({
						('ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 1, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. \
									non-fastq entries will be discarded.', ],\
						('skipDoneAlignment', 0, int):[0, 'K', 0, 'skip alignment whose path is a valid file'],\
						("alignmentPerLibrary", 0, int): [0, '', 0, 'toggle to run alignment for each library of one individual_sequence'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_alignment and/or individual_alignment.path'],\
						})

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
	
	def addAllAlignmentJobs(self, db_vervet, individualSequenceID2FilePairLs=None, \
					data_dir=None, \
					isqLs=None,\
					refSequence=None, refFastaFList=None, refIndexJob=None,
					workflow=None, bwa=None, additionalArguments=None, samtools=None, mkdirWrap=None, mv=None,\
					java=None, mergeSamFilesJava=None, mergeSamFilesJar=None, \
					MarkDuplicatesJava=None, MarkDuplicatesJar=None, tmpDir='/tmp',\
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					SortSamFilesJava=None, SortSamFilesJar=None, \
					addOrReplaceReadGroupsJava=None, addOrReplaceReadGroupsJar=None,\
					alignment_method_name='bwa-short-read', alignment_format='bam',\
					namespace='workflow', version='1.0', transferOutput=False,\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None, \
					no_of_aln_threads=3, stampy=None, skipDoneAlignment=False, \
					alignmentPerLibrary=False, outputDirPrefix="", **keywords):
		"""
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
		
		
		for individual_sequence  in isqLs:
			if individual_sequence is not None and individual_sequence.format=='fastq':
				library2Data = individual_sequence.library2Data
				AlignmentJobAndOutputLs = []
				# get or add alignment method
				if individual_sequence.sequencer=='454':
					alignment_method = db_vervet.getAlignmentMethod("bwa-long-read")
				else:
					alignment_method = db_vervet.getAlignmentMethod(alignment_method_name)
				"""
				#2012.4.5 have all this commented out
				elif individual_sequence.sequencer=='GA':
					if individual_sequence.sequence_type=='SR':	#single-end
						alignment_method = db_vervet.getAlignmentMethod("bwa-short-read-SR")
					else:	#default is PE
						alignment_method = db_vervet.getAlignmentMethod("bwa-short-read")
				"""
				#alignment for the whole individual_sequence
				individual_alignment = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
												individual_sequence_id=individual_sequence.id,\
									path_to_original_alignment=None, sequencer=individual_sequence.sequencer, \
									sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
									ref_individual_sequence_id=refSequence.id, \
									alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
									individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
									data_dir=data_dir)
				skipIndividualAlignment = False
				if individual_alignment.path:
					alignmentAbsPath= os.path.join(data_dir, individual_alignment.path)
					#2012.3.29	check if the alignment exists or not. if it already exists, no alignment jobs.
					if skipDoneAlignment and os.path.isfile(alignmentAbsPath):
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
					if alignmentPerLibrary:
						#alignment for this library of the individual_sequence
						oneLibraryAlignmentEntry = db_vervet.getAlignment(individual_code=individual_sequence.individual.code, \
												individual_sequence_id=individual_sequence.id,\
									path_to_original_alignment=None, sequencer=individual_sequence.sequencer, \
									sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
									ref_individual_sequence_id=refSequence.id, \
									alignment_method_name=alignment_method.short_name, alignment_format=alignment_format,\
									individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
									data_dir=data_dir, individual_sequence_file_raw_id=minIsqFileRawID)
						skipLibraryAlignment = False
						if oneLibraryAlignmentEntry.path:
							alignmentAbsPath= os.path.join(data_dir, oneLibraryAlignmentEntry.path)
							#2012.3.29	check if the alignment exists or not. if it already exists, no alignment jobs.
							if skipDoneAlignment and os.path.isfile(alignmentAbsPath):
								skipLibraryAlignment = True
					else:
						skipLibraryAlignment = True
					if skipIndividualAlignment and skipLibraryAlignment:	#2012.9.19 if both skip flags are true, then yes
						continue
					for splitOrder in splitOrderLs:
						splitOrderIndex = splitOrder2Index[splitOrder]
						fileObjectLs = fileObjectPairLs[splitOrderIndex]

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
							java=java, SortSamFilesJava=SortSamFilesJava, SortSamFilesJar=SortSamFilesJar,\
							addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
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
												addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
												parentJobLs=[alignmentJob, mkdirJob], extraDependentInputLs=None, \
												extraArguments=None, job_max_memory = 2500, transferOutput=False)
							AlignmentJobAndOutputLs.append([addRGJob, addRGJob.output])
						if not skipLibraryAlignment:
							#2012.9.19 add a AddReadGroup job for the library bam file
							outputRGBAM = File(os.path.join(tmpOutputDir, "%s.isq_library_%s_RG.bam"%(fname_prefix, library)))
							addRGJob = self.addReadGroupInsertionJob(workflow=workflow, individual_alignment=oneLibraryAlignmentEntry, \
												inputBamFile=alignmentJob.output, \
												outputBamFile=outputRGBAM,\
												addOrReplaceReadGroupsJava=addOrReplaceReadGroupsJava, \
												addOrReplaceReadGroupsJar=addOrReplaceReadGroupsJar,\
												parentJobLs=[alignmentJob, mkdirJob], extraDependentInputLs=None, \
												extraArguments=None, job_max_memory = 2500, transferOutput=False)
							oneLibraryAlignmentJobAndOutputLs.append([addRGJob, addRGJob.output])
					if alignmentPerLibrary and not skipLibraryAlignment and oneLibraryAlignmentJobAndOutputLs:	#2012.9.19
						fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
														os.path.basename(oneLibraryAlignmentEntry.constructRelativePath()))[0]
						mergedBamFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_merged.bam'%(fname_prefix, library)))
						alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, \
											AlignmentJobAndOutputLs=oneLibraryAlignmentJobAndOutputLs, \
											outputBamFile=mergedBamFile, \
											samtools=samtools, java=java, \
											mergeSamFilesJava=mergeSamFilesJava, mergeSamFilesJar=mergeSamFilesJar, \
											BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
											mv=mv, namespace=namespace, version=version, \
											transferOutput=False, parentJobLs=[oneLibraryAlignmentFolderJob])
						
						finalBamFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_dupMarked.bam'%(fname_prefix, library)))
						markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[alignmentMergeJob, bamIndexJob], \
								inputBamF=alignmentMergeJob.output, \
								inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
								MarkDuplicatesJava=MarkDuplicatesJava, MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
								BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
								namespace=namespace, version=version, transferOutput=False)
						no_of_merging_jobs += 1
						
						#2012.9.19 add/copy the alignment file to db-affliated storage
						#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
						logFile = File(os.path.join(oneLibraryAlignmentFolder, '%s_%s_2db.log'%(fname_prefix, library)))
						alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
											inputFile=markDupJob.output, otherInputFileList=[markDupJob.MarkDupOutputMetricF], \
											individual_alignment_id=oneLibraryAlignmentEntry.id, \
											individual_sequence_file_raw_id=minIsqFileRawID,\
											logFile=logFile, data_dir=data_dir, \
											parentJobLs=[markDupJob, markDupBamIndexJob], \
											extraDependentInputLs=[markDupBamIndexJob.output,], \
											extraArguments=None, transferOutput=transferOutput, \
											job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
						
				if AlignmentJobAndOutputLs and not skipIndividualAlignment:
					#2012.3.29	merge alignment output only when there is something to merge!
					fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(\
													os.path.basename(individual_alignment.constructRelativePath()))[0]
					mergedBamFile = File(os.path.join(alignmentFolder, '%s_merged.bam'%(fname_prefix)))
					alignmentMergeJob, bamIndexJob = self.addAlignmentMergeJob(workflow, AlignmentJobAndOutputLs=AlignmentJobAndOutputLs, \
										outputBamFile=mergedBamFile, \
										samtools=samtools, java=java, \
										mergeSamFilesJava=mergeSamFilesJava, mergeSamFilesJar=mergeSamFilesJar, \
										BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
										mv=mv, namespace=namespace, version=version, \
										transferOutput=False, parentJobLs=[alignmentFolderJob])
					#relative path in the scratch
					finalBamFile = File(os.path.join(alignmentFolder, '%s_dupMarked.bam'%(fname_prefix)))
					
					markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[alignmentMergeJob, bamIndexJob],
							inputBamF=alignmentMergeJob.output, \
							inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
							MarkDuplicatesJava=MarkDuplicatesJava, MarkDuplicatesJar=MarkDuplicatesJar, tmpDir=tmpDir,\
							BuildBamIndexFilesJava=workflow.IndexMergedBamIndexJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar, \
							namespace=namespace, version=version, transferOutput=False)
					no_of_merging_jobs += 1
					
					#2012.9.19 add/copy the alignment file to db-affliated storage
					#add the metric file to AddAlignmentFile2DB.py as well (to be moved into db-affiliated storage)
					logFile = File(os.path.join(alignmentFolder, '%s_2db.log'%(fname_prefix)))
					alignment2DBJob = self.addAddAlignmentFile2DBJob(workflow=workflow, executable=self.AddAlignmentFile2DB, \
										inputFile=markDupJob.output, otherInputFileList=[markDupJob.MarkDupOutputMetricF],\
										individual_alignment_id=individual_alignment.id, \
										logFile=logFile, data_dir=data_dir, \
										parentJobLs=[markDupJob, markDupBamIndexJob], \
										extraDependentInputLs=[markDupBamIndexJob.output], \
										extraArguments=None, transferOutput=transferOutput, \
										job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel, commit=True)
					
		sys.stderr.write("%s alignment jobs; %s merge alignment jobs.\n"%(no_of_alignment_jobs, no_of_merging_jobs))
	
	def registerFileToWorkflow(self, filePair, workflow, data_dir=None):
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
	
	def addAlignmentJobPipe(self, workflow, filePair, data_dir=None, refFastaFile=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3, **keywords):
		"""
		2011-9-7
			version that uses shell(bash) pipes to skip files.
		"""
		#in MB, 2.5GB for one  aln, 5G for sampe (input is 4.5G gzipped fastq, versus 480K contigs (total size~3G)), 2G for sort
		# paired end
		pe_job_max_memory = 11000
		#pe_job_max_memory = 7000	#lie about the usage. on usc cluster, this is translated into min memory requirement.
		## single end
		se_job_max_memory = 9000
		#se_job_max_memory = 6500	#lie about the usage. on usc cluster, this is translated into min memory requirement.
		if len(filePair)==1:	#single end
			fileRecord = filePair[0]
			fastqF, format, sequence_type, sequencer = fileRecord[:4]
			relativePath = fastqF.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
			sortBamF = File('%s.bam'%(bam_output_fname_prefix))
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				alignmentJob = Job(namespace=namespace, name=ShortSEAlignmentByBWA.name, version=version)
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				alignmentJob = Job(namespace=namespace, name=LongSEAlignmentByBWA.name, version=version)
			
			alignmentJob.addArguments(refFastaFile, fastqF, bam_output_fname_prefix)
			alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
			alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%se_job_max_memory))
			alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%se_job_max_memory))
			#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
		
		elif len(filePair)==2:	#paired end
			fastqF1, format, sequence_type = filePair[0][:3]
			fastqF2, format, sequence_type = filePair[1][:3]
			relativePath = fastqF1.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			fname_prefix = fname_prefix[:-2]
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
			sortBamF = File('%s.bam'%(bam_output_fname_prefix))
			alignmentJob = Job(namespace=namespace, name=PEAlignmentByBWA.name, version=version)
			alignmentJob.addArguments(refFastaFile, fastqF1, fastqF2, bam_output_fname_prefix)
			alignmentJob.uses(fastqF1, transfer=True, register=True, link=Link.INPUT)
			alignmentJob.uses(fastqF2, transfer=True, register=True, link=Link.INPUT)
			
			alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%pe_job_max_memory))
			alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%pe_job_max_memory))
			#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
			
		alignmentJob.uses(sortBamF, transfer=False, register=True, link=Link.OUTPUT)
		workflow.addJob(alignmentJob)
		if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
			alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		if mkdirJob:
			workflow.depends(parent=mkdirJob, child=alignmentJob)
		return alignmentJob, sortBamF
	
	def addAlignmentJobNoPipe(self, workflow, filePair, data_dir=None, refFastaFile=None, bwa=None, additionalArguments=None, \
					samtools=None, refIndexJob=None, mkdirJob=None,\
					alignment_method=None, outputDir=None, namespace='workflow', version='1.0',\
					PEAlignmentByBWA=None, ShortSEAlignmentByBWA=None, LongSEAlignmentByBWA=None,\
					no_of_aln_threads=3, **keywords):
		"""
		2011-8-28
		"""
		aln_job_max_memory = 3000	#in MB, 2.5GB is enough for 4.5G gzipped fastq versus 480K contigs (total size~3G)
		if len(filePair)==1:	#single end
			fileRecord = filePair[0]
			fastqF, format, sequence_type, sequencer = fileRecord[:4]
			relativePath = fastqF.name
			
			if alignment_method.command=='aln' and sequencer!='454':	#short single-end read
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sai'%fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments,"-t %s"%no_of_aln_threads, \
										"-f", saiOutput, refFastaFile, \
										fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				
				sai2samJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				sai2samJob.addArguments("samse", "-f", outputFname, refFastaFile, saiOutput, fastqF)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
				job_max_memory = 2500	#in MB
				sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
				
				
				workflow.addJob(sai2samJob)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
				
			elif alignment_method.command=='bwasw' or sequencer=='454':	#long single-end read
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				outputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
				alignmentSamF = File(outputFname)
				alignmentJob.addArguments('bwasw', additionalArguments, "-t %s"%(no_of_aln_threads), "-f", alignmentSamF, refFastaFile, \
									fastqF)
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				workflow.addJob(alignmentJob)
				#fake a sai2samJob
				sai2samJob = alignmentJob
				
			if refIndexJob:
				workflow.depends(parent=refIndexJob, child=alignmentJob)
			if mkdirJob:
				workflow.depends(parent=mkdirJob, child=alignmentJob)
			
		elif len(filePair)==2:	#paired end
			### run sampe to combine two paired-end results into one sam file
			"bwa sampe -a 1000 -P hsref.fa ga1.sai ga2.sai ga1.fq ga2.fq | gzip > ga.sam.gz"
			sai2samJob = Job(namespace=namespace, name=bwa.name, version=version)
			# -P of sampe speeds things up but requires 4-5G memory for a human-size genome
			# "-a 1000" means maximum insert size is 1000.
			sai2samJob.addArguments('sampe', "-P")
			job_max_memory = 6000	#in MB
			sai2samJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
			sai2samJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
			sai2samJob.addProfile(Profile(Namespace.CONDOR, key="request_memory", value="%s"%job_max_memory))
			fileRecord, fileRecord2 = filePair[:2]
			fastqF, format, sequence_type = fileRecord[:3]
			relativePath = fastqF.name
			fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
			# remove the _1 or _2 in the end of the two paired-end filenames.
			fname_prefix = fname_prefix[:-2]
			samOutputFname = os.path.join(outputDir, '%s.sam'%fname_prefix)
			alignmentSamF = File(samOutputFname)
			sai2samJob.addArguments("-f", alignmentSamF)
			sai2samJob.addArguments(refFastaFile)
			if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
				sai2samJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
			sai2samJob.uses(alignmentSamF, transfer=False, register=True, link=Link.OUTPUT)
			workflow.addJob(sai2samJob)
			for fileRecord in filePair:
				fastqF, format, sequence_type = fileRecord[:3]
				relativePath = fastqF.name
				#relativePath, format, sequence_type = fileRecord[:3]
				tmp_fname_prefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(relativePath))[0]
				
				alignmentJob = Job(namespace=namespace, name=bwa.name, version=version)
				outputFname = os.path.join(outputDir, '%s.sai'%tmp_fname_prefix)
				saiOutput = File(outputFname)
				alignmentJob.addArguments(alignment_method.command, additionalArguments, "-t %s"%(no_of_aln_threads), \
										"-f", saiOutput, refFastaFile, fastqF)
				alignmentJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%aln_job_max_memory))
				alignmentJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%aln_job_max_memory))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="RequestCpus", value="%s"%no_of_aln_threads))
				#alignmentJob.addProfile(Profile(Namespace.CONDOR, key="request_cpus", value="%s"%no_of_aln_threads))
				alignmentJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				if refIndexJob:	#if this job is present, yes transfer it. otherwise assume it's there already
					alignmentJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
				alignmentJob.uses(saiOutput, transfer=False, register=True, link=Link.OUTPUT)
				workflow.addJob(alignmentJob)
				
				if refIndexJob:
					workflow.depends(parent=refIndexJob, child=alignmentJob)
				if mkdirJob:
					workflow.depends(parent=mkdirJob, child=alignmentJob)
				
				sai2samJob.addArguments(saiOutput)
				sai2samJob.uses(fastqF, transfer=True, register=True, link=Link.INPUT)
				sai2samJob.uses(saiOutput, transfer=True, register=True, link=Link.INPUT)
				workflow.depends(parent=alignmentJob, child=sai2samJob)
				
			
			#add a pair of fastq files to sampe in the end
			for fileRecord in filePair:
				fastqF = fileRecord[0]
				sai2samJob.addArguments(fastqF)
		
		
		## convert sam into bam and remove unmapped reads (or queries)
		sam_convert_job = Job(namespace=namespace, name=samtools.name, version=version)
		bamOutputF = File(os.path.join(outputDir, "%s.bam"%(fname_prefix)))
		sam_convert_job.addArguments('view',  '-F', '4', '-bSh', '-o', bamOutputF, alignmentSamF)
		sam_convert_job.uses(alignmentSamF, transfer=True, register=True, link=Link.INPUT)
		sam_convert_job.uses(bamOutputF, transfer=False, register=True, link=Link.OUTPUT)
		workflow.addJob(sam_convert_job)
		workflow.depends(parent=sai2samJob, child=sam_convert_job)
		
		"""
		# 2010-2-4
			sort it so that it could be used for merge
		"""
		bam_output_fname_prefix = '%s.sorted'%(os.path.join(outputDir, fname_prefix))
		sam_sort_job = Job(namespace=namespace, name=samtools.name, version=version)
		sam_sort_job.addArguments('sort', '-m', '2000000000', bamOutputF, bam_output_fname_prefix)	#maxMemory is down to 2G
		sortBamF = File('%s.bam'%(bam_output_fname_prefix))
		sam_sort_job.uses(bamOutputF, transfer=True, register=True, link=Link.INPUT)
		sam_sort_job.uses(sortBamF, transfer=False, register=True, link=Link.OUTPUT)
		job_max_memory = 2000	#in MB
		sam_sort_job.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
		sam_sort_job.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
		workflow.addJob(sam_sort_job)
		workflow.depends(parent=sam_convert_job, child=sam_sort_job)
		
		return sam_sort_job, sortBamF

	
	@classmethod
	def addLocalRealignmentJob(cls, workflow, samtoolsCalmd=None, inputBamF=None, \
					refFastaFList=None, outputBamF=None, \
					parentJob=None, \
					BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None, \
					namespace='workflow', version='1.0', transferOutput=True,\
					**keywords):
		"""
		2011-11-20
			create a sub-workflow to do local re-alignment for one inputBAM
			for each 2million bp interval
				split a small bam out of input bam (add 1kb on both ends of interval)
				index the small bam
				run RealignerTargetCreator
				run IndelRealigner
				(extract the exact interval out of bam if 1kb is added to both ends of input interval)
			merge all intervals back into one bam
			index the merged bam
			
		"""

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
		
		
		#individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.ind_seq_id_ls, data_dir=self.local_data_dir)
		isqLs = db_vervet.getISQDBEntryLsForAlignment(self.ind_seq_id_ls, data_dir=self.data_dir, \
												filtered=None, ignoreEmptyReadFile=self.ignoreEmptyReadFile)
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		
		#2011-11-16 new way of registering reference fasta file. but still dont' want to trasnfer 7Gb of data
		refFastaFname = os.path.join(self.data_dir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaFile = refFastaFList[0]

		if self.needRefIndexJob:
			if self.alignment_method_name.find('bwa')>=0:
				refIndexJob = self.addBWAReferenceIndexJob(workflow, refFastaFList=refFastaFList, \
												refSequenceBaseCount=refSequence.base_count, bwa=workflow.bwa)
			elif self.alignment_method_name.find('stampy')>=0:
				refIndexJob = self.addStampyGenomeIndexHashJob(workflow, executable=workflow.stampy, refFastaFList=refFastaFList, \
						parentJobLs=[], job_max_memory=3000, job_max_walltime = 200, \
						extraDependentInputLs=[], \
						transferOutput=True)
			else:
				refIndexJob = None
		else:
			refIndexJob = None
		
		self.addAllAlignmentJobs(db_vervet, individualSequenceID2FilePairLs=None, \
					isqLs = isqLs,\
					data_dir=self.data_dir,\
					refSequence=refSequence, refFastaFList=refFastaFList, refIndexJob=refIndexJob,
					workflow=workflow, bwa=workflow.bwa, additionalArguments=self.additionalArguments, \
					samtools=workflow.samtools, \
					mkdirWrap=workflow.mkdirWrap, mv=workflow.cp, \
					java=workflow.java, \
					mergeSamFilesJava=workflow.mergeSamFilesJava, mergeSamFilesJar=workflow.mergeSamFilesJar, \
					MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=self.tmpDir,\
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					SortSamFilesJava=workflow.SortSamFilesJava, SortSamFilesJar=workflow.SortSamFilesJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar,\
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
