#!/usr/bin/env python
"""
Examples:
	#2012.1.10 downsample two trio to 1X,1X,4X (father,mother,child)
	%s  -u yh -a 524 -s 2 -z localhost -o DownsampleAln558To1_649To1_618To4_557To1_615To1_626To4_5sampling_top5Contigs.xml 
		-j hcondor -l hcondor -N 5 -F pedigreeVRCMerlin.2012.1.10.txt  -O 1 -C 30 -S 447 -e /u/home/eeskin/polyacti/ 
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-n 558:1,649:1,618:4,557:1,615:1,626:4 -p secret -f 5
	
	# 2012.1.10 down-sample a  VRC (-S 447) trio to 1X coverage, top 5 contigs, 20 samplings (-f 20)
	# "-O 1" controls clustering of calling programs.
	%s -u yh -a 524 -s 2 -z localhost -o DownsampleAln558To1_618To1_649To1_20sampling_top5Contigs.xml -j hcondor -l hcondor
		-N 5 -F pedigreeVRCMerlin.2012.1.9.txt -O 1 -C 30 -S 447 -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-f 20 -n 558:1,618:1,649:1
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "-O 3" controls clustering of calling programs.
	# "-C 30" controls clustering for other programs.
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
import VervetDB, csv, random
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from AlignmentToTrioCallPipeline import AlignmentToTrioCallPipeline
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class DownsampleAlignmentToTrioCallWorkflow(AlignmentToTrioCallPipeline):
	__doc__ = __doc__
	option_default_dict = AlignmentToTrioCallPipeline.option_default_dict.copy()
	option_default_dict.update({
						("probToSample", 1, float): [0.1, 'T', 1, 'probability for a read in a bam to be included in down-sampled bam, overridden by alnId2targetDepth'],\
						("no_of_sampling", 1, int): [50, 'f', 1, 'how many samplings to run'],\
						("alnId2targetDepth", 1, ): [None, 'n', 1, 'a coma-separated list, in the form of alignment_id:targetDepth. 620:1,632:1,648:1'],\
						('minDepth', 1, float): [1, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AlignmentToTrioCallPipeline.__init__(self, **keywords)
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

	
	def registerCustomExecutables(self, workflow):
		"""
		2011-11-28
		"""
		AlignmentToTrioCallPipeline.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		downsampleSamJava = Executable(namespace=namespace, name="DownsampleSamJava", version=version, os=operatingSystem, \
									arch=architecture, installed=True)
		downsampleSamJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		#splitReadFileJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(downsampleSamJava)
		workflow.downsampleSamJava = downsampleSamJava
		
	def registerJars(self, workflow, ):
		"""
		2012.1.6
			some custom jars
		"""
		AlignmentToTrioCallPipeline.registerJars(self, workflow)
		
		site_handler = self.site_handler
		
		abs_path = os.path.join(self.picard_path, 'DownsampleSam.jar')
		downsampleSamJar = File(abs_path)	#using abs_path avoids add this jar to every job as Link.INPUT
		downsampleSamJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(downsampleSamJar)
		workflow.downsampleSamJar = downsampleSamJar
	
	def addDownsampleSamJob(self, workflow, downsampleSamJava=None, downsampleSamJar=None, \
							inputF=None, probToSample=0.1, outputF=None, \
							parentJobLs=[], job_max_memory=3000, extraDependentInputLs=[], \
							transferOutput=False, **keywords):
		"""
		2012.1.6
		"""
		javaMemRequirement = "-Xms128m -Xmx%sm -XX:-UseGCOverheadLimit"%job_max_memory
		job = Job(namespace=workflow.namespace, name=downsampleSamJava.name, version=workflow.version)
		job.addArguments(javaMemRequirement, "-jar", downsampleSamJar, "VALIDATION_STRINGENCY=LENIENT")
		#2012.1.10 Changing the RANDOM_SEED argument is important as otherwise it is identical series of random integers in java program
		#	given the same default RANDOM_SEED (=1). another option is to set it 'null'.
		job.addArguments("I=", inputF,  "P=%s"%(probToSample), "O=", outputF, "RANDOM_SEED=%s"%(int(random.random()*10000)))
		#
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory, max_walltime=600)	#10 hour walltime
		workflow.addJob(job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		for parentJob in parentJobLs:
			if parentJob:
				workflow.depends(parent=parentJob, child=job)
		return job
	
	def addDownsampleJobToSelectedAlignment(self, workflow, alignmentDataLs=[], alnId2targetDepth={}, sampleNo=1):
		"""
		2012.1.9
		"""
		sys.stderr.write("Adding downsample jobs to %s selected alignments out of %s total alignments ..."%(len(alnId2targetDepth), len(alignmentDataLs)))
		returnData = []
		
		downsampleDir = "downsampleBam_%s"%(sampleNo)
		downsampleDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=downsampleDir)
		no_of_jobs = 1
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
					outputBamF = File(os.path.join(downsampleDir, os.path.basename(alignment.path)))
					if alignment.median_depth:
						downsamplerMaxMemory= int(alignment.median_depth*800)	#in Mb
					else:
						downsamplerMaxMemory = defaultDownsamplerMaxMemory
					downsampleJob = self.addDownsampleSamJob(workflow, workflow.downsampleSamJava, workflow.downsampleSamJar, inputF=bamF, \
									probToSample=probToSample, outputF=outputBamF, parentJobLs=parentJobLs+[downsampleDirJob], \
									job_max_memory=downsamplerMaxMemory, extraDependentInputLs=[baiF], transferOutput=False)
					index_sam_job = self.addBAMIndexJob(workflow, BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
						inputBamF=outputBamF, parentJobLs=[downsampleJob], stageOutFinalOutput=False, javaMaxMemory=2500)
					newAlignmentData = PassingData(alignment=alignment)
					#don't modify the old alignmentData as it will affect the original alignmentDataLs, which should remain same across different samplings
					newAlignmentData.jobLs = [downsampleJob, index_sam_job]	#downsampleJob has to be included otherwise its output (bamF) will be wiped out after index_sam_job is done
					newAlignmentData.bamF = index_sam_job.bamFile
					newAlignmentData.baiF = index_sam_job.baiFile
					no_of_jobs += 2
				else:
					newAlignmentData = alignmentData
			else:
				newAlignmentData = alignmentData
			returnData.append(newAlignmentData)
		sys.stderr.write("%s jobs added.\n"%(no_of_jobs))
		return returnData
	
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
		
		from CalculateTrioInconsistencyPipeline import CalculateTrioInconsistencyPipeline
		calculateTrioInconsistencyPipeline_ins = CalculateTrioInconsistencyPipeline(drivername=self.drivername, hostname=self.hostname, dbname=self.dbname, \
							schema=self.schema, db_user=self.db_user, db_passwd=self.db_passwd, ref_ind_seq_id=self.ref_ind_seq_id, \
							samtools_path=self.samtools_path, picard_path=self.picard_path, gatk_path=self.gatk_path,\
							vervetSrcPath=self.vervetSrcPath, home_path=self.home_path, tabixPath=self.tabixPath, javaPath=self.javaPath,\
							dataDir=self.dataDir, localDataDir=self.localDataDir, site_handler=self.site_handler,\
							input_site_handler=self.input_site_handler, clusters_size=self.clusters_size,
							outputFname=self.outputFname, checkEmptyVCFByReading=False,\
							debug=self.debug, report=self.report, inputDir="random", maxContigID=self.topNumberOfContigs+10,\
							minDepth=self.minDepth)
							#checkEmptyVCFByReading is not used in addJobs()
		calculateTrioInconsistencyPipeline_ins.registerCustomExecutables(workflow)
		
		origAlignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
		
		#reduce the trio consistency by the same trio across different samplings 
		allTrioInconsistencyFile = File('trio_inconsistency_avg_all_samples.tsv')
		allTrioInconsistencyJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByAverageColumnsWithSameKey, \
						outputF=allTrioInconsistencyFile, extraArguments='-k 0 -v 1,2,3', parentJobLs=[], \
						extraDependentInputLs=[], transferOutput=True)
		
		trioLs = calculateTrioInconsistencyPipeline_ins.getDuoTrioFromAlignmentLs(db_vervet, alignmentLs)
		
		for i in xrange(self.no_of_sampling):
			alignmentDataLs = self.addDownsampleJobToSelectedAlignment(workflow, alignmentDataLs=origAlignmentDataLs, \
											alnId2targetDepth=self.alnId2targetDepth, sampleNo=i+1)
			
			genotypeCallJobData = self.addGenotypeCallJobs(workflow, alignmentDataLs, refName2size, samtools=workflow.samtools, \
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
						intervalSize=2000000, intervalOverlapSize=100000, site_type=self.site_type, dataDir=self.dataDir,\
						outputDirPrefix="downsampleBam_%s/"%(i+1))
			
			
			trioInconsistencJobData = calculateTrioInconsistencyPipeline_ins.addJobs(workflow, inputData=genotypeCallJobData, trioLs=trioLs, db_vervet=db_vervet,\
											addTrioContigSpecificPlotJobs=False,\
											outputDirPrefix="downsampleBam_%s/"%(i+1))
			#add trio inconsistency summary output to reduction job
			self.addInputToStatMergeJob(workflow, statMergeJob=allTrioInconsistencyJob, \
								inputF=trioInconsistencJobData.trioInconsistencySummaryJob.output, \
								parentJobLs=[trioInconsistencJobData.trioInconsistencySummaryJob])
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = DownsampleAlignmentToTrioCallWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
