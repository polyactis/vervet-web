#!/usr/bin/env python
"""
Examples:
	
Description:
	2012.9.17 abstract workflow that takes both alignment and VCF as input. re-factored out of AlignmentToCallPipeline.py 
"""
import sys, os, math
__doc__ = __doc__%()


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from AbstractVervetWorkflow import AbstractVervetWorkflow
from Pegasus.DAX3 import *
from vervet.src import VervetDB


class AbstractAlignmentAndVCFWorkflow(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	#change the short option to nothing to avoid conflict
	option_default_dict[('inputDir', 0, )] = ['', 'L', 1, 'input folder that contains vcf or vcf.gz files', ]

	#option_default_dict.pop(('inputDir', 0, ))
	commonAlignmentWorkflowOptionDict = {
						('ind_seq_id_ls', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('ind_aln_id_ls', 0, ): ['', 'I', 1, 'a comma/dash-separated list of IndividualAlignment.id. This overrides ind_seq_id_ls.', ],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'To filter alignments. None: whatever; 0: unfiltered sequences, 1: filtered sequences: 2: ...'],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						("site_id_ls", 0, ): ["", 'S', 1, 'comma/dash-separated list of site IDs. individuals must come from these sites.'],\
						("country_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						("tax_id_ls", 0, ): ["", '', 1, 'comma/dash-separated list of country IDs. individuals must come from these countries.'],\
						('defaultSampleAlignmentDepth', 1, int): [10, '', 1, "when database doesn't have median_depth info for one alignment, use this number instead.", ],\
						('individual_sequence_file_raw_id_type', 1, int): [1, '', 1, "1: only all-library-fused libraries,\n\
		2: only library-specific alignments,\n\
		3: both all-library-fused and library-specific alignments", ],\
						}
	partitionWorkflowOptionDict= {
						("needFastaIndexJob", 0, int): [0, 'A', 0, 'need to add a reference index job by samtools?'],\
						("needFastaDictJob", 0, int): [0, 'B', 0, 'need to add a reference dict job by picard CreateSequenceDictionary.jar?'],\
						("selectedRegionFname", 0, ): ["", 'R', 1, 'the file is in bed format, tab-delimited, chr start stop.\
		used to restrict SAMtools/GATK to only make calls at this region. \
		start and stop are 0-based. i.e. start=0, stop=100 means bases from 0-99.\
		This overrides the contig/chromosome selection approach defined by --contigMaxRankBySize and --contigMinRankBySize. \
		This file would be split into maxNoOfRegionsPerJob lines.'],\
						('maxNoOfRegionsPerJob', 1, int): [5000, 'K', 1, 'Given selectedRegionFname, this dictates the maximum number of regions each job would handle,\
		The actual number could be lower because the regions are first grouped into chromosomes. If one chromosome has <maxNoOfRegionsPerJob, then that job handles less.', ],\
						}
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		#if self.option_default_dict.get('ligateVcfPerlPath'):
		#	self.pathToInsertHomePathList.append('ligateVcfPerlPath')
		AbstractVervetWorkflow.__init__(self, **keywords)
		
		listArgumentName_data_type_ls = [('ind_seq_id_ls', int), ("ind_aln_id_ls", int), \
								("site_id_ls", int), ('country_id_ls', int), ('tax_id_ls', int)]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
	
	def getAlignments(self, ref_ind_seq_id=None, ind_seq_id_ls=[], ind_aln_id_ls=[], alignment_method_id=2, \
					dataDir=None, sequence_type=None):
		"""
		2012.4.13
			moved to VervetDB.py
		2012.4.5
			select alignment using AND between all input arguments
		2011-11-27
			add argument sequence_type
		2011-9-16
			order each alignment by id. It is important because this is the order that gatk&samtools take input bams.
			#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
			# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		2011-9-13
			add argument dataDir, to filter out alignments that don't exist on file storage
		2011-8-31
			add argument aln_method_id
		2011-7-12
		
		"""
		return self.db_vervet.getAlignments(ref_ind_seq_id=ref_ind_seq_id, ind_seq_id_ls=ind_seq_id_ls, ind_aln_id_ls=ind_aln_id_ls, alignment_method_id=alignment_method_id, \
					dataDir=dataDir, sequence_type=sequence_type)
	
	def filterAlignments(self, alignmentLs, max_coverage=None, individual_site_id=447, sequence_filtered=0):
		"""
		2012.4.13
			moved to VervetDB.py
		2012.4.2
			add argument sequence_filtered
		2011-11-22
			447 in "individual_site_id=447" is VRC.
		"""
		return self.db_vervet.filterAlignments(alignmentLs, max_coverage=max_coverage, individual_site_id=individual_site_id, \
								sequence_filtered=sequence_filtered)
	
	def preReduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachChromosome(self, workflow=None, alignmentData=None, chromosome=None,\
				VCFFile=None, passingData=None, reduceBeforeEachAlignmentData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def map(self, workflow=None, alignmentData=None, intervalData=None,\
		VCFFile=None, passingData=None, mapEachChromosomeData=None, transferOutput=True, **keywords):
		"""
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def mapEachInterval(self, **keywords):
		"""
		2012.9.22 link to map()
		"""
		return self.map(**keywords)
	
	
	def linkMapToReduce(self, workflow=None, mapEachIntervalData=None, preReduceReturnData=None, passingData=None, transferOutput=True, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData

	def mapEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9.22
			similar to reduceBeforeEachAlignmentData() but for mapping programs that run on one alignment each.
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduceAfterEachChromosome(self, workflow=None, chromosome=None, passingData=None, transferOutput=True, \
								mapEachIntervalDataLs=None, **keywords):
		"""
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		return returnData
	
	def reduceBeforeEachAlignment(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.9 setup some reduce jobs before loop over all intervals of one alignment begins.
			these reduce jobs will collect stuff from each map() job.
			the link will be established in linkMapToReduce().
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
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
		2012.9.17
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.reduceAfterEachAlignmentDataLs = reduceAfterEachAlignmentDataLs
		return returnData
	
	def addAllJobs(self, workflow=None, inputVCFData=None, alignmentDataLs=None, chr2IntervalDataLs=None, samtools=None, \
				genomeAnalysisTKJar=None, \
				mergeSamFilesJar=None, \
				createSequenceDictionaryJava=None, createSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, \
				refFastaFList=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				dataDir=None, no_of_gatk_threads = 1, \
				outputDirPrefix="", transferOutput=True, **keywords):
		"""
		2012.7.26
		"""
				#2012.8.26 so that each recalibration will pick up the right vcf
		chr2VCFFile = {}
		for jobData in inputVCFData.jobDataLs:
			inputF = jobData.file
			chr = self.getChrFromFname(os.path.basename(inputF.name))
			chr2VCFFile[chr] = inputF
		chrIDSet = set(chr2VCFFile.keys())&set(chr2IntervalDataLs.keys())
		
		sys.stderr.write("Adding alignment+VCF-related jobs for %s chromosomes/contigs ..."%(len(chrIDSet)))
		refFastaF = refFastaFList[0]
		
		topOutputDir = "%sMap"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		if needFastaDictJob:	# the .dict file is required for GATK
			fastaDictJob = self.addRefFastaDictJob(workflow, createSequenceDictionaryJava=createSequenceDictionaryJava, \
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
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		#2012.9.22 AlignmentJobAndOutputLs is a relic.
		#	but it's similar to mapEachIntervalDataLs but designed for addAlignmentMergeJob(),
		#	so AlignmentJobAndOutputLs gets re-set for every alignment.
		# 	mapEachAlignmentDataLs is never reset.
		#	mapEachChromosomeDataLs is reset right after a new alignment is chosen.
		#	mapEachIntervalDataLs is reset right after each chromosome is chosen.
		#	all reduce dataLs never gets reset.
		passingData = PassingData(AlignmentJobAndOutputLs=[], bamFnamePrefix=None, topOutputDirJob=topOutputDirJob,\
					outputDirPrefix=outputDirPrefix, refFastaFList=refFastaFList, \
					
					mapEachAlignmentData = None,\
					mapEachChromosomeData=None, \
					mapEachIntervalData=None,\
					reduceBeforeEachAlignmentData = None, \
					reduceAfterEachAlignmentData=None,\
					reduceAfterEachChromosomeData=None,\
					
					mapEachAlignmentDataLs = [],\
					mapEachChromosomeDataLs=[], \
					mapEachIntervalDataLs=[],\
					reduceBeforeEachAlignmentDataLs = [], \
					reduceAfterEachAlignmentDataLs=[],\
					reduceAfterEachChromosomeDataLs=[],\
					
					gzipReduceAfterEachChromosomeFolderJob=None,\
					gzipReduceBeforeEachAlignmentFolderJob = None,\
					gzipReduceAfterEachAlignmentFolderJob = None,\
					gzipPreReduceFolderJob = None,\
					gzipReduceFolderJob=None,\
					)
		preReduceReturnData = self.preReduce(workflow=workflow, passingData=passingData, transferOutput=False,\
											**keywords)
		passingData.preReduceReturnData = preReduceReturnData
		
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			parentJobLs = alignmentData.jobLs
			bamF = alignmentData.bamF
			baiF = alignmentData.baiF
			
			bamFnamePrefix = alignment.getReadGroup()
			
			passingData.AlignmentJobAndOutputLs = []
			passingData.bamFnamePrefix = bamFnamePrefix
			passingData.individual_alignment = alignment
			
			mapEachAlignmentData = self.mapEachAlignment(workflow=workflow, passingData=passingData, transferOutput=False, \
												preReduceReturnData=preReduceReturnData, **keywords)
			mapEachAlignmentDataLs.append(mapEachAlignmentData)
			passingData.mapEachAlignmentData = mapEachAlignmentData
			
			reduceBeforeEachAlignmentData = self.reduceBeforeEachAlignment(workflow=workflow, passingData=passingData, \
													preReduceReturnData=preReduceReturnData, transferOutput=False, \
													**keywords)
			passingData.reduceBeforeEachAlignmentData = reduceBeforeEachAlignmentData
			passingData.reduceBeforeEachAlignmentDataLs.append(reduceBeforeEachAlignmentData)
			
			mapEachChromosomeDataLs = passingData.mapEachChromosomeDataLs
			mapEachChromosomeDataLs = []
			
			for chr in chrIDSet:
				intervalDataLs = chr2IntervalDataLs.get(chr)
				VCFFile = chr2VCFFile.get(chr)
				if VCFFile is None:
					sys.stderr.write("WARNING: no VCFFile for chromosome %s. skipped.\n"%(chr))
					continue
				mapEachChromosomeData = self.mapEachChromosome(workflow=workflow, alignmentData=alignmentData, chromosome=chr, \
									VCFFile=VCFFile, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
									mapEachAlignmentData=mapEachAlignmentData,\
									transferOutput=False, **keywords)
				passingData.mapEachChromosomeData = mapEachChromosomeData
				mapEachChromosomeDataLs.append(mapEachChromosomeData)
				
				mapEachIntervalDataLs = passingData.mapEachIntervalDataLs
				mapEachIntervalDataLs = []
				
				for intervalData in intervalDataLs:
					if intervalData.file:
						mpileupInterval = intervalData.interval
						bcftoolsInterval = intervalData.file
					else:
						mpileupInterval = intervalData.interval
						bcftoolsInterval = intervalData.interval
					intervalFnameSignature = intervalData.intervalFnameSignature
					overlapInterval = intervalData.overlapInterval
					overlapFilenameSignature = intervalData.overlapIntervalFnameSignature
					
					mapEachIntervalData = self.mapEachInterval(workflow=workflow, alignmentData=alignmentData, intervalData=intervalData,\
										VCFFile=VCFFile, passingData=passingData, reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
										mapEachAlignmentData=mapEachAlignmentData,\
										mapEachChromosomeData=mapEachChromosomeData, transferOutput=False, **keywords)
					passingData.mapEachIntervalData = mapEachIntervalData
					mapEachIntervalDataLs.append(mapEachIntervalData)
					
					linkMapToReduceData = self.linkMapToReduce(workflow=workflow, mapEachIntervalData=mapEachIntervalData, preReduceReturnData=preReduceReturnData, \
										reduceBeforeEachAlignmentData=reduceBeforeEachAlignmentData,\
										mapEachAlignmentData=mapEachAlignmentData,\
										passingData=passingData, \
										**keywords)
				
				reduceAfterEachChromosomeData = self.reduceAfterEachChromosome(workflow=workflow, chromosome=chr, passingData=passingData, \
									mapEachIntervalDataLs=passingData.mapEachIntervalDataLs,\
									transferOutput=False, dataDir=dataDir, \
									**keywords)
				passingData.reduceAfterEachChromosomeData = reduceAfterEachChromosomeData
				passingData.reduceAfterEachChromosomeDataLs.append(reduceAfterEachChromosomeData)
				
				gzipReduceAfterEachChromosomeData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceAfterEachChromosomeData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduceAfterEachChromosome"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceAfterEachChromosomeFolderJob, report=False)
				passingData.gzipReduceAfterEachChromosomeFolderJob = gzipReduceAfterEachChromosomeData.topOutputDirJob
				
			reduceAfterEachAlignmentData = self.reduceAfterEachAlignment(workflow=workflow, \
													mapEachAlignmentData=mapEachAlignmentData,\
													mapEachChromosomeDataLs=mapEachChromosomeDataLs,\
													reduceAfterEachChromosomeDataLs=reduceAfterEachChromosomeDataLs,\
													passingData=passingData, \
													transferOutput=False, dataDir=dataDir, **keywords)
			passingData.reduceAfterEachAlignmentData = reduceAfterEachAlignmentData
			passingData.reduceAfterEachAlignmentDataLs.append(reduceAfterEachAlignmentData)
			
			gzipReduceBeforeEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceBeforeEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduceBeforeEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceBeforeEachAlignmentFolderJob, report=False)
			passingData.gzipReduceBeforeEachAlignmentFolderJob = gzipReduceBeforeEachAlignmentData.topOutputDirJob
			
			gzipReduceAfterEachAlignmentData = self.addGzipSubWorkflow(workflow=workflow, \
						inputData=reduceAfterEachAlignmentData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduceAfterEachAlignment"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceAfterEachAlignmentFolderJob, \
						report=False)
			passingData.gzipReduceAfterEachAlignmentFolderJob = gzipReduceAfterEachAlignmentData.topOutputDirJob
		reduceReturnData = self.reduce(workflow=workflow, passingData=passingData, \
							mapEachAlignmentData=mapEachAlignmentData, \
							reduceAfterEachAlignmentDataLs=passingData.reduceAfterEachAlignmentDataLs,\
							**keywords)
		passingData.reduceReturnData = reduceReturnData
		
		
		#2012.9.18 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=preReduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%spreReduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipPreReduceFolderJob, \
						report=False)
		passingData.gzipPreReduceFolderJob = newReturnData.topOutputDirJob
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=reduceReturnData, transferOutput=transferOutput,\
						outputDirPrefix="%sreduce"%(outputDirPrefix), \
						topOutputDirJob=passingData.gzipReduceFolderJob, \
						report=False)
		passingData.gzipReduceFolderJob = newReturnData.topOutputDirJob
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def registerAlignmentAndItsIndexFile(self, workflow=None, alignmentLs=None, dataDir=None, checkFileExistence=True):
		"""
		2012.9.18 copied from AlignmentToCallPipeline.py
		2012.6.12
			add argument checkFileExistence
		2012.1.9
			register the input alignments and return in a data structure usd by several other functions
		"""
		sys.stderr.write("Registering %s alignments ..."%(len(alignmentLs)))
		returnData = []
		for alignment in alignmentLs:
			inputFname = os.path.join(dataDir, alignment.path)
			input = File(alignment.path)	#relative path, induces symlinking or stage-in
			baiFilepath = '%s.bai'%(inputFname)
			if checkFileExistence and (not os.path.isfile(inputFname) or not os.path.isfile(baiFilepath)):
				continue
			input.addPFN(PFN("file://" + inputFname, workflow.input_site_handler))
			workflow.addFile(input)
			baiF = File('%s.bai'%alignment.path)
			baiF.addPFN(PFN("file://" + baiFilepath, workflow.input_site_handler))
			workflow.addFile(baiF)
			alignmentData = PassingData(alignment=alignment, jobLs = [], bamF=input, baiF=baiF)
			returnData.append(alignmentData)
		sys.stderr.write("Done.\n")
		return returnData
	
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		2011-11-28
		"""
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
	def run(self):
		"""
		2011-9-28
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
		
		chr2size = self.getTopNumberOfContigs(contigMaxRankBySize=self.contigMaxRankBySize, contigMinRankBySize=self.contigMinRankBySize)
		#chr2size = set(['Contig149'])	#temporary when testing Contig149
		#chr2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		chrLs = chr2size.keys()
		chr2IntervalDataLs = self.getChr2IntervalDataLsBySplitChrSize(chr2size=chr2size, \
													intervalSize=self.intervalSize, \
													intervalOverlapSize=self.intervalOverlapSize)
		
		alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, dataDir=self.localDataDir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type,\
										country_id_ls=self.country_id_ls, tax_id_ls=self.tax_id_ls)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
									individual_site_id_set=set(self.site_id_ls),\
									mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls))
		
		workflow = self.initiateWorkflow()
		
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
		
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
		
		self.addAllJobs(workflow=workflow, inputVCFData=inputData, alignmentDataLs=alignmentDataLs, \
					chr2IntervalDataLs=chr2IntervalDataLs, samtools=workflow.samtools, \
				genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
				mergeSamFilesJar=workflow.mergeSamFilesJar, \
				createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar,\
				mv=workflow.mv, \
				refFastaFList=refFastaFList,\
				needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
				dataDir=self.dataDir, no_of_gatk_threads = 1, transferOutput=True,\
				outputDirPrefix=self.pegasusFolderName)
		
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)