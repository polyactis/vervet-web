#!/usr/bin/env python
"""
Examples:
	%s
	
	dirPrefix=AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_single_sample_condor_20111105T0143/556_
	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o dags/CheckTwoVCF/4HighCovVRC_isq_15_18_vs_524_top804Contigs_gatk_vs_samtools_overlap_stat.xml -z uclaOffice -u yh -k genome
		-C 100
	
	#2012.8.3 compare two VCF on hcondor, do per-sample checking
	%s -i AlignmentToCall_ISQ643_646_vs_524_Contig731.2012.8.2T1530/samtools/
		-I AlignmentToCall_ISQ643_646_vs_524_method7Contig731Sites.2012.8.2T1610/samtools/
		-l hcondor -j hcondor -o dags/CheckTwoVCF/CheckTwoVCF_ISQ643_646_vs_524_MultiSample_vs_Method7ROICalling_Contig731.xml
		-z localhost -u yh -k genome -C 1
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--perSampleMatchFraction
	
Description:
	2011-11-7 pegasus workflow that compares overlap between two vcf files (mapper/CheckTwoVCFOverlap.py),
			calculate mismatch rate, pi statistics based on the intersection
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

from Pegasus.DAX3 import File, Executable, PFN
from pymodule import ProcessOptions, getListOutOfStr, PassingData, NextGenSeq, utils
from pymodule.pegasus import yh_pegasus
from vervet.src import VervetDB, AbstractVervetWorkflow

parentClass = AbstractVervetWorkflow
class CheckTwoVCFOverlapPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('vcf1Dir', 1, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 1, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('perSampleMatchFraction', 0, ): [0, '', 0, 'whether calculating per-sample mismatch fraction or not.', ],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
	def addCheckingVCFOverlapSubWorkflow(self, workflow=None, chr2size=None, inputVCFData1=None, inputVCFData2=None, \
					registerReferenceData=None, outputDirPrefix="", **keywords):
		"""
		2013.09.05
		"""
		if workflow is None:
			workflow = self
		if registerReferenceData is None:
			registerReferenceData = self.registerReferenceData
		
		sys.stderr.write("Adding Check-VCF overlap jobs between %s (patch 1) and %s (patch 2), job count=%s..."%
						(len(inputVCFData1.jobDataLs), len(inputVCFData2.jobDataLs), self.no_of_jobs))
		returnData = PassingData()
		
		mapDirJob = self.addMkDirJob(outputDir="%sMap"%(outputDirPrefix))
		reduceDirJob = self.addMkDirJob(outputDir="%sReduce"%(outputDirPrefix))
		plotOutputDirJob = self.addMkDirJob(outputDir="%sPlot"%(outputDirPrefix))
		
		overlapStatF = File(os.path.join(reduceDirJob.output, 'overlapSites.perChromosome.stat.tsv.gz'))
		overlapSitesByChromosomeMergeJob=self.addStatMergeJob(statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
					outputF=overlapStatF, parentJobLs=[reduceDirJob], \
					extraDependentInputLs=None, transferOutput=True, extraArguments=None)
		
		overlapSitesMergeJob=self.addStatMergeJob(statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
					outputF=File(os.path.join(reduceDirJob.output, "overlapSites.tsv.gz")), parentJobLs=[reduceDirJob], \
					extraDependentInputLs=None, transferOutput=True, extraArguments=None)
		
		perSampleMatchFractionFile = File(os.path.join(reduceDirJob.output, 'perSampleMatchFraction.tsv.gz'))
		perSampleMatchFractionReduceJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixBySumSameKeyColsAndThenDivide, \
					outputF=perSampleMatchFractionFile, parentJobLs=[reduceDirJob], extraDependentInputLs=[], transferOutput=True, \
					extraArguments='-k 0 -v 1-2')
		returnData.perSampleMatchFractionReduceJob = perSampleMatchFractionReduceJob
		
		outputFile = File( os.path.join(plotOutputDirJob.output, 'perSampleMatchFraction_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[perSampleMatchFractionFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="no_of_matches_by_no_of_non_NA_pairs", whichColumnPlotLabel="matchFraction", \
					logY=None, logCount=True, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, perSampleMatchFractionReduceJob ], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		overlapStatSumF = File(os.path.join(reduceDirJob.output, 'overlapSites.wholeGenome.stat.tsv'))
		overlapStatSumJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
						outputF=overlapStatSumF, parentJobLs=[reduceDirJob], extraDependentInputLs=[], transferOutput=True, \
						extraArguments='-k 1000000 -v 1-25000')	#The key column (-k 1000000) doesn't exist.
						# essentially merging every rows into one 
						##25000 is a random big upper limit. 100 monkeys => 101*3 + 9 => 312 columns
						#2012.8.17 the number of columns no longer expand as the number of samples because it's split into perSampleMatchFractionFile.
		self.addInputToStatMergeJob(statMergeJob=overlapStatSumJob, inputF=overlapStatF, \
							parentJobLs=[overlapSitesByChromosomeMergeJob])
		
		vcfJobDataRBTree1 = self.constructGenomeFileRBTreeByFilenameInterval(jobDataStructure=inputVCFData1, chr2size=chr2size)
		vcfJobDataRBTree2 = self.constructGenomeFileRBTreeByFilenameInterval(jobDataStructure=inputVCFData2, chr2size=chr2size)
		
		noOfPairs=0
		for vcfJobDataNode1 in vcfJobDataRBTree1:
			chromosome = vcfJobDataNode1.key.chromosome
			chrLength = chr2size.get(chromosome)
			if chrLength is None:
				sys.stderr.write("Warning: size for chromosome %s is unknown. set it to 1000.\n"%(chromosome))
				chrLength = 1000
			jobData1 = vcfJobDataNode1.value
			
			vcfJobDataNodeListInTree2 = []
			vcfJobDataRBTree2.findNodes(key=vcfJobDataNode1.key, node_ls=vcfJobDataNodeListInTree2)
			for vcfJobDataNode2 in vcfJobDataNodeListInTree2:
				noOfPairs += 1
				jobData2 = vcfJobDataNode2.value
				
				#narrow down either VCF file based on the interval info
				overlap_start = max(vcfJobDataNode1.key.start, vcfJobDataNode2.key.start)
				overlap_stop = min(vcfJobDataNode1.key.stop, vcfJobDataNode2.key.stop)
				if overlap_start!=vcfJobDataNode1.key.start or overlap_stop!=vcfJobDataNode1.key.stop:
					fileBasenamePrefix = "%s"%(utils.getFileBasenamePrefixFromPath(jobData1.file.name))
					outputF = File(os.path.join(mapDirJob.output, "%s_%s_%s_%s.vcf"%(fileBasenamePrefix, chromosome, overlap_start, overlap_stop)))
					selectVCF1Job = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
										inputF=jobData1.file, outputF=outputF, \
										interval="%s:%s-%s"%(chromosome, overlap_start, overlap_stop),\
										refFastaFList=registerReferenceData.refFastaFList, \
										parentJobLs=[mapDirJob] + jobData1.jobLs, extraDependentInputLs=jobData1.fileLs[1:], transferOutput=False, \
										extraArguments=None, extraArgumentList=None, job_max_memory=2000, walltime=None)
					jobData1 = self.constructJobDataFromJob(selectVCF1Job)
				if overlap_start!=vcfJobDataNode2.key.start or overlap_stop!=vcfJobDataNode2.key.stop:
					fileBasenamePrefix = "%s"%(utils.getFileBasenamePrefixFromPath(jobData2.file.name))
					outputF = File(os.path.join(mapDirJob.output, "%s_%s_%s_%s.vcf"%(fileBasenamePrefix, chromosome, overlap_start, overlap_stop)))
					selectVCF2Job = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
										inputF=jobData2.file, outputF=outputF, \
										interval="%s:%s-%s"%(chromosome, overlap_start, overlap_stop),\
										refFastaFList=registerReferenceData.refFastaFList, \
										parentJobLs=[mapDirJob] + jobData2.jobLs, extraDependentInputLs=jobData2.fileLs[1:], transferOutput=False, \
										extraArguments=None, extraArgumentList=None, job_max_memory=2000, walltime=None)
					jobData2 = self.constructJobDataFromJob(selectVCF2Job)
				
				fileBasenamePrefix = "%s_vs_%s"%(utils.getFileBasenamePrefixFromPath(jobData1.file.name), 
												utils.getFileBasenamePrefixFromPath(jobData2.file.name))
				
				outputFnamePrefix = os.path.join(mapDirJob.output, fileBasenamePrefix)
				outputFile = File("%s.tsv.gz"%(outputFnamePrefix))
				perSampleConcordanceOutputFile = File("%s_perSample.tsv.gz"%(outputFnamePrefix))
				overlapSiteOutputFile = File("%s_overlapSitePos.tsv.gz"%(outputFnamePrefix))
				checkTwoVCFOverlapJob = self.addCheckTwoVCFOverlapJob(executable=workflow.CheckTwoVCFOverlapCC, \
						vcf1=jobData1.file, vcf2=jobData2.file, chromosome=chromosome, chrLength=chrLength, \
						outputFile=outputFile, perSampleConcordanceOutputFile=perSampleConcordanceOutputFile, \
						overlapSiteOutputFile=overlapSiteOutputFile,\
						parentJobLs=[mapDirJob] + jobData1.jobLs + jobData2.jobLs, \
						extraDependentInputLs=jobData1.fileLs[1:] + jobData2.fileLs[1:], \
						transferOutput=False, extraArguments=None,\
						#"--minDepth %s "%(self.minDepth),\
						job_max_memory=1000)
				
				self.addInputToStatMergeJob(statMergeJob=overlapSitesByChromosomeMergeJob, \
							inputF=checkTwoVCFOverlapJob.output, \
							parentJobLs=[checkTwoVCFOverlapJob], extraDependentInputLs=[])
				self.addInputToStatMergeJob(statMergeJob=overlapSitesMergeJob, \
							inputF=checkTwoVCFOverlapJob.overlapSitePosFile, \
							parentJobLs=[checkTwoVCFOverlapJob], extraDependentInputLs=[])
				self.addInputToStatMergeJob(statMergeJob=perSampleMatchFractionReduceJob, \
							inputF=checkTwoVCFOverlapJob.perSampleFile, \
							parentJobLs=[checkTwoVCFOverlapJob], extraDependentInputLs=[])
		sys.stderr.write("%s pairs of VCF files, %s jobs.\n"%(noOfPairs, self.no_of_jobs))
		return returnData
	
	
	addAllJobs = addCheckingVCFOverlapSubWorkflow
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		#db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, db_user=self.db_user,
		#				db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema="genome")
		#db_genome.setup(create_tables=False)
		#db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		pdata = self.setup_run()
		workflow = pdata.workflow
		
		
		inputVCFData1 = self.registerAllInputFiles(inputDir=self.vcf1Dir, input_site_handler=None, \
					checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName='%s_vcf1'%(self.pegasusFolderName),\
					maxContigID=self.maxContigID, minContigID=self.minContigID, \
					db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLociInVCF=None, includeIndelVCF=True)
		#split the VCFs
		outputDirJob = self.addMkDirJob(outputDir='%s_vcf1_split'%(self.pegasusFolderName))
		inputVCFData1 = self.addSplitVCFSubWorkflow(inputVCFData=inputVCFData1, intervalSize=self.intervalSize, chr2size=self.chr2size,
							registerReferenceData=self.registerReferenceData, outputDirJob=outputDirJob)
		inputVCFData2 = self.registerAllInputFiles(inputDir=self.vcf2Dir, input_site_handler=None, \
					checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName='%s_vcf2'%(self.pegasusFolderName),\
					maxContigID=self.maxContigID, minContigID=self.minContigID, \
					db_vervet=None, needToKnowNoOfLoci=False,
					minNoOfLociInVCF=None, includeIndelVCF=True)
		#split the VCFs
		outputDirJob = self.addMkDirJob(outputDir='%s_vcf2_split'%(self.pegasusFolderName))
		inputVCFData2 = self.addSplitVCFSubWorkflow(inputVCFData=inputVCFData2, intervalSize=self.intervalSize, chr2size=self.chr2size,
							registerReferenceData=self.registerReferenceData, outputDirJob=outputDirJob)
		self.addCheckingVCFOverlapSubWorkflow(workflow=workflow, chr2size = self.chr2size, inputVCFData1=inputVCFData1, inputVCFData2=inputVCFData2, \
							registerReferenceData=pdata.registerReferenceData,\
							outputDirPrefix="")

		"""
		vcfFileID2object_1 = self.getVCFFileID2object(self.vcf1Dir)
		vcfFileID2object_2 = self.getVCFFileID2object(self.vcf2Dir)
		sharedVCFFileIDSet = set(vcfFileID2object_1.keys())&set(vcfFileID2object_2.keys())
		sys.stderr.write("%s shared vcf files.\n"%(len(sharedVCFFileIDSet)))
		
		for vcfFileID in sharedVCFFileIDSet:
			gatkVCFAbsPath = vcfFileID2object_1.get(vcfFileID).vcfFilePath
			samtoolsVCFAbsPath = vcfFileID2object_2.get(vcfFileID).vcfFilePath
		"""
		
		
		self.end_run()
		


	
if __name__ == '__main__':
	main_class = CheckTwoVCFOverlapPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
