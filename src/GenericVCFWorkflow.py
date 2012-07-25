#!/usr/bin/env python
"""
Examples:
	#2012.5.11 convert alignment read group (sample id) into UCLAID
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/ 
		-o workflow/SampleIDInUCLAID_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml 
		-u yh -y4 -l hcondor -j hcondor  -z localhost
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
	
	# 2012.5.10 run on hoffman2 condor, minMAC=1 (-n 1), minMAF=0.1 (-M 0.1), maxSNPMissingRate=0 (-N 0)   (turn on checkEmptyVCFByReading, -E)
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/
		-o workflow/SubsetTo36RNASamplesAndPlink_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml
		-i ~/script/vervet/data/RNADevelopment_eQTL/36monkeys.phenotypes.txt
		-w ~/script/vervet/data/RNADevelopment_eQTL/36monkeys.inAlignmentReadGroup.tsv
		-n1 -M 0.1 -N 0 -y3 -E
		-l hcondor -j hcondor  -u yh -z localhost
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
	# 2012.7.16 convert a folder of VCF files into plink
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719/trioCaller_vcftoolsFilter/ 
		-o workflow/ToPlinkFilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719.xml
		-y 2 -x 200 -E
		-l condorpool -j condorpool
		-u yh -z uclaOffice  -C 4
	
	# 2012.7.25 calculate haplotype distance & majority call support stats
	%s -I AlignmentToTrioCall_VRC_FilteredSeq.2012.7.21T0248_VCFWithReplicates/
		-o workflow/GetReplicateHaplotypeStat_TrioCall_VRC_FilteredSeq.2012.7.21T0248_VCFWithReplicats.xml
		-y 5 -E -l condorpool -j condorpool -u yh -z uclaOffice  -C 1 -a 524
	
Description:
	#2012.5.9
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.VCFFile import VCFFile

class GenericVCFWorkflow(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('individualUCLAIDFname', 0, ): [None, 'i', 1, 'a file containing individual ucla_id in each row. one column with header UCLAID. ', ],\
						('vcfSampleIDFname', 0, ): [None, 'w', 1, 'a file containing the sample ID (a composite ID including ucla_id) each row. \
						if not present, infer it from individualUCLAIDFname + first VCF file header.', ],\
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', '', 1, 'the 2nd input folder that contains vcf or vcf.gz files.', ],\
						('maxContigID', 0, int): [None, 'x', 1, 'if contig ID is beyond this number, it will not be included. If None or 0, no restriction.', ],\
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \
							1: subset VCF based on input file containing sample IDs;\
							2: convert to plink format; \
							3: subset + convert-2-plink. MAC & MAF & maxSNPMissingRate applied in the convert-to-plink step.\
							4: ConvertAlignmentReadGroup2UCLAIDInVCF jobs.\
							5: addMergeVCFReplicateHaplotypesJobs to get haplotype distance & majority call support stats.\
							?: Combine VCF files from two input folder, chr by chr. (not done yet. check CheckTwoVCFOverlapPipeline.py for howto)', ],\
						("minMAC", 0, int): [2, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
						("minMAF", 0, float): [None, 'M', 1, 'minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRate", 0, float): [0, 'N', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		self.inputDir = os.path.abspath(self.inputDir)
	
	def addVCF2PlinkJobs(self, workflow, inputData=None, db_vervet=None, minMAC=None, minMAF=None,\
						maxSNPMissingRate=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix=""):
		"""
		2012.7.19 
			add a modifyTPEDJob that modify 2nd column (snp-id) of tped output from default 0 to chr_pos.
			added a GzipSubworkflow in the end to gzip the final merged tped file
			all previous intermediate files are not transferred.
		2012.5.9
		"""
		sys.stderr.write("Adding VCF2plink jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sVCF2Plink"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		#each contig in each trio gets a summary.
		mergedTPEDFile = File(os.path.join(topOutputDir, 'merged.tped'))
		#each input has no header
		tpedFileMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=mergedTPEDFile, transferOutput=False, parentJobLs=[topOutputDirJob], \
							extraArguments='-n')
		no_of_jobs += 1
		returnData = PassingData()
		returnData.jobDataLs = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputF = jobData.vcfFile
			if maxContigID:
				contig_id = self.getContigIDFromFname(inputF.name)
				try:
					contig_id = int(contig_id)
					if contig_id>maxContigID:	#skip the small contigs
						continue
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
			else:
				transferOneContigPlinkOutput = False
			vcf2plinkJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=inputF, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[topOutputDirJob]+jobData.jobLs, \
						snpMisMatchStatFile=None, \
						minMAC=minMAC, minMAF=minMAF, \
						maxSNPMissingRate=maxSNPMissingRate,\
						extraDependentInputLs=[jobData.tbi_F], outputFormat='--plink-tped', transferOutput=transferOneContigPlinkOutput)
			#2012.7.20 modify the TPED 2nd column, to become chr_pos (rather than 0)
			outputF = os.path.join(topOutputDir, '%s.tped.gz'%(commonPrefix))
			modifyTPEDJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.modifyTPED, \
						inputF=vcf2plinkJob.outputF, outputF=outputF, \
						parentJobLs=[vcf2plinkJob], transferOutput=transferOutput, job_max_memory=200,\
						extraArguments=None, extraDependentInputLs=[])
			
			#add output to some reduce job
			self.addInputToStatMergeJob(workflow, statMergeJob=tpedFileMergeJob, \
								inputF=modifyTPEDJob.output, \
								parentJobLs=[modifyTPEDJob])
			no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		#include the tfam (outputList[1]) into the fileList
		returnData.jobDataLs.append(PassingData(jobLs=[tpedFileMergeJob], file=mergedTPEDFile, \
											fileList=[mergedTPEDFile, vcf2plinkJob.outputList[1]]))
		#2012.7.21 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
						outputDirPrefix="")
		return newReturnData
	
	def addVCFSubsetJobs(self, workflow, inputData=None, db_vervet=None, sampleIDFile=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix=""):
		"""
		2012.5.9
		"""
		sys.stderr.write("Adding vcf-subset jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sVCFSubset"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			if maxContigID:
				contig_id = self.getContigIDFromFname(inputF.name)
				try:
					contig_id = int(contig_id)
					if contig_id>maxContigID:	#skip the small contigs
						continue
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputVCF = File(os.path.join(topOutputDir, '%s.subset.vcf'%(commonPrefix)))
			vcfSubsetJob = self.addVCFSubsetJob(workflow, executable=workflow.vcfSubset, vcfSubsetPath=workflow.vcfSubsetPath, \
						sampleIDFile=sampleIDFile,\
						inputVCF=inputF, outputF=outputVCF, \
						parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
						extraArguments=None, extraDependentInputLs=[])
			VCFGzipOutputF = File("%s.gz"%outputVCF.name)
			VCFGzipOutput_tbi_F = File("%s.gz.tbi"%outputVCF.name)
			bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
					parentJobLs=[vcfSubsetJob], inputF=vcfSubsetJob.output, outputF=VCFGzipOutputF, \
					transferOutput=transferOutput)
			
			
			returnData.jobDataLs.append(PassingData(jobLs=[bgzip_tabix_VCF_job], vcfFile=VCFGzipOutputF, \
									tbi_F=VCFGzipOutput_tbi_F, \
									fileList=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
			
			no_of_jobs += 2
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	

	def addSubsetAndVCF2PlinkJobs(self, workflow, inputData=None, db_vervet=None, minMAC=None, minMAF=None,\
						maxSNPMissingRate=None, sampleIDFile=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix=""):
		"""
		2012.5.9
		"""
		vcfSubsetJobData = self.addVCFSubsetJobs(workflow, inputData=inputData, db_vervet=db_vervet, sampleIDFile=sampleIDFile, \
							transferOutput=True, maxContigID=maxContigID, outputDirPrefix="")
		vcf2plinkJobData = self.addVCF2PlinkJobs(workflow, inputData=vcfSubsetJobData, db_vervet=db_vervet, \
						minMAC=minMAC, minMAF=minMAF, maxSNPMissingRate=maxSNPMissingRate, transferOutput=transferOutput,\
						maxContigID=maxContigID, outputDirPrefix="")
	
	def addAlignmentReadGroup2UCLAIDJobs(self, workflow, inputData=None, db_vervet=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix=""):
		"""
		2012.5.9
		"""
		sys.stderr.write("Adding alignment read-group -> UCLAID jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sSampleInUCLAID"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			if maxContigID:
				contig_id = self.getContigIDFromFname(inputF.name)
				try:
					contig_id = int(contig_id)
					if contig_id>maxContigID:	#skip the small contigs
						continue
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputVCF = File(os.path.join(topOutputDir, '%s.UCLAID.vcf'%(commonPrefix)))
			abstractMapperJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.ConvertAlignmentReadGroup2UCLAIDInVCF, \
					inputVCF=inputF, outputF=outputVCF, \
					parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=[])
			
			VCFGzipOutputF = File("%s.gz"%outputVCF.name)
			VCFGzipOutput_tbi_F = File("%s.gz.tbi"%outputVCF.name)
			bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
					parentJobLs=[abstractMapperJob], inputF=abstractMapperJob.output, outputF=VCFGzipOutputF, \
					transferOutput=transferOutput)
			
			returnData.jobDataLs.append(PassingData(jobLs=[bgzip_tabix_VCF_job], vcfFile=VCFGzipOutputF, \
									tbi_F=VCFGzipOutput_tbi_F, \
									fileList=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
			
			no_of_jobs += 2
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	
	def addSplitNamVCFJobs(self, workflow, inputData=None, db_vervet=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix=""):
		"""
		2012.5.11
			not functional. don't know what to do the fact that SplitNamVCFIntoMultipleSingleChrVCF outputs into a folder
				multiple VCF files (one per chromosome)
		"""
		sys.stderr.write("Adding split Nam VCF-file jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sSampleInUCLAID"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			if maxContigID:
				contig_id = self.getContigIDFromFname(inputF.name)
				try:
					contig_id = int(contig_id)
					if contig_id>maxContigID:	#skip the small contigs
						continue
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputVCF = File(os.path.join(topOutputDir, '%s.vcf'%(commonPrefix)))
			abstractMapperJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.SplitNamVCFIntoMultipleSingleChrVCF, \
					inputVCF=inputF, outputF=outputVCF, \
					parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=False, job_max_memory=200,\
					extraArguments=None, extraDependentInputLs=[])
			
			VCFGzipOutputF = File("%s.gz"%outputVCF.name)
			VCFGzipOutput_tbi_F = File("%s.gz.tbi"%outputVCF.name)
			bgzip_tabix_VCF_job = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
					parentJobLs=[abstractMapperJob], inputF=abstractMapperJob.output, outputF=VCFGzipOutputF, \
					transferOutput=transferOutput)
			
			returnData.jobDataLs.append(PassingData(jobLs=[bgzip_tabix_VCF_job], vcfFile=VCFGzipOutputF, \
									tbi_F=VCFGzipOutput_tbi_F, \
									fileList=[VCFGzipOutputF, VCFGzipOutput_tbi_F]))
			
			no_of_jobs += 2
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	
	def addMergeVCFReplicateHaplotypesJobs(self, workflow, inputData=None, db_vervet=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix="",replicateIndividualTag='copy', refFastaFList=None ):
		"""
		2012.7.25
			input vcf is output of TrioCaller with replicates.
			this workflow outputs extra debug statistics
				1. replicate haplotype distance to the consensus haplotype
				2. majority support for the consensus haplotype
		"""
		sys.stderr.write("Adding MergeVCFReplicateHaplotype jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sMergeVCFReplicateHaplotypeStat"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		
		haplotypeDistanceMergeFile = File(os.path.join(topOutputDir, 'haplotypeDistanceMerge.tsv'))
		haplotypeDistanceMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=haplotypeDistanceMergeFile, transferOutput=False, parentJobLs=[topOutputDirJob])
		majoritySupportMergeFile = File(os.path.join(topOutputDir, 'majoritySupportMerge.tsv'))
		majoritySupportMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=majoritySupportMergeFile, transferOutput=False, parentJobLs=[topOutputDirJob])
		no_of_jobs += 2
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputVCF = File(os.path.join(topOutputDir, '%s.vcf'%(commonPrefix)))
			debugHaplotypeDistanceFile = File(os.path.join(topOutputDir, '%s.haplotypeDistance.tsv'%(commonPrefix)))
			debugMajoritySupportFile = File(os.path.join(topOutputDir, '%s.majoritySupport.tsv'%(commonPrefix)))
			#2012.4.2
			fileSize = utils.getFileOrFolderSize(yh_pegasus.getAbsPathOutOfFile(inputF))
			memoryRequest = 45000
			memoryRequest = min(42000, max(4000, int(33000*(fileSize/950452059.0))) )
				#extrapolates (33,000Mb memory for a ungzipped VCF file with size=950,452,059)
				#upper bound is 42g. lower bound is 4g.
			#mergeReplicateOutputF = File(os.path.join(trioCallerOutputDirJob.folder, '%s.noReplicate.vcf'%vcfBaseFname))
			#noOfAlignments= len(alignmentDataLs)
			#entireLength = stopPos - startPos + 1	#could be very small for shorter reference contigs
			#memoryRequest = min(42000, max(4000, int(20000*(noOfAlignments/323.0)*(entireLength/2600000.0))) )
				#extrapolates (20000Mb memory for a 323-sample + 2.6Mbase reference length/26K loci)
				#upper bound is 42g. lower bound is 4g.
			mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(workflow, \
								executable=workflow.MergeVCFReplicateHaplotypesJava,\
								genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
								inputF=inputF, outputF=outputVCF, \
								replicateIndividualTag=replicateIndividualTag, \
								refFastaFList=refFastaFList, \
								debugHaplotypeDistanceFile=debugHaplotypeDistanceFile, \
								debugMajoritySupportFile=debugMajoritySupportFile,\
								parentJobLs=[topOutputDirJob]+jobData.jobLs, \
								extraDependentInputLs=[], transferOutput=False, \
								extraArguments=None, job_max_memory=memoryRequest)
			
			#add output to some reduce job
			self.addInputToStatMergeJob(workflow, statMergeJob=haplotypeDistanceMergeJob, \
								inputF=mergeVCFReplicateColumnsJob.outputLs[1] , \
								parentJobLs=[mergeVCFReplicateColumnsJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=majoritySupportMergeJob, \
								inputF=mergeVCFReplicateColumnsJob.outputLs[2] , \
								parentJobLs=[mergeVCFReplicateColumnsJob])
			no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		
		returnData.jobDataLs.append(PassingData(jobLs=[haplotypeDistanceMergeJob], file=haplotypeDistanceMergeFile, \
											fileList=[haplotypeDistanceMergeFile]))
		returnData.jobDataLs.append(PassingData(jobLs=[majoritySupportMergeJob], file=majoritySupportMergeFile, \
											fileList=[majoritySupportMergeFile]))
		#2012.7.21 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
						outputDirPrefix="")
		return newReturnData
	
	def generateVCFSampleIDFilenameFromIndividualUCLAIDFname(self, db_vervet=None, individualUCLAIDFname=None, \
													vcfSampleIDFname=None, oneSampleVCFFname=None):
		"""
		2012.5.9
			
		"""
		sys.stderr.write("Generating vcfSampleIDFname %s from individualUCLAIDFname %s ..."%(vcfSampleIDFname, individualUCLAIDFname))
		reader = csv.reader(open(individualUCLAIDFname), delimiter=figureOutDelimiter(individualUCLAIDFname))
		header = reader.next()
		colName2Index = getColName2IndexFromHeader(header)
		UCLAID_col_index = colName2Index.get('UCLAID')
		individualUCLAIDSet = set()
		for row in reader:
			individualUCLAID=row[UCLAID_col_index].strip()
			individualUCLAIDSet.add(individualUCLAID)
		sys.stderr.write(" %s uclaIDs. "%(len(individualUCLAIDSet)))
		del reader
		
		writer = csv.writer(open(vcfSampleIDFname, 'w'), delimiter='\t')
		from pymodule.VCFFile import VCFFile
		vcfFile = VCFFile(inputFname=oneSampleVCFFname, minDepth=0)
		no_of_samples = 0
		for sample_id in vcfFile.sample_id_ls:
			individual_code = db_vervet.parseAlignmentReadGroup(sample_id).individual_code
			if individual_code in individualUCLAIDSet:
				no_of_samples += 1
				writer.writerow([sample_id])
		del writer, vcfFile
		sys.stderr.write("%s vcf samples selected.\n"%(no_of_samples))
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableList = []
		
		ConvertAlignmentReadGroup2UCLAIDInVCF = Executable(namespace=namespace, name="ConvertAlignmentReadGroup2UCLAIDInVCF", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		ConvertAlignmentReadGroup2UCLAIDInVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ConvertAlignmentReadGroup2UCLAIDInVCF.py"), \
														site_handler))
		executableList.append(ConvertAlignmentReadGroup2UCLAIDInVCF)
	
		SplitNamVCFIntoMultipleSingleChrVCF = Executable(namespace=namespace, name="SplitNamVCFIntoMultipleSingleChrVCF", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		SplitNamVCFIntoMultipleSingleChrVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/SplitNamVCFIntoMultipleSingleChrVCF.py"), \
														site_handler))
		executableList.append(SplitNamVCFIntoMultipleSingleChrVCF)
	
		ModifyTPED = Executable(namespace=namespace, name="ModifyTPED", \
							version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		ModifyTPED.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/ModifyTPED.py"), \
							site_handler))
		executableList.append(ModifyTPED)
		
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			workflow.addExecutable(executable)
			setattr(workflow, executable.name, executable)
	
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
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
		
		if self.individualUCLAIDFname and os.path.isfile(self.individualUCLAIDFname):
			self.generateVCFSampleIDFilenameFromIndividualUCLAIDFname(db_vervet=db_vervet, individualUCLAIDFname=self.individualUCLAIDFname, \
												vcfSampleIDFname=self.vcfSampleIDFname,\
												oneSampleVCFFname=inputData.jobDataLs[0].vcfFile.abspath)
			sampleIDFile = self.registerOneInputFile(workflow, self.vcfSampleIDFname)
		elif self.vcfSampleIDFname and os.path.isfile(self.vcfSampleIDFname):
			sampleIDFile = self.registerOneInputFile(workflow, self.vcfSampleIDFname)
		else:
			sampleIDFile = None
		
		if self.run_type==1:
			if sampleIDFile is None:
				sys.stderr.write("sampleIDFile is None.\n")
				sys.exit(0)
			self.addVCFSubsetJobs(workflow, inputData=inputData, db_vervet=db_vervet, sampleIDFile=sampleIDFile, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="")
		elif self.run_type==2:
			self.addVCF2PlinkJobs(workflow, inputData=inputData, db_vervet=db_vervet, minMAC=self.minMAC, minMAF=self.minMAF,\
						maxSNPMissingRate=self.maxSNPMissingRate, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="")
		elif self.run_type==3:
			if sampleIDFile is None:
				sys.stderr.write("sampleIDFile is None.\n")
				sys.exit(0)
			self.addSubsetAndVCF2PlinkJobs(workflow, inputData=inputData, db_vervet=db_vervet, minMAC=self.minMAC, \
							minMAF=self.minMAF,\
							maxSNPMissingRate=self.maxSNPMissingRate, sampleIDFile=sampleIDFile, transferOutput=True,\
							maxContigID=self.maxContigID, outputDirPrefix="")
		elif self.run_type==4:
			self.addAlignmentReadGroup2UCLAIDJobs(workflow, inputData=inputData, db_vervet=db_vervet, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="")
		elif self.run_type==5:
			refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
			refFastaFname = os.path.join(self.dataDir, refSequence.path)
			refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
								input_site_handler=self.input_site_handler,\
								checkAffiliateFileExistence=True)
			self.addMergeVCFReplicateHaplotypesJobs(workflow, inputData=inputData, db_vervet=db_vervet, transferOutput=True,\
						maxContigID=None, outputDirPrefix="",replicateIndividualTag='copy', refFastaFList=refFastaFList )
		else:
			sys.stderr.write("run_type %s not supported.\n"%(self.run_type))
			sys.exit(0)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


if __name__ == '__main__':
	main_class = GenericVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
