#!/usr/bin/env python
"""
Examples:
	#2012.5.11 convert alignment read group (sample id) into UCLAID
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/ 
		-o dags/SampleIDInUCLAID_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml 
		--db_user yh -y4 --site_handler hcondor --input_site_handler hcondor  --hostname localhost
		--home_path /u/home/eeskin/polyacti/ --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
	
	# 2012.8.10 IBD check
	# add --kinshipFname kinshipFile if you want comparison (table&figures) between IBD pi-hat and kinship 
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_14/ -o dags/PlinkIBDCheck/PlinkIBDCheck_Method14.xml --clusters_size 1 
		--needSSHDBTunnel --site_handler hcondor --input_site_handler hcondor  --db_user yh --hostname localhost
		--localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--hostname localhost  -y3 --mergeListFname ./aux/Method14_LDPrune_merge_list.2012.8.10T0441.txt
		#--kinshipFname ~/NetworkData/vervet/Kinx2Apr2012.txt
	
	# 2012.8.9 LD-prune a folder of VCF files into plink, need the db tunnel (--needSSHDBTunnel) for output pedigree in tfam
	# "--minContigID 90 --maxContigID 100" are used to restrict contig IDs between 90 and 100.
	# LDPruneMinR2=0.3 (--LDPruneMinR2), LDPruneWindowSize=500 (--LDPruneWindowSize), 
	# LDPruneWindowShiftSize=100 (--LDPruneWindowShiftSize)
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719/trioCaller_vcftoolsFilter/ 
		-o dags/ToPlinkFilterVCF/ToPlinkFilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719.xml
		-y 2  --checkEmptyVCFByReading
		--site_handler condorpool --input_site_handler condorpool
		--db_user yh --hostname uclaOffice  --clusters_size 4 --needSSHDBTunnel --LDPruneMinR2 0.3 --LDPruneWindowSize 500
		--LDPruneWindowShiftSize 100
		#--minContigID 90 --maxContigID 100
	
	# 2012.8.10 Plink Mendel Error
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_14/ -o dags/PlinkMendelError/PlinkMendelError_Method14.xml
		--checkEmptyVCFByReading --clusters_size 4  --needSSHDBTunnel --site_handler hcondor --input_site_handler hcondor
		--db_user yh --hostname localhost
		--localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-y1
	
	# 2012.8.13 sex check using the top 195 contigs (--maxContigID 195) (Contig 83, 149,193 are sex chromosomes)
	# no clustering (--clusters_size 1)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-o dags/PlinkSexCheck/PlinkSexCheck_Method14_W100Z10R0.4_maxContigID195.xml
		--LDPruneWindowSize 100 --LDPruneWindowShiftSize 10 --LDPruneMinR2 0.4 --clusters_size 1  --needSSHDBTunnel
		--site_handler hcondor --input_site_handler hcondor  --db_user yh --hostname localhost
		--localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--hostname localhost  -y4 --mergeListFname ./aux/Method14_LDPrune_merge_list.2012.8.13T1702.txt --maxContigID 195
	
	#2013.2.1 mark mendel error-calls missing (-y5)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_38/
		-o dags/MarkMendelErrorCallMissing/MarkMendelErrorCallMissing_method38.xml
		-y5 --clusters_size 1  --needSSHDBTunnel --site_handler hcondor --input_site_handler hcondor
		--db_user yh --hostname localhost
		--localDataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --dataDir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--mergeListFname ./aux/Method38_MarkMendelErrorCallMissing_merge_list.20130201.txt --minContigID 90 --maxContigID 100
	
Description:
	2012.8.14 a plink workflow on folder of VCF files.
		1: mendel errors,\
		2: LD pruning, \
		3: IBD check,\
		4: sex check,
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
from pymodule import GenomeDB
from pymodule import AbstractVCFWorkflow
from vervet.src.pegasus.GenericVCFWorkflow import GenericVCFWorkflow

class PlinkOnVCFWorkflow(GenericVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVCFWorkflow.option_default_dict)
	plinkWorkflowOptionDict= {
						('LDPruneMinR2', 0, float): [0.4, 'R', 1, 'minimum r2 for LD pruning', ],\
						('locusSamplingRate', 0, float): [0.0001, 'c', 1, 'how many loci to sample (unused)', ],\
						('mergeListFname', 0, ): [None, 'g', 1, 'the file to contain the merge-list for plink, required for run_type>1', ],\
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \n\
	1: mendel errors,\n\
	2: LD pruning, \n\
	3: IBD check,\n\
	4: sex check,\n\
	5: mark genotype calls associated with mendel inconsistency as missing', ],\
						('LDPruneWindowSize', 1, int): [50, 'W', 1, ' window size (in the number of SNPs, not bp) for plink LD pruning'],\
						('LDPruneWindowShiftSize', 1, int): [20, '', 1, 'adjacent window shift (in the number of SNPs), not bp '],\
						('kinshipFname', 0, ): ["", 's', 1, 'the kinship file from Sue (in turn from Solar, based on pedigree )\
		if given, y3 (IBD check) workflow will compare IBD result with kinship'],\
							}
	option_default_dict.update(plinkWorkflowOptionDict)
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		GenericVCFWorkflow.__init__(self, **keywords)

	def addPlinkMendelErrorJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix="", returnMode=2, **keywords):
		"""
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink mendel error jobs for %s  files ... "%(len(inputData.jobDataLs)))
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		topOutputDir = "%sMendelError"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		
		mendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_mendel.tsv'))
		#each input has no header
		mendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=mendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob])
		returnData.jobDataLs.append(PassingData(jobLs=[mendelMergeJob], file=mendelMergeJob.output, \
											fileList=mendelMergeJob.outputLs))
		
		imendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_imendel.tsv'))
		#each input has no header
		imendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=imendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob], \
							extraArguments='--keyColumnLs 1 --valueColumnLs 2')
		
		returnData.jobDataLs.append(PassingData(jobLs=[imendelMergeJob], file=imendelMergeJob.output, \
											fileList=imendelMergeJob.outputLs))
		
		outputFile = File( os.path.join(plotOutputDir, 'individualMendelErrorHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(executable=workflow.DrawHistogram, inputFileList=[imendelMergeFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="N", whichColumnPlotLabel="log_NoOfMendelErrors", \
					logY=1, logCount=True, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, imendelMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		
		fmendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_fmendel.tsv'))
		fmendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=fmendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob])
		returnData.jobDataLs.append(PassingData(jobLs=[fmendelMergeJob], file=fmendelMergeJob.output, \
											fileList=fmendelMergeJob.outputLs))
		
		lmendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_lmendel.tsv'))
		lmendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=lmendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob], \
							extraArguments='')
		returnData.jobDataLs.append(PassingData(jobLs=[lmendelMergeJob], file=lmendelMergeJob.output, \
											fileList=lmendelMergeJob.outputLs))
		
		outputFile = File( os.path.join(plotOutputDir, 'locusMendelErrorHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		# samplingRate=1 because plink doesn't output zero-mendel-error sites.
		self.addDrawHistogramJob(executable=workflow.DrawHistogram, inputFileList=[lmendelMergeFile], \
							outputFile=outputFile, \
					whichColumn=2, whichColumnHeader=None, whichColumnPlotLabel="NoOfMendelErrors", \
					logY=0, logCount=True, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, lmendelMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		
		#2013.1.8
		meanMedianModePerLocusMendelErrorFile = File(os.path.join(mergedOutputDirJob.output, "meanMedianModePerLocusMendelError.tsv"))
		meanMedianModePerLocusMendelErrorJob = self.addCalculateDepthMeanMedianModeJob(\
							executable=workflow.CalculateMedianMeanOfInputColumn, \
							inputFile=lmendelMergeJob.output, outputFile=meanMedianModePerLocusMendelErrorFile, alignmentID=0, \
							fractionToSample=1, \
							noOfLinesInHeader=1, whichColumn=2, maxNumberOfSamplings=1E8, inputStatName="MendelError",\
							parentJobLs=[mergedOutputDirJob, lmendelMergeJob], job_max_memory = 500, extraArguments=None, \
							transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[meanMedianModePerLocusMendelErrorJob], \
									file=meanMedianModePerLocusMendelErrorJob.output, \
									fileList=meanMedianModePerLocusMendelErrorJob.outputLs))
		
		#2013.1.8 next add a job that calculates the number of nuclear families (plink definition) and divide above meanMendelError with that
		# => meanMendelErrorRate per site
		# the number of nuclear families => the number of unique parents pairs in the tfam file (both parents are not missing)
		# ...
		
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
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
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
			else:
				transferOneContigPlinkOutput = False
			i += 1
			mendelFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkMendelJob = self.addPlinkJob(executable=self.plink, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=mendelFnamePrefix, outputOption='--out',\
					calculateMendelError=True, \
					extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ jobData.jobLs)
			
			#add output to some reduce job
			self.addInputToStatMergeJob(workflow, statMergeJob=mendelMergeJob, \
								inputF=plinkMendelJob.mendelFile, \
								parentJobLs=[plinkMendelJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=imendelMergeJob, \
								inputF=plinkMendelJob.imendelFile, \
								parentJobLs=[plinkMendelJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=fmendelMergeJob, \
								inputF=plinkMendelJob.fmendelFile, \
								parentJobLs=[plinkMendelJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=lmendelMergeJob, \
								inputF=plinkMendelJob.lmendelFile, \
								parentJobLs=[plinkMendelJob])
		
		##2012.7.21 gzip the final output
		returnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
					outputDirPrefix=outputDirPrefix)
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def writePlinkMergeListFile(self, outputFname=None, extractedFilenameTupleList=[]):
		"""
		2012.9.12
			exit if outputFname exists already
		"""
		if os.path.isfile(outputFname):
			sys.stderr.write("Error: plink merge file %s already exists. Rename it to something else or delete it (if its workflow has finished.)\n"%\
							(outputFname))
			sys.exit(2)
		outf = open(outputFname, 'w')
		for data_tuple in extractedFilenameTupleList:
			outf.write('%s %s %s\n'%(data_tuple[0], data_tuple[1], data_tuple[2]))
		del outf
	
	def markMendelErrorLociMissingSubWorkflow(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, tfamJob=None, outputDirPrefix="", **keywords):
		"""
		2013.1.29
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding jobs that mark mendel-error loci missing for %s  files ... "%(len(inputData.jobDataLs)))
		
		topOutputDir = "%sMarkMendelLocusMissing"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
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
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
				transferOneContigModifyTPEDOutput = True
			else:
				transferOneContigPlinkOutput = False
				transferOneContigModifyTPEDOutput = transferOutput
				
			mendelFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkMendelJob = self.addPlinkJob(executable=self.plink, \
									tpedFile=inputJob.output, tfamFile=tfamJob.tfamFile,\
					outputFnamePrefix=mendelFnamePrefix, outputOption='--out',\
					calculateMendelError=True, \
					extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob, tfamJob]+ jobData.jobLs)
			
			
			#2012.7.20 modify the TPED 2nd column, to become chr_pos (rather than 0)
			modifyTPEDFnamePrefix = os.path.join(topOutputDir, '%s_markMissing'%(commonPrefix))
			outputF = File('%s.tped'%(modifyTPEDFnamePrefix))
			modifyTPEDJobExtraArguments = "--run_type 4"
			extraArgumentList = [" --mendelErrorFname", plinkMendelJob.mendelFile, \
								"--tfamFname ",tfamJob.tfamFile ]
			modifyTPEDJob = self.addAbstractMapperLikeJob(workflow, executable=workflow.ModifyTPED, \
						inputF=inputJob.output, outputF=outputF, \
						parentJobLs=jobData.jobLs + [plinkMendelJob, tfamJob], \
						transferOutput=transferOneContigModifyTPEDOutput, job_max_memory=200,\
						extraArguments=modifyTPEDJobExtraArguments, extraArgumentList=extraArgumentList,\
						extraDependentInputLs=[plinkMendelJob.mendelFile, tfamJob.tfamFile])
			returnData.jobDataLs.append(PassingData(jobLs=[modifyTPEDJob], file=modifyTPEDJob.output, \
											fileList=modifyTPEDJob.outputLs))
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def addPlinkBinaryConversionJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, tfamJob=None, outputDirPrefix="", returnMode=1, \
						mergeListFile=None, **keywords):
		"""
		2013.1.29	returnMode:
			1: individual contig/chr binary file/job
			2: final merged binary file/job 
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink binary converting jobs for %s  files ... "%(len(inputData.jobDataLs)))
		
		topOutputDir = "%sPlinkBinary"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		mergedOutputDir = "%sPlinkMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		
		returnData = PassingData()
		returnData.jobDataLs = []
		plinkTPED2BEDJobList = []
		extractedFilenameTupleList = []
		plinkMergeExtraDependentInputList = []
		
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
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
			
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
			else:
				transferOneContigPlinkOutput = False
			
			#convert single plink tped file into binary bed file
			bedFnamePrefix = os.path.join(topOutputDirJob.output, '%s_bed'%(commonPrefix))
			convertSingleTPED2BEDJob = self.addPlinkJob(executable=self.plinkConvert, inputFileList=[], 
								tpedFile=inputJob.output, tfamFile=tfamJob.tfamFile,\
								inputFnamePrefix=None, inputOption=None, \
						outputFnamePrefix=bedFnamePrefix, outputOption='--out',\
						makeBED=True, \
						extraDependentInputLs=None, transferOutput=transferOutput, \
						extraArguments=None, job_max_memory=2000,\
						parentJobLs = [inputJob, topOutputDirJob, tfamJob])
			plinkTPED2BEDJobList.append(convertSingleTPED2BEDJob)
			
			#returnData.jobDataLs.append(PassingData(jobLs=[convertSingleTPED2BEDJob], file=convertSingleTPED2BEDJob.bedFile, \
			#							fileList=convertSingleTPED2BEDJob.outputLs))
		
			if i>0:	#the 1st one is to be directly added to the plink merge job 
				extractedFilenameTupleList.append([convertSingleTPED2BEDJob.bedFile.name, convertSingleTPED2BEDJob.bimFile.name, convertSingleTPED2BEDJob.famFile.name])
				plinkMergeExtraDependentInputList.extend([convertSingleTPED2BEDJob.bedFile, convertSingleTPED2BEDJob.bimFile, convertSingleTPED2BEDJob.famFile])
			i += 1
			if returnMode==1:
				returnData.jobDataLs.append(PassingData(jobLs=[convertSingleTPED2BEDJob], file=convertSingleTPED2BEDJob.bedFile, \
											fileList=convertSingleTPED2BEDJob.outputLs))
		
		#fill up the mergeListFile
		self.writePlinkMergeListFile(outputFname=yh_pegasus.getAbsPathOutOfFile(mergeListFile),\
									extractedFilenameTupleList=extractedFilenameTupleList)
		
		plinkMergeFnamePrefix = os.path.join(mergedOutputDir, 'LDPrunedMerged')
		firstPlinkTPED2BEDJob = plinkTPED2BEDJobList[0]
		plinkMergeJob = self.addPlinkJob(executable=self.plinkMerge, \
						bedFile=firstPlinkTPED2BEDJob.bedFile, famFile=firstPlinkTPED2BEDJob.famFile, \
						bimFile=firstPlinkTPED2BEDJob.bimFile, \
						outputFnamePrefix=plinkMergeFnamePrefix, outputOption='--out',\
						mergeListFile=mergeListFile, makeBED=True, \
						extraDependentInputLs=plinkMergeExtraDependentInputList, transferOutput=transferOutput, \
						extraArguments=None, job_max_memory=2000,\
						parentJobLs =[mergedOutputDirJob] + plinkTPED2BEDJobList)
		
		if returnMode==2:
			returnData.jobDataLs.append(PassingData(jobLs=[plinkMergeJob], file=plinkMergeJob.bedFile, \
											fileList=plinkMergeJob.outputLs))
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	
	def addPlinkLDPruneJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, LDPruneMinR2=0.1, outputDirPrefix="", returnMode=1, \
						LDPruneWindowSize=100, LDPruneWindowShiftSize=5, mergeListFile=None, **keywords):
		"""
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink LD pruning jobs for %s  files ... "%(len(inputData.jobDataLs)))
		
		topOutputDir = "%sLDPrune"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		plinkExtractJobList = []
		extractedFilenameTupleList = []
		plinkMergeExtraDependentInputList = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
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
			LDPruneFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkLDPruneJob = self.addPlinkJob(executable=self.plinkLDPrune, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=LDPruneFnamePrefix, outputOption='--out',\
					LDPruneWindowSize=LDPruneWindowSize, LDPruneWindowShiftSize=LDPruneWindowShiftSize, \
					LDPruneByPairwiseR2=True, LDPruneMinR2=LDPruneMinR2,\
					extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ jobData.jobLs)
			
			extractFnamePrefix = os.path.join(topOutputDir, '%s_extract'%(commonPrefix))
			plinkExtractJob = self.addPlinkJob(executable=self.plinkExtract, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=extractFnamePrefix, outputOption='--out',\
					extractSNPFile = plinkLDPruneJob.prune_inFile, makeBED=True, \
					extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob, plinkLDPruneJob ]+ jobData.jobLs)
			plinkExtractJobList.append(plinkExtractJob)
			if i>0:	#the 1st one is to be directly added to the plink merge job 
				extractedFilenameTupleList.append([plinkExtractJob.bedFile.name, plinkExtractJob.bimFile.name, plinkExtractJob.famFile.name])
				plinkMergeExtraDependentInputList.extend([plinkExtractJob.bedFile, plinkExtractJob.bimFile, plinkExtractJob.famFile])
			i += 1
			if returnMode==2:
				returnData.jobDataLs.append(PassingData(jobLs=[plinkExtractJob], file=plinkExtractJob.bedFile, \
											fileList=plinkExtractJob.outputLs))
		
		#fill up the mergeListFile
		self.writePlinkMergeListFile(outputFname=yh_pegasus.getAbsPathOutOfFile(mergeListFile),\
									extractedFilenameTupleList=extractedFilenameTupleList)
		
		plinkMergeFnamePrefix = os.path.join(mergedOutputDir, 'LDPrunedMerged')
		firstPlinkExtractJob = plinkExtractJobList[0]
		plinkMergeJob = self.addPlinkJob(executable=self.plinkMerge, \
						bedFile=firstPlinkExtractJob.bedFile, famFile=firstPlinkExtractJob.famFile, \
						bimFile=firstPlinkExtractJob.bimFile, \
						outputFnamePrefix=plinkMergeFnamePrefix, outputOption='--out',\
						mergeListFile=mergeListFile, makeBED=True, \
						extraDependentInputLs=plinkMergeExtraDependentInputList, transferOutput=transferOutput, \
						extraArguments=None, job_max_memory=2000,\
						parentJobLs =[mergedOutputDirJob] + plinkExtractJobList)
		
		
		if returnMode==1 or returnMode==3:
			returnData.jobDataLs.append(PassingData(jobLs=[plinkMergeJob], file=plinkMergeJob.bedFile, \
											fileList=plinkMergeJob.outputLs))
			
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def addDetectWrongLabelByCompKinshipVsIBDJob(self, workflow=None, executable=None, inputFile=None, \
							plinkIBDCheckOutputFile=None, outputFile=None, outputFnamePrefix=None, \
							iterativeAlgorithm=False, minAbsDeltaForOutlier=0, \
							parentJobLs=None, \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.8.28
			inputFile is the kinship file from Sue.
		
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = []
		key2ObjectForJob = {}
		extraOutputLs = []
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		if outputFnamePrefix is None and outputFile:
			outputFnamePrefix = os.path.splitext(outputFile.name)[0]
		if outputFnamePrefix:
			extraArgumentList.append("--outputFnamePrefix %s"%(outputFnamePrefix))
			suffixAndNameTupleList.extend([['_PCAOnAbsKinshipIBDDelta.tsv', 'PCAOnAbsKinshipIBDDelta'], ['_orderByPC1.tsv', 'orderByPC1']])
			
		if iterativeAlgorithm:
			extraArgumentList.append("-A")
		if plinkIBDCheckOutputFile:
			extraArgumentList.extend(["--plinkIBDCheckOutputFname", plinkIBDCheckOutputFile])
			extraDependentInputLs.append(plinkIBDCheckOutputFile)
		if minAbsDeltaForOutlier:
			extraArgumentList.extend(['--minAbsDeltaForOutlier %s'%(minAbsDeltaForOutlier)])
			if outputFnamePrefix:
				suffixAndNameTupleList.extend([['_minAbsDelta%s.tsv'%(minAbsDeltaForOutlier), 'minAbsDeltaOutlier'], \
									['_monkeyFrequencyInMinAbsDelta%sPairs.tsv'%(minAbsDeltaForOutlier), 'monkeyFrequencyInMinAbsDeltaOutliers']])
			
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if outputFnamePrefix:
			self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		return job
	
	
	def addPlotPedigreeKinshipVsGeneticIBDJob(self, workflow=None, executable=None, inputFile=None, \
							plinkIBDCheckOutputFile=None, outputFile=None, outputFnamePrefix=None, \
							whichColumnPlotLabel=None, xColumnPlotLabel=None, \
							minNoOfTotal=100,\
							figureDPI=300, samplingRate=1, doPairwiseLabelCheck=False, \
							parentJobLs=None, \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=True, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2012.8.21
			inputFile is the kinship file from Sue.
		
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		extraArgumentList = ['--minNoOfTotal %s'%(minNoOfTotal), \
							'--figureDPI %s'%(figureDPI), '--samplingRate %s'%(samplingRate)]
		key2ObjectForJob = {}
		extraOutputLs = []
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		if outputFnamePrefix is None and outputFile:
			outputFnamePrefix = os.path.splitext(outputFile.name)[0]
		if outputFnamePrefix:
			extraArgumentList.append("--outputFnamePrefix %s"%(outputFnamePrefix))
			suffixAndNameTupleList.extend([['_scatter.png', 'scatter'], ['_table.tsv', 'table']])		#the file with kinship vs. PI_HAT
			if doPairwiseLabelCheck:
				extraArgumentList.append("--doPairwiseLabelCheck ")
				suffixAndNameTupleList.extend([['_SumAbsDelta.tsv', 'sumAbsDelta'], \
											['_pairwiseCorOfKinshipIBDDelta.tsv', 'pairwiseCorOfKinshipIBDDelta']])
			
		if whichColumnPlotLabel:
			extraArgumentList.append("--whichColumnPlotLabel %s"%(whichColumnPlotLabel))
		if xColumnPlotLabel:
			extraArgumentList.append("--xColumnPlotLabel %s"%(xColumnPlotLabel))
		if plinkIBDCheckOutputFile:
			extraArgumentList.extend(["--plinkIBDCheckOutputFname", plinkIBDCheckOutputFile])
			extraDependentInputLs.append(plinkIBDCheckOutputFile)
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if outputFnamePrefix:
			self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, \
				sshDBTunnel=sshDBTunnel, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		return job
	
	def addPlinkIBDCheckJobs(self, workflow=None, inputData=None, kinshipFile=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix="", returnMode=1, useAlleleFrequencyFromNonFounders=False,\
						fontFile=None, no_of_individuals=None):
		"""
		2012.9.12 add no_of_individuals, for plotting kinship-IBD matrix
		2012.8.28
			add argument useAlleleFrequencyFromNonFounders
				fontFile, for DrawMatrix job
		2012.8.21
			add plinkIBDCheckOutputFile
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink IBD checking jobs for %s  files ... "%(len(inputData.jobDataLs)))
		
		
		topOutputDir = "%sIBDCheck"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
			else:
				transferOneContigPlinkOutput = False
			if useAlleleFrequencyFromNonFounders:	#2012.8.28
				freqFnamePrefix = os.path.join(topOutputDir, '%s_nonFounders_frequency'%(commonPrefix))
				plinkFrqJob = self.addPlinkJob(executable=self.plinkIBD, \
										bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
						outputFnamePrefix=freqFnamePrefix, outputOption='--out',\
						estimateAlleFrequency=True, \
						extraDependentInputLs=None, transferOutput=transferOutput, \
						extraArguments="--nonfounders", job_max_memory=2000,\
						parentJobLs =[topOutputDirJob]+ jobData.jobLs)
				frqFile = plinkFrqJob.frqFile
				plinkIBDCheckParentJobLs = [plinkFrqJob] + jobData.jobLs
			else:
				frqFile = None
				plinkIBDCheckParentJobLs = jobData.jobLs
				
			IBDCheckFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkIBDCheckJob = self.addPlinkJob(executable=self.plinkIBD, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=IBDCheckFnamePrefix, outputOption='--out',\
					estimatePairwiseGenomeWideIBD=True, estimatePairwiseGenomeWideIBDFreqFile=frqFile, \
					extraDependentInputLs=None, transferOutput=transferOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ plinkIBDCheckParentJobLs)
			#2012.8.15 transform it to tab-delimited matrix
			outputFile = File(os.path.join(topOutputDir, "%s_ibdCheck.tsv"%(commonPrefix)))
			toTsvMatrixJob = self.addAbstractMatrixFileWalkerJob(executable=self.AbstractMatrixFileWalker, \
					inputFile=plinkIBDCheckJob.genomeFile, outputFile=outputFile, \
					outputFnamePrefix=None, whichColumn=None, whichColumnHeader=None, \
					logY=None, valueForNonPositiveYValue=-1, \
					minNoOfTotal=10,\
					samplingRate=1, \
					parentJobLs=[topOutputDirJob, plinkIBDCheckJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput, job_max_memory=2000)
			
			outputFile = File( os.path.join(plotOutputDir, '%s_ibdCheck_PI_HAT_hist.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFile=toTsvMatrixJob.output, \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="PI_HAT", whichColumnPlotLabel="proportionIBD", \
						logY=None, logCount=True, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						figureDPI=100, samplingRate=1, \
						parentJobLs=[plotOutputDirJob, toTsvMatrixJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			
			returnData.jobDataLs.append(PassingData(jobLs=[toTsvMatrixJob], file=toTsvMatrixJob.output, \
											fileList=toTsvMatrixJob.outputLs))
			if kinshipFile:
				wrongLabelOutputFnamePrefix = 'wrongLabelChiSq_greedy'
				outputFile = File(os.path.join(topOutputDir, '%s.tsv'%(wrongLabelOutputFnamePrefix)))
				detectWrongLabelJob = self.addDetectWrongLabelByCompKinshipVsIBDJob(executable=self.DetectWrongLabelByCompKinshipVsIBD, \
																	inputFile=kinshipFile, \
							plinkIBDCheckOutputFile=toTsvMatrixJob.output, outputFile=outputFile, outputFnamePrefix=None, \
							iterativeAlgorithm=True, minAbsDeltaForOutlier=0.2, \
							parentJobLs=[toTsvMatrixJob, topOutputDirJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
				
				if fontFile:
					#maxBlockSize is the max number of individuals to be put in one block of image. \
					#(if matrix is too large, DrawMatrix will split them).
					if no_of_individuals:
						maxBlockSize = no_of_individuals + 2
					else:
						maxBlockSize = 2000
					outputFile = File(os.path.join(plotOutputDir, '%s_orderByPC1.png'%(wrongLabelOutputFnamePrefix)))
					drawMatrixJob = self.addGenericJob(executable=self.DrawMatrix, inputFile=detectWrongLabelJob.orderByPC1File, \
													outputFile=outputFile, \
							parentJobLs=[detectWrongLabelJob, plotOutputDirJob], extraDependentInputLs=[fontFile], \
							extraOutputLs=None,\
							transferOutput=transferOutput, \
							extraArgumentList=['--font_size 10', '--font_path', fontFile, '--no_grid --blockColUnit %s --blockRowUnit %s'%(maxBlockSize, maxBlockSize)], job_max_memory=4000)
				
				outputFile = File(os.path.join(plotOutputDir, '%s_chiSqPvalue_hist.png'%(wrongLabelOutputFnamePrefix)))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, \
							inputFile=detectWrongLabelJob.output, \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="chiSqPvalue", whichColumnPlotLabel="log10chiSqPvalue", \
							logY=1, logCount=True, valueForNonPositiveYValue=-1,\
							minNoOfTotal=10,\
							figureDPI=150, samplingRate=1, \
							parentJobLs=[plotOutputDirJob, detectWrongLabelJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFnamePrefix = os.path.join(plotOutputDir, '%s_KinshipVsIBDPI'%(commonPrefix))
				plotPedigreeKinshipVsIBDJob = self.addPlotPedigreeKinshipVsGeneticIBDJob(workflow=workflow, executable=self.PlotPedigreeKinshipVsGeneticIBD, \
								inputFile=kinshipFile, \
								plinkIBDCheckOutputFile=toTsvMatrixJob.output, outputFile=None, outputFnamePrefix=outputFnamePrefix, \
								whichColumnPlotLabel="IBDCheckPIHat", xColumnPlotLabel="Kinship", \
								minNoOfTotal=10,\
								figureDPI=200, samplingRate=1, doPairwiseLabelCheck=True, \
								parentJobLs=[toTsvMatrixJob, plotOutputDirJob], \
								extraDependentInputLs=None, \
								extraArguments=None, transferOutput=transferOutput, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
				
				outputFile = File('%s_sumAbsDelta_hist.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, \
							inputFile=plotPedigreeKinshipVsIBDJob.sumAbsDeltaFile, \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="sumAbsDelta", whichColumnPlotLabel="sumAbsDelta", \
							logY=None, logCount=True, valueForNonPositiveYValue=-1,\
							minNoOfTotal=10,\
							figureDPI=150, samplingRate=1, \
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_avgAbsDelta_hist.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, \
							inputFile=plotPedigreeKinshipVsIBDJob.sumAbsDeltaFile, \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="avgAbsDelta", whichColumnPlotLabel="avgAbsDelta", \
							logY=False, logCount=True, valueForNonPositiveYValue=-1,\
							minNoOfTotal=10,\
							figureDPI=150, samplingRate=1, \
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_medianAbsDelta_hist.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, \
							inputFile=plotPedigreeKinshipVsIBDJob.sumAbsDeltaFile, \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="medianAbsDelta", whichColumnPlotLabel="medianAbsDelta", \
							logY=False, logCount=True, valueForNonPositiveYValue=-1,\
							minNoOfTotal=10,\
							figureDPI=150, samplingRate=1, \
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_sumAbsDelta_vs_medianAbsDelta.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
							inputFileList=[plotPedigreeKinshipVsIBDJob.sumAbsDeltaFile], \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="sumAbsDelta", whichColumnPlotLabel="sumAbsDelta", \
							logY=False, valueForNonPositiveYValue=50,\
							xColumnHeader="medianAbsDelta", xColumnPlotLabel="medianAbsDelta", \
							minNoOfTotal=20,\
							figureDPI=150, samplingRate=1,\
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_noOfNonMissing_vs_medianAbsDelta.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
							inputFileList=[plotPedigreeKinshipVsIBDJob.sumAbsDeltaFile], \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="medianAbsDelta", whichColumnPlotLabel="medianAbsDelta", \
							logY=False, valueForNonPositiveYValue=50,\
							xColumnHeader="noOfNonMissing", xColumnPlotLabel="noOfNonMissing", \
							minNoOfTotal=20,\
							figureDPI=150, samplingRate=1,\
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_pairwiseCorOfKinshipIBDDelta_hist.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, \
							inputFile=plotPedigreeKinshipVsIBDJob.pairwiseCorOfKinshipIBDDeltaFile, \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="corr", whichColumnPlotLabel="pairwiseCorOfKinshipIBDDelta", \
							logY=False, logCount=True, valueForNonPositiveYValue=-1,\
							minNoOfTotal=10,\
							figureDPI=150, samplingRate=1, \
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_pairwiseCorOfKinshipIBDDelta_vs_monkey1_medianAbsDelta.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
							inputFileList=[plotPedigreeKinshipVsIBDJob.pairwiseCorOfKinshipIBDDeltaFile], \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="corr", whichColumnPlotLabel="pairwiseCorOfKinshipIBDDelta", \
							logY=False, valueForNonPositiveYValue=50,\
							xColumnHeader="monkey1_medianAbsDelta", xColumnPlotLabel="monkey1_medianAbsDelta", \
							minNoOfTotal=20,\
							figureDPI=150, samplingRate=1,\
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
				
				outputFile = File('%s_pairwiseCorOfKinshipIBDDelta_vs_monkey2_medianAbsDelta.png'%(outputFnamePrefix))
				#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
				self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
							inputFileList=[plotPedigreeKinshipVsIBDJob.pairwiseCorOfKinshipIBDDeltaFile], \
							outputFile=outputFile, \
							whichColumn=None, whichColumnHeader="corr", whichColumnPlotLabel="pairwiseCorOfKinshipIBDDelta", \
							logY=False, valueForNonPositiveYValue=50,\
							xColumnHeader="monkey2_medianAbsDelta", xColumnPlotLabel="monkey2_medianAbsDelta", \
							minNoOfTotal=20,\
							figureDPI=150, samplingRate=1,\
							parentJobLs=[plotOutputDirJob, plotPedigreeKinshipVsIBDJob], \
							extraDependentInputLs=None, \
							extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnData
	
	def addPlinkSexCheckJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix="", returnMode=1,):
		"""
		2012.8.13
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink Sex checking jobs for %s  files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sSexCheck"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		no_of_jobs += 1
		
		
		returnData = PassingData()
		returnData.jobDataLs = []
		for i in xrange(len(inputData.jobDataLs)):
			jobData = inputData.jobDataLs[i]
			inputJob = jobData.jobLs[0]
			inputF = jobData.file
			inputFBaseName = os.path.basename(inputF.name)
			commonPrefix = inputFBaseName.split('.')[0]
			outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			if i ==0:	#need at least one tfam file. 
				transferOneContigPlinkOutput = True
			else:
				transferOneContigPlinkOutput = False
			outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkJob = self.addPlinkJob(executable=self.plink, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=outputFnamePrefix, outputOption='--out',\
					checkSex=True,\
					extraDependentInputLs=None, transferOutput=transferOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ jobData.jobLs)
			returnData.jobDataLs.append(PassingData(jobLs=[plinkJob], file=plinkJob.sexcheckFile, \
											fileList=plinkJob.outputLs))
			
			outputFile = File( os.path.join(plotOutputDir, '%s_inbreedCoeffByChrX_vs_sexByInput.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, inputFileList=[plinkJob.sexcheckFile], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="F", whichColumnPlotLabel="inbreedCoeffByChrX", \
						logY=False, valueForNonPositiveYValue=50,\
						xColumnHeader="PEDSEX", xColumnPlotLabel="sexByInput", \
						minNoOfTotal=2,\
						figureDPI=100, samplingRate=1,\
						parentJobLs=[plotOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			no_of_jobs +=2
			
			outputFile = File(os.path.join(topOutputDir, "%s.male.sexCheck"%(commonPrefix)))
			selectMaleDataJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
								inputFileList=[plinkJob.sexcheckFile], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="PEDSEX", \
						logY=False, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						samplingRate=1,\
						parentJobLs=[topOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments="--minWhichColumnValue 1 --maxWhichColumnValue 1", transferOutput=False,  job_max_memory=2000)
			outputFile = File( os.path.join(plotOutputDir, '%s.maleInbreedCoeffByChrX_Hist.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selectMaleDataJob.output], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="F", whichColumnPlotLabel="inbreedCoeffByChrX", \
						logY=False, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						figureDPI=100, samplingRate=1,\
						parentJobLs=[plotOutputDirJob, selectMaleDataJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			no_of_jobs +=2

			outputFile = File(os.path.join(topOutputDir, "%s.female.sexCheck"%(commonPrefix)))
			selectFemaleDataJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
								inputFileList=[plinkJob.sexcheckFile], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="PEDSEX", \
						logY=False, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						samplingRate=1,\
						parentJobLs=[topOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments="--minWhichColumnValue 2 --maxWhichColumnValue 2", transferOutput=False,  job_max_memory=2000)
			outputFile = File( os.path.join(plotOutputDir, '%s.femaleInbreedCoeffByChrX_Hist.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selectFemaleDataJob.output], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="F", whichColumnPlotLabel="inbreedCoeffByChrX", \
						logY=False, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						figureDPI=100, samplingRate=1,\
						parentJobLs=[plotOutputDirJob, selectFemaleDataJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			no_of_jobs +=2
		
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		#GenericVCFWorkflow.registerCustomExecutables(self, workflow=workflow)
		#2013.1.6 equivalent to above
		GenericVCFWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		
		PlotPedigreeKinshipVsGeneticIBD = Executable(namespace=namespace, name="PlotPedigreeKinshipVsGeneticIBD", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		PlotPedigreeKinshipVsGeneticIBD.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "plot/PlotPedigreeKinshipVsGeneticIBD.py"), \
							site_handler))
		executableClusterSizeMultiplierList.append((PlotPedigreeKinshipVsGeneticIBD, 0))
		
		DetectWrongLabelByCompKinshipVsIBD = Executable(namespace=namespace, name="DetectWrongLabelByCompKinshipVsIBD", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		DetectWrongLabelByCompKinshipVsIBD.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "pedigree/DetectWrongLabelByCompKinshipVsIBD.py"), \
							site_handler))
		executableClusterSizeMultiplierList.append((DetectWrongLabelByCompKinshipVsIBD, 0))
		
		ConvertPlinkBIM = Executable(namespace=namespace, name="ConvertPlinkBIM", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		ConvertPlinkBIM.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/converter/ConvertPlinkBIM.py"), \
							site_handler))
		executableClusterSizeMultiplierList.append((ConvertPlinkBIM, 0))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.run_type not in [1,5]:	#the run_types that involve LD-pruning
			#without commenting out db_vervet connection code. schema "genome" wont' be default path.
			db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, db_user=self.db_user,
							db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema="genome")
			db_genome.setup(create_tables=False)
			#chrOrder=2 means chromosomes are not ordered alphabetically but by their sizes (descendingly)
			oneGenomeData = db_genome.getOneGenomeData(tax_id=60711, chr_gap=0, chrOrder=2, sequence_type_id=9,\
													maxPseudoChrSize=1000000000)	#plink could not handle a chromosome of >2^31 bp length
			chr_id2cumu_chr_start = oneGenomeData.chr_id2cumu_chr_start
			if self.run_type==4:
				#fake selected contigs as X chromosomes and ignore others.
				sexContigList = ['Contig83', 'Contig149', 'Contig193']
				for contig in sexContigList:
					if contig in chr_id2cumu_chr_start:
						newChr, cumuStart = chr_id2cumu_chr_start.get(contig)
						chr_id2cumu_chr_start[contig] = ['X', cumuStart]	#assign some contigs to X
				#2012.9.12 contigs beyond 193 is a waste of time.
				self.maxContigID=195
				sys.stderr.write("Warning: maxContigID is manually set to 195 since contigs with ID bigger than 193 are useless.\n")
			ModifyTPEDRunType = 3	#plink LD-prune will skip non-recognizable chromosomes (1-24, human), so fake the chromosome
			if not self.mergeListFname:
				sys.stderr.write("Error: mergeListFname %s is nothing. required for this run_type %s.\n"%(self.mergeListFname,\
																										self.run_type))
				sys.exit(3)
			
			#for LD pruning, need to get rid files with too few loci
			needToKnowNoOfLoci = True
			#files with too few loci (like 1) cause problem by becoming empty after LD-prune (why? at least one SNP).
			minNoOfLoci = 2
			
			#only plink mendel job needs full VRC pedigree
			treatEveryOneIndependent = True
		else:
			#only plink mendel-error (-y1) and mark mendel-error-call missing (-y5) job needs full VRC pedigree
			treatEveryOneIndependent = False
			
			chr_id2cumu_chr_start = None
			ModifyTPEDRunType = 1	#plink mendel doesn't skip non-human chromosomes
			
			needToKnowNoOfLoci = False
			minNoOfLoci = None	#2012.10.19 bugfix
		
		db_vervet = self.db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		# Create a abstract dag
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
											checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
											pegasusFolderName=self.pegasusFolderName,\
											maxContigID=self.maxContigID, \
											minContigID=self.minContigID,\
											db_vervet=db_vervet, \
											needToKnowNoOfLoci=needToKnowNoOfLoci,\
											minNoOfLoci=minNoOfLoci)	#files with too few loci cause problem by becoming empty after LD-prune
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
		
		if self.run_type==5:
			vcf2PlinkReturnMode = 4
		else:
			vcf2PlinkReturnMode = 2
		vcf2PlinkJobData = self.addVCF2PlinkJobs(workflow, inputData=inputData, db_vervet=db_vervet, minMAC=None, minMAF=None,\
						maxSNPMissingRate=None, transferOutput=False,\
						maxContigID=self.maxContigID, outputDirPrefix="vcf2plink", outputPedigreeAsTFAM=True,\
						treatEveryOneIndependent=treatEveryOneIndependent,\
						returnMode=vcf2PlinkReturnMode, ModifyTPEDRunType=ModifyTPEDRunType, \
						chr_id2cumu_chr_start=chr_id2cumu_chr_start)
		
		if self.run_type==1:	#plink mendel
			mendelJobData = self.addPlinkMendelErrorJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="mendel", \
						returnMode=2)
		elif self.run_type==2:	#LD-prune
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
		elif self.run_type==3:	#IBD check
			no_of_individuals = inputData.jobDataLs[0].file.no_of_individuals
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
			if self.kinshipFname:
				kinshipFile = self.registerOneInputFile(inputFname=self.kinshipFname, folderName=self.pegasusFolderName)
			else:
				kinshipFile = None
			fontFile = self.registerOneInputFile(inputFname=os.path.expanduser("~/FreeSerif.ttf"), folderName=self.pegasusFolderName)
			self.addPlinkIBDCheckJobs(workflow=None, inputData=LDPruneJobData, kinshipFile=kinshipFile, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="ibdCheck", fontFile=fontFile, no_of_individuals=no_of_individuals)
			self.addPlinkIBDCheckJobs(workflow=None, inputData=LDPruneJobData, kinshipFile=kinshipFile, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="ibdCheckWithNonFounderFreq", useAlleleFrequencyFromNonFounders=True,\
						fontFile=fontFile, no_of_individuals=no_of_individuals)
		elif self.run_type==4:	#sex check
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
			self.addPlinkSexCheckJobs(workflow=None, inputData=LDPruneJobData, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="sexCheck")
		elif self.run_type==5:
			missingJobData = self.markMendelErrorLociMissingSubWorkflow(inputData=vcf2PlinkJobData, \
													transferOutput=False, \
													maxContigID=self.maxContigID, tfamJob=vcf2PlinkJobData.tfamJob, \
													outputDirPrefix='markMendelErrorLociMissing')
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			self.addPlinkBinaryConversionJobs(inputData=missingJobData, transferOutput=True, maxContigID=self.maxContigID, \
											tfamJob=vcf2PlinkJobData.tfamJob, \
											outputDirPrefix='plinkBinary', returnMode=1, mergeListFile=mergeListFile)
		else:
			sys.stderr.write("run_type %s not supported.\n"%(self.run_type))
			sys.exit(0)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


if __name__ == '__main__':
	main_class = PlinkOnVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()