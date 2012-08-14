#!/usr/bin/env python
"""
Examples:
	#2012.5.11 convert alignment read group (sample id) into UCLAID
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/ 
		-o workflow/SampleIDInUCLAID_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml 
		-u yh -y4 -l hcondor -j hcondor  -z localhost
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
	
	# 2012.8.10 IBD check
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_14/ -o workflow/PlinkIBDCheck_Method14.xml -C 1 
		-H -l hcondor -j hcondor  -u yh -z localhost
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		 -z localhost  -y3 -c 0.01 -g ./aux/Method14_LDPrune_merge_list.2012.8.10T0441.txt
	
	# 2012.8.9 LD-prune a folder of VCF files into plink, need the db tunnel (-H) for output pedigree in tfam
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	# -c 1 is to sample all locus data (samplingRate). LDPruneMinR2=0.3 (-R), LDPruneWindowSize=500 (-W), 
	# LDPruneWindowShiftSize=100 (-Z)
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719/trioCaller_vcftoolsFilter/ 
		-o workflow/ToPlinkFilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719.xml
		-y 2  -E
		-l condorpool -j condorpool
		-u yh -z uclaOffice  -C 4 -H -R 0.3 -W 500 -Z 100
		#-V 90 -x 100 -c 1
	
	# 2012.8.10 LD pruning
	%s  -I ~/NetworkData/vervet/db/genotype_file/method_14/ -o workflow/PlinkMendelError_Method14.xml
		-E -C 4  -H -l hcondor -j hcondor  -u yh -z localhost
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
		-z localhost  -y2  -c 0.01 -g ./aux/Method14_LDPrune_merge_list.txt
	
	# 2012.8.13 sex check using the top 195 contigs (-x 195) (Contig 83, 149,193 are sex chromosomes)
	# no clustering (-C 1)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-o workflow/PlinkSexCheck_Method14_W100Z10R0.4_maxContigID195.xml
		-W 100 -Z 10 -R 0.4 -C 1  -H -l hcondor -j hcondor  -u yh -z localhost
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-z localhost  -y4 -c 0.01 -g ./aux/Method14_LDPrune_merge_list.2012.8.13T1702.txt -x 195
		
Description:
	2012.8.14 a plink workflow on folder of VCF files.
		1: mendel errors,\
		2: LD pruning, \
		3: IBD check,\
		4: sex check,
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq, \
	figureOutDelimiter, getColName2IndexFromHeader, utils
from pymodule import GenomeDB

from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.VCFFile import VCFFile
from GenericVCFWorkflow import GenericVCFWorkflow

class PlinkOnVCFWorkflow(GenericVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('LDPruneMinR2', 0, float): [0.2, 'R', 1, 'minimum r2 for LD pruning', ],\
						('locusSamplingRate', 0, float): [0.0001, 'c', 1, 'how many loci to sample', ],\
						('mergeListFname', 0, ): [None, 'g', 1, 'the file to contain the merge-list for plink, required for run_type>1', ],\
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \
							1: mendel errors,\
							2: LD pruning, \
							3: IBD check,\
							4: sex check,', ],\
						('LDPruneWindowSize', 1, int): [500, 'W', 1, ' window size for plink LD pruning'],\
						('LDPruneWindowShiftSize', 1, int): [100, 'Z', 1, 'adjacent window shift'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		GenericVCFWorkflow.__init__(self, **keywords)

	def addPlinkMendelErrorJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, locusSamplingRate=0.0001, outputDirPrefix="", returnMode=2):
		"""
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink mendel error jobs for %s  files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sMendelError"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		no_of_jobs += 1
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		no_of_jobs += 1
		
		mendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_mendel.tsv'))
		#each input has no header
		mendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=mendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob])
		no_of_jobs += 1
		
		imendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_imendel.tsv'))
		#each input has no header
		imendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=imendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob], \
							extraArguments='-k 1 -v 2')
		outputFile = File( os.path.join(plotOutputDir, 'individualMendelErrorHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(executable=workflow.DrawHistogram, inputFileList=[imendelMergeFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="N", whichColumnPlotLabel="NoOfMendelErrors", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, imendelMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2
		
		fmendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_fmendel.tsv'))
		fmendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=fmendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob])
		no_of_jobs += 1
		
		lmendelMergeFile = File(os.path.join(mergedOutputDir, 'merged_lmendel.tsv'))
		lmendelMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=lmendelMergeFile, transferOutput=False, parentJobLs=[mergedOutputDirJob], \
							extraArguments='')
		outputFile = File( os.path.join(plotOutputDir, 'locusMendelErrorHist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(executable=workflow.DrawHistogram, inputFileList=[lmendelMergeFile], \
							outputFile=outputFile, \
					whichColumn=2, whichColumnHeader=None, whichColumnPlotLabel="NoOfMendelErrors", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=locusSamplingRate,\
					parentJobLs=[plotOutputDirJob, lmendelMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 2
		
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
			outputFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
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
			
			no_of_jobs += 1
			
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
		
		returnData.jobDataLs.append(PassingData(jobLs=[mendelMergeJob], file=mendelMergeJob.output, \
											fileList=mendelMergeJob.outputLs))
		returnData.jobDataLs.append(PassingData(jobLs=[imendelMergeJob], file=imendelMergeJob.output, \
											fileList=imendelMergeJob.outputLs))
		returnData.jobDataLs.append(PassingData(jobLs=[fmendelMergeJob], file=fmendelMergeJob.output, \
											fileList=fmendelMergeJob.outputLs))
		returnData.jobDataLs.append(PassingData(jobLs=[lmendelMergeJob], file=lmendelMergeJob.output, \
											fileList=lmendelMergeJob.outputLs))
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		##2012.7.21 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
						outputDirPrefix=outputDirPrefix)
		return newReturnData
	
	def writePlinkMergeListFile(self, outputFname=None, extractedFilenameTupleList=[]):
		"""
		"""
		outf = open(outputFname, 'w')
		for tuple in extractedFilenameTupleList:
			outf.write('%s %s %s\n'%(tuple[0], tuple[1], tuple[2]))
		del outf
		
	def addPlinkLDPruneJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, LDPruneMinR2=0.1, outputDirPrefix="", returnMode=1, \
						LDPruneWindowSize=100, LDPruneWindowShiftSize=5, mergeListFile=None):
		"""
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink LD pruning jobs for %s  files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sLDPrune"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		no_of_jobs += 1
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		no_of_jobs += 1
		
		
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
			plinkLDPruneJob = self.addPlinkJob(executable=self.plink, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=LDPruneFnamePrefix, outputOption='--out',\
					LDPruneWindowSize=LDPruneWindowSize, LDPruneWindowShiftSize=LDPruneWindowShiftSize, \
					LDPruneByPairwiseR2=True, LDPruneMinR2=LDPruneMinR2,\
					extraDependentInputLs=None, transferOutput=transferOneContigPlinkOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ jobData.jobLs)
			
			no_of_jobs += 1
			
			
			extractFnamePrefix = os.path.join(topOutputDir, '%s_extract'%(commonPrefix))
			plinkExtractJob = self.addPlinkJob(executable=self.plinkConvert, \
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
			no_of_jobs += 1
			i += 1
			
		
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
		no_of_jobs += 1
		
		
		if returnMode==1 or returnMode==3:
			returnData.jobDataLs.append(PassingData(jobLs=[plinkMergeJob], file=plinkMergeJob.bedFile, \
											fileList=plinkMergeJob.outputLs))
			
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		return returnData
	
	
	def addPlinkIBDCheckJobs(self, workflow=None, inputData=None, transferOutput=True,\
						maxContigID=None, outputDirPrefix="", returnMode=1,):
		"""
		2012.8.9
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding plink IBD checking jobs for %s  files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sIBDCheck"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		mergedOutputDir = "%sMerged"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
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
			IBDCheckFnamePrefix = os.path.join(topOutputDir, '%s'%(commonPrefix))
			plinkIBDCheckJob = self.addPlinkJob(executable=self.plinkIBD, \
									bedFile=inputJob.bedFile, famFile=inputJob.famFile, bimFile=inputJob.bimFile, \
					outputFnamePrefix=IBDCheckFnamePrefix, outputOption='--out',\
					estimatePairwiseGenomeWideIBD=True,\
					extraDependentInputLs=None, transferOutput=transferOutput, \
					extraArguments=None, job_max_memory=2000,\
					parentJobLs =[topOutputDirJob]+ jobData.jobLs)
			
			
		
			returnData.jobDataLs.append(PassingData(jobLs=[plinkIBDCheckJob], file=plinkIBDCheckJob.genomeFile, \
											fileList=plinkIBDCheckJob.outputLs))
			
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
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
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
						xColumnHeader="PEDSEX", xColumnPlotLabel="sexByInput", \
						minNoOfTotal=2,\
						figureDPI=100, samplingRate=1,\
						parentJobLs=[plotOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			no_of_jobs +=2
			
			outputFile = File(os.path.join(topOutputDir, "%s.male.sexCheck"%(commonPrefix)))
			selectMaleDataJob = self.addDrawHistogramJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
								inputFileList=[plinkJob.sexcheckFile], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="PEDSEX", whichColumnPlotLabel="PEDSEX", \
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						samplingRate=1,\
						parentJobLs=[topOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments="-V 1 -x 1", transferOutput=False,  job_max_memory=2000)
			outputFile = File( os.path.join(plotOutputDir, '%s.maleInbreedCoeffByChrX_Hist.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selectMaleDataJob.output], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="F", whichColumnPlotLabel="inbreedCoeffByChrX", \
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						figureDPI=100, samplingRate=1,\
						parentJobLs=[plotOutputDirJob, selectMaleDataJob], \
						extraDependentInputLs=None, \
						extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
			no_of_jobs +=2

			outputFile = File(os.path.join(topOutputDir, "%s.female.sexCheck"%(commonPrefix)))
			selectFemaleDataJob = self.addDrawHistogramJob(workflow=workflow, executable=workflow.SelectRowsFromMatrix, \
								inputFileList=[plinkJob.sexcheckFile], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="PEDSEX", whichColumnPlotLabel="PEDSEX", \
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
						minNoOfTotal=10,\
						samplingRate=1,\
						parentJobLs=[topOutputDirJob, plinkJob], \
						extraDependentInputLs=None, \
						extraArguments="-V 2 -x 2", transferOutput=False,  job_max_memory=2000)
			outputFile = File( os.path.join(plotOutputDir, '%s.femaleInbreedCoeffByChrX_Hist.png'%(commonPrefix)))
			#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
			self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selectFemaleDataJob.output], \
								outputFile=outputFile, \
						whichColumn=None, whichColumnHeader="F", whichColumnPlotLabel="inbreedCoeffByChrX", \
						logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
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
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.run_type!=1:
			#without commenting out db_vervet connection code. schema "genome" wont' be default path.
			db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
							password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
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
						
		else:
			chr_id2cumu_chr_start = None
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
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
											minContigID=self.minContigID)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("No VCF files in this folder , %s.\n"%self.inputDir)
			sys.exit(0)
		
		sampleIDFile = None
		
		if self.run_type!=1:
			ModifyTPEDRunType = 3	#plink will skip non-recognizable chromosomes (1-24, human), so fake the chromosome
			if not self.mergeListFname:
				sys.stderr.write("Error: mergeListFname %s is nothing. required for this run_type %s.\n"%(self.mergeListFname,\
																										self.run_type))
				sys.exit(3)
		else:
			ModifyTPEDRunType = 1
		vcf2PlinkJobData = self.addVCF2PlinkJobs(workflow, inputData=inputData, db_vervet=db_vervet, minMAC=None, minMAF=None,\
						maxSNPMissingRate=None, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="vcf2plink", outputPedigreeAsTFAM=True,\
						returnMode=2, ModifyTPEDRunType=ModifyTPEDRunType, chr_id2cumu_chr_start=chr_id2cumu_chr_start)
		if self.run_type==1:
			mendelJobData = self.addPlinkMendelErrorJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="mendel", locusSamplingRate=self.locusSamplingRate, 
						returnMode=2)
		elif self.run_type==2:
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
		elif self.run_type==3:
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
			self.addPlinkIBDCheckJobs(workflow=None, inputData=LDPruneJobData, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="ibdCheck")
		elif self.run_type==4:
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName)
			LDPruneJobData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
						maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
						LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
						outputDirPrefix="ldPrune", returnMode=1, \
						mergeListFile=mergeListFile)
			self.addPlinkSexCheckJobs(workflow=None, inputData=LDPruneJobData, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="sexCheck")
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