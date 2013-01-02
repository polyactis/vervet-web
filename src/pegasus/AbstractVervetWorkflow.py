#!/usr/bin/env python
"""
2012.6.12
	a NGS-workflow that derives from AbstractVCFWorkflow and specific for vervet repository
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from vervet.src import VervetDB

class AbstractVervetWorkflow(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update(AbstractVCFWorkflow.db_option_dict.copy())
	
	option_default_dict.update({
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.6.12
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
	
	def outputAlignmentDepthAndOthersForFilter(self, db_vervet=None, outputFname=None, ref_ind_seq_id=524, \
											foldChange=2, minGQ=30):
		"""
		2012.6.12
			added argument db_vervet, moved from FilterVCFPipeline.py
		2011-9-2
		"""
		sys.stderr.write("Outputting sequence coverage to %s ..."%outputFname)
		import csv
		counter = 0
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['aln.id', 'minDepth', 'maxDepth', 'minGQ'])
		TableClass = VervetDB.IndividualAlignment
		query = TableClass.query.filter(TableClass.median_depth!=None)
		if ref_ind_seq_id:
			query = query.filter(TableClass.ref_ind_seq_id==ref_ind_seq_id)
		query = query.order_by(TableClass.id)
		for row in query:
			minDepth = row.median_depth/float(foldChange)
			if abs(minDepth-0)<=0.001:	#if it's too close to 0, regard it as 0.
				minDepth = 0
			writer.writerow([row.getReadGroup(), minDepth, \
							row.median_depth*float(foldChange), minGQ])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
	def addTranslateIDInIBDCheckResultJob(self, workflow=None, plinkIBDCheckOutputFile=None, pop_country_id_ls_str=None, \
										pop_site_id_ls_str=None, popHeader=None,\
										readGroupFile=None, parentJobLs=None, sampleIDHeader='sampleID',\
										transferOutput=False):
		"""
		2012.10.24
			moved from popgen/CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.py
		2012.10.15
		"""
		job = None
		if plinkIBDCheckOutputFile:
			pop_country_id_set = set(getListOutOfStr(pop_country_id_ls_str))
			pop_site_id_set = set(getListOutOfStr(pop_site_id_ls_str))
			if 447 in pop_site_id_set or 135 in pop_country_id_set:	#either site = VRC or country = USA
				commonPrefix = os.path.splitext(plinkIBDCheckOutputFile.name)[0]
				outputFile = File('%s_%s_withReadGroup.tsv'%(commonPrefix, popHeader))
				extraArgumentList = [" --readGroupFname", readGroupFile, "--readGroupHeader %s"%(sampleIDHeader), \
									'--replaceColumnHeaderLs IID1,IID2']
				job = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.ReplaceIndividualIDInMatrixFileWithReadGroup, \
							inputFileList=None, inputFile=plinkIBDCheckOutputFile, outputFile=outputFile, \
							outputFnamePrefix=None, whichColumn=None, whichColumnHeader=sampleIDHeader, \
							minNoOfTotal=0,\
							samplingRate=1.0, \
							parentJob=None, parentJobLs=parentJobLs, \
							extraDependentInputLs=[readGroupFile], \
							extraArgumentList=extraArgumentList, transferOutput=transferOutput,  job_max_memory=1000,\
							sshDBTunnel=self.needSSHDBTunnel, \
							objectWithDBArguments=self,)
		
		return job
	
	def addExtractSampleIDJob(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=False, \
							pop_tax_id_ls_str=None, pop_site_id_ls_str=None, pop_country_id_ls_str=None, popHeader=None,\
							pop_sampleSize=None, returnData=None, **keywords):
		"""
		2012.10.24
			moved from popgen/CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.py
		2012.10.15
		"""
		if workflow is None:
			workflow = self
		#use a random VCF file as input
		jobData = passingData.chr2jobDataLs.values()[0][0]
		VCFFile = jobData.vcfFile
		inputFBaseName = os.path.basename(VCFFile.name)
		commonPrefix = inputFBaseName.split('.')[0]
		reduceOutputDirJob = passingData.reduceOutputDirJob
		#ExtractSamplesFromVCF for the 1st population
		outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID.tsv'%(commonPrefix, popHeader)))
		extraArgumentList = ['--outputFormat 2']
		if pop_tax_id_ls_str:
			extraArgumentList.append("--tax_id_ls %s"%(pop_tax_id_ls_str))
		if pop_site_id_ls_str:
			extraArgumentList.append("--site_id_ls %s"%(pop_site_id_ls_str))
		if pop_country_id_ls_str:
			extraArgumentList.append("--country_id_ls %s"%(pop_country_id_ls_str))
		
		extractPopSampleIDJob = self.addGenericDBJob(workflow=workflow, executable=self.ExtractSamplesFromVCF, inputFile=VCFFile, \
					outputFile=outputFile, inputFileList=None, \
					parentJobLs=[reduceOutputDirJob]+jobData.jobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=extraArgumentList, job_max_memory=2000, \
					sshDBTunnel=self.needSSHDBTunnel, \
					key2ObjectForJob=None)
		returnData.jobDataLs.append(PassingData(jobLs=[extractPopSampleIDJob], file=extractPopSampleIDJob.output, \
											fileList=[extractPopSampleIDJob.output]))
		if pop_sampleSize and pop_sampleSize>1:
			sampleIDHeader='sampleID'	#determined by extractPopSampleIDJob
			#. SelectRowsWithinCoverageRange
			minMedianDepth = 2
			maxMedianDepth = 15
			extraArguments = " --minMedianDepth %s --maxMedianDepth %s "%(minMedianDepth, maxMedianDepth)
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleID_depth%s_%s.tsv'%(commonPrefix, popHeader,\
																			minMedianDepth, maxMedianDepth)))
			selectRowsWithinCoverageRangeJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.SelectRowsWithinCoverageRange, \
						inputFileList=None, inputFile=extractPopSampleIDJob.output, outputFile=outputFile, \
						outputFnamePrefix=None, whichColumn=None, whichColumnHeader=sampleIDHeader, \
						logY=False, valueForNonPositiveYValue=-1, \
						minNoOfTotal=0,\
						samplingRate=1.0, \
						parentJobLs=[extractPopSampleIDJob], \
						extraDependentInputLs=None, \
						extraArguments=extraArguments, transferOutput=transferOutput, job_max_memory=100,\
						sshDBTunnel=self.needSSHDBTunnel, \
						objectWithDBArguments=self,)
			returnData.jobDataLs.append(PassingData(jobLs=[selectRowsWithinCoverageRangeJob], file=selectRowsWithinCoverageRangeJob.output, \
											fileList=[selectRowsWithinCoverageRangeJob.output]))
			#. optional, a ReplaceIndividualIDInMatrixFileWithReadGroup job (for VRC) on the IBD check result
			translateIDInIBDResultJob = self.addTranslateIDInIBDCheckResultJob(workflow=workflow, plinkIBDCheckOutputFile=self.plinkIBDCheckOutputFile, \
										pop_country_id_ls_str=pop_country_id_ls_str, \
										pop_site_id_ls_str=pop_site_id_ls_str, popHeader=popHeader,\
										readGroupFile=selectRowsWithinCoverageRangeJob.output, parentJobLs=[selectRowsWithinCoverageRangeJob], \
										sampleIDHeader=sampleIDHeader, transferOutput=transferOutput)
			#. SampleRows job
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_pop%s_sampleSize%s.tsv'%(commonPrefix, popHeader, pop_sampleSize)))
			extraArgumentList = [" --sampleSize %s "%(pop_sampleSize), "--maxIBDSharing %s"%(self.maxIBDSharing)]
			if translateIDInIBDResultJob:
				extraArgumentList.extend(["--plinkIBDCheckOutputFname", translateIDInIBDResultJob.output])
				extraDependentInputLs = [translateIDInIBDResultJob.output]
				returnData.jobDataLs.append(PassingData(jobLs=[translateIDInIBDResultJob], file=translateIDInIBDResultJob.output, \
											fileList=[translateIDInIBDResultJob.output]))
			else:
				extraDependentInputLs = None
			sampleRowsJob = self.addAbstractMatrixFileWalkerJob(workflow=workflow, executable=self.SampleRows, \
						inputFileList=None, inputFile=selectRowsWithinCoverageRangeJob.output, outputFile=outputFile, \
						outputFnamePrefix=None, whichColumn=None, whichColumnHeader=sampleIDHeader, \
						logY=False, valueForNonPositiveYValue=-1, \
						minNoOfTotal=0, \
						samplingRate=1.0, \
						parentJob=translateIDInIBDResultJob, parentJobLs=[selectRowsWithinCoverageRangeJob], \
						extraDependentInputLs=extraDependentInputLs, \
						extraArgumentList=extraArgumentList, transferOutput=transferOutput,  job_max_memory=1000)
			
			returnData.jobDataLs.append(PassingData(jobLs=[sampleRowsJob], file=sampleRowsJob.output, \
											fileList=[sampleRowsJob.output]))
			#rename the extractPopSampleIDJob
			extractPopSampleIDJob = sampleRowsJob
		return extractPopSampleIDJob
	
	def connectDB(self):
		"""
		2012.9.24
			establish db connection for all derivative classes
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
			
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(self, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		self.refFastaFList = refFastaFList
	
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		
		AbstractVCFWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		#mergeSameHeaderTablesIntoOne is used here on per chromosome basis, so allow clustering
		executableClusterSizeMultiplierList.append((workflow.mergeSameHeaderTablesIntoOne, 1))
		
		
		ReplaceIndividualIDInMatrixFileWithReadGroup = Executable(namespace=namespace, name="ReplaceIndividualIDInMatrixFileWithReadGroup", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		ReplaceIndividualIDInMatrixFileWithReadGroup.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/ReplaceIndividualIDInMatrixFileWithReadGroup.py"), site_handler))
		executableClusterSizeMultiplierList.append((ReplaceIndividualIDInMatrixFileWithReadGroup, 0.5))
		
		SelectRowsWithinCoverageRange = Executable(namespace=namespace, name="SelectRowsWithinCoverageRange", version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		SelectRowsWithinCoverageRange.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/SelectRowsWithinCoverageRange.py"), site_handler))
		executableClusterSizeMultiplierList.append((SelectRowsWithinCoverageRange, 0.5))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		