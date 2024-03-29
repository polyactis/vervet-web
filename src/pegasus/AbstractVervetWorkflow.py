#!/usr/bin/env python
"""
2012.6.12 a NGS-workflow that derives from AbstractVCFWorkflow and specific for vervet repository
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from Pegasus.DAX3 import Executable, PFN, File
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from vervet.src import VervetDB

parentClass = AbstractVCFWorkflow

class AbstractVervetWorkflow(parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(parentClass.option_default_dict)
	option_default_dict.update(parentClass.db_option_dict.copy())
	
	option_default_dict.update({
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.6.12
		"""
		parentClass.__init__(self, **keywords)
		
		self.db = self.db_vervet	#2013.1.25 main db

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
	
	def addExtractSampleIDJob(self, workflow=None, passingData=None, inputFile=None, \
							outputFile=None,\
							min_coverage=None, max_coverage=None,\
							pop_tax_id_ls_str=None, pop_site_id_ls_str=None, pop_country_id_ls_str=None, popHeader=None,\
							pop_sampleSize=None, is_contaminated=None, outputFormat=2, returnData=None,\
							transferOutput=False, extraArguments=None, extraArgumentList=None, job_max_memory=1000, \
							parentJobLs=None,\
							**keywords):
		"""
		
		argument popHeader and pop_sampleSize are required if you want to uniformly sample a few from extracted ones.
		pass a PassingData object to returnData if you want extractPopSampleIDJob to be included in final transfer-out:
				returnData.jobDataLs.append(PassingData(jobLs=[extractPopSampleIDJob], file=extractPopSampleIDJob.output, \
											fileLs=[extractPopSampleIDJob.output]))
		argument inputFile is a vcf file.
			if inputFile is None, this function will try to derive it from passingData
		if outputFile or parentJobLs is not given, this function will try to derive it from passingData.
		
		2013.07.03 added argument is_contaminated (whether to fetch contaminated samples or not)
		2013.06.24 expose argument outputFormat (default=2)
			1: a subset VCF file; \n\
			2: file with a list of sample IDs (one per line), with header
			3: file with a list of sample IDs (one per line), without header
		2013.06.11
			added argument inputFile and outputFile (previously these are derived from passingData)
			added argument min_coverage, max_coverage
			removed argument outputDirPrefix
		2012.10.24
			moved from popgen/CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.py
		2012.10.15
		"""
		if workflow is None:
			workflow = self
		if extraArgumentList is None:
			extraArgumentList = []
		if parentJobLs is None:
			parentJobLs = []
		
		#2013.06.11 handle inputFile
		if inputFile is None and passingData is not None:
			#use a random VCF file as input
			jobData = passingData.chr2jobDataLs.values()[0][0]
			inputFile = jobData.vcfFile
			parentJobLs.extend(jobData.jobLs)
			
			reduceOutputDirJob = passingData.reduceOutputDirJob
			parentJobLs.append(reduceOutputDirJob)
		
		inputFBaseName = os.path.basename(inputFile.name)
		commonPrefix = inputFBaseName.split('.')[0]
		#ExtractSamplesFromVCF for the 1st population
		if popHeader is not None:	#2013.06.11
			commonPrefix = '%s_pop%s'%(commonPrefix, popHeader)
		
		#2013.06.11 handle outputFile
		if outputFile is None and passingData is not None:
			reduceOutputDirJob = passingData.reduceOutputDirJob
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_sampleID.tsv'%(commonPrefix)))
		
		if outputFormat is not None:
			extraArgumentList.append('--outputFormat %s'%(outputFormat))	#output a list of sample IDs
		if min_coverage is not None:
			extraArgumentList.append("--min_coverage %s"%(min_coverage))
		if max_coverage is not None:
			extraArgumentList.append("--max_coverage %s"%(max_coverage))
		if is_contaminated is not None:
			extraArgumentList.append("--is_contaminated %s"%(is_contaminated))
		if pop_tax_id_ls_str:
			extraArgumentList.append("--tax_id_ls %s"%(pop_tax_id_ls_str))
		if pop_site_id_ls_str:
			extraArgumentList.append("--site_id_ls %s"%(pop_site_id_ls_str))
		if pop_country_id_ls_str:
			extraArgumentList.append("--country_id_ls %s"%(pop_country_id_ls_str))
		
		extractPopSampleIDJob = self.addGenericDBJob(workflow=workflow, executable=self.ExtractSamplesFromVCF, inputFile=inputFile, \
					outputFile=outputFile, inputFileList=None, \
					parentJobLs=parentJobLs, extraDependentInputLs=None, \
					extraOutputLs=None, transferOutput=transferOutput, \
					extraArguments=extraArguments, extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, \
					sshDBTunnel=self.needSSHDBTunnel, \
					key2ObjectForJob=None)
		if returnData:
			returnData.jobDataLs.append(PassingData(jobLs=[extractPopSampleIDJob], file=extractPopSampleIDJob.output, \
											fileLs=[extractPopSampleIDJob.output]))
		if pop_sampleSize and pop_sampleSize>1:	#uniform random sampling + some other filtering
			sampleIDHeader='sampleID'	#determined by extractPopSampleIDJob
			#. SelectRowsWithinCoverageRange
			minMedianDepth = 2
			maxMedianDepth = 15
			extraArguments = " --minMedianDepth %s --maxMedianDepth %s "%(minMedianDepth, maxMedianDepth)
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_sampleID_depth%s_%s.tsv'%(commonPrefix,\
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
			if returnData:
				returnData.jobDataLs.append(PassingData(jobLs=[selectRowsWithinCoverageRangeJob], file=selectRowsWithinCoverageRangeJob.output, \
											fileLs=[selectRowsWithinCoverageRangeJob.output]))
			#. optional, a ReplaceIndividualIDInMatrixFileWithReadGroup job (for VRC) on the IBD check result
			translateIDInIBDResultJob = self.addTranslateIDInIBDCheckResultJob(workflow=workflow, plinkIBDCheckOutputFile=self.plinkIBDCheckOutputFile, \
										pop_country_id_ls_str=pop_country_id_ls_str, \
										pop_site_id_ls_str=pop_site_id_ls_str, popHeader=popHeader,\
										readGroupFile=selectRowsWithinCoverageRangeJob.output, parentJobLs=[selectRowsWithinCoverageRangeJob], \
										sampleIDHeader=sampleIDHeader, transferOutput=transferOutput)
			#. SampleRows job
			outputFile = File(os.path.join(reduceOutputDirJob.output, '%s_sampleSize%s.tsv'%(commonPrefix, pop_sampleSize)))
			extraArgumentList = [" --sampleSize %s "%(pop_sampleSize), "--maxIBDSharing %s"%(self.maxIBDSharing)]
			if translateIDInIBDResultJob:
				extraArgumentList.extend(["--plinkIBDCheckOutputFname", translateIDInIBDResultJob.output])
				extraDependentInputLs = [translateIDInIBDResultJob.output]
				if returnData:
					returnData.jobDataLs.append(PassingData(jobLs=[translateIDInIBDResultJob], file=translateIDInIBDResultJob.output, \
											fileLs=[translateIDInIBDResultJob.output]))
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
			
			if returnData:
				returnData.jobDataLs.append(PassingData(jobLs=[sampleRowsJob], file=sampleRowsJob.output, \
											fileLs=[sampleRowsJob.output]))
			#rename the extractPopSampleIDJob
			extractPopSampleIDJob = sampleRowsJob
		return extractPopSampleIDJob
	
	def addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(self, executable=None, inputFile=None, outputFile=None,\
						sampleID2FamilyCountF=None, polymuttDatFile=None, treatEveryOneIndependent=False, \
						outputFileFormat=None,\
						replicateIndividualTag=None, addUngenotypedDuoParents=False, \
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						extraArguments=None, job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2013.07.19 added argument addUngenotypedDuoParents
		2012.9.13
			add argument treatEveryOneIndependent for OutputVRCPedigreeInTFAMGivenOrderFromFile.
		2012.8.9
		"""
		extraArgumentList = []
		extraOutputLs = []
		if extraArguments:
			extraArgumentList.append(extraArguments)
		if outputFileFormat is not None:
			extraArgumentList.append("--outputFileFormat %s"%(outputFileFormat))
		if treatEveryOneIndependent:
			extraArgumentList.append("--treatEveryOneIndependent")
		if replicateIndividualTag is not None:
			extraArgumentList.append("--replicateIndividualTag %s"%(replicateIndividualTag))
		if addUngenotypedDuoParents:
			extraArgumentList.append("--addUngenotypedDuoParents 1")
		
		if sampleID2FamilyCountF:
			extraArgumentList.extend(["--sampleID2FamilyCountFname", sampleID2FamilyCountF])
			extraOutputLs.append(sampleID2FamilyCountF)
		if polymuttDatFile:
			extraArgumentList.extend(["--polymuttDatFname", polymuttDatFile])
			extraOutputLs.append(polymuttDatFile)
		
		job= self.addGenericDBJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, \
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, objectWithDBArguments=self, **keywords)
		job.sampleID2FamilyCountF = sampleID2FamilyCountF
		job.polymuttDatFile = polymuttDatFile
		return job
	
	def addOutputVCFAlignmentDepthRangeJob(self, executable=None, inputFile=None, \
						ref_ind_seq_id=None, depthFoldChange=None, minGQ=None,\
						outputFile=None, outputFileFormat=None,\
						data_dir=None,
						sequence_filtered=None, excludeContaminant=None, \
						local_realigned=None, reduce_reads=None, completedAlignment=None, \
						alignment_method_id=None, alignment_outdated_index=None,\
						extraArgumentList=None,\
						parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
						job_max_memory=2000, sshDBTunnel=False, **keywords):
		"""
		2013.06.13
		"""
		extraOutputLs = []
		if extraDependentInputLs is None:
			extraArgumentList = []
		if data_dir is not None:
			extraArgumentList.append("--data_dir %s"%(data_dir))
		if outputFileFormat is not None:
			extraArgumentList.append("--outputFileFormat %s"%(outputFileFormat))
		if ref_ind_seq_id is not None:
			extraArgumentList.append("--ref_ind_seq_id %s"%(ref_ind_seq_id))
		if depthFoldChange is not None:
			extraArgumentList.append("--depthFoldChange %s"%(depthFoldChange))
		if minGQ is not None:
			extraArgumentList.append("--minGQ %s"%(minGQ))
		if sequence_filtered is not None:
			extraArgumentList.append("--sequence_filtered %s"%(sequence_filtered))
		if excludeContaminant is not None:
			extraArgumentList.append("--excludeContaminant")
		if local_realigned is not None:
			extraArgumentList.append("--local_realigned %s"%(local_realigned))
		if reduce_reads is not None:
			extraArgumentList.append("--reduce_reads %s"%(reduce_reads))
		if completedAlignment is not None:
			extraArgumentList.append("--completedAlignment %s"%(completedAlignment))
		if alignment_method_id is not None:
			extraArgumentList.append("--alignment_method_id %s"%(alignment_method_id))
		if alignment_outdated_index is not None:
			extraArgumentList.append("--alignment_outdated_index %s"%(alignment_outdated_index))
		job= self.addGenericDBJob(executable=executable, inputFile=inputFile, outputFile=outputFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArguments=None, extraArgumentList=extraArgumentList, \
						job_max_memory=job_max_memory, sshDBTunnel=sshDBTunnel, objectWithDBArguments=self, **keywords)
		return job
	
	def connectDB(self):
		"""
		2012.9.24
			establish db connection for all derivative classes
		"""
		parentClass.connectDB(self)	#2013.11.25
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, db_user=self.db_user,
					db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		self.db = db_vervet	#2013.04.09
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		#self.refFastaFList = self.getReferenceSequence(workflow=self)	#2013.1.25 done in run()
	
	def getReferenceSequence(self, workflow=None, **keywords):
		"""
		2013.3.20 yh_pegasus.registerRefFastaFile() returns a PassingData
		2013.1.25
		"""
		sys.stderr.write("Getting reference sequences ...")
		if workflow is None:
			workflow = self
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.data_dir, refSequence.path)
		registerReferenceData = yh_pegasus.registerRefFastaFile(workflow=workflow, refFastaFname=refFastaFname, \
							registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		
		sys.stderr.write(" %s files.\n"%(len(registerReferenceData.refFastaFList)))
		return registerReferenceData
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
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
		
		#2013.06.13
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "db/output/OutputVCFAlignmentDepthRange.py"), \
											name='OutputVCFAlignmentDepthRange', \
											clusterSizeMultipler=1)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "db/output/OutputVRCPedigreeInTFAMGivenOrderFromFile.py"), \
											name='OutputVRCPedigreeInTFAMGivenOrderFromFile', \
											clusterSizeMultipler=0.8)

if __name__ == '__main__':
	main_class = AbstractVervetWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()