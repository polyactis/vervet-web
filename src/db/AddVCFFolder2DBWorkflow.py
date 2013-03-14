#!/usr/bin/env python
"""
Examples:
	#2012.5.11 
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/ 
		-o workflow/2DB/AddVCF2DB_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml 
		-s ... -u yh -l hcondor -j hcondor  -z localhost 
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ 
	
	# 2012.5.10 run on hoffman2 condor, turn on checkEmptyVCFByReading (-E), require db connection (--needSSHDBTunnel), no clustering (-C1)
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8T21.42/trioCaller_vcftoolsFilter/
		-o workflow/2DB/AddVCF2DB_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.2012.5.6_trioCaller.2012.5.8.xml
		-s ... -E --needSSHDBTunnel -C1
		-l hcondor -j hcondor  -u yh -z localhost 
		-e /u/home/eeskin/polyacti/ 
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
	# 2012.7.16 add a folder of VCF files to DB without checking zero-loci VCF
	%s -I FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719/trioCaller_vcftoolsFilter/ 
		-o workflow/2DB/AddVCF2DB_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs.MAC10.MAF.05_trioCaller.2012.5.21T1719.xml
		-s ... -l condorpool -j condorpool
		-u yh -z uclaOffice -C1
	
Description:
	#2012.5.9
		the usual -c (commit) is not here. All DB jobs are run with commit=True.
	2012.8.3 if such a workflow with clustering on (several AddVCFFolder2DB jobs crammed into one) fails halfway,
		you can safely re-run it. Already-imported files would be checked and not be imported again (MD5SUM).
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, NextGenSeq, figureOutDelimiter, getColName2IndexFromHeader
from pymodule import VCFFile
from pymodule.pegasus import yh_pegasus
from vervet.src.pegasus.GenericVCFWorkflow import GenericVCFWorkflow
from vervet.src import VervetDB, AbstractVervetWorkflow

class AddVCFFolder2DBWorkflow(GenericVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('maxContigID', 0, int): [None, 'x', 1, 'if contig ID is beyond this number, it will not be included. If None or 0, no restriction.', ],\
						('genotypeMethodShortName', 1, ):[None, 's', 1, 'column short_name of GenotypeMethod table,\
			will be created if not present in db.'],\
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. '],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		self.inputDir = os.path.abspath(self.inputDir)
	
	def registerCustomJars(self, workflow, ):
		"""
		2012.2.10
		"""
		site_handler = self.site_handler
	
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
		
		AddGenotypeMethod2DB = Executable(namespace=namespace, name="AddGenotypeMethod2DB", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		AddGenotypeMethod2DB.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/input/AddGenotypeMethod2DB.py"), site_handler))
		executableList.append(AddGenotypeMethod2DB)
		
		UpdateGenotypeMethodNoOfLoci = Executable(namespace=namespace, name="UpdateGenotypeMethodNoOfLoci", \
											version=version, \
											os=operatingSystem, arch=architecture, installed=True)
		UpdateGenotypeMethodNoOfLoci.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "db/UpdateGenotypeMethodNoOfLoci.py"), site_handler))
		executableList.append(UpdateGenotypeMethodNoOfLoci)
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			workflow.addExecutable(executable)
			setattr(workflow, executable.name, executable)
	
	
	def addAddGenotypeMethod2DBJob(self, executable=None, inputFile=None, genotypeMethodShortName=None,\
								logFile=None, data_dir=None, commit=False, parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
								extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.6.27
		"""
		extraArgumentList = ['--genotypeMethodShortName', genotypeMethodShortName]
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
		if data_dir:
			extraArgumentList.extend(['--data_dir', data_dir])
		if commit:
			extraArgumentList.append('--commit')
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericJob(executable=executable, inputFile=inputFile, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[logFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		return job
	
	def addUpdateGenotypeMethodNoOfLociJob(self, executable=None, genotypeMethodShortName=None, genotypeMethodID=None,\
								logFile=None, data_dir=None, commit=False, parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
								extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.6.27
		"""
		extraArgumentList = []
		if logFile:
			extraArgumentList.extend(["--logFilename", logFile])
		if data_dir:
			extraArgumentList.extend(['--data_dir', data_dir])
		if commit:
			extraArgumentList.append('--commit')
		if genotypeMethodShortName:
			extraArgumentList.extend(['--genotypeMethodShortName', genotypeMethodShortName, ])
		if genotypeMethodID:
			extraArgumentList.extend(['--genotypeMethodID', genotypeMethodID, ])
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=[logFile],\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory, **keywords)
		self.addDBArgumentsToOneJob(job=job, objectWithDBArguments=self)
		return job
	
	def addJobs(self, workflow=None, inputData=None, db_vervet=None, genotypeMethodShortName=None, commit=None,\
			data_dir=None, checkEmptyVCFByReading=False, transferOutput=True,\
			maxContigID=None, outputDirPrefix="", needSSHDBTunnel=False):
		"""
		2012.5.9
		"""
		sys.stderr.write("Adding VCF2DB jobs for %s vcf files ... "%(len(inputData.jobDataLs)))
		no_of_jobs= 0
		
		
		topOutputDir = "%sVCF2DB"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		firstVCFFile = inputData.jobDataLs[0].vcfFile
		logFile = File(os.path.join(topOutputDir, 'AddGenotypeMethod2DB.log'))
		addGM2DBJob = self.addAddGenotypeMethod2DBJob(executable=self.AddGenotypeMethod2DB, inputFile=firstVCFFile, \
								genotypeMethodShortName=genotypeMethodShortName,\
								logFile=logFile, data_dir=data_dir, commit=commit, parentJobLs=[], extraDependentInputLs=[], transferOutput=True, \
								extraArguments=None, job_max_memory=10, sshDBTunnel=needSSHDBTunnel)
		updateGMlogFile = File(os.path.join(topOutputDir, 'updateGM.log'))
		updateGMNoOfLociJob = self.addUpdateGenotypeMethodNoOfLociJob(executable=self.UpdateGenotypeMethodNoOfLoci, \
																	genotypeMethodShortName=genotypeMethodShortName,\
								logFile=updateGMlogFile, data_dir=data_dir, commit=commit, parentJobLs=[topOutputDirJob], \
								extraDependentInputLs=[], transferOutput=True, \
								extraArguments=None, job_max_memory=20, sshDBTunnel=needSSHDBTunnel)
		
		no_of_jobs += 2
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
			logFile = File(os.path.join(topOutputDir, 'AddVCFFile2DB_%s.log'%(self.getChrFromFname(inputF.name))))
			addVCFJob = self.addAddVCFFile2DBJob(executable=self.AddVCFFile2DB, inputFile=inputF, genotypeMethodShortName=genotypeMethodShortName,\
						logFile=logFile, format="VCF", data_dir=data_dir, checkEmptyVCFByReading=checkEmptyVCFByReading, commit=commit, \
						parentJobLs=[addGM2DBJob]+jobData.jobLs, extraDependentInputLs=[], transferOutput=True, \
						extraArguments=None, job_max_memory=1000, sshDBTunnel=needSSHDBTunnel)
			workflow.depends(parent=addVCFJob, child=updateGMNoOfLociJob)
			no_of_jobs += 1
		sys.stderr.write("%s jobs. Done.\n"%(no_of_jobs))
		#include the tfam (outputList[1]) into the fileList
		returnData.jobDataLs.append(PassingData(jobLs=[updateGMNoOfLociJob], file=updateGMlogFile, \
											fileList=[updateGMlogFile]))
		return returnData
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = self.db_vervet
		
		# Create a abstract dag
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		inputData = self.registerAllInputFiles(workflow, inputDir=self.inputDir, input_site_handler=self.input_site_handler, \
										checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
										pegasusFolderName=self.pegasusFolderName)
		if len(inputData.jobDataLs)<=0:
			sys.stderr.write("Warning: Number of VCF files in this folder, %s, <=0.\n"%self.inputDir)
			sys.exit(0)
		
		if self.run_type==1:
			self.addJobs(workflow=workflow, inputData=inputData, db_vervet=db_vervet, genotypeMethodShortName=self.genotypeMethodShortName, \
						commit=True,\
						data_dir=self.data_dir, checkEmptyVCFByReading=self.checkEmptyVCFByReading, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="", needSSHDBTunnel=self.needSSHDBTunnel)
		else:
			sys.stderr.write("run_type %s not supported.\n"%(self.run_type))
			sys.exit(0)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


if __name__ == '__main__':
	main_class = AddVCFFolder2DBWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()