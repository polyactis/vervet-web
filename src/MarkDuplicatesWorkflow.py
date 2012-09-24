#!/usr/bin/env python
"""
Examples:
	
	#2012.3.21 run on condorpool
	%s -j condorpool -l condorpool -u yh -z uclaOffice -o MarkDup.xml -i ...
	
	#2012.3.21 run on uschpc (input data is on uschpc), for each top contig as well
	%s -j uschpc -l uschpc -u yh -z uclaOffice -o MarkDup....xml
		-i ... -e /home/cmb-03/mn/yuhuang/ -m /home/cmb-03/mn/yuhuang/tmp/
		-J /usr/usc/jdk/default/bin/java -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/ -D /Network/Data/vervet/db/
	
	#2012-3-21 on hoffman2's condor pool
	%s -j hcondor -l hcondor -u yh -i ShortRead2Alignment_Isq_631_700_vs_524_hcondor.2012.2.24T1902/
		-o MarkDuplicatesShortRead2Alignment_Isq_631_700_vs_524_hcondor.2012.2.24T1902.xml
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -J ~/bin/jdk/bin/java -m /work/
	 
Description:
	2012.3.21
		This workflow is a patch upon ShortRead2AlignmentPipeline.py because the latter
		fails to run jobs after MergeSam in proper parallel fashion.
		
		outputFname: if inputFname is like xxx_merged.bam, outputFname is xxx.bam. If not, it'll become xxx.markDup.bam.
		bwa and stampy are not used but added because it inherits ShortRead2AlignmentPipeline's registerCustomExecutables().
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow
from ShortRead2AlignmentPipeline import ShortRead2AlignmentPipeline

class MarkDuplicatesWorkflow(ShortRead2AlignmentPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						("bwa_path", 1, ): ["%s/bin/bwa", '', 1, 'bwa binary'],\
						("stampy_path", 1, ): ["%s/bin/stampy.py", '', 1, 'path to stampy.py'],\
						('inputFolder', 0, ): ['', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						("tmpDir", 1, ): ["/tmp/", 'm', 1, 'for MarkDuplicates.jar, default is /tmp/ but sometimes too small'],\
						})

	def __init__(self, **keywords):
		"""
		2012.3.21
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		#bwa and stampy are not used but added because it inherits ShortRead2AlignmentPipeline's registerCustomExecutables().
		self.bwa_path =  self.insertHomePath(self.bwa_path, self.home_path)
		self.stampy_path =  self.insertHomePath(self.stampy_path, self.home_path)
	
	def getMarkDupOutputFnameBasedOnInputFname(self, inputFname=None):
		"""
		2012.3.21
			outputFname: if inputFname is like xxx_merged.bam, outputFname is xxx.bam. If not, it'll become xxx.markDup.bam.
			
		"""
		baseFilename = os.path.basename(inputFname)
		filePrefix, fileSuffix = os.path.splitext(baseFilename)
		if filePrefix[-7:]=='_merged':
			outputFname = filePrefix[:-7] + fileSuffix
		else:
			outputFname = '%s.markDup%s'%(filePrefix, fileSuffix)
		return outputFname
	
	def addJobs(self, workflow, inputData=None, pegasusFolderName="", tmpDir="/tmp"):
		"""
		2012.3.21
		"""
		sys.stderr.write("Adding MarkDuplicates jobs on %s input datasets ..."%(len(inputData.jobDataLs)))
		returnJobData = PassingData(jobDataLs = [])
		
		topOutputDir = pegasusFolderName
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		no_of_jobs = 1
		for jobData in inputData.jobDataLs:
			inputFile = jobData.output
			
			bamIndexJob = self.addBAMIndexJob(workflow, BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, \
											BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
											inputBamF=inputFile,\
											parentJobLs=[topOutputDirJob]+jobData.jobLs, transferOutput=False)
			
			outputFname = self.getMarkDupOutputFnameBasedOnInputFname(inputFile.abspath)
			
			finalBamFileName = os.path.join(topOutputDir, outputFname)
			finalBamFile = File(finalBamFileName)
			
			markDupJob, markDupBamIndexJob = self.addMarkDupJob(workflow, parentJobLs=[bamIndexJob]+jobData.jobLs, \
						inputBamF=bamIndexJob.bamFile, \
						inputBaiF=bamIndexJob.output, outputBamFile=finalBamFile,\
						MarkDuplicatesJava=workflow.MarkDuplicatesJava, MarkDuplicatesJar=workflow.MarkDuplicatesJar, tmpDir=tmpDir,\
						BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
						transferOutput=True)
			
			no_of_jobs += 3
			returnJobData.jobDataLs.append(PassingData(output=finalBamFile, jobLs=[markDupJob, markDupBamIndexJob]))
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return returnJobData
	
	def run(self):
		"""
		2012.3.21
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		"""
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		#find all hdf5 correlation files
		inputFnameLs = self.getFilesWithProperSuffixFromFolder(self.inputFolder, suffix='.bam')
		inputData = self.registerAllInputFiles(workflow, inputFnameLs=inputFnameLs, input_site_handler=self.input_site_handler, \
								pegasusFolderName=self.pegasusFolderName)
		
		self.addJobs(workflow, inputData=inputData, pegasusFolderName=self.pegasusFolderName, tmpDir=self.tmpDir)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


	
if __name__ == '__main__':
	main_class = MarkDuplicatesWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
