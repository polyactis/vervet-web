#!/usr/bin/env python
"""
Examples:
	%s
	
	# 2011-9-29 (including depth=0 sites, -m 0), all files should stay in "pca" folder (-F pca)
	# clustering jobs in unit of 5
	%s -I ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o workflow/PCAOnVCFWorkflow_VRC105_Top1000Contigs.xml  -P ./PCAOnVCFWorkflow_VRC105_Top1000Contigs.par
		-l condorpool -j condorpool  -u yh -z uclaOffice -C 10 -m 0 -F pca
	
	# 2011.12.16 run on hoffman2 condor  (site depth >=1, -m 1) (turn on checkEmptyVCFByReading, -E)
	%s -I .. -o ... -P ... -l hcondor -j hcondor  -u yh -z uclaOffice
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/  -C 5 -E -m 1
	
	
Description:
	2012.2.27
		a program which generates a pegasus workflow dag (xml file) to 
		
		1. convert each vcf file into EIGENStrat format (3 files: geno, individual, locus)
		2. concatenate the geno and locus EIGENStrat file into one giant file.
		3. run smartpca.pl
			binary "smartpca", "ploteig" (which uses gnuplot, ps2pdf, etc.) must be in PATH
			OR write the parfile by the workflow and run smartpca directly
		4. run twstats
		
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from pymodule.VCFFile import VCFFile

class PCAOnVCFWorkflow(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('smartpca_path', 1, ): ['%s/script/polyactis/EIG3.0/bin/smartpca', '', 1, 'path to smartpca binary', ],\
						('pegasusFolderName', 0, ): ['', 'F', 1, 'the folder relative to pegasus workflow root to contain input & output.\
								It will be created during the pegasus staging process. It is useful to separate multiple workflows.\
								If empty, everything is in the pegasus root.', ],\
						('maxContigID', 1, int): [1000, 'x', 1, 'if contig ID is beyond this number, it will not be included', ],\
						('smartpcaParameterFname', 1, ): ['', 'P', 1, 'file to store the smartpca parameters', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
		self.inputDir = os.path.abspath(self.inputDir)
		self.smartpca_path = self.insertHomePath(self.smartpca_path, self.home_path)
	
	
	def registerCustomExecutables(self, workflow):
		"""
		2012.2.27
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		ConvertVCF2EigenStrat = Executable(namespace=namespace, name="ConvertVCF2EigenStrat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		ConvertVCF2EigenStrat.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ConvertVCF2EigenStrat.py"), site_handler))
		ConvertVCF2EigenStrat.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(ConvertVCF2EigenStrat)
		workflow.ConvertVCF2EigenStrat = ConvertVCF2EigenStrat
		
		
		smartpca = Executable(namespace=namespace, name="smartpca", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		smartpca.addPFN(PFN("file://" + self.smartpca_path, site_handler))
		smartpca.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(smartpca)
		workflow.smartpca = smartpca
	
	
	def addConvertVCF2EigenStratJob(self, workflow, executable=None, inputF=None, outputFnamePrefix=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.3.1
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-i', inputF, "-O", outputFnamePrefix)
		if extraArguments:
			job.addArguments(extraArguments)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		genoOutputF = File('%s.geno'%(outputFnamePrefix))
		locusOutputF = File('%s.snp'%(outputFnamePrefix))
		indOutputF = File('%s.ind'%(outputFnamePrefix))
		outputFLs = [genoOutputF, locusOutputF, indOutputF]
		for outputF in outputFLs:
			job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.genoOutputF = genoOutputF
		job.indOutputF = indOutputF
		job.locusOutputF = locusOutputF
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job

	def addSmartpcaJob(self, workflow, executable=None, smartpcaParameterFile=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					outputFileList=[], \
					job_max_memory=100, **keywords):
		"""
		2012.3.1
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-p', smartpcaParameterFile)
		job.uses(smartpcaParameterFile, transfer=True, register=True, link=Link.INPUT)
		if extraArguments:
			job.addArguments(extraArguments)
		for outputF in outputFileList:
			job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		workflow.addJob(job)
		yh_pegasus.setJobProperRequirement(job, job_max_memory=job_max_memory)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		
		return job
	
	def outputSmartpcaParameters(self, smartpcaParameterFname=None, smartpcaGenotypeInputFile=None, smartpcaLocusInputFile=None,\
						smartpcaIndFile=None, smartpcaCorFile=None, smartpcaEvecFile=None, smartpcaEvalFile=None):
		"""
		2012.3.1
			smartpcaCorFile stores the correlation output of smartpca.
			smartpcaEvecFile stores different principle component values for each individual.
			smartpcaEvalFile stores the eigen value for each principle component.
			numoutevec is the number of principle components to output.
			outliersigmathresh is the sigma outside which an outlier is considered.
			
		"""
		sys.stderr.write("Outputting smartpca parameters ...")
		outf = open(smartpcaParameterFname, 'w')
		parameterList = ["genotypename: %s"%(smartpcaGenotypeInputFile.name),\
						"snpname: %s"%(smartpcaLocusInputFile.name),\
						"indivname: %s"%(smartpcaIndFile.name),\
						"evecoutname: %s"%(smartpcaEvecFile.name),\
						"evaloutname: %s"%(smartpcaEvalFile.name),\
						"cor_outfilename: %s"%(smartpcaCorFile.name),\
						'altnormstyle: NO', 'numoutevec: 40', 'numoutlieriter: 5', 'numoutlierevec: 2', \
						'outliersigmathresh: 15.0', 'qtmode: 0']
		for parameter in parameterList:
			outf.write("%s\n"%parameter)
		del outf
		return smartpcaParameterFname
	
	def addJobs(self, workflow, inputData=None, db_vervet=None,\
			smartpcaParameterFname="", pegasusFolderName="", maxContigID=None):
		"""
		2011.1.8
			add outputDirPrefix to differentiate one run from another if multiple trio call workflows are run simultaneously
			outputDirPrefix could contain "/" to denote sub-folders.
		"""
		
		sys.stderr.write("Adding smartpca jobs on %s VCFs (contig_id<=%s) ..."%(len(inputData.jobDataLs), maxContigID))
		returnJobData = PassingData()
		
		no_of_jobs = 0
		
		topOutputDir = pegasusFolderName
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		smartpcaGenotypeInputFile = File(os.path.join(topOutputDir, 'smartpca.geno'))
		smartpcaLocusInputFile = File(os.path.join(topOutputDir, 'smartpca.snp'))
		
		smartpcaGenotypeMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeFiles, \
						outputF=smartpcaGenotypeInputFile, transferOutput=True, extraArguments='', parentJobLs=[topOutputDirJob])
		smartpcaLocusMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeFiles, \
						outputF=smartpcaLocusInputFile, transferOutput=True, extraArguments='', parentJobLs=[topOutputDirJob])
		
		smartpcaIndFile = None
		smartpcaIndJob = None
		no_of_jobs += 3
		
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			contig_id = self.getContigIDFromFname(inputF.name)
			try:
				if maxContigID:
					contig_id = int(contig_id)
					if contig_id>maxContigID:	#skip the small contigs
						continue
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			outputFnamePrefix = os.path.join(topOutputDir, 'Contig%s'%(contig_id))
			convertJob = self.addConvertVCF2EigenStratJob(workflow, executable=workflow.ConvertVCF2EigenStrat, inputF=inputF, \
								outputFnamePrefix=outputFnamePrefix, parentJobLs=[topOutputDirJob] + jobData.jobLs, \
								extraDependentInputLs=[], transferOutput=False, extraArguments=None, \
								job_max_memory=100)
			if smartpcaIndFile is None:	#every VCF has the same order of individuals
				smartpcaIndFile = convertJob.indOutputF
				smartpcaIndJob = convertJob
			self.addInputToStatMergeJob(workflow, statMergeJob=smartpcaGenotypeMergeJob, \
								inputF=convertJob.genoOutputF, parentJobLs=[convertJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=smartpcaLocusMergeJob, \
								inputF=convertJob.locusOutputF, parentJobLs=[convertJob])
			no_of_jobs += 1
		
		smartpcaCorFile = File(os.path.join(topOutputDir, 'smartpca.cor'))
		smartpcaEvecFile = File(os.path.join(topOutputDir, 'smartpca.evec'))
		smartpcaEvalFile = File(os.path.join(topOutputDir, 'smartpca.eval'))
		self.outputSmartpcaParameters(smartpcaParameterFname=smartpcaParameterFname, smartpcaGenotypeInputFile=smartpcaGenotypeInputFile, \
								smartpcaLocusInputFile=smartpcaLocusInputFile,\
								smartpcaIndFile=smartpcaIndFile, smartpcaCorFile=smartpcaCorFile, \
								smartpcaEvecFile=smartpcaEvecFile, smartpcaEvalFile=smartpcaEvalFile)
		smartpcaParameterFile = self.registerOneInputFile(workflow, smartpcaParameterFname, folderName=pegasusFolderName)
		
		smartpcaJob = self.addSmartpcaJob(workflow, executable=workflow.smartpca, smartpcaParameterFile=smartpcaParameterFile, \
					parentJobLs=[smartpcaGenotypeMergeJob, smartpcaLocusMergeJob, smartpcaIndJob], \
					extraDependentInputLs=[smartpcaGenotypeInputFile, smartpcaLocusInputFile, smartpcaIndFile], \
					transferOutput=True, extraArguments=None, \
					outputFileList=[smartpcaCorFile, smartpcaEvecFile, smartpcaEvalFile], \
					job_max_memory=1000)
		
		no_of_jobs += 1
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return smartpcaJob
	
	def run(self):
		"""
		2012.2.27
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		# Create a abstract dag
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
							checkEmptyVCFByReading=self.checkEmptyVCFByReading, pegasusFolderName=self.pegasusFolderName)
		
		self.addJobs(workflow, inputData, db_vervet=db_vervet, smartpcaParameterFname=self.smartpcaParameterFname, \
					pegasusFolderName=self.pegasusFolderName, maxContigID=self.maxContigID)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


if __name__ == '__main__':
	main_class = PCAOnVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
