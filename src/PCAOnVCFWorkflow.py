#!/usr/bin/env python
"""
Examples:
	%s
	
	# 2011-9-29 (including depth=0 sites, -m 0), all files should stay in "pca" folder (-F pca)
	# clustering jobs in unit of 5
	%s -I ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o dags/PCAOnVCF/PCAOnVCFWorkflow_VRC105_Top1000Contigs.xml  -P ./PCAOnVCFWorkflow_VRC105_Top1000Contigs.par
		-l condorpool -j condorpool  -u yh -z uclaOffice -C 10 -m 0 -F pca
	
	# 2011.12.16 run on hoffman2 condor  (site depth >=1, -m 1) (turn on checkEmptyVCFByReading, -E)
	# add "--missingCallAsRefBase" to assign missing call as reference base
	%s -I .. -o ... -P ... -l hcondor -j hcondor  -u yh -z localhost
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/  -C 5 -E -H -m 1
		## --missingCallAsRefBase
	
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
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from pymodule import VCFFile
from vervet.src import VervetDB, AbstractVervetWorkflow

class PCAOnVCFWorkflow(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('smartpca_path', 1, ): ['%s/script/polyactis/EIG3.0/bin/smartpca', '', 1, 'path to smartpca binary', ],\
						('smartpcaParameterFname', 1, ): ['', 'P', 1, 'file to store the smartpca parameters', ],\
						('missingCallAsRefBase', 0, int):[0, '', 0, 'toggle to regard all missing calls as homozygous reference'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
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
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		#noClusteringExecutableSet = set()	#2012.8.2 you don't want to cluster for some jobs.
		
		
		ConvertVCF2EigenStrat = Executable(namespace=namespace, name="ConvertVCF2EigenStrat", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		ConvertVCF2EigenStrat.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/ConvertVCF2EigenStrat.py"), site_handler))
		executableClusterSizeMultiplierList.append((ConvertVCF2EigenStrat, 1))
		
		
		
		smartpca = Executable(namespace=namespace, name="smartpca", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		smartpca.addPFN(PFN("file://" + self.smartpca_path, site_handler))
		executableClusterSizeMultiplierList.append((smartpca, 0))
		#smartpca.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
	def addConvertVCF2EigenStratJob(self, workflow, executable=None, inputF=None, outputFnamePrefix=None, \
								missingCallAsRefBase=None, \
					parentJobLs=[], extraDependentInputLs=[], transferOutput=True, extraArguments=None, \
					job_max_memory=100, **keywords):
		"""
		2012.9.11 add argument missingCallAsRefBase
		2012.3.1
		"""
		job = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		job.addArguments('-i', inputF, "-O", outputFnamePrefix)
		if missingCallAsRefBase:
			job.addArguments("--missingCallAsRefBase")
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
		if parentJobLs:
			for parentJob in parentJobLs:
				if parentJob:
					workflow.depends(parent=parentJob, child=job)
		if extraDependentInputLs:
			for input in extraDependentInputLs:
				if input:
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
			if outputF:
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
		2012.9.5
			smartpcaCorFile is optional.
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
						'altnormstyle: NO', 'numoutevec: 40', 'numoutlieriter: 5', 'numoutlierevec: 2', \
						'outliersigmathresh: 15.0', 'qtmode: 0']
		if smartpcaCorFile:
			parameterList.append("cor_outfilename: %s"%(smartpcaCorFile.name))
		for parameter in parameterList:
			outf.write("%s\n"%parameter)
		del outf
		return smartpcaParameterFname
	
	def addJobs(self, workflow, inputData=None, db_vervet=None,\
			smartpcaParameterFname="", pegasusFolderName="", maxContigID=None, missingCallAsRefBase=0, transferOutput=True):
		"""
		2012.9.11
			add argument missingCallAsRefBase
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
						outputF=smartpcaGenotypeInputFile, transferOutput=transferOutput, extraArguments='', parentJobLs=[topOutputDirJob])
		smartpcaLocusMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeFiles, \
						outputF=smartpcaLocusInputFile, transferOutput=transferOutput, extraArguments='', parentJobLs=[topOutputDirJob])
		
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
								outputFnamePrefix=outputFnamePrefix, missingCallAsRefBase=missingCallAsRefBase,\
								parentJobLs=[topOutputDirJob] + jobData.jobLs, \
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
		
		#smartpcaCorFile = File(os.path.join(topOutputDir, 'smartpca.cor'))
		smartpcaEvecFile = File(os.path.join(topOutputDir, 'smartpca.evec'))
		smartpcaEvalFile = File(os.path.join(topOutputDir, 'smartpca.eval'))
		self.outputSmartpcaParameters(smartpcaParameterFname=smartpcaParameterFname, smartpcaGenotypeInputFile=smartpcaGenotypeInputFile, \
								smartpcaLocusInputFile=smartpcaLocusInputFile,\
								smartpcaIndFile=smartpcaIndFile, smartpcaCorFile=None, \
								smartpcaEvecFile=smartpcaEvecFile, smartpcaEvalFile=smartpcaEvalFile)
		smartpcaParameterFile = self.registerOneInputFile(workflow, smartpcaParameterFname, folderName=pegasusFolderName)
		
		smartpcaJob = self.addSmartpcaJob(workflow, executable=workflow.smartpca, smartpcaParameterFile=smartpcaParameterFile, \
					parentJobLs=[smartpcaGenotypeMergeJob, smartpcaLocusMergeJob, smartpcaIndJob], \
					extraDependentInputLs=[smartpcaGenotypeInputFile, smartpcaLocusInputFile, smartpcaIndFile], \
					transferOutput=transferOutput, extraArguments=None, \
					outputFileList=[None, smartpcaEvecFile, smartpcaEvalFile], \
					job_max_memory=18000)
		
		#2012.9.5 add the job to append meta info (country, sex, latitude, etc. of each monkey)
		outputF = File(os.path.join(topOutputDir, 'smartpca_evec_withMetaInfo.tsv'))
		appendInfo2SmartPCAOutputJob = self.addGenericJob(executable=self.AppendInfo2SmartPCAOutput, inputFile=smartpcaEvecFile, \
				outputFile=outputF, \
				parentJobLs=[smartpcaJob], extraDependentInputLs=None, \
				extraOutputLs=None,\
				transferOutput=transferOutput, \
				extraArgumentList=None, extraArguments='--inversePCValue', key2ObjectForJob=None, job_max_memory=2000)
		self.addDBArgumentsToOneJob(job=appendInfo2SmartPCAOutputJob, objectWithDBArguments=self)
		
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
		
		db_vervet = self.db_vervet
		
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
					pegasusFolderName=self.pegasusFolderName, maxContigID=self.maxContigID,\
					missingCallAsRefBase=self.missingCallAsRefBase)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)


if __name__ == '__main__':
	main_class = PCAOnVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
