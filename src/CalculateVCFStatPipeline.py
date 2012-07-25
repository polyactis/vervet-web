#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s  -a 524 -i ./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/gatk/
		-o VCFStat_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_gatk.xml -l condorpool -j condorpool  -u yh -z uclaOffice
	
	#different window size for pi and TiTv
	%s  -a 524 -i ./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/gatk/
		-o VCFStat_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_gatk_window20Mb.xml -l condorpool -j condorpool  -u yh -z uclaOffice
		-w 20000000
	
	#2012.5.11 run on hoffman2 condor
	%s ...
		-j hcondor -l hcondor -u yh -z localhost  -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
	
Description:
	2011-10-6
		a program which generates a pegasus workflow dag (xml file) that for each vcf file,
		generates these statistics and plots and more 
		
		1. pi (by two ways) => theta
		2. snp density
		3. AAF over position
		4. LD (r^2) over distance
		5. TsTv
		6. heterozygosity
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow

class CalculateVCFStatPipeline(AbstractNGSWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('windowSize', 1, int): [1000000, 'w', 1, 'window size for TiTv, pi, snpDensity calculation by vcftools', ],\
						('vcf1Dir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('includeIndelVCF', 1, int): [0, '', 1, 'toggle this to include indel VCF, filename with "indel"'],\
						})

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractNGSWorkflow.__init__(self, **keywords)
		
		self.vcf1Dir = os.path.abspath(self.vcf1Dir)
	
	def addChrLengthAppendingJob(self, workflow, AddChromosomeLengthToTSVFile=None, inputF=None, outputF=None, \
							divideByLength=True, chrLength=1000, divideStartingColumn=2, parentJobLs=[], \
							extraDependentInputLs=[], transferOutput=False, **keywords):
		job = Job(namespace=workflow.namespace, name=AddChromosomeLengthToTSVFile.name, version=workflow.version)
		job.addArguments('-i', inputF,'-o', outputF, "-l %s"%chrLength, "-s %s"%divideStartingColumn)
		if divideByLength:
			job.addArguments("-v")
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		job.output = outputF
		workflow.addJob(job)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=job)
		for input in extraDependentInputLs:
			job.uses(input, transfer=True, register=True, link=Link.INPUT)
		return job
	
	
	def addTallyAACFromVCFJob(self, workflow, TallyAACFromVCF=None, genomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputF=None, parentJob=None, \
							job_max_memory=1000, extraDependentInputLs=[],\
							input_site_handler=None, transferOutput=False, **keywords):
		# Add a mkdir job for any directory.
		TallyAACFromVCFJob = Job(namespace=workflow.namespace, name=TallyAACFromVCF.name, version=workflow.version)
		refFastaF = refFastaFList[0]
		TallyAACFromVCFJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", genomeAnalysisTKJar, "-R", refFastaF, 
						"-T TallyAACFromVCF", "--variant", inputVCFF, "-o", outputF)
		for refFastaFile in refFastaFList:
			TallyAACFromVCFJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		TallyAACFromVCFJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)
		vcfTabixF = File('%s.tbi'%(inputVCFF.name))	#relative path
		vcfTabixF.absPath = '%s.tbi'%inputVCFF.absPath
		vcfTabixF.addPFN(PFN("file://" + vcfTabixF.absPath, input_site_handler))
		workflow.addFile(vcfTabixF)
		TallyAACFromVCFJob.uses(vcfTabixF, transfer=True, register=True, link=Link.INPUT)
		
		for input in extraDependentInputLs:
			TallyAACFromVCFJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		TallyAACFromVCFJob.uses(outputF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		TallyAACFromVCFJob.output = outputF
		workflow.addJob(TallyAACFromVCFJob)
		if parentJob:
			workflow.depends(parent=parentJob, child=TallyAACFromVCFJob)
		yh_pegasus.setJobProperRequirement(TallyAACFromVCFJob, job_max_memory=job_max_memory)
		return TallyAACFromVCFJob
	
		
	def addVCFToolsJob(self, workflow, executable=None, \
							refFastaFList=None, inputVCFF=None, outputF=None, parentJob=None, \
							job_max_memory=1000, extraDependentInputLs=[],\
							input_site_handler=None, transferOutput=False, **keywords):
		"""
		2011.11.28
			not finished yet
		"""
		vcftoolsJob = Job(namespace=workflow.namespace, name=executable.name, version=workflow.version)
		outputFnamePrefix = os.path.join(statOutputDir, "%s.vcftools"%(commonPrefix))
		piFile = File('%s.windowed.pi'%(outputFnamePrefix))
		TsTvFile = File('%s.TsTv'%(outputFnamePrefix))
		TsTvSummaryFile = File('%s.TsTv.summary'%(outputFnamePrefix))
		snpDensityFile = File('%s.snpden'%(outputFnamePrefix))
		vcftoolsLogFile = File('%s.log'%(outputFnamePrefix))
		vcftoolsJob.addArguments(workflow.vcftoolsPath)
		if vcf1.name[-2:]=='gz':
			vcftoolsJob.addArguments('--gzvcf', vcf1)
		else:
			vcftoolsJob.addArguments('--vcf', vcf1)
		vcftoolsJob.addArguments('--out', outputFnamePrefix, "--depth --het ", \
						" --TsTv %s --window-pi %s --SNPdensity %s"%(self.windowSize, self.windowSize, self.windowSize))
		#"--freq --counts --freq2 --counts2 --hardy",
		vcftoolsJob.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
		vcftoolsJob.uses(vcftoolsLogFile, transfer=True, register=True, link=Link.OUTPUT)	#transfer out the log file
		for outputFile in [piFile, TsTvFile, snpDensityFile, TsTvSummaryFile]:	#no transfer for these
			vcftoolsJob.uses(outputFile, transfer=False, register=True, link=Link.OUTPUT)
		yh_pegasus.setJobProperRequirement(vcftoolsJob, job_max_memory=500)
		workflow.addJob(vcftoolsJob)
		workflow.depends(parent=statOutputDirJob, child=vcftoolsJob)
		return vcftoolsJob
	
	def registerCustomExecutables(self, workflow):
		"""
		2011-11-28
		"""
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		countHomoHetInOneVCF = Executable(namespace=namespace, name="CountHomoHetInOneVCF", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		countHomoHetInOneVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CountHomoHetInOneVCF.py"), site_handler))
		countHomoHetInOneVCF.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(countHomoHetInOneVCF)
		workflow.countHomoHetInOneVCF = countHomoHetInOneVCF
		
		AddChromosomeLengthToTSVFile = Executable(namespace=namespace, name="AddChromosomeLengthToTSVFile", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		AddChromosomeLengthToTSVFile.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mapper/AddChromosomeLengthToTSVFile.py"), site_handler))
		AddChromosomeLengthToTSVFile.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(AddChromosomeLengthToTSVFile)
		workflow.AddChromosomeLengthToTSVFile = AddChromosomeLengthToTSVFile
		
		TallyAACFromVCF = Executable(namespace=namespace, name="TallyAACFromVCF", version=version, os=operatingSystem, arch=architecture, installed=True)
		TallyAACFromVCF.addPFN(PFN("file://" + self.javaPath, site_handler))
		TallyAACFromVCF.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(TallyAACFromVCF)
		workflow.TallyAACFromVCF = TallyAACFromVCF
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
											sequence_type_id=9)
		
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir)
		
		statOutputDir = "statDir"
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=statOutputDir, \
												namespace=workflow.namespace, version=workflow.version)
		
		import re
		chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		#2011.11.16 initiate vervet db connection after genome db connection
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
						input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		homoHetCountFinalOutputF = File('homoHetCountPerSamplePerContig.tsv')
		homoHetCountMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=homoHetCountFinalOutputF)
		
		homoHetCountReduceOutputF = File('homoHetCountPerSample.tsv')
		homoHetCountReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=homoHetCountReduceOutputF, \
							extraArguments='-k 0 -v 2-5,7')	#reduce by sampleId and aggregate only certain columns
		
		TiTvFinalOutputF = File('TiTvWindowSize%s.tsv'%(self.windowSize))
		TiTvMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=TiTvFinalOutputF)
		windowedPiFinalOutputF = File('PiWindowSize%s.tsv'%(self.windowSize))
		windowedPiMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=windowedPiFinalOutputF)
		
		snpDensityOutputF = File('SNPDensityByWindowSize%s.tsv'%(self.windowSize))
		snpDensityMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=snpDensityOutputF)
		
		AACTallyReduceOutputF = File('AAC_tally.tsv')
		AACTallyReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=AACTallyReduceOutputF)
		
		TiTvReduceOutputF = File('TiTv.summary.tsv')
		TiTvReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
						outputF=TiTvReduceOutputF, extraArguments='-k 0 -v 1')
		
		counter = 0
		no_of_vcf_files = 0
		no_of_vcf_non_empty_files = 0
		for inputFname in os.listdir(self.vcf1Dir):
			vcfAbsPath = os.path.join(os.path.abspath(self.vcf1Dir), inputFname)
			no_of_vcf_files += 1
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=self.includeIndelVCF) and \
					not NextGenSeq.isVCFFileEmpty(vcfAbsPath, checkContent=self.checkEmptyVCFByReading):
				no_of_vcf_non_empty_files += 1
				commonPrefix = inputFname.split('.')[0]
				chr = chr_pattern.search(inputFname).group(1)
				chr_size = chr2size.get(chr)
				if chr_size is None:
					sys.stderr.write("size for chr %s is unknown. set it to 1000.\n"%(chr))
					chr_size = 1000
				
				vcf1 = File(os.path.join(vcf1Name, inputFname))	#relative path
				vcf1.addPFN(PFN("file://" + vcfAbsPath, input_site_handler))
				vcf1.absPath = vcfAbsPath
				workflow.addFile(vcf1)
				
				countHomoHetInOneVCFJob = Job(namespace=workflow.namespace, name=workflow.countHomoHetInOneVCF.name, version=workflow.version)
				outputF = File(os.path.join(statOutputDir, "%s.homoHetCountPerSample.tsv"%(os.path.splitext(inputFname)[0])))
				countHomoHetInOneVCFJob.addArguments("-i", vcf1, "-c", chr, "-l %s"%(chr_size), '-o', outputF)
				countHomoHetInOneVCFJob.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
				countHomoHetInOneVCFJob.uses(outputF, transfer=False, register=True, link=Link.OUTPUT)
				countHomoHetInOneVCFJob.output = outputF
				yh_pegasus.setJobProperRequirement(countHomoHetInOneVCFJob, job_max_memory=1000)
				workflow.addJob(countHomoHetInOneVCFJob)
				workflow.depends(parent=statOutputDirJob, child=countHomoHetInOneVCFJob)
				
				self.addInputToStatMergeJob(workflow, statMergeJob=homoHetCountMergeJob, inputF=outputF, \
							parentJobLs=[countHomoHetInOneVCFJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=homoHetCountReduceJob, inputF=countHomoHetInOneVCFJob.output, \
							parentJobLs=[countHomoHetInOneVCFJob])
				
				### vcftools job
				vcftoolsJob = Job(namespace=workflow.namespace, name=workflow.vcftoolsWrapper.name, version=workflow.version)
				outputFnamePrefix = os.path.join(statOutputDir, "%s.vcftools"%(commonPrefix))
				piFile = File('%s.windowed.pi'%(outputFnamePrefix))
				TsTvFile = File('%s.TsTv'%(outputFnamePrefix))
				TsTvSummaryFile = File('%s.TsTv.summary'%(outputFnamePrefix))
				snpDensityFile = File('%s.snpden'%(outputFnamePrefix))
				vcftoolsLogFile = File('%s.log'%(outputFnamePrefix))
				vcftoolsJob.addArguments(workflow.vcftoolsPath)
				if vcf1.name[-2:]=='gz':
					vcftoolsJob.addArguments('--gzvcf', vcf1)
				else:
					vcftoolsJob.addArguments('--vcf', vcf1)
				vcftoolsJob.addArguments('--out', outputFnamePrefix, "--depth --het ", \
								" --TsTv %s --window-pi %s --SNPdensity %s"%(self.windowSize, self.windowSize, self.windowSize))
				#"--freq --counts --freq2 --counts2 --hardy",
				vcftoolsJob.uses(vcf1, transfer=True, register=True, link=Link.INPUT)
				vcftoolsJob.uses(vcftoolsLogFile, transfer=True, register=True, link=Link.OUTPUT)	#transfer out the log file
				for outputFile in [piFile, TsTvFile, snpDensityFile, TsTvSummaryFile]:	#no transfer for these
					vcftoolsJob.uses(outputFile, transfer=False, register=True, link=Link.OUTPUT)
				yh_pegasus.setJobProperRequirement(vcftoolsJob, job_max_memory=500)
				workflow.addJob(vcftoolsJob)
				workflow.depends(parent=statOutputDirJob, child=vcftoolsJob)
				
				addChrLengthToTsTvFileJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=TsTvFile, outputF=File('%s.divByLength'%(TsTvFile.name)), divideByLength=True, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[vcftoolsJob])
				addChrLengthToPiFileJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=piFile, outputF=File('%s.divByLength'%(piFile.name)), divideByLength=True, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[vcftoolsJob])
				addChrLengthToSNPDensityFileJob = self.addChrLengthAppendingJob(workflow, AddChromosomeLengthToTSVFile=workflow.AddChromosomeLengthToTSVFile, \
								inputF=snpDensityFile, outputF=File('%s.divByLength'%(snpDensityFile.name)), divideByLength=True, \
								chrLength=chr_size, divideStartingColumn=2, parentJobLs=[vcftoolsJob])
				
				self.addInputToStatMergeJob(workflow, statMergeJob=TiTvMergeJob, inputF=addChrLengthToTsTvFileJob.output, \
							parentJobLs=[addChrLengthToTsTvFileJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=windowedPiMergeJob, inputF=addChrLengthToPiFileJob.output, \
							parentJobLs=[addChrLengthToPiFileJob])
				self.addInputToStatMergeJob(workflow, statMergeJob=snpDensityMergeJob, inputF=addChrLengthToSNPDensityFileJob.output, \
							parentJobLs=[addChrLengthToSNPDensityFileJob])
				
				self.addInputToStatMergeJob(workflow, statMergeJob=TiTvReduceJob, inputF=TsTvSummaryFile, \
							parentJobLs=[vcftoolsJob])
				
				### TallyAACFromVCF
				outputF = File(os.path.join(statOutputDir, "%s_AAC_tally.tsv"%(commonPrefix)))
				TallyAACFromVCFJob = self.addTallyAACFromVCFJob(workflow, TallyAACFromVCF=workflow.TallyAACFromVCF, \
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf1, outputF=outputF, parentJob=statOutputDirJob, \
							job_max_memory=500, extraDependentInputLs=[],\
							input_site_handler=input_site_handler, transferOutput=False)
				
				self.addInputToStatMergeJob(workflow, statMergeJob=AACTallyReduceJob, inputF=TallyAACFromVCFJob.output, \
							parentJobLs=[TallyAACFromVCFJob])
				counter += 6
		sys.stderr.write("%s jobs from %s non-empty vcf files (%s total files).\n"%(counter+1, \
																	no_of_vcf_non_empty_files, no_of_vcf_files))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
	
if __name__ == '__main__':
	main_class = CalculateVCFStatPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
