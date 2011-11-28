#!/usr/bin/env python
"""
Examples:
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor_20111106T1554/;

	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor.xml  -z uclaOffice -u yh -q alnStatForFilter.tsv
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;
	
	%s ~/script/vervet/src/FilterVCFPipeline.py  -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor.xml -z uclaOffice -u yh
		-q alnStatForFilter.tsv
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;
	
	%s ~/script/vervet/src/FilterVCFPipeline.py  -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_minGQ30_multi_sample_condor.xml -z uclaOffice -u yh
		-q alnStatForFilter.tsv	-G 30

Description:
	2011-11-7 pipeline that filters VCF by depth, GQ, MAC, SNP missing rate, mismatch rate between two input VCFs
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from AlignmentToCallPipeline import AlignmentToCallPipeline

class FilterVCFPipeline(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('ref_ind_seq_id', 1, int): [524, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype', ],\
						('depthFoldChange', 1, float): [2.0, 'F', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						("maxSNPMismatchRate", 0, float): [0, 'N', 1, 'maximum SNP mismatch rate between two vcf calls'],\
						("minMAC", 0, int): [2, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
						("minMAF", 0, int): [None, 'm', 1, 'minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRate", 0, float): [0, 'x', 1, 'minimum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('vcf1Dir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", 'S', 1, 'vervet source code folder'],\
						("picard_path", 1, ): ["%s/script/picard/dist", '', 1, 'picard folder containing its jar binaries'],\
						("gatk_path", 1, ): ["%s/script/gatk/dist", '', 1, 'GATK folder containing its jar binaries'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("javaPath", 1, ): ["/usr/bin/java", 'J', 1, 'java interpreter binary'],\
						('alnStatForFilterFname', 1, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
						("localDataDir", 0, ): ["", 'D', 1, 'localDataDir should contain same files as dataDir but accessible locally.\
									If not given, use the default stored in db. This argument is used to find all input files available.'],\
						
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.picard_path = self.picard_path%self.home_path
		self.gatk_path = self.gatk_path%self.home_path
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
		
	def outputAlignmentDepthAndOthersForFilter(self, outputFname, ref_ind_seq_id=524, foldChange=2, minGQ=30):
		"""
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
			writer.writerow([row.getReadGroup(), row.median_depth/float(foldChange), \
							row.median_depth*float(foldChange), minGQ])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
	
	def addFilterVCFByDepthJob(self, workflow, FilterVCFByDepthJava=None, genomeAnalysisTKJar=None, \
							refFastaFList=None, inputVCFF=None, outputVCFF=None, parentJob=None, alnStatForFilterF=None, \
							namespace=None, version=None, job_max_memory=1000, extraDependentInputLs=[]):
		# Add a mkdir job for any directory.
		filterByDepthJob = Job(namespace=namespace, name=FilterVCFByDepthJava.name, version=version)
		refFastaF = refFastaFList[0]
		filterByDepthJob.addArguments("-Xmx%sm"%(job_max_memory), "-jar", genomeAnalysisTKJar, "-R", refFastaF, "-T FilterVCFByDepth", \
							"--variant", inputVCFF, "-o", outputVCFF, "-depthFile", alnStatForFilterF)
		for refFastaFile in refFastaFList:
			filterByDepthJob.uses(refFastaFile, transfer=True, register=True, link=Link.INPUT)
		filterByDepthJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)		
		for input in extraDependentInputLs:
			filterByDepthJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		filterByDepthJob.uses(alnStatForFilterF, transfer=True, register=True, link=Link.INPUT)
		filterByDepthJob.uses(outputVCFF, transfer=False, register=True, link=Link.OUTPUT)
		filterByDepthJob.output = outputVCFF
		workflow.addJob(filterByDepthJob)
		if parentJob:
			workflow.depends(parent=parentJob, child=filterByDepthJob)
		yh_pegasus.setJobProperRequirement(filterByDepthJob, job_max_memory=job_max_memory)
		return filterByDepthJob
	
	def addFilterJobByvcftools(self, workflow, vcftoolsWrapper=None, inputVCFF=None, outputFnamePrefix=None, \
							parentJobLs=None, snpMisMatchStatFile=None, minMAC=2, minMAF=None, maxSNPMissingRate=0.9,\
							namespace=None, version=None, extraDependentInputLs=[], transferOutput=False):
		"""
		2011-11-21
			argument vcftools is replaced with a wrapper, which takes vcftools path as 1st argument
		"""
		# Add a mkdir job for any directory.
		vcftoolsJob = Job(namespace=namespace, name=vcftoolsWrapper.name, version=version)
		vcftoolsJob.addArguments(vcftoolsWrapper.vcftoolsPath)	#2011-11-21
		if inputVCFF.name[-2:]=='gz':
			vcftoolsJob.addArguments("--gzvcf", inputVCFF)
		else:
			vcftoolsJob.addArguments("--vcf", inputVCFF)
		vcftoolsJob.addArguments("--out", outputFnamePrefix, "--recode", "--positions", snpMisMatchStatFile)

		if maxSNPMissingRate is not None:
			vcftoolsJob.addArguments("--geno %s"%(1-maxSNPMissingRate))
		if minMAF is not None:
			vcftoolsJob.addArguments("--maf %s"%(minMAF))
		if minMAC is not None:
			vcftoolsJob.addArguments("--mac %s"%(minMAC))
		
		vcftoolsJob.uses(inputVCFF, transfer=True, register=True, link=Link.INPUT)
		vcftoolsJob.uses(snpMisMatchStatFile, transfer=True, register=True, link=Link.INPUT)
		for input in extraDependentInputLs:
			vcftoolsJob.uses(input, transfer=True, register=True, link=Link.INPUT)
		outputVCFF = File("%s.recode.vcf"%(outputFnamePrefix))
		vcftoolsJob.uses(outputVCFF, transfer=transferOutput, register=True, link=Link.OUTPUT)
		vcftoolsJob.output = outputVCFF
		workflow.addJob(vcftoolsJob)
		for parentJob in parentJobLs:
			workflow.depends(parent=parentJob, child=vcftoolsJob)
		return vcftoolsJob
	
	def registerVCFAndItsTabixIndex(self, workflow, vcfF, input_site_handler='local'):
		"""
		2011-11-11
			vcfF.absPath is path to its physical path
			register both vcf and its tabix
		"""
		vcfF.addPFN(PFN("file://" + vcfF.absPath, input_site_handler))
		workflow.addFile(vcfF)
		vcfF.tbi_F = File("%s.tbi"%vcfF.name)
		vcfF.tbi_F.addPFN(PFN("file://" + "%s.tbi"%vcfF.absPath, input_site_handler))
		workflow.addFile(vcfF.tbi_F)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		"""
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = ADAG(workflowName)
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		clusters_size = 20
		
		abs_path = os.path.join(self.gatk_path, 'GenomeAnalysisTK.jar')
		genomeAnalysisTKJar = File(abs_path)
		genomeAnalysisTKJar.addPFN(PFN("file://" + abs_path, site_handler))
		workflow.addFile(genomeAnalysisTKJar)
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mkdirWrap)
		
		"""
		vcftools = Executable(namespace=namespace, name="vcftools", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcftools.addPFN(PFN("file://" + os.path.join(self.home_path, 'bin/vcftools/vcftools'), site_handler))
		vcftools.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(vcftools)
		"""
		
		vcftoolsPath = os.path.join(self.home_path, "bin/vcftools/vcftools")
		vcftoolsWrapper = Executable(namespace=namespace, name="vcftoolsWrapper", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		vcftoolsWrapper.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "shell/vcftoolsWrapper.sh"), site_handler))
		vcftoolsWrapper.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		vcftoolsWrapper.vcftoolsPath = vcftoolsPath
		workflow.addExecutable(vcftoolsWrapper)
		
		
		SelectVariantsJava = Executable(namespace=namespace, name="SelectVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		SelectVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		SelectVariantsJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(SelectVariantsJava)
		
		CombineVariantsJava = Executable(namespace=namespace, name="CombineVariantsJava", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		CombineVariantsJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		CombineVariantsJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CombineVariantsJava)
		
		FilterVCFByDepthJava = Executable(namespace=namespace, name="FilterVCFByDepth", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		FilterVCFByDepthJava.addPFN(PFN("file://" + self.javaPath, site_handler))
		FilterVCFByDepthJava.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(FilterVCFByDepthJava)
		
		CalculateSNPMismatchRateOfTwoVCF = Executable(namespace=namespace, name="CalculateSNPMismatchRateOfTwoVCF", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		CalculateSNPMismatchRateOfTwoVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CalculateSNPMismatchRateOfTwoVCF.py"), site_handler))
		CalculateSNPMismatchRateOfTwoVCF.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CalculateSNPMismatchRateOfTwoVCF)
		
		bgzip_tabix = Executable(namespace=namespace, name="bgzip_tabix", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		bgzip_tabix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "bgzip_tabix.sh"), site_handler))
		bgzip_tabix.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(bgzip_tabix)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.outputAlignmentDepthAndOthersForFilter(self.alnStatForFilterFname, ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
		alnStatForFilterF = File(os.path.basename(self.alnStatForFilterFname))
		alnStatForFilterF.addPFN(PFN("file://" + os.path.abspath(self.alnStatForFilterFname), \
											self.input_site_handler))
		workflow.addFile(alnStatForFilterF)
		
		#name to distinguish between vcf1Dir, and vcf2Dir
		vcf1Name = os.path.basename(os.path.dirname(self.vcf1Dir))	#dirname is used because it removes the trailing "/" which makes basename return ""
		if not vcf1Name:
			vcf1Name = "vcf1"
		vcf2Name = os.path.basename(os.path.dirname(self.vcf2Dir))
		if vcf2Name==vcf1Name or not vcf2Name:
			vcf2Name = "vcf2"
		
		vcf1DepthFilterDir = "%s_DepthFilter"%(vcf1Name)
		vcf1DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=vcf1DepthFilterDir, namespace=namespace, version=version)
		vcf2DepthFilterDir = "%s_DepthFilter"%(vcf2Name)
		vcf2DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=vcf2DepthFilterDir, namespace=namespace, version=version)
		
		SNPMismatchStatDir = "SNPMismatchStat"
		SNPMismatchStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=SNPMismatchStatDir, namespace=namespace, version=version)
		
		vcf1_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf1Name)
		vcf1_vcftoolsFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=vcf1_vcftoolsFilterDir, namespace=namespace, version=version)
		vcf2_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf2Name)
		vcf2_vcftoolsFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=vcf2_vcftoolsFilterDir, namespace=namespace, version=version)
		
		
		#import re
		#chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		counter = 0
		for inputFname in os.listdir(self.vcf1Dir):
			vcf1AbsPath = os.path.join(os.path.abspath(self.vcf1Dir), inputFname)
			vcf2AbsPath = os.path.join(os.path.abspath(self.vcf2Dir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(vcf1AbsPath):
				if not NextGenSeq.isVCFFileEmpty(vcf2AbsPath):	#make sure the samtools vcf exists
					#chr = chr_pattern.search(inputFname).group(1)
					commonPrefix = inputFname.split('.')[0]
					vcf1 = File(os.path.join(vcf1Name, inputFname))	#relative path
					vcf1.absPath = vcf1AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf1, input_site_handler)
					
					vcf1AfterDepthFilter = File(os.path.join(vcf1DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=FilterVCFByDepthJava, \
							genomeAnalysisTKJar=genomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf1, outputVCFF=vcf1AfterDepthFilter, parentJob=vcf1DepthFilterDirJob, \
							alnStatForFilterF=alnStatForFilterF, \
							namespace=namespace, version=version, extraDependentInputLs=[vcf1.tbi_F])
				
					vcf1AfterDepthFilterGzip = File("%s.gz"%vcf1AfterDepthFilter.name)
					vcf1AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf1AfterDepthFilter.name)
					vcf1FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
							parentJob=vcf1FilterByDepthJob, inputF=vcf1AfterDepthFilter, outputF=vcf1AfterDepthFilterGzip, \
							namespace=namespace, version=version, transferOutput=False)
					
					vcf2 = File(os.path.join(vcf2Name, inputFname))	#relative path
					vcf2.absPath = vcf2AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf2, input_site_handler)
					
					vcf2AfterDepthFilter = File(os.path.join(vcf2DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf2FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=FilterVCFByDepthJava, \
							genomeAnalysisTKJar=genomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf2, outputVCFF=vcf2AfterDepthFilter, parentJob=vcf2DepthFilterDirJob, \
							alnStatForFilterF=alnStatForFilterF, \
							namespace=namespace, version=version, extraDependentInputLs=[vcf2.tbi_F])
					
					vcf2AfterDepthFilterGzip = File("%s.gz"%vcf2AfterDepthFilter.name)
					vcf2AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf2AfterDepthFilter.name)
					vcf2FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
							parentJob=vcf2FilterByDepthJob, inputF=vcf2AfterDepthFilter, outputF=vcf2AfterDepthFilterGzip, \
							namespace=namespace, version=version, transferOutput=False)
					
					snpMisMatchStatFile = File(os.path.join(SNPMismatchStatDir, '%s_snpMismatchStat.tsv'%(os.path.splitext(commonPrefix)[0])))
					calculateSNPMismatchRateOfTwoVCFJob = Job(namespace=namespace, name=CalculateSNPMismatchRateOfTwoVCF.name, version=version)
					calculateSNPMismatchRateOfTwoVCFJob.addArguments("-i", vcf1AfterDepthFilterGzip, "-j", vcf2AfterDepthFilterGzip, \
									"-m %s"%self.maxSNPMismatchRate, '-o', snpMisMatchStatFile)
					calculateSNPMismatchRateOfTwoVCFJob.uses(vcf1AfterDepthFilterGzip, transfer=True, register=True, link=Link.INPUT)
					calculateSNPMismatchRateOfTwoVCFJob.uses(vcf2AfterDepthFilterGzip, transfer=True, register=True, link=Link.INPUT)
					calculateSNPMismatchRateOfTwoVCFJob.uses(snpMisMatchStatFile, transfer=True, register=True, link=Link.OUTPUT)
					yh_pegasus.setJobProperRequirement(calculateSNPMismatchRateOfTwoVCFJob, job_max_memory=1000)
					workflow.addJob(calculateSNPMismatchRateOfTwoVCFJob)
					workflow.depends(parent=vcf1FilterByDepthBGZipTabixJob, child=calculateSNPMismatchRateOfTwoVCFJob)
					workflow.depends(parent=vcf2FilterByDepthBGZipTabixJob, child=calculateSNPMismatchRateOfTwoVCFJob)
					workflow.depends(parent=SNPMismatchStatDirJob, child=calculateSNPMismatchRateOfTwoVCFJob)
					
					
					outputFnamePrefix = os.path.join(vcf1_vcftoolsFilterDir, '%s.filter_by_vcftools'%(commonPrefix))
					vcf1FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftools=vcftools, inputVCFF=vcf1AfterDepthFilterGzip, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[vcf1FilterByDepthBGZipTabixJob, vcf1_vcftoolsFilterDirJob, calculateSNPMismatchRateOfTwoVCFJob], \
							snpMisMatchStatFile=snpMisMatchStatFile, \
							minMAC=self.minMAC, minMAF=self.minMAF, \
							maxSNPMissingRate=self.maxSNPMissingRate,\
							namespace=namespace, version=version, extraDependentInputLs=[vcf1FilterByDepthBGZipTabixJob.tbi_F])
					vcf1FilterByvcftoolsGzip = File("%s.gz"%vcf1FilterByvcftoolsJob.output.name)
					vcf1FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
							parentJob=vcf1FilterByvcftoolsJob, inputF=vcf1FilterByvcftoolsJob.output, outputF=vcf1FilterByvcftoolsGzip, \
							namespace=namespace, version=version, transferOutput=True)
					
					outputFnamePrefix = os.path.join(vcf2_vcftoolsFilterDir, '%s.filter_by_vcftools'%(commonPrefix))
					vcf2FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftools=vcftools, inputVCFF=vcf2AfterDepthFilterGzip, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[vcf2FilterByDepthBGZipTabixJob, vcf2_vcftoolsFilterDirJob, calculateSNPMismatchRateOfTwoVCFJob],
							snpMisMatchStatFile=snpMisMatchStatFile, \
							minMAC=self.minMAC, minMAF=self.minMAF, \
							maxSNPMissingRate=self.maxSNPMissingRate,\
							namespace=namespace, version=version, extraDependentInputLs=[vcf2FilterByDepthBGZipTabixJob.tbi_F])
					
					vcf2FilterByvcftoolsGzip = File("%s.gz"%vcf2FilterByvcftoolsJob.output.name)
					vcf2FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=bgzip_tabix, \
							parentJob=vcf2FilterByvcftoolsJob, inputF=vcf2FilterByvcftoolsJob.output, outputF=vcf2FilterByvcftoolsGzip, \
							namespace=namespace, version=version, transferOutput=True)
					
					counter += 9
		sys.stderr.write("%s jobs.\n"%(counter+1))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = FilterVCFPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
