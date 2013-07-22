#!/usr/bin/env python
"""
Examples:
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor_20111106T1554/;

	%s -I $dirPrefix\gatk/ -i $dirPrefix\samtools/ -l condorpool -j condorpool
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor.xml  -z uclaOffice -u yh
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;
	
	%s -I $dirPrefix\gatk/ -i $dirPrefix\samtools/ -l condorpool -j condorpool -a 524
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor.xml -z uclaOffice -u yh
		#--alnStatForFilterFname alnStatForFilter.tsv
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;	
	%s  -I $dirPrefix\gatk/ -i $dirPrefix\samtools/ -l condorpool -j condorpool -a 524
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_minGQ30_multi_sample_condor.xml -z uclaOffice -u yh
		--minGQ 30
	
	#2011.11.28 (--onlyKeepBiAllelicSNP = keep only bi-allelic SNPs), --checkEmptyVCFByReading
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s  -I $dirPrefix\gatk -i $dirPrefix\call/ -l condorpool -j condorpool
		-o FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25.xml
		-z uclaOffice -u yh
		--minGQ 1 --onlyKeepBiAllelicSNP --maxSNPMismatchRate 0.1
		--minMAC 5 --maxSNPMissingRate 0.25 -a 524 -C 50 --checkEmptyVCFByReading
	
	#2011.12.9 to remove SNPs that are not in a file. no other filters.
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s -I $dirPrefix\gatk -i $dirPrefix\call/ -l condorpool -j condorpool
		-o Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4.xml
		-z uclaOfficeTemp -u yh --minGQ 0 --maxSNPMismatchRate 1 --minMAC 0
		--maxSNPMissingRate 1 --depthFoldChange 100000 -a 524 -C 50
		--keepSNPPosFname ./4HighCovVRC_inter_minMAC4_vs_LowPass_top7559Contigs_no12eVarFilter_inter.2011.12.9T0107/overlapPos.tsv
	
	#2013.07.04 remove contaminants in method 101
	%s -I ~/NetworkData/vervet/db/genotype_file/method_101/
		--ref_ind_seq_id 3280 --checkEmptyVCFByReading -H
		-o dags/FilterVariants/FilterGenotypeMethod101_RemoveContaminant.xml
		-l hcondor -j hcondor
		-t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/
		--intervalOverlapSize 0 --intervalSize 10000
		-u yh -C 10 --contigMaxRankBySize 2000 --siteTotalCoverageINFOFieldName DP --is_contaminated 0
		
	#2012.5.1 filter trioCaller output with total depth (no minGQ filter anymore), minMAC=10 (--minMAC 10),
	# maxSNPMismatchRate=1 (--maxSNPMissingRate 1.0)
	# minMAF=0.05 (--minMAF 0.05), no depth-band filter (--depthFoldChange 0)
	%s -I AlignmentToTrioCallPipeline_VRC_top7559Contigs.2011.12.15T0057/trioCaller/ -l condorpool -j condorpool
		-z uclaOffice -u yh --minMAC 10 --maxSNPMissingRate 1.0
		--minMAF 0.05 --depthFoldChange 0 -a 524 -C 50 --checkEmptyVCFByReading -H
		-o FilterVCF_trioCallerTop7559Contigs.xml
	
	#2012.8.1 FilterGenotypeMethod5_ByMethod7Sites (--keepSNPPosFname ...) NoDepthFilter (--depthFoldChange 0) MaxSNPMissing0.5 (--maxSNPMissingRate 0.5)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_5/
		--maxSNPMissingRate 0.5 -a 524  --checkEmptyVCFByReading -H
		-o dags/FilterVariants/FilterGenotypeMethod5_ByMethod7Sites_NoDepthFilter_MaxSNPMissing0.5.xml
		-l hcondor -j hcondor -t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/
		-u yh -C 5 --keepSNPPosFname ./method7_sites.tsv --depthFoldChange 0
	
	#2012.8.1 FilterGenotypeMethod6_ByMaskingZeroDPSite (--minDepth 1) 2FoldDepthFilter (--depthFoldChange 2) MaxSNPMissing1.0 (--maxSNPMissingRate 1.0)
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	# "--minNeighborDistance 5" to the minimum distance between neighboring SNPs, "--minMAF 0.1", minMAF=0.1
	%s -I ~/NetworkData/vervet/db/genotype_file/method_6/
		--maxSNPMissingRate 1.0 -a 524 --checkEmptyVCFByReading -H
		-o dags/FilterVariants/FilterGenotypeMethod6_MinDP1_2FoldDepthFilter_MaxSNPMissing1.0MinNeighborDistance5MinMAF0.1.xml
		-l hcondor -j hcondor
		-t ~/NetworkData/vervet/db/ -D ~/NetworkData/vervet/db/
		-u yh -C 5 --minDepth 1 --depthFoldChange 2 --minNeighborDistance 5 --minMAF 0.1
		#-V 90 -x 100 
		--excludeFilteredSites 2
		--siteTotalCoverageINFOFieldName DP
	
Description:
	2012.9.12 pipeline that runs VCF2plink, plink Mendel, then filter VCF by max mendel error on top of filters by depth, GQ, MAC, 
		SNP missing rate.
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0]\
				, sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, PassingData, yh_pegasus, utils, NextGenSeq
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
parentClass = AbstractVervetWorkflow
class FilterVCFPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	common_filter_option_dict = {
						("onlyKeepBiAllelicSNP", 0, int): [0, 'K', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('keepSNPPosFname', 0, ): ['', '', 1, 'a tab-delimited file with 1st two columns as chr and pos.\
		should have a header.\
		to filter out SNPs that are absent in this file.\
		this step will be applied before all other filter jobs.\
		extra columns are fine.'],\
						}
	option_default_dict.update(common_filter_option_dict)
	option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype. 2012.5.1 no longer enforced in FilterVCFByDepth.java', ],\
						('depthFoldChange', 0, float): [0, '', 1, 'a variant is retained if its depth within this fold change of meanDepth,\
		set this to 0 or below to eliminate this step of filtering.', ],\
						("maxSNPMismatchRate", 0, float): [0, '', 1, 'maximum SNP mismatch rate between two vcf calls'],\
						("minMAC", 0, int): [None, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
						("minMAF", 0, float): [None, 'f', 1, 'minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('minNeighborDistance', 0, int): [None, 'g', 1, 'minimum distance between two adjacent SNPs'],\
						('excludeFilteredSites', 0, int): [0, '', 1, '0: no such filter, 1: remove sites whose FILTER!=PASS, 2: remove sites whose FILTER!=PASS and is not a SNP (indels+MNP)', ],\
						("is_contaminated", 0, int): [None, '', 1, 'if not None, this will trigger a step-0 filter:\n\
	choose samples whose individual_sequence.is_contaminated is equal as this argument value (0=non-contaminant, 1=contaminant)'],\
						('vcf1Dir', 1, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files. If not provided, filter vcf1Dir without applying maxSNPMismatchRate filter.', ],\
						('siteTotalCoverageINFOFieldName', 1, ): ['DP', '', 1, 'used in the depthFoldChange filter step,  by GATK SelectVariants to parse the depth of entire site.\n\
		SAMtools, GATK output uses DP, Platypus output uses TC', ],\
						})
	#set no overlap between adjacent intervals 
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	#('alnStatForFilterFname', 0, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
	#	2013.06.13 alnStatForFilterFname is no longer used. 
	#("minDepthPerGenotype", 0, int): [0, 'Z', 1, 'mask genotype with below this depth as ./. (other fields retained), \
	#	esp. necessary for SAMtools, which output homozygous reference if no read for one sample.'],\
	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
		# 2012.8.3 relative path causes stage-in problem as the pegasus program can't find the files.
		if getattr(self, 'vcf1Dir', None):
			self.vcf1Dir = os.path.abspath(self.vcf1Dir)
		if getattr(self, 'vcf2Dir', None):
			self.vcf2Dir = os.path.abspath(self.vcf2Dir)
		"""
		if getattr(self, 'depthFoldChange', None) and self.depthFoldChange>0 and not self.alnStatForFilterFname:
			sys.stderr.write("Error: alnStatForFilterFname (%s) is nothing while depthFoldChange=%s.\n"%\
							(self.alnStatForFilterFname, self.depthFoldChange))
			sys.exit(3)
		"""
		self.minDepthPerGenotype = self.minDepth
		
		#2013.07.18 offer child classes option to turn it off
		self.needGzipPreReduceReturnData = False	#all stat data in preReduce() are marked with transferOutput=True
		self.needGzipReduceReturnData = False	#no need to gzip reduce return data as they are already marked as transferOutput=True and bgzipped in reduceEachVCF() 
	
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
		
		#in order for two different input's FilterVCFByDepth to be merged into different-named clustered jobs
		FilterVCFByDepth2Java = Executable(namespace=namespace, name="FilterVCFByDepth2", version=version, os=operatingSystem,\
											arch=architecture, installed=True)
		FilterVCFByDepth2Java.addPFN(PFN("file://" + self.javaPath, site_handler))
		FilterVCFByDepth2Java.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(FilterVCFByDepth2Java)
		workflow.FilterVCFByDepth2Java = FilterVCFByDepth2Java
		
		CalculateSNPMismatchRateOfTwoVCF = Executable(namespace=namespace, name="CalculateSNPMismatchRateOfTwoVCF", \
							version=version, os=operatingSystem, arch=architecture, installed=True)
		CalculateSNPMismatchRateOfTwoVCF.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CalculateSNPMismatchRateOfTwoVCF.py"), site_handler))
		CalculateSNPMismatchRateOfTwoVCF.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CalculateSNPMismatchRateOfTwoVCF)
		workflow.CalculateSNPMismatchRateOfTwoVCF = CalculateSNPMismatchRateOfTwoVCF
		
		#2013.05.20
		self.setOrChangeExecutableClusterSize(executable=self.SelectVariantsJava, clusterSizeMultipler=1.0, \
									defaultClustersSize=self.clusters_size)
	
	def addJobsToFilterTwoVCFDir(self, workflow=None, vcf1Dir=None, vcf2Dir=None, registerReferenceData=None, \
							alnStatForFilterF=None, keepSNPPosF=None, onlyKeepBiAllelicSNP=True,\
							minMAC=None, minMAF=None, maxSNPMissingRate=None, maxSNPMismatchRate=None):
		"""
		2013.04.08 added argument maxSNPMismatchRate
		2012.5.10
			add extraArguments="--recode-INFO-all" to addFilterJobByvcftools()
		2012.1.14
		"""
		if workflow is None:
			workflow = self
		refFastaFList = registerReferenceData.refFastaFList
		refFastaF = refFastaFList[0]
		
		#name to distinguish between vcf1Dir, and vcf2Dir
		vcf1Name = self.findProperVCFDirIdentifier(vcf1Dir, defaultName='vcf1')
		vcf2Name = self.findProperVCFDirIdentifier(vcf2Dir, defaultName='vcf2')
		if vcf2Name==vcf1Name or not vcf2Name:
			vcf2Name = "vcf2"
		
		vcf1DepthFilterDir = "%s_DepthFilter"%(vcf1Name)
		vcf1DepthFilterDirJob = self.addMkDirJob(outputDir=vcf1DepthFilterDir)
		vcf2DepthFilterDir = "%s_DepthFilter"%(vcf2Name)
		vcf2DepthFilterDirJob = self.addMkDirJob(outputDir=vcf2DepthFilterDir)
		
		SNPMismatchStatDir = "SNPMismatchStat"
		SNPMismatchStatDirJob = self.addMkDirJob(outputDir=SNPMismatchStatDir)
		
		vcf1_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf1Name)
		vcf1_vcftoolsFilterDirJob = self.addMkDirJob(outputDir=vcf1_vcftoolsFilterDir)
		vcf2_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf2Name)
		vcf2_vcftoolsFilterDirJob = self.addMkDirJob(outputDir=vcf2_vcftoolsFilterDir)
		
		#import re
		#chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		counter = 0
		for inputFname in os.listdir(vcf1Dir):
			vcf1AbsPath = os.path.join(os.path.abspath(vcf1Dir), inputFname)
			vcf2AbsPath = os.path.join(os.path.abspath(vcf2Dir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(vcf1AbsPath):
				if not NextGenSeq.isVCFFileEmpty(vcf2AbsPath, checkContent=self.checkEmptyVCFByReading):	#make sure the samtools vcf exists
					#chr = chr_pattern.search(inputFname).group(1)
					commonPrefix = inputFname.split('.')[0]
					vcf1 = File(os.path.join(vcf1Name, inputFname))	#relative path
					vcf1.absPath = vcf1AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf1, input_site_handler)
					vcf2 = File(os.path.join(vcf2Name, inputFname))	#relative path
					vcf2.absPath = vcf2AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf2, input_site_handler)
					
					if keepSNPPosF:
						#toss SNPs that are not in this keepSNPPosFname file
						outputFnamePrefix = os.path.join(vcf1_vcftoolsFilterDir, '%s.keepGivenSNP'%(commonPrefix))
						vcf1KeepGivenSNPByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
								inputVCFF=vcf1, \
								outputFnamePrefix=outputFnamePrefix, \
								parentJobLs=[vcf1_vcftoolsFilterDirJob], \
								snpMisMatchStatFile=keepSNPPosF, \
								minMAC=None, minMAF=None, \
								maxSNPMissingRate=None,\
								extraDependentInputLs=[vcf1.tbi_F], extraArguments="--recode-INFO-all")
						
						vcf1KeepGivenSNPByvcftoolsGzip = File("%s.gz"%vcf1KeepGivenSNPByvcftoolsJob.output.name)
						vcf1KeepGivenSNPByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf1KeepGivenSNPByvcftoolsJob, inputF=vcf1KeepGivenSNPByvcftoolsJob.output, \
							outputF=vcf1KeepGivenSNPByvcftoolsGzip, \
							transferOutput=True)
					
						outputFnamePrefix = os.path.join(vcf2_vcftoolsFilterDir, '%s.keepGivenSNP'%(commonPrefix))
						vcf2KeepGivenSNPByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
								inputVCFF=vcf2, \
								outputFnamePrefix=outputFnamePrefix, \
								parentJobLs=[vcf2_vcftoolsFilterDirJob],
								snpMisMatchStatFile=keepSNPPosF, \
								minMAC=None, minMAF=None, \
								maxSNPMissingRate=None,\
								extraDependentInputLs=[vcf2.tbi_F], extraArguments="--recode-INFO-all")
						vcf2KeepGivenSNPByvcftoolsGzip = File("%s.gz"%vcf2KeepGivenSNPByvcftoolsJob.output.name)
						vcf2KeepGivenSNPByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf2KeepGivenSNPByvcftoolsJob, inputF=vcf2KeepGivenSNPByvcftoolsJob.output, \
							outputF=vcf2KeepGivenSNPByvcftoolsGzip, \
							transferOutput=True)
					
						vcf1filterByDepthInput=vcf1KeepGivenSNPByvcftoolsJob.output
						lastRoundJobLs=[vcf1DepthFilterDirJob, vcf1KeepGivenSNPByvcftoolsJob]
						lastRoundExtraDependentInputLs=[]
						vcf2filterByDepthInput=vcf2KeepGivenSNPByvcftoolsJob.output
						vcf2filterByDepthParentJobLs=[vcf2DepthFilterDirJob, vcf2KeepGivenSNPByvcftoolsJob]
						vcf2filterByDepthExtraDependentInputLs=[]
						counter += 4
						continue	#skip the rest
					else:
						vcf1filterByDepthInput=vcf1
						lastRoundJobLs=[vcf1DepthFilterDirJob]
						lastRoundExtraDependentInputLs=[vcf1.tbi_F]
						vcf2filterByDepthInput=vcf2
						vcf2filterByDepthParentJobLs=[vcf2DepthFilterDirJob]
						vcf2filterByDepthExtraDependentInputLs=[vcf2.tbi_F]
					vcf1AfterDepthFilter = File(os.path.join(vcf1DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
							GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf1filterByDepthInput, outputVCFF=vcf1AfterDepthFilter, \
							parentJobLs=lastRoundJobLs, \
							alnStatForFilterF=alnStatForFilterF, \
							extraDependentInputLs=lastRoundExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
					
					
					vcf1AfterDepthFilterGzip = File("%s.gz"%vcf1AfterDepthFilter.name)
					vcf1AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf1AfterDepthFilter.name)
					vcf1FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf1FilterByDepthJob, inputF=vcf1AfterDepthFilter, outputF=vcf1AfterDepthFilterGzip, \
							transferOutput=False)
					
					vcf2AfterDepthFilter = File(os.path.join(vcf2DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf2FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepth2Java, \
							GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf2filterByDepthInput, outputVCFF=vcf2AfterDepthFilter, \
							parentJobLs=vcf2filterByDepthParentJobLs, \
							alnStatForFilterF=alnStatForFilterF, \
							extraDependentInputLs=vcf2filterByDepthExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
					
					vcf2AfterDepthFilterGzip = File("%s.gz"%vcf2AfterDepthFilter.name)
					vcf2AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf2AfterDepthFilter.name)
					vcf2FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf2FilterByDepthJob, inputF=vcf2AfterDepthFilter, outputF=vcf2AfterDepthFilterGzip, \
							transferOutput=False)
					
					snpMisMatchStatFile = File(os.path.join(SNPMismatchStatDir, '%s_snpMismatchStat.tsv'%(os.path.splitext(commonPrefix)[0])))
					calculateSNPMismatchRateOfTwoVCFJob = self.addCalculateTwoVCFSNPMismatchRateJob(workflow, \
							executable=workflow.CalculateSNPMismatchRateOfTwoVCF, \
							vcf1=vcf1AfterDepthFilterGzip, vcf2=vcf2AfterDepthFilterGzip, snpMisMatchStatFile=snpMisMatchStatFile, \
							maxSNPMismatchRate=maxSNPMismatchRate, \
							parentJobLs=[vcf1FilterByDepthBGZipTabixJob, vcf2FilterByDepthBGZipTabixJob, SNPMismatchStatDirJob], \
							job_max_memory=1000, extraDependentInputLs=[], \
							transferOutput=False)
					
					
					outputFnamePrefix = os.path.join(vcf1_vcftoolsFilterDir, '%s.filter_by_vcftools'%(commonPrefix))
					vcf1FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
							inputVCFF=vcf1AfterDepthFilterGzip, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[vcf1FilterByDepthBGZipTabixJob, vcf1_vcftoolsFilterDirJob, calculateSNPMismatchRateOfTwoVCFJob], \
							snpMisMatchStatFile=snpMisMatchStatFile, \
							minMAC=self.minMAC, minMAF=self.minMAF, \
							maxSNPMissingRate=maxSNPMissingRate,\
							extraDependentInputLs=[vcf1FilterByDepthBGZipTabixJob.tbi_F], extraArguments="--recode-INFO-all")
					vcf1FilterByvcftoolsGzip = File("%s.gz"%vcf1FilterByvcftoolsJob.output.name)
					vcf1FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf1FilterByvcftoolsJob, inputF=vcf1FilterByvcftoolsJob.output, outputF=vcf1FilterByvcftoolsGzip, \
							transferOutput=True)
					
					outputFnamePrefix = os.path.join(vcf2_vcftoolsFilterDir, '%s.filter_by_vcftools'%(commonPrefix))
					vcf2FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
							inputVCFF=vcf2AfterDepthFilterGzip, \
							outputFnamePrefix=outputFnamePrefix, \
							parentJobLs=[vcf2FilterByDepthBGZipTabixJob, vcf2_vcftoolsFilterDirJob, calculateSNPMismatchRateOfTwoVCFJob],
							snpMisMatchStatFile=snpMisMatchStatFile, \
							minMAC=self.minMAC, minMAF=self.minMAF, \
							maxSNPMissingRate=maxSNPMissingRate,\
							extraDependentInputLs=[vcf2FilterByDepthBGZipTabixJob.tbi_F], extraArguments="--recode-INFO-all")
					
					vcf2FilterByvcftoolsGzip = File("%s.gz"%vcf2FilterByvcftoolsJob.output.name)
					vcf2FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf2FilterByvcftoolsJob, inputF=vcf2FilterByvcftoolsJob.output, outputF=vcf2FilterByvcftoolsGzip, \
							transferOutput=True)
					
					counter += 9
		sys.stderr.write("%s jobs.\n"%(counter+1))
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.06.30
		"""
		returnData = parentClass.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix, \
										passingData=passingData, transferOutput=transferOutput, **keywords)
		
		self.auxDirJob = self.addMkDirJob(outputDir="%sAuxilliary"%(outputDirPrefix))
		
		self.chooseContaminantOrNotDirJob = self.addMkDirJob(outputDir="%sChoose_is_contaminated%s"%\
														(outputDirPrefix, self.is_contaminated))
		
		self.filterPASSDirJob = self.addMkDirJob(outputDir="%sFILTER_PASS"%(outputDirPrefix))
		
		self.vcf1DepthFilterDirJob = self.addMkDirJob(outputDir= "%sDepthFilter"%(outputDirPrefix))
		
		self.SNPMismatchStatDirJob = self.addMkDirJob(outputDir="%sSNPMismatchStat"%(outputDirPrefix))
		
		self.filterDirJob = self.addMkDirJob(outputDir="%sVCFtoolsFilter"%(outputDirPrefix))
		
		self.filterStatDirJob = self.addMkDirJob(outputDir="%sFilterStat"%(outputDirPrefix))
		
		if self.is_contaminated is not None:
			#a job outputs contaminant samples's alignment read-groups
			inputFileBasenamePrefix = utils.getFileBasenamePrefixFromPath(self.firstVCFJobData.file.name)
			outputFile = File(os.path.join(self.auxDirJob.output, '%s.is_contaminated%s.sampleIDList.tsv'%\
										(inputFileBasenamePrefix, self.is_contaminated)))
			self.outputContaminantOrNotSampleIDListJob = self.addExtractSampleIDJob(inputFile=self.firstVCFJobData.file, \
								outputFile=outputFile,\
								is_contaminated=self.is_contaminated, outputFormat=3,\
								returnData=returnData,\
								transferOutput=True, \
								parentJobLs=[self.firstVCFJobData.jobLs, self.auxDirJob])
			
			filterByChooseContaminantOrNotStatMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByChooseContaminantOrNot_s0.tsv'))
			self.filterByChooseContaminantOrNotStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByChooseContaminantOrNotStatMergeFile, transferOutput=True, \
								parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")
			
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByChooseContaminantOrNotStatMergeJob))
		
		if self.keepSNPPosF:
			filterByGivenSitesStatMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByGivenSitesStat_s1.tsv'))
			self.filterByGivenSitesStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByGivenSitesStatMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByGivenSitesStatMergeJob))
		
		if self.excludeFilteredSites:
			filterPassStatMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByFILTER_PASS_s1.5.tsv'))
			self.filterPassStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterPassStatMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterPassStatMergeJob))
		
		if self.minDepthPerGenotype:
			filterByMinDP1MergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByMinDP%s_s2.tsv'%(self.minDepthPerGenotype)))
			self.filterByMinDP1MergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinDP1MergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByMinDP1MergeJob))
		
		if self.cumulativeMedianDepth and self.depthFoldChange:
			filterByDepthStatMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByDepthStat_%sFoldMedianDepth%s_s3.tsv'%\
														(self.depthFoldChange, self.cumulativeMedianDepth)))
			self.filterByDepthStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByDepthStatMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByDepthStatMergeJob))
		
		if self.minMAC is not None:
			filterByMinMACMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByMinMAC%s_s4.tsv'%(self.minMAC)))
			self.filterByMinMACMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinMACMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByMinMACMergeJob))
		
		if self.minMAF is not None:
			filterByMinMAFMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByMinMAFS%s_s5.tsv'%(self.minMAF)))
			self.filterByMinMAFMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinMAFMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
								extraArguments="-k 1 -v 2-4")
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByMinMAFMergeJob))
		
		if self.maxSNPMissingRate is not None and self.maxSNPMissingRate>=0 and self.maxSNPMissingRate<1:
			filterByMaxSNPMissingRateMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByMaxSNPMissingRate%s_s6.tsv'%(self.maxSNPMissingRate)))
			self.filterByMaxSNPMissingRateMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=filterByMaxSNPMissingRateMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
							extraArguments="-k 1 -v 2-4")
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByMaxSNPMissingRateMergeJob))
		
		if self.minNeighborDistance is not None and self.minNeighborDistance>=0:
			filterByMinNeighborDistanceMergeFile = File(os.path.join(self.filterStatDirJob.output, 'filterByMinNeighborDistance%s_s7.tsv'%(self.minNeighborDistance)))
			self.filterByMinNeighborDistanceMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=filterByMinNeighborDistanceMergeFile, transferOutput=True, parentJobLs=[self.filterStatDirJob],\
							extraArguments="-k 1 -v 2-4")
			returnData.jobDataLs.append(self.constructJobDataFromJob(job=self.filterByMinNeighborDistanceMergeJob))
		
		return returnData
	
	def mapEachInterval(self, workflow=None, VCFJobData=None, chromosome=None,intervalData=None,\
					mapEachChromosomeData=None, passingData=None, transferOutput=False, \
					**keywords):
		"""
		2013.06.30
		
			argument VCFJobData looks like PassingData(file=splitVCFFile, vcfFile=splitVCFFile, fileLs=[splitVCFFile], \
																		job=splitVCFJob, jobLs=[splitVCFJob], tbi_F=None)
		"""
		
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		passingData.intervalFileBasenamePrefix
		passingData.splitVCFFile
		passingData.unitNumber
		"""
		## 2013.06.19 structures available from passingData, specific to the interval
		passingData.splitVCFFile = splitVCFFile
		passingData.unitNumber = unitNumber
		passingData.intervalFileBasenamePrefix = '%s_%s_splitVCF_u%s'%(chromosome, commonPrefix, unitNumber)
		passingData.noOfIndividuals = jobData.file.noOfIndividuals
		passingData.span = self.intervalSize + self.intervalOverlapSize*2 	#2013.06.19 for memory/walltime gauging
		"""
		commonPrefix = passingData.intervalFileBasenamePrefix
			
		lastRoundJobLs= VCFJobData.jobLs
		if lastRoundJobLs:
			lastVCFJob = lastRoundJobLs[0]
			if not hasattr(lastVCFJob, 'tbi_F'):
				lastVCFJob.tbi_F = None
			lastVCFJob.output = VCFJobData.file	#2013.07.18 otherwise it's always unit1
		else:
			lastVCFJob = PassingData(output=VCFJobData.file, tbi_F=None)	#2012.8.3 fake, not a job. only useful when all filtering jobs are skipped.
		lastRoundExtraDependentInputLs =[]
		
		noTransferFlagJobSet = set()
		
		if self.is_contaminated is not None:
			########## 
			### step 0 remove, contaminants (marked by non-0 value of is_contaminated and outdated_index in individual_sequence table) 
			# job to exclude contamnant	
			outputF = File(os.path.join(self.chooseContaminantOrNotDirJob.output, '%s.is_contaminated%s.vcf'%\
									(commonPrefix, self.is_contaminated)))
			selectVariantsContaminantOrNotJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
						inputF=lastVCFJob.output, outputF=outputF, \
						interval=None,\
						refFastaFList=self.registerReferenceData.refFastaFList, \
						sampleIDKeepFile=self.outputContaminantOrNotSampleIDListJob.output, snpIDKeepFile=None, \
						sampleIDExcludeFile=None, \
						parentJobLs=lastRoundJobLs + [self.outputContaminantOrNotSampleIDListJob, self.chooseContaminantOrNotDirJob], \
						extraDependentInputLs=lastRoundExtraDependentInputLs, transferOutput=False, \
						extraArguments="--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES", \
						extraArgumentList=None, job_max_memory=6000, walltime=None)
			currentVCFJob = selectVariantsContaminantOrNotJob
			
			#check how much sites got filtered by this filter
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByChooseContaminantOrNotStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastVCFJob.output, currentVCFJob=currentVCFJob, \
									statMergeJob=self.filterByChooseContaminantOrNotStatMergeJob, \
									parentJobLs=lastRoundJobLs + [self.filterStatDirJob])
			
			noTransferFlagJobSet.add(currentVCFJob)
			
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
			lastRoundExtraDependentInputLs=[]
		
		if self.keepSNPPosF:
			#toss SNPs that are not in this keepSNPPosFname file
			outputFnamePrefix = os.path.join(self.filterDirJob.output, '%s.keepGivenSNP'%(commonPrefix))
			parentJobLs = lastRoundJobLs + [self.filterDirJob]
			if self.keepSNPPosParentJobLs:
				parentJobLs += self.keepSNPPosParentJobLs
			vcf1KeepGivenSNPByvcftoolsJob = self.addFilterJobByvcftools(vcftoolsWrapper=workflow.vcftoolsWrapper, \
					inputVCFF=lastVCFJob.output,\
					outputFnamePrefix=outputFnamePrefix, \
					parentJobLs = parentJobLs, \
					snpMisMatchStatFile=self.keepSNPPosF, \
					minMAC=None, minMAF=None, \
					maxSNPMissingRate=None,\
					extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
			
			currentVCFJob = vcf1KeepGivenSNPByvcftoolsJob
			#check how much sites got filtered
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByGivenSitesStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=lastVCFJob.output, \
								currentVCFJob=currentVCFJob,\
								statMergeJob=self.filterByGivenSitesStatMergeJob, \
								parentJobLs=[self.filterStatDirJob]+lastRoundJobLs)
		
			lastVCFJob = currentVCFJob
			noTransferFlagJobSet.add(currentVCFJob)
			lastRoundJobLs=[currentVCFJob]
			lastRoundExtraDependentInputLs=[]
		else:
			lastRoundJobLs= lastRoundJobLs
			lastRoundExtraDependentInputLs =lastRoundExtraDependentInputLs
			lastVCFJob = lastVCFJob	#faking it
		
		if self.excludeFilteredSites:
			#2013.05.20
			if self.excludeFilteredSites==1:
				selectExpression="vc.isNotFiltered()"
			else:
				selectExpression = "vc.isNotFiltered() && vc.isSNP()"
			vcfAfterFILTERPASS = File(os.path.join(self.filterPASSDirJob.output, '%s_filterByFILTERPASS.vcf'%(commonPrefix)))
			FILTER_PASS_Job = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=lastVCFJob.output, outputF=vcfAfterFILTERPASS, \
				refFastaFList=self.registerReferenceData.refFastaFList, \
				sampleIDKeepFile=None, snpIDKeepFile=None, sampleIDExcludeFile=None, \
				interval=None,\
				parentJobLs=lastRoundJobLs + self.filterPASSDirJob, extraDependentInputLs=lastRoundExtraDependentInputLs, transferOutput=False, \
				extraArguments="""--select_expressions "%s" """%(selectExpression), job_max_memory=2000, walltime=None)
			
			currentVCFJob = FILTER_PASS_Job
			#check how much sites got filtered
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s_FILTER_PASS_Stat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=lastVCFJob.output, \
								currentVCFJob=currentVCFJob,\
								statMergeJob=self.filterPassStatMergeJob, \
								parentJobLs=lastRoundJobLs + [self.filterStatDirJob])
			
			lastVCFJob = currentVCFJob
			noTransferFlagJobSet.add(currentVCFJob)
			lastRoundJobLs=[currentVCFJob]
			lastRoundExtraDependentInputLs = []
		#2012.8.1 mask zero-depth sites
		if self.minDepthPerGenotype:
			outputFnamePrefix = os.path.join(self.vcf1DepthFilterDirJob.output, '%s.minDP%s'%(commonPrefix, self.minDepthPerGenotype))
			maskZeroDepthGenotypeAsMissingJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
					inputVCFF=lastVCFJob.output, \
					outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=lastRoundJobLs + [self.vcf1DepthFilterDirJob], \
					minMAC=None, minMAF=None, \
					maxSNPMissingRate=None,\
					extraDependentInputLs=lastRoundExtraDependentInputLs, outputFormat='--recode', \
					extraArguments="--recode-INFO-all --minDP %s"%(self.minDepthPerGenotype))
			
			#check how much sites got filtered
			currentVCFJob = maskZeroDepthGenotypeAsMissingJob
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByminDP%s.tsv'%(commonPrefix, self.minDepthPerGenotype)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=lastVCFJob.output, \
								currentVCFJob=currentVCFJob,\
								statMergeJob=self.filterByMinDP1MergeJob, parentJobLs=lastRoundJobLs + [self.filterStatDirJob])
		
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs=[currentVCFJob]
			
		if self.cumulativeMedianDepth and self.depthFoldChange:
			vcf1AfterDepthFilter = File(os.path.join(self.vcf1DepthFilterDirJob.output, '%s.filterByDepth.vcf'%(commonPrefix)))
			#2013.05.20 TC stands for total coverage in platypus output
			selectExpression = "%s>=%s && %s<=%s"%(self.siteTotalCoverageINFOFieldName, \
										self.cumulativeMedianDepth/self.depthFoldChange, \
										self.siteTotalCoverageINFOFieldName, self.cumulativeMedianDepth*self.depthFoldChange)
			vcf1FilterByDepthJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=lastVCFJob.output, outputF=vcf1AfterDepthFilter, \
				refFastaFList=self.registerReferenceData.refFastaFList, \
				sampleIDKeepFile=None, snpIDKeepFile=None, sampleIDExcludeFile=None, \
				interval=None,\
				parentJobLs=lastRoundJobLs + [self.vcf1DepthFilterDirJob], \
				extraDependentInputLs=lastRoundExtraDependentInputLs, transferOutput=False, \
				extraArguments="""--select_expressions "%s" """%(selectExpression), job_max_memory=2000, walltime=None)
			"""
			vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
					GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
					refFastaFList=refFastaFList, inputVCFF=lastVCFJob.output, outputVCFF=vcf1AfterDepthFilter, \
					parentJobLs=lastRoundJobLs, \
					alnStatForFilterF=alnStatForFilterF, \
					extraDependentInputLs=lastRoundExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
			"""
		
			currentVCFJob = vcf1FilterByDepthJob
			#check how much sites got filtered by depth filter
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByDepthStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastVCFJob.output, parentJobLs=lastRoundJobLs + [self.filterStatDirJob], \
									currentVCFJob=currentVCFJob,\
									statMergeJob=self.filterByDepthStatMergeJob)
			
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
			
		
			
		if self.minMAC is not None:
			outputFnamePrefix = os.path.join(self.filterDirJob.output, '%s.filterByMinMAC'%(commonPrefix))
			vcf1FilterByMinMACJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
					inputVCFF=lastVCFJob.output, \
					outputFnamePrefix=outputFnamePrefix, \
					parentJobLs=lastRoundJobLs+ [self.filterDirJob], \
					snpMisMatchStatFile=None, \
					minMAC=self.minMAC, minMAF=None, \
					maxSNPMissingRate=None,\
					extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
			
			#check how much sites got filtered by maxSNPMissingRate filter
			currentVCFJob = vcf1FilterByMinMACJob
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByMinMACStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=lastVCFJob.output, currentVCFJob=currentVCFJob,\
								parentJobLs=lastRoundJobLs + [self.filterStatDirJob], statMergeJob=self.filterByMinMACMergeJob)
			
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
		
		if self.minMAF is not None:
			outputFnamePrefix = os.path.join(self.filterDirJob.output, '%s.filterByMinMAF'%(commonPrefix))
			vcf1FilterByMinMAFJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
					inputVCFF=lastVCFJob.output, \
					outputFnamePrefix=outputFnamePrefix, \
					parentJobLs= lastRoundJobLs + [self.filterDirJob], \
					snpMisMatchStatFile=None, \
					minMAC=None, minMAF=self.minMAF, \
					maxSNPMissingRate=None,\
					extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
			
			#check how much sites got filtered by maxSNPMissingRate filter
			currentVCFJob = vcf1FilterByMinMAFJob
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByMinMAFStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
								vcf1=lastVCFJob.output, currentVCFJob=currentVCFJob,\
								statMergeJob=self.filterByMinMAFMergeJob, \
								parentJobLs=lastRoundJobLs + [self.filterStatDirJob])
			
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
		
		if self.maxSNPMissingRate is not None and self.maxSNPMissingRate>=0 and self.maxSNPMissingRate<1:
			outputFnamePrefix = os.path.join(self.filterDirJob.output, '%s.filterByMaxSNPMissingRate'%(commonPrefix))
			vcf1FilterByMaxSNPMissingRateJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
					inputVCFF=lastVCFJob.output, \
					outputFnamePrefix=outputFnamePrefix, \
					parentJobLs= lastRoundJobLs + [self.filterDirJob], \
					snpMisMatchStatFile=None, \
					minMAC=None, minMAF=None, \
					maxSNPMissingRate=self.maxSNPMissingRate,\
					extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
			
			currentVCFJob = vcf1FilterByMaxSNPMissingRateJob
			#check how much sites got filtered by maxSNPMissingRate filter
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByMaxSNPMissingRateStat.tsv'%(commonPrefix)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastVCFJob.output, currentVCFJob=currentVCFJob, \
									statMergeJob=self.filterByMaxSNPMissingRateMergeJob, \
									parentJobLs=lastRoundJobLs)
			
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
			lastRoundExtraDependentInputLs=None
		
		if self.minNeighborDistance is not None and self.minNeighborDistance>=0:
			outputFile = File(os.path.join(self.filterDirJob.output, '%s.filterByMinNeighborDistance%s.vcf'%(commonPrefix, self.minNeighborDistance)))
			filterByMinNeighborDistanceJob = self.addGenericJob(executable=self.FilterVCFSNPCluster, inputFile=lastVCFJob.output, \
				outputFile=outputFile, \
				parentJobLs=lastRoundJobLs + [self.filterDirJob], extraDependentInputLs=lastRoundExtraDependentInputLs, \
				extraOutputLs=None,\
				transferOutput=False, \
				extraArgumentList=None, extraArguments="--minNeighborDistance %s"%(self.minNeighborDistance), \
				key2ObjectForJob=None, job_max_memory=2000)
			
			
			currentVCFJob = filterByMinNeighborDistanceJob
			#check how much sites got filtered by this filter
			outputF = File(os.path.join(self.filterStatDirJob.output, '%s.filterByMinNeighborDistanceJob%s.tsv'%(commonPrefix, self.minNeighborDistance)))
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastVCFJob.output, currentVCFJob=currentVCFJob, \
									statMergeJob=self.filterByMinNeighborDistanceMergeJob, parentJobLs=lastRoundJobLs)
			
			noTransferFlagJobSet.add(currentVCFJob)
			lastVCFJob = currentVCFJob
			lastRoundJobLs = [currentVCFJob]
			lastRoundExtraDependentInputLs=[]
		
		#bgzip last VCF job's output
		bgzipFile = File("%s.gz"%lastVCFJob.output.name)
		bgzipJob = self.addBGZIP_tabix_Job(bgzip_tabix=workflow.bgzip_tabix, \
			parentJobLs=lastRoundJobLs, inputF=lastVCFJob.output, \
			outputF=bgzipFile, \
			transferOutput=False)
		
		lastVCFJob = bgzipJob
		lastRoundJobLs = [bgzipJob]
		lastRoundExtraDependentInputLs = [bgzipJob.tbi_F]
		
		lastBGZipTabixJobOutputFile = getattr(lastVCFJob, 'output', None)	#could be None if all filter jobs are skipped
		lastBGZipTabixJobTbiF = getattr(lastVCFJob, 'tbi_F', None)
		
		returnData.lastVCFJob = lastVCFJob
		returnData.jobLs=lastRoundJobLs
		returnData.vcfFile=lastBGZipTabixJobOutputFile
		returnData.tbi_F=lastBGZipTabixJobTbiF
		returnData.fileLs=[lastBGZipTabixJobOutputFile, lastBGZipTabixJobTbiF]
		
		return returnData
	
	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		2013.06.30 make sure something like this is the final return structure out of addJobsToFilterOneVCFDir()
		  
			returnData.jobDataLs.append(PassingData(jobLs=lastRoundJobLs, vcfFile=lastBGZipTabixJobOutputFile, \
								tbi_F=lastBGZipTabixJobTbiF, \
								fileLs=[lastBGZipTabixJobOutputFile, lastBGZipTabixJobTbiF]))
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		
		
		baseInputVolume = 200*20000
		noOfIndividuals = getattr(passingData.jobData.file, "noOfIndividuals", 200)
		noOfLoci = getattr(passingData.jobData.file, "noOfLoci", 20000)
		realInputVolume = noOfIndividuals * noOfLoci
		
		#base is 4X coverage in 20Mb region => 120 minutes
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
							minJobPropertyValue=2000, maxJobPropertyValue=8000).value
		if self.intervalOverlapSize==0:
			returnData.concatenateIntoOneVCFJobData = self.concatenateIntervalsIntoOneVCFSubWorkflow(chromosome=chromosome, \
						refFastaFList=self.registerReferenceData.refFastaFList,\
						passingData=passingData, \
						intervalJobLs=[pdata.lastVCFJob for pdata in mapEachIntervalDataLs],\
						outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		else:
			returnData.concatenateIntoOneVCFJobData = self.concatenateOverlapIntervalsIntoOneVCFSubWorkflow(chromosome=chromosome, \
						passingData=passingData, \
						intervalJobLs=[pdata.lastVCFJob for pdata in mapEachIntervalDataLs],\
						outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		#returnData.concatenateIntoOneVCFJobData.noOfIndividuals = passingData.jobData.file.noOfIndividuals
		#returnData.concatenateIntoOneVCFJobData.noOfLoci = passingData.jobData.file.noOfLoci
		return returnData
	
	def reduce(self, workflow=None, reduceEachChromosomeDataLs=None, \
			mapEachChromosomeDataLs=None, passingData=None, transferOutput=True, \
			**keywords):
		"""
		2013.07.03 
			the data returned by this function would be the eventual return of AbstractVCFWorkflow.addAllJobs()
			in this case, another workflow BootstrapFilterCalculateVCFWorkflow.py needs this return data to further
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachChromosomeDataLs = mapEachChromosomeDataLs
		returnData.reduceEachChromosomeDataLs = reduceEachChromosomeDataLs
		for reduceEachVCFDataLs in passingData.reduceEachVCFDataLsLs:
			for reduceEachVCFData in reduceEachVCFDataLs:
				returnData.jobDataLs.append(reduceEachVCFData.concatenateIntoOneVCFJobData)
			"""
			returnData.jobDataLs.append(PassingData(jobLs=[concatenateIntoOneVCFJob], \
										vcfFile=concatenateIntoOneVCFJob.output, \
								tbi_F=getattr(concatenateIntoOneVCFJob, 'tbi_F', None), \
								fileLs=concatenateIntoOneVCFJob.outputLs))
			"""
		return returnData
	
	def addJobsToFilterOneVCFDir(self, workflow=None, inputData=None, registerReferenceData=None, \
								cumulativeMedianDepth=None, depthFoldChange=None, keepSNPPosF=None, \
						onlyKeepBiAllelicSNP=True,\
						minDepthPerGenotype=False, minMAC=None, minMAF=None, maxSNPMissingRate=None, outputDirPrefix="",\
						minNeighborDistance=None, transferOutput=True, keepSNPPosParentJobLs=None, excludeFilteredSites=0,\
						**keywords):
		"""
		2013.05.20 add argument excludeFilteredSites, cumulativeMedianDepth, depthFoldChange
		2012.9.11 add argument keepSNPPosParentJobLs
		2012.9.6 add argument minNeighborDistance
		2012.7.30 add stat collecting jobs
		2012.5.10
			add extraArguments="--recode-INFO-all" to addFilterJobByvcftools()
		2012.1.14
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding filter-VCF jobs for %s vcf files ... \n"%(len(inputData.jobDataLs)))
		
		#2013.06.30 pass keepSNPPosParentJobLs to self so that jobs could see it
		self.keepSNPPosParentJobLs = keepSNPPosParentJobLs
		self.keepSNPPosF = keepSNPPosF
		
		return self.addAllJobs(workflow=workflow, inputVCFData=inputData, \
				chr2IntervalDataLs=self.chr2IntervalDataLs, samtools=workflow.samtools, \
				GenomeAnalysisTKJar=workflow.GenomeAnalysisTKJar, \
				CreateSequenceDictionaryJava=workflow.CreateSequenceDictionaryJava, \
				CreateSequenceDictionaryJar=workflow.CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexJar=workflow.BuildBamIndexJar,\
				mv=workflow.mv, \
				registerReferenceData=registerReferenceData,\
				needFastaIndexJob=getattr(self, 'needFastaIndexJob',False), \
				needFastaDictJob=getattr(self, 'needFastaDictJob', False), \
				data_dir=self.data_dir, no_of_gatk_threads = 1,\
				intervalSize=self.intervalSize, intervalOverlapSize=self.intervalOverlapSize, \
				outputDirPrefix="%s%s"%(outputDirPrefix, self.pegasusFolderName), \
				transferOutput=transferOutput,)
	
	
	def setup_run(self):
		"""
		wrap all standard pre-run() related functions into this function.
		setting up for run(), called by run()
		
		2013.06.11 assign all returned data to self, rather than pdata (pdata has become self)
		2013.05.20 parse cumulativeMedianDepth from a sample VCF file.
		
		"""
		pdata = parentClass.setup_run(self)
		
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		"""
		firstVCFFile = None
		firstVCFJobData = None
		inputData = None
		cumulativeMedianDepth = None
		if self.vcf1Dir and  not getattr(self, 'vcf2Dir', None):
			#a relative-path name for vcf1Dir
			vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir, defaultName='vcf1')
			inputData = self.registerAllInputFiles(workflow=pdata.workflow, inputDir=self.vcf1Dir, \
										input_site_handler=self.input_site_handler, \
										checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
										pegasusFolderName="%s_%s"%(self.pegasusFolderName, vcf1Name), \
										maxContigID=self.maxContigID, \
										minContigID=self.minContigID, \
										db_vervet=getattr(self, 'db_vervet', None), \
										needToKnowNoOfLoci=getattr(self, 'needToKnowNoOfLoci', True),\
										minNoOfLociInVCF=getattr(self, 'minNoOfLociInVCF', 10))
			
			if inputData.jobDataLs:
				#self.firstVCFJobData is not set in parentClass.setup_run() because self.inputDir is missing
				firstVCFJobData = inputData.jobDataLs[0]
				firstVCFFile = inputData.jobDataLs[0].file
		sys.stderr.write("One sample VCF file is %s (used to get alignments).\n"%(firstVCFFile))
		if self.depthFoldChange and self.depthFoldChange>0 and firstVCFFile:
			alignmentLs = self.db.getAlignmentsFromVCFFile(inputFname=yh_pegasus.getAbsPathOutOfFile(firstVCFFile))
			cumulativeMedianDepth = self.db.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
										defaultSampleAlignmentDepth=8)
			"""
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
			#alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
		#else:
			#alnStatForFilterF = None
			"""
		
		if self.keepSNPPosFname:
			keepSNPPosF = self.registerOneInputFile(inputFname=os.path.abspath(self.keepSNPPosFname),\
														folderName=self.pegasusFolderName)
		else:
			keepSNPPosF = None
		self.inputData = inputData
		self.keepSNPPosF = keepSNPPosF
		self.cumulativeMedianDepth = cumulativeMedianDepth
		self.firstVCFJobData = firstVCFJobData
		self.firstVCFFile = firstVCFFile
		sys.stderr.write("cumulativeMedianDepth for all samples is %s.\n"%(cumulativeMedianDepth))
		return self
	
	def run(self):
		"""
		"""
		
		pdata = self.setup_run()
		
		
		if self.vcf1Dir and self.vcf2Dir:
			self.addJobsToFilterTwoVCFDir(vcf1Dir=self.vcf1Dir, vcf2Dir=self.vcf2Dir, \
							registerReferenceData=pdata.registerReferenceData, \
							keepSNPPosF=self.keepSNPPosF, \
							onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, minMAC=self.minMAC, minMAF=self.minMAF, \
							maxSNPMissingRate=self.maxSNPMissingRate, maxSNPMismatchRate=self.maxSNPMismatchRate)
		elif self.vcf1Dir:
			# 2012.5.1 filter only on the 1st vcf folder
			
			self.addJobsToFilterOneVCFDir(inputData=pdata.inputData, registerReferenceData=pdata.registerReferenceData, \
									cumulativeMedianDepth=pdata.cumulativeMedianDepth, depthFoldChange=self.depthFoldChange, \
									keepSNPPosF=self.keepSNPPosF, \
									onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
									minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate,\
									minDepthPerGenotype=self.minDepthPerGenotype, outputDirPrefix="",\
									minNeighborDistance=self.minNeighborDistance, keepSNPPosParentJobLs=None,\
									excludeFilteredSites=self.excludeFilteredSites)
		self.end_run()
		


	
if __name__ == '__main__':
	main_class = FilterVCFPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
