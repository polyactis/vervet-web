#!/usr/bin/env python
"""
Examples:
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor_20111106T1554/;

	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top4Contigs_multi_sample_condor.xml  -z uclaOffice -u yh -q alnStatForFilter.tsv
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;
	
	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool -a 524
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor.xml -z uclaOffice -u yh
		-q alnStatForFilter.tsv
	
	#2011-11-11
	dirPrefix=./AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_multi_sample_condor_20111106T1554/;	
	%s  -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool -a 524
		-o FilterVCF_4HighCovVRC_isq_15_18_vs_524_top804Contigs_minGQ30_multi_sample_condor.xml -z uclaOffice -u yh
		-q alnStatForFilter.tsv	-G 30
	
	#2011.11.28 (-K = keep only bi-allelic SNPs), -E = checkEmptyVCFByReading
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s  -i $dirPrefix\gatk -I $dirPrefix\call/ -l condorpool -j condorpool
		-o FilterVCF_LowPass_top7559Contigs_no12eVarFilter_minGQ1_maxSNPMisMatch0.1_minMAC5_maxSNPMissing0.25.xml
		-z uclaOffice -u yh -q ./alnStatForFilter.2011.12.9T0207.tsv
		-G1 -K -N 0.1 -n 5 -x 0.25 -a 524 -C 50 -E
	
	#2011.12.9 to remove SNPs that are not in a file. no other filters.
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s -i $dirPrefix\gatk -I $dirPrefix\call/ -l condorpool -j condorpool
		-o Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4.xml
		-z uclaOfficeTemp -u yh -q ./alnStatForFilter.2011.12.9T0207.tsv -G0 -N 1 -n 0 -x 1 -A 100000 -a 524 -C 50
		-S ./4HighCovVRC_inter_minMAC4_vs_LowPass_top7559Contigs_no12eVarFilter_inter.2011.12.9T0107/overlapPos.tsv
	
	#2011.12.19 run on hoffman2's condorpool
	%s ....
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		
	#2012.5.1 filter trioCaller output with total depth (no minGQ filter anymore), minMAC=10 (-n 10), maxSNPMismatchRate=1 (-x 1.0)
	#
	%s -i AlignmentToTrioCallPipeline_VRC_top7559Contigs.2011.12.15T0057/trioCaller/ -l condorpool -j condorpool
		-z uclaOffice -u yh -q ./alnStatForFilter.2012.5.1T1430.tsv  -n 10 -x 1.0 -a 524 -C 50 -E
		-o FilterVCF_trioCallerTop7559Contigs.xml

	
Description:
	2011-11-7 pipeline that filters VCF by depth, GQ, MAC, SNP missing rate, mismatch rate between two input VCFs
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from AbstractVervetWorkflow import AbstractVervetWorkflow

class FilterVCFPipeline(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype. 2012.5.1 no longer enforced in FilterVCFByDepth.java', ],\
						('depthFoldChange', 1, float): [2.0, 'A', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						("maxSNPMismatchRate", 0, float): [0, 'N', 1, 'maximum SNP mismatch rate between two vcf calls'],\
						("minMAC", 0, int): [2, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
						("minMAF", 0, float): [None, 'm', 1, 'minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRate", 0, float): [0, 'x', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('vcf1Dir', 1, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files. If not provided, filter vcf1Dir without applying maxSNPMismatchRate filter.', ],\
						("onlyKeepBiAllelicSNP", 0, int): [0, 'K', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 1, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						('keepSNPPosFname', 0, ): ['', '', 1, 'a tab-delimited file (optional), chr pos 0, to pre-filter SNPs that are absent in this file.'],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
	
	
	
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
	
	def addJobsToFilterTwoVCFDir(self, workflow, vcf1Dir, vcf2Dir, refFastaFList, alnStatForFilterF, keepSNPPosF, onlyKeepBiAllelicSNP=True,\
							minMAC=None, minMAF=None, maxSNPMissingRate=None):
		"""
		2012.5.10
			add extraArguments="--recode-INFO-all" to addFilterJobByvcftools()
		2012.1.14
		"""
		refFastaF = refFastaFList[0]
		
		#name to distinguish between vcf1Dir, and vcf2Dir
		vcf1Name = self.findProperVCFDirIdentifier(vcf1Dir, defaultName='vcf1')
		vcf2Name = self.findProperVCFDirIdentifier(vcf2Dir, defaultName='vcf2')
		if vcf2Name==vcf1Name or not vcf2Name:
			vcf2Name = "vcf2"
		
		vcf1DepthFilterDir = "%s_DepthFilter"%(vcf1Name)
		vcf1DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1DepthFilterDir)
		vcf2DepthFilterDir = "%s_DepthFilter"%(vcf2Name)
		vcf2DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf2DepthFilterDir)
		
		SNPMismatchStatDir = "SNPMismatchStat"
		SNPMismatchStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=SNPMismatchStatDir)
		
		vcf1_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf1Name)
		vcf1_vcftoolsFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1_vcftoolsFilterDir)
		vcf2_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf2Name)
		vcf2_vcftoolsFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf2_vcftoolsFilterDir)
		
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
						vcf1filterByDepthParentJobLs=[vcf1DepthFilterDirJob, vcf1KeepGivenSNPByvcftoolsJob]
						vcf1filterByDepthExtraDependentInputLs=[]
						vcf2filterByDepthInput=vcf2KeepGivenSNPByvcftoolsJob.output
						vcf2filterByDepthParentJobLs=[vcf2DepthFilterDirJob, vcf2KeepGivenSNPByvcftoolsJob]
						vcf2filterByDepthExtraDependentInputLs=[]
						counter += 4
						continue	#skip the rest
					else:
						vcf1filterByDepthInput=vcf1
						vcf1filterByDepthParentJobLs=[vcf1DepthFilterDirJob]
						vcf1filterByDepthExtraDependentInputLs=[vcf1.tbi_F]
						vcf2filterByDepthInput=vcf2
						vcf2filterByDepthParentJobLs=[vcf2DepthFilterDirJob]
						vcf2filterByDepthExtraDependentInputLs=[vcf2.tbi_F]
					vcf1AfterDepthFilter = File(os.path.join(vcf1DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf1filterByDepthInput, outputVCFF=vcf1AfterDepthFilter, \
							parentJobLs=vcf1filterByDepthParentJobLs, \
							alnStatForFilterF=alnStatForFilterF, \
							extraDependentInputLs=vcf1filterByDepthExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
				
					vcf1AfterDepthFilterGzip = File("%s.gz"%vcf1AfterDepthFilter.name)
					vcf1AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf1AfterDepthFilter.name)
					vcf1FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf1FilterByDepthJob, inputF=vcf1AfterDepthFilter, outputF=vcf1AfterDepthFilterGzip, \
							transferOutput=False)
					
					vcf2AfterDepthFilter = File(os.path.join(vcf2DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
					vcf2FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepth2Java, \
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
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
							minMAC=minMAC, minMAF=minMAF, \
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
							minMAC=minMAC, minMAF=minMAF, \
							maxSNPMissingRate=maxSNPMissingRate,\
							extraDependentInputLs=[vcf2FilterByDepthBGZipTabixJob.tbi_F], extraArguments="--recode-INFO-all")
					
					vcf2FilterByvcftoolsGzip = File("%s.gz"%vcf2FilterByvcftoolsJob.output.name)
					vcf2FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=vcf2FilterByvcftoolsJob, inputF=vcf2FilterByvcftoolsJob.output, outputF=vcf2FilterByvcftoolsGzip, \
							transferOutput=True)
					
					counter += 9
		sys.stderr.write("%s jobs.\n"%(counter+1))
	
	def addJobsToFilterOneVCFDir(self, workflow, vcf1Dir, refFastaFList, alnStatForFilterF, keepSNPPosF, onlyKeepBiAllelicSNP=True,\
							minMAC=None, minMAF=None, maxSNPMissingRate=None):
		"""
		2012.5.10
			add extraArguments="--recode-INFO-all" to addFilterJobByvcftools()
		2012.1.14
		"""
		
		refFastaF = refFastaFList[0]
		
		#a relative-path name for vcf1Dir
		vcf1Name = self.findProperVCFDirIdentifier(vcf1Dir, defaultName='vcf1')
		
		vcf1DepthFilterDir = "%s_DepthFilter"%(vcf1Name)
		vcf1DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1DepthFilterDir)
		
		SNPMismatchStatDir = "SNPMismatchStat"
		SNPMismatchStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=SNPMismatchStatDir)
		
		vcf1_vcftoolsFilterDir = "%s_vcftoolsFilter"%(vcf1Name)
		vcf1_vcftoolsFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1_vcftoolsFilterDir)
		
		#import re
		#chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		counter = 0
		no_of_vcf_files = 0
		for inputFname in os.listdir(vcf1Dir):
			vcf1AbsPath = os.path.join(os.path.abspath(vcf1Dir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and \
					not NextGenSeq.isVCFFileEmpty(vcf1AbsPath, checkContent=self.checkEmptyVCFByReading):
				#chr = chr_pattern.search(inputFname).group(1)
				no_of_vcf_files += 1
				if no_of_vcf_files%100==0:
					sys.stderr.write("%s%s VCFs. "%('\x08'*40, no_of_vcf_files))
				commonPrefix = inputFname.split('.')[0]
				vcf1 = File(os.path.join(vcf1Name, inputFname))	#relative path
				vcf1.absPath = vcf1AbsPath
				self.registerVCFAndItsTabixIndex(workflow, vcf1, input_site_handler)
				
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
				
				
					vcf1filterByDepthInput=vcf1KeepGivenSNPByvcftoolsJob.output
					vcf1filterByDepthParentJobLs=[vcf1DepthFilterDirJob, vcf1KeepGivenSNPByvcftoolsJob]
					vcf1filterByDepthExtraDependentInputLs=[]
					counter += 2
					continue	#skip the rest
				else:
					vcf1filterByDepthInput=vcf1
					vcf1filterByDepthParentJobLs=[vcf1DepthFilterDirJob]
					vcf1filterByDepthExtraDependentInputLs=[vcf1.tbi_F]
				vcf1AfterDepthFilter = File(os.path.join(vcf1DepthFilterDir, '%s.depthFiltered.vcf'%(commonPrefix)))
				vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						refFastaFList=refFastaFList, inputVCFF=vcf1filterByDepthInput, outputVCFF=vcf1AfterDepthFilter, \
						parentJobLs=vcf1filterByDepthParentJobLs, \
						alnStatForFilterF=alnStatForFilterF, \
						extraDependentInputLs=vcf1filterByDepthExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
			
				vcf1AfterDepthFilterGzip = File("%s.gz"%vcf1AfterDepthFilter.name)
				vcf1AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf1AfterDepthFilter.name)
				vcf1FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByDepthJob, inputF=vcf1AfterDepthFilter, outputF=vcf1AfterDepthFilterGzip, \
						transferOutput=False)
				
				outputFnamePrefix = os.path.join(vcf1_vcftoolsFilterDir, '%s.filter_by_vcftools'%(commonPrefix))
				vcf1FilterByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=vcf1AfterDepthFilterGzip, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[vcf1FilterByDepthBGZipTabixJob, vcf1_vcftoolsFilterDirJob], \
						snpMisMatchStatFile=None, \
						minMAC=minMAC, minMAF=minMAF, \
						maxSNPMissingRate=maxSNPMissingRate,\
						extraDependentInputLs=[vcf1FilterByDepthBGZipTabixJob.tbi_F], extraArguments="--recode-INFO-all")
				
				vcf1FilterByvcftoolsGzip = File("%s.gz"%vcf1FilterByvcftoolsJob.output.name)
				vcf1FilterByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByvcftoolsJob, inputF=vcf1FilterByvcftoolsJob.output, outputF=vcf1FilterByvcftoolsGzip, \
						transferOutput=True)
				
				counter += 9
		sys.stderr.write("%s jobs.\n"%(counter+1))
	
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
		
		workflow = self.initiateWorkflow()
		
		self.registerJars(workflow)
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
						input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		
		self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
		alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
		
		if self.keepSNPPosFname:
			keepSNPPosF = File(os.path.basename(self.keepSNPPosFname))	#relative path
			keepSNPPosF.absPath = os.path.abspath(self.keepSNPPosFname)
			keepSNPPosF.addPFN(PFN("file://" + keepSNPPosF.absPath, self.input_site_handler))
			workflow.addFile(keepSNPPosF)
		else:
			keepSNPPosF = None
		
		if self.vcf1Dir and self.vcf2Dir:
			self.addJobsToFilterTwoVCFDir(workflow, self.vcf1Dir, self.vcf2Dir, refFastaFList, alnStatForFilterF, keepSNPPosF, \
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, minMAC=self.minMAC, minMAF=self.minMAF, \
						maxSNPMissingRate=self.maxSNPMissingRate)
		elif self.vcf1Dir:
			# 2012.5.1 filter only on the 1st vcf folder
			self.addJobsToFilterOneVCFDir(workflow, self.vcf1Dir, refFastaFList, alnStatForFilterF, keepSNPPosF, \
									onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
									minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate)
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = FilterVCFPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
