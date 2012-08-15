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
		-G1 -K -R 0.1 -n 5 -L 0.25 -a 524 -C 50 -E
	
	#2011.12.9 to remove SNPs that are not in a file. no other filters.
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s -i $dirPrefix\gatk -I $dirPrefix\call/ -l condorpool -j condorpool
		-o Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4.xml
		-z uclaOfficeTemp -u yh -q ./alnStatForFilter.2011.12.9T0207.tsv -G0 -R 1 -n 0 -L 1 -A 100000 -a 524 -C 50
		-S ./4HighCovVRC_inter_minMAC4_vs_LowPass_top7559Contigs_no12eVarFilter_inter.2011.12.9T0107/overlapPos.tsv
	
	#2011.12.19 run on hoffman2's condorpool
	%s ....
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		
	#2012.5.1 filter trioCaller output with total depth (no minGQ filter anymore), minMAC=10 (-n 10), maxSNPMismatchRate=1 (-L 1.0)
	# minMAF=0.05 (-f 0.05), no depth-band filter (-A 0)
	%s -i AlignmentToTrioCallPipeline_VRC_top7559Contigs.2011.12.15T0057/trioCaller/ -l condorpool -j condorpool
		-z uclaOffice -u yh -q ./alnStatForFilter.2012.5.1T1430.tsv  -n 10 -L 1.0 -f 0.05 -A 0 -a 524 -C 50 -E
		-o FilterVCF_trioCallerTop7559Contigs.xml
	
	#2012.8.1 FilterGenotypeMethod5_ByMethod7Sites (-S ...) NoDepthFilter (-A 0) MaxSNPMissing0.5 (-L 0.5)
	%s -i ~/NetworkData/vervet/db/genotype_file/method_5/ -q ./aux/alnStatForFilter.2012.7.30T1542.tsv
		-L 0.5 -a 524  -E -o workflow/FilterGenotypeMethod5_ByMethod7Sites_NoDepthFilter_MaxSNPMissing0.5.xml  -l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh -C 5 -S ./method7_sites.tsv -A 0
	
	#2012.8.1 FilterGenotypeMethod6_ByMaskingZeroDPSite (-Z 1) 2FoldDepthFilter (-A 2) MaxSNPMissing1.0 (-L 1.0)
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	%s -i ~/NetworkData/vervet/db/genotype_file/method_6/ -q ./aux/alnStatForFilter.2012.8.1T1805.tsv
		-L 1.0 -a 524  -E -o workflow/FilterGenotypeMethod6_ByMaskingZeroDPSite_2FoldDepthFilter_MaxSNPMissing1.0.xml 
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -C 5 -Z 1 -A 2
		#-V 90 -x 100 
	
Description:
	2012.7.30 set depthFoldChange=0 to skip the filter by depth.
	2011-11-7 pipeline that filters VCF by depth, GQ, MAC, SNP missing rate, mismatch rate between two input VCFs
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0]\
				, sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from AbstractVervetWorkflow import AbstractVervetWorkflow
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow

class FilterVCFPipeline(AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVervetWorkflow.option_default_dict.copy()
	common_filter_option_dict = {
						("onlyKeepBiAllelicSNP", 0, int): [0, 'K', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 0, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						('keepSNPPosFname', 0, ): ['', 'S', 1, 'a tab-delimited file (optional) with 1st two columns as chr and pos. \
			to filter out SNPs that are absent in this file.\
			this step will be applied before all other filter jobs.'],\
						}
	option_default_dict.update(common_filter_option_dict)
	option_default_dict.update({
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype. 2012.5.1 no longer enforced in FilterVCFByDepth.java', ],\
						('depthFoldChange', 1, float): [2.0, 'A', 1, 'a variant is retained if its depth within this fold change of meanDepth,\
				set this to 0 or below to eliminate this step of filtering.', ],\
						("maxSNPMismatchRate", 0, float): [0, 'R', 1, 'maximum SNP mismatch rate between two vcf calls'],\
						("minDepthPerGenotype", 0, int): [0, 'Z', 1, 'mask genotype with below this depth as ./. (other fields retained), \
	esp. necessary for SAMtools, which output homozygous reference if no read for one sample.'],\
						("minMAC", 0, int): [None, 'n', 1, 'minimum MinorAlleleCount (by chromosome)'],\
						("minMAF", 0, float): [None, 'f', 1, 'minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRate", 0, float): [1.0, 'L', 1, 'maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('vcf1Dir', 1, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files. If not provided, filter vcf1Dir without applying maxSNPMismatchRate filter.', ],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVervetWorkflow.__init__(self, **keywords)
		# 2012.8.3 relative path causes stage-in problem as the pegasus program can't find the files.
		if getattr(self, 'vcf1Dir', None):
			self.vcf1Dir = os.path.abspath(self.vcf1Dir)
		if getattr(self, 'vcf2Dir', None):
			self.vcf2Dir = os.path.abspath(self.vcf2Dir)
		if getattr(self, 'depthFoldChange', None) and self.depthFoldChange>0 and not self.alnStatForFilterFname:
			sys.stderr.write("Error: alnStatForFilterFname (%s) is nothing while depthFoldChange=%s.\n"%\
							(self.alnStatForFilterFname, self.depthFoldChange))
			sys.exit(3)
	
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
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
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
	
	def addJobsToFilterOneVCFDir(self, workflow=None, inputData=None, refFastaFList=None, alnStatForFilterF=None, keepSNPPosF=None, \
						onlyKeepBiAllelicSNP=True,\
						minDepthPerGenotype=False, minMAC=None, minMAF=None, maxSNPMissingRate=None, outputDirPrefix="",\
						transferOutput=True):
		"""
		2012.7.30 add stat collecting jobs
		2012.5.10
			add extraArguments="--recode-INFO-all" to addFilterJobByvcftools()
		2012.1.14
		"""
		sys.stderr.write("Adding filter-VCF jobs for %s vcf files ... \n"%(len(inputData.jobDataLs)))
		refFastaF = refFastaFList[0]
		no_of_jobs = 0

		vcf1DepthFilterDir = "%s_DepthFilter"%(outputDirPrefix)
		vcf1DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1DepthFilterDir)
		
		SNPMismatchStatDir = "%s_SNPMismatchStat"%(outputDirPrefix)
		SNPMismatchStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=SNPMismatchStatDir)
		
		filterDir = "%s_vcftoolsFilter"%(outputDirPrefix)
		filterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=filterDir)
		
		topOutputDir = "%sFilterStat"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 4
		
		if keepSNPPosF:
			filterByGivenSitesStatMergeFile = File(os.path.join(topOutputDir, 'filterByGivenSitesStat_s1.tsv'))
			filterByGivenSitesStatMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByGivenSitesStatMergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			no_of_jobs += 1
		if minDepthPerGenotype:
			filterByMinDP1MergeFile = File(os.path.join(topOutputDir, 'filterByMinDP1.tsv_s2'))
			filterByMinDP1MergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinDP1MergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			no_of_jobs += 1
		if alnStatForFilterF:
			filterByDepthStatMergeFile = File(os.path.join(topOutputDir, 'filterByDepthStat_s3.tsv'))
			filterByDepthStatMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByDepthStatMergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
			no_of_jobs += 1
		
		
		if minMAC is not None:
			filterByMinMACMergeFile = File(os.path.join(topOutputDir, 'filterByMinMACStat_s4.tsv'))
			filterByMinMACMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinMACMergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
								extraArguments="-k 1 -v 2-4")
			no_of_jobs += 1
		
		if minMAF is not None:
			filterByMinMAFMergeFile = File(os.path.join(topOutputDir, 'filterByMinMAFStat_s5.tsv'))
			filterByMinMAFMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByMinMAFMergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
								extraArguments="-k 1 -v 2-4")
			no_of_jobs += 1
		if maxSNPMissingRate is not None and maxSNPMissingRate>=0 and maxSNPMissingRate<1:
			filterByMaxSNPMissingRateMergeFile = File(os.path.join(topOutputDir, 'filterByMaxSNPMissingRateStat_s6.tsv'))
			filterByMaxSNPMissingRateMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
							outputF=filterByMaxSNPMissingRateMergeFile, transferOutput=transferOutput, parentJobLs=[topOutputDirJob],\
							extraArguments="-k 1 -v 2-4")
			no_of_jobs += 1
		
		input_site_handler = self.input_site_handler
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		counter = 0
		no_of_vcf_files = 0
		for jobData in inputData.jobDataLs:
			inputF = jobData.vcfFile
			tbi_F = jobData.tbi_F
			inputJobLs = jobData.jobLs
			chromosome = self.getChrFromFname(os.path.basename(inputF.name))
			
			no_of_vcf_files += 1
			if no_of_vcf_files%100==0:
				sys.stderr.write("%s%s VCFs. "%('\x08'*40, no_of_vcf_files))
			commonPrefix = os.path.basename(inputF.name).split('.')[0]
			
			lastRoundJobLs= inputJobLs
			lastBGZipTabixJob = PassingData(output=inputF, tbi_F=tbi_F)	#2012.8.3 fake, not a job. only useful when all filtering jobs are skipped.
			lastRoundExtraDependentInputLs =[]
			
			if keepSNPPosF:
				#toss SNPs that are not in this keepSNPPosFname file
				outputFnamePrefix = os.path.join(filterDir, '%s.keepGivenSNP'%(commonPrefix))
				vcf1KeepGivenSNPByvcftoolsJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=inputF, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[filterDirJob] + inputJobLs, \
						snpMisMatchStatFile=keepSNPPosF, \
						minMAC=None, minMAF=None, \
						maxSNPMissingRate=None,\
						extraDependentInputLs=[tbi_F], extraArguments="--recode-INFO-all")
				
				vcf1KeepGivenSNPByvcftoolsGzip = File("%s.gz"%vcf1KeepGivenSNPByvcftoolsJob.output.name)
				vcf1KeepGivenSNPByvcftoolsBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
					parentJob=vcf1KeepGivenSNPByvcftoolsJob, inputF=vcf1KeepGivenSNPByvcftoolsJob.output, \
					outputF=vcf1KeepGivenSNPByvcftoolsGzip, \
					transferOutput=False)
				
				currentBGZipTabixJob = vcf1KeepGivenSNPByvcftoolsBGZipTabixJob
				#check how much sites got filtered
				outputF = File(os.path.join(topOutputDir, '%s.filterByGivenSitesStat.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=inputF, \
									currentBGZipTabixJob=currentBGZipTabixJob,\
									statMergeJob=filterByGivenSitesStatMergeJob, parentJobLs=[topOutputDirJob]+inputJobLs)
			
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs=[vcf1DepthFilterDirJob, vcf1KeepGivenSNPByvcftoolsBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				no_of_jobs += 3
			else:
				lastRoundJobLs=[ vcf1DepthFilterDirJob] + inputJobLs
				lastRoundExtraDependentInputLs =[tbi_F]
				lastBGZipTabixJob = vcf1DepthFilterDirJob
				lastBGZipTabixJob.output = inputF	#faking it
			
			#2012.8.1 mask zero-depth sites
			if minDepthPerGenotype:
				outputFnamePrefix = os.path.join(vcf1DepthFilterDir, '%s.minDP%s'%(commonPrefix, minDepthPerGenotype))
				maskZeroDepthGenotypeAsMissingJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=lastBGZipTabixJob.output, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=lastRoundJobLs, \
						minMAC=None, minMAF=None, \
						maxSNPMissingRate=None,\
						extraDependentInputLs=lastRoundExtraDependentInputLs, outputFormat='--recode', \
						extraArguments="--recode-INFO-all --minDP %s"%(minDepthPerGenotype))
				
				maskVCFGzipFile = File("%s.gz"%maskZeroDepthGenotypeAsMissingJob.output.name)
				maskVCFGzipJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
					parentJob=maskZeroDepthGenotypeAsMissingJob, inputF=maskZeroDepthGenotypeAsMissingJob.output, \
					outputF=maskVCFGzipFile, \
					transferOutput=False)
				
				currentBGZipTabixJob = maskVCFGzipJob
				#check how much sites got filtered
				outputF = File(os.path.join(vcf1DepthFilterDir, '%s.filterByminDP1.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastBGZipTabixJob.output, \
									currentBGZipTabixJob=currentBGZipTabixJob,\
									statMergeJob=filterByMinDP1MergeJob, parentJobLs=lastRoundJobLs)
			
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs=[currentBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				no_of_jobs += 3
				
			if alnStatForFilterF:
				vcf1AfterDepthFilter = File(os.path.join(vcf1DepthFilterDir, '%s.filterByDepth.vcf'%(commonPrefix)))
				vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
						genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
						refFastaFList=refFastaFList, inputVCFF=lastBGZipTabixJob.output, outputVCFF=vcf1AfterDepthFilter, \
						parentJobLs=lastRoundJobLs, \
						alnStatForFilterF=alnStatForFilterF, \
						extraDependentInputLs=lastRoundExtraDependentInputLs, onlyKeepBiAllelicSNP=onlyKeepBiAllelicSNP)
			
			
				vcf1AfterDepthFilterGzip = File("%s.gz"%vcf1AfterDepthFilter.name)
				vcf1AfterDepthFilterGzip_tbi_F = File("%s.gz.tbi"%vcf1AfterDepthFilter.name)
				vcf1FilterByDepthBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByDepthJob, inputF=vcf1FilterByDepthJob.output, outputF=vcf1AfterDepthFilterGzip, \
						transferOutput=False)
				currentBGZipTabixJob = vcf1FilterByDepthBGZipTabixJob
				#check how much sites got filtered by depth filter
				outputF = File(os.path.join(vcf1DepthFilterDir, '%s.filterByDepthStat.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
										vcf1=lastBGZipTabixJob.output, parentJobLs=lastRoundJobLs, \
										currentBGZipTabixJob=currentBGZipTabixJob,\
										statMergeJob=filterByDepthStatMergeJob)
				
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs = [currentBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				
				no_of_jobs += 3
			
				
			if minMAC is not None:
				outputFnamePrefix = os.path.join(filterDir, '%s.filterByMinMAC'%(commonPrefix))
				vcf1FilterByMinMACJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=lastBGZipTabixJob.output, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[filterDirJob] + lastRoundJobLs, \
						snpMisMatchStatFile=None, \
						minMAC=minMAC, minMAF=None, \
						maxSNPMissingRate=None,\
						extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
				
				#check how much sites got filtered by maxSNPMissingRate filter
				
				vcf1FilterByMinMACGzip = File("%s.gz"%vcf1FilterByMinMACJob.output.name)
				vcf1FilterByMinMACBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByMinMACJob, inputF=vcf1FilterByMinMACJob.output, outputF=vcf1FilterByMinMACGzip, \
						transferOutput=False)
				currentBGZipTabixJob = vcf1FilterByMinMACBGZipTabixJob
				outputF = File(os.path.join(filterDir, '%s.filterByMinMACStat.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastBGZipTabixJob.output, currentBGZipTabixJob=currentBGZipTabixJob,\
									parentJobLs=lastRoundJobLs, statMergeJob=filterByMinMACMergeJob)
				
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs = [currentBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				no_of_jobs += 3
			
			if minMAF is not None:
				outputFnamePrefix = os.path.join(filterDir, '%s.filterByMinMAF'%(commonPrefix))
				vcf1FilterByMinMAFJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=lastBGZipTabixJob.output, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[filterDirJob] + lastRoundJobLs, \
						snpMisMatchStatFile=None, \
						minMAC=None, minMAF=minMAF, \
						maxSNPMissingRate=None,\
						extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
				
				#check how much sites got filtered by maxSNPMissingRate filter
				
				vcf1FilterByMinMAFGzip = File("%s.gz"%vcf1FilterByMinMAFJob.output.name)
				vcf1FilterByMinMAFBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByMinMAFJob, inputF=vcf1FilterByMinMAFJob.output, outputF=vcf1FilterByMinMAFGzip, \
						transferOutput=False)
				
				currentBGZipTabixJob = vcf1FilterByMinMAFBGZipTabixJob
				outputF = File(os.path.join(filterDir, '%s.filterByMinMAFStat.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									vcf1=lastBGZipTabixJob.output, currentBGZipTabixJob=currentBGZipTabixJob,\
									statMergeJob=filterByMinMAFMergeJob, parentJobLs=lastRoundJobLs)
				
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs = [currentBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				no_of_jobs += 3
			
			if maxSNPMissingRate is not None and maxSNPMissingRate>=0 and maxSNPMissingRate<1:
				outputFnamePrefix = os.path.join(filterDir, '%s.filterByMaxSNPMissingRate'%(commonPrefix))
				vcf1FilterByMaxSNPMissingRateJob = self.addFilterJobByvcftools(workflow, vcftoolsWrapper=workflow.vcftoolsWrapper, \
						inputVCFF=lastBGZipTabixJob.output, \
						outputFnamePrefix=outputFnamePrefix, \
						parentJobLs=[filterDirJob] + lastRoundJobLs, \
						snpMisMatchStatFile=None, \
						minMAC=None, minMAF=None, \
						maxSNPMissingRate=maxSNPMissingRate,\
						extraDependentInputLs=lastRoundExtraDependentInputLs, extraArguments="--recode-INFO-all")
				
				#check how much sites got filtered by maxSNPMissingRate filter
				
				vcf1FilterByMaxSNPMissingRateGzip = File("%s.gz"%vcf1FilterByMaxSNPMissingRateJob.output.name)
				vcf1FilterByMaxSNPMissingRateGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
						parentJob=vcf1FilterByMaxSNPMissingRateJob, inputF=vcf1FilterByMaxSNPMissingRateJob.output, outputF=vcf1FilterByMaxSNPMissingRateGzip, \
						transferOutput=transferOutput)
				currentBGZipTabixJob = vcf1FilterByMaxSNPMissingRateGZipTabixJob
				outputF = File(os.path.join(filterDir, '%s.filterByMaxSNPMissingRateStat.tsv'%(commonPrefix)))
				self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
										vcf1=lastBGZipTabixJob.output, currentBGZipTabixJob=currentBGZipTabixJob, \
										statMergeJob=filterByMaxSNPMissingRateMergeJob, parentJobLs=lastRoundJobLs)
				
				lastBGZipTabixJob = currentBGZipTabixJob
				lastRoundJobLs = [currentBGZipTabixJob]
				lastRoundExtraDependentInputLs=[currentBGZipTabixJob.tbi_F]
				no_of_jobs += 3
			
			lastBGZipTabixJobOutputFile = getattr(lastBGZipTabixJob, 'output', None)	#could be None if all filter jobs are skipped
			lastBGZipTabixJobTbiF = getattr(lastBGZipTabixJob, 'tbi_F', None) 
			returnData.jobDataLs.append(PassingData(jobLs=lastRoundJobLs, vcfFile=lastBGZipTabixJobOutputFile, \
									tbi_F=lastBGZipTabixJobTbiF, \
									fileList=[lastBGZipTabixJobOutputFile, lastBGZipTabixJobTbiF]))
		
		sys.stderr.write("%s%s VCFs. "%('\x08'*40, no_of_vcf_files))
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		return returnData
	
	def addVCFBeforeAfterFilterStatJob(self, chromosome=None, outputF=None, vcf1=None, vcf2=None,\
									lastBGZipTabixJob=None, currentBGZipTabixJob=None,\
									statMergeJob=None, parentJobLs=None):
		"""
		2012.7.30
			examples:
			
			self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
									currentBGZipTabixJob=currentBGZipTabixJob, lastBGZipTabixJob=lastBGZipTabixJob,\
									statMergeJob=filterByMaxSNPMissingRateMergeJob)
		"""
		if vcf1 is None and lastBGZipTabixJob:
			vcf1 = lastBGZipTabixJob.output
		if vcf2 is None and currentBGZipTabixJob:
			vcf2 = currentBGZipTabixJob.output
		if parentJobLs is None:
			parentJobLs = []
		if lastBGZipTabixJob:
			parentJobLs.append(lastBGZipTabixJob)
		if currentBGZipTabixJob:
			parentJobLs.append(currentBGZipTabixJob)
		vcfFilterStatJob = self.addCheckTwoVCFOverlapJob(executable=self.CheckTwoVCFOverlap, \
										vcf1=vcf1, \
						vcf2=vcf2, chromosome=chromosome, chrLength=None, \
				outputF=outputF, parentJobLs=parentJobLs, \
				extraDependentInputLs=None, transferOutput=False, extraArguments=None, job_max_memory=1000, \
				perSampleMismatchFraction=False)
		self.addInputToStatMergeJob(statMergeJob=statMergeJob, \
							inputF=vcfFilterStatJob.output , \
							parentJobLs=[vcfFilterStatJob])
		return vcfFilterStatJob
			
	
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
		if self.depthFoldChange>0:
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
			alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
		else:
			alnStatForFilterF = None
		
		if self.keepSNPPosFname:
			keepSNPPosF = self.registerOneInputFile(inputFname=os.path.abspath(self.keepSNPPosFname),\
														folderName=self.pegasusFolderName)
		else:
			keepSNPPosF = None
		
		if self.vcf1Dir and self.vcf2Dir:
			self.addJobsToFilterTwoVCFDir(workflow, self.vcf1Dir, self.vcf2Dir, refFastaFList, alnStatForFilterF, keepSNPPosF, \
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP, minMAC=self.minMAC, minMAF=self.minMAF, \
						maxSNPMissingRate=self.maxSNPMissingRate)
		elif self.vcf1Dir:
			# 2012.5.1 filter only on the 1st vcf folder
			#a relative-path name for vcf1Dir
			vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir, defaultName='vcf1')
			inputData = self.registerAllInputFiles(workflow, self.vcf1Dir, input_site_handler=self.input_site_handler, \
												checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
												pegasusFolderName="%s_%s"%(self.pegasusFolderName, vcf1Name), \
												maxContigID=self.maxContigID, \
												minContigID=self.minContigID)
			self.addJobsToFilterOneVCFDir(workflow, inputData=inputData, refFastaFList=refFastaFList, \
									alnStatForFilterF=alnStatForFilterF, keepSNPPosF=keepSNPPosF, \
									onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
									minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate,\
									minDepthPerGenotype=self.minDepthPerGenotype, outputDirPrefix="")
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = FilterVCFPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
