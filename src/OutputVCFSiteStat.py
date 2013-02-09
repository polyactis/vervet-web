#!/usr/bin/env python
"""
Examples:
	
	#2011.11.28 (-K = keep only bi-allelic SNPs), -E = checkEmptyVCFByReading
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s -B ./TrioInconsistency_4HC101LCToCallLowPass_top7559Contigs_inter.2011.12.16T0428/avgTrioInconsistencyByPosition.tsv.gz
		-i $dirPrefix\call/ -I $dirPrefix\gatk/ -o OutputSiteStat_4HC101LC_maxContigID1000_minGQ1_2foldMedianDepth.xml
		-u yh -q ./alnStatForFilter.2011.12.21T0315.tsv -G1 -A 2 -K -a 524 -C 10
		-l condorpool -j condorpool
	
	#2011.12.9 to remove SNPs that are not in a file. no other filters.
	dirPrefix=./AlignmentToCallLowPass_top7559Contigs_no12eVarFilter_2011.11.23T1620/
	%s -i $dirPrefix\call/ -I $dirPrefix\gatk/ -l condorpool -j condorpool
		-o Keep_LowPass_top7559Contigs_no12eVarFilter_SNPs_PresentIn4HC_inter_minMAC4.xml
		-z uclaOfficeTemp -u yh -q ./alnStatForFilter.2011.12.9T0207.tsv -G0 -K -A 100000 -a 524 -C 10
		-S ./4HighCovVRC_inter_minMAC4_vs_LowPass_top7559Contigs_no12eVarFilter_inter.2011.12.9T0107/overlapPos.tsv
	
	#2011.12.19 run on hoffman2's condorpool
	%s ....
		-z localhost -l hcondor -j hcondor -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
	
Description:
	2011-12-21 a pegasus workflow that outputs variant statistics, mismatch rate with variants in another vcf folder,
		and trio inconsistency rate
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
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from pymodule.pegasus.AbstractNGSWorkflow import AbstractNGSWorkflow
from pymodule.utils import runLocalCommand
from vervet.src import VervetDB, AbstractVervetWorkflow
from vervet.src.qc.FilterVCFPipeline import FilterVCFPipeline

class OutputVCFSiteStat(FilterVCFPipeline):
	__doc__ = __doc__
	option_default_dict = AbstractNGSWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype', ],\
						('depthFoldChange', 1, float): [2.0, 'A', 1, 'a variant is retained if its depth within this fold change of meanDepth,', ],\
						('vcf1Dir', 0, ): ['', 'i', 1, 'main input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz for SNP mismatch rate calculation', ],\
						("onlyKeepBiAllelicSNP", 0, int): [0, 'K', 0, 'toggle this to remove all SNPs with >=3 alleles?'],\
						('alnStatForFilterFname', 1, ): ['', 'q', 1, 'The alignment stat file for FilterVCFByDepthJava. tab-delimited: individual_alignment.id minCoverage maxCoverage minGQ'],\
						('keepSNPPosFname', 0, ): ['', '', 1, 'a tab-delimited file (optional), chr pos 0, to pre-filter SNPs that are absent in this file.'],\
						('maxContigID', 1, int): [1000, 'm', 1, 'if contig ID is beyond this number, it will not be included', ],\
						('trioInconsistencyByPosistionFname', 1, ): ['', 'B', 1, 'This file contains the average trio/duo inconsistency rate for each SNP. bgzip&tabix-indexed', ],\
						
						})
	
	
	def __init__(self,  **keywords):
		"""
		"""
		FilterVCFPipeline.__init__(self, **keywords)
	
	def registerCustomExecutables(self, workflow):
		"""
		2011-11-28
		"""
		FilterVCFPipeline.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		ReduceMatrixByMergeColumnsWithSameKey = Executable(namespace=namespace, name="ReduceMatrixByMergeColumnsWithSameKey", \
								version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		ReduceMatrixByMergeColumnsWithSameKey.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "reducer/ReduceMatrixByMergeColumnsWithSameKey.py"), \
												site_handler))
		ReduceMatrixByMergeColumnsWithSameKey.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(ReduceMatrixByMergeColumnsWithSameKey)
		workflow.ReduceMatrixByMergeColumnsWithSameKey = ReduceMatrixByMergeColumnsWithSameKey
	
	def getChrListInTrioInconsistencyFile(self, tabixPath, trioInconsistencyByPosistionFname=None):
		"""
		2011.12.21
		"""
		sys.stderr.write("Getting list of chromosomes out of %s ..."%(trioInconsistencyByPosistionFname))
		chr_id_ls = []
		commandline = "%s -l %s"%(tabixPath, trioInconsistencyByPosistionFname)
		return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=False)
		for chr in return_data.output_stdout:
			chr_id_ls.append(chr.strip())
		sys.stderr.write(" %s chromosomes.\n"%(len(chr_id_ls)))
		return chr_id_ls
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = self.db_vervet
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		
		refFastaFname = os.path.join(self.data_dir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
						input_site_handler=self.input_site_handler,\
						checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.outputAlignmentDepthAndOthersForFilter(self.alnStatForFilterFname, ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
		alnStatForFilterF = self.registerOneInputFile(workflow, self.alnStatForFilterFname)
		
		#name to distinguish between vcf1Dir, and vcf2Dir
		vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir, defaultName='vcf1')
		vcf2Name = self.findProperVCFDirIdentifier(self.vcf2Dir, defaultName='vcf2')
		if vcf2Name==vcf1Name or not vcf2Name:
			vcf2Name = "vcf2"
		
		no_of_jobs = 0
		vcf1DepthFilterDir = "%s_DepthFilter"%(vcf1Name)
		vcf1DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf1DepthFilterDir)
		#vcf2DepthFilterDir = "%s_DepthFilter"%(vcf2Name)
		#vcf2DepthFilterDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=vcf2DepthFilterDir)
		
		trioInconsistencyDir = "trioInconsistency"
		trioInconsistencyDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=trioInconsistencyDir)
		
		
		SNPMismatchStatDir = "SNPMismatchStat"
		SNPMismatchStatDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=SNPMismatchStatDir)
		
		input_site_handler = self.input_site_handler
		
		
		#whole genome reduction job.
		wholeGenomeSiteStatFile = File('siteStatAndTrioInconsistency.tsv')
		wholeGenomeSiteStatMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=wholeGenomeSiteStatFile,transferOutput=False)
		
		wholeGenomeSiteStatBGzipFile = File("%s.gz"%wholeGenomeSiteStatFile.name)
		wholeGenomeSiteStatBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
							parentJob=wholeGenomeSiteStatMergeJob, inputF=wholeGenomeSiteStatFile, \
							outputF=wholeGenomeSiteStatBGzipFile, \
							transferOutput=True, tabixArguments="-s 1 -b 2 -e 2")
		no_of_jobs += 5
		
		#read the trioInconsistencyByPosistionFname and figure out how many contigs in it and add an extraction job for each contig
		chrLs = self.getChrListInTrioInconsistencyFile(self.tabixPath, self.trioInconsistencyByPosistionFname)
		chr2tabixRetrieveJob = {}
		trioInconsistencyByPosistionF = self.registerOneInputFile(workflow, self.trioInconsistencyByPosistionFname)
		trioInconsistencyByPosistion_tbi_Fname = '%s.tbi'%(self.trioInconsistencyByPosistionFname)
		trioInconsistencyByPosistion_tbi_F = self.registerOneInputFile(workflow, trioInconsistencyByPosistion_tbi_Fname)
		
		for chr in chrLs:
			outputF = File(os.path.join(trioInconsistencyDir, '%s.trioInconsistency.tsv'%chr))
			tabixRetrieveJob = self.addTabixRetrieveJob(workflow, executable=workflow.tabixRetrieve, tabixPath=self.tabixPath, \
							inputF=trioInconsistencyByPosistionF, outputF=outputF, regionOfInterest=chr, includeHeader=True,\
							parentJobLs=[trioInconsistencyDirJob], job_max_memory=100, extraDependentInputLs=[trioInconsistencyByPosistion_tbi_F], \
							transferOutput=False)
			chr2tabixRetrieveJob[chr] = tabixRetrieveJob
			no_of_jobs += 1
		
		counter = 0
		no_of_vcf = 0
		no_of_good_vcf = 0
		for inputFname in os.listdir(self.vcf1Dir):
			counter += 1
			if counter%500==0:
				sys.stderr.write("%s %s jobs %s good vcf, %s total vcf, %s total files"%('\x08'*180, no_of_jobs, \
														no_of_good_vcf, no_of_vcf, counter))
			
			vcf1AbsPath = os.path.join(os.path.abspath(self.vcf1Dir), inputFname)
			vcf2AbsPath = os.path.join(os.path.abspath(self.vcf2Dir), inputFname)
			if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False) and not NextGenSeq.isVCFFileEmpty(vcf1AbsPath):
				if not NextGenSeq.isVCFFileEmpty(vcf2AbsPath, checkContent=self.checkEmptyVCFByReading):	#make sure the samtools vcf exists
					no_of_vcf += 1
					chr = self.getChrFromFname(inputFname)
					if not chr or chr not in chr2tabixRetrieveJob:
						continue
					no_of_good_vcf += 1
					#find the contig id and  the matching tabix job
					commonPrefix = inputFname.split('.')[0]
					vcf1 = File(os.path.join(vcf1Name, inputFname))	#relative path
					vcf1.absPath = vcf1AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf1, input_site_handler)
					vcf2 = File(os.path.join(vcf2Name, inputFname))	#relative path
					vcf2.absPath = vcf2AbsPath
					self.registerVCFAndItsTabixIndex(workflow, vcf2, input_site_handler)
					
					outputSiteStatF = File(os.path.join(vcf1DepthFilterDir, '%s.siteStat.tsv'%(commonPrefix)))
					vcf1FilterByDepthJob = self.addFilterVCFByDepthJob(workflow, FilterVCFByDepthJava=workflow.FilterVCFByDepthJava, \
							genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
							refFastaFList=refFastaFList, inputVCFF=vcf1, outputVCFF=None, outputSiteStatF=outputSiteStatF,\
							parentJobLs=[vcf1DepthFilterDirJob], \
							alnStatForFilterF=alnStatForFilterF, \
							extraDependentInputLs=[vcf1.tbi_F], onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP)
					
					snpMisMatchStatFile = File(os.path.join(SNPMismatchStatDir, '%s_snpMismatchStat.tsv'%(os.path.splitext(commonPrefix)[0])))
					calculateSNPMismatchRateOfTwoVCFJob = self.addCalculateTwoVCFSNPMismatchRateJob(workflow, \
							executable=workflow.CalculateSNPMismatchRateOfTwoVCF, \
							vcf1=vcf1, vcf2=vcf2, snpMisMatchStatFile=snpMisMatchStatFile, \
							maxSNPMismatchRate=1.0, parentJobLs=[SNPMismatchStatDirJob], \
							job_max_memory=1000, extraDependentInputLs=[], \
							transferOutput=False)
					
					#add a ReduceMatrixByMergeColumnsWithSameKey job
					chrMergingStatF = File('%s_variantSiteStatAndTrioInconsistencyRate.tsv'%(chr))
					chrMergingStatJob = self.addStatMergeJob(workflow, \
									statMergeProgram=workflow.ReduceMatrixByMergeColumnsWithSameKey, \
									outputF=chrMergingStatF, extraArguments='-k 0,1', transferOutput=False)
					tabixRetrieveJob = chr2tabixRetrieveJob[chr]
					self.addInputToStatMergeJob(workflow, statMergeJob=chrMergingStatJob, \
								inputF=tabixRetrieveJob.output, \
								parentJobLs=[tabixRetrieveJob])
					
					self.addInputToStatMergeJob(workflow, statMergeJob=chrMergingStatJob, \
								inputF=outputSiteStatF, \
								parentJobLs=[vcf1FilterByDepthJob])
					self.addInputToStatMergeJob(workflow, statMergeJob=chrMergingStatJob, \
								inputF=snpMisMatchStatFile, \
								parentJobLs=[calculateSNPMismatchRateOfTwoVCFJob])
					
					#add to the whole genome reduction job
					self.addInputToStatMergeJob(workflow, statMergeJob=wholeGenomeSiteStatMergeJob, \
								inputF=chrMergingStatJob.output, \
								parentJobLs=[chrMergingStatJob])
					no_of_jobs += 3
		
		sys.stderr.write("%s %s jobs %s good vcf, %s total vcf, %s total files.\n"%('\x08'*180, no_of_jobs, \
														no_of_good_vcf, no_of_vcf, counter))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = OutputVCFSiteStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
