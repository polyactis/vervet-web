#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 ...  -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -I ... -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	# set pop1ss/pop2ss to 0 (or negative) if u want to take all samples 
	m=64; country1=135; country2=151; pop1=VRC; pop2=Gambia; maxCID=100; pop1ss=10; pop2ss=10; n=10;
	%s  -I /u/home/eeskin/polyacti/NetworkData/vervet/db/genotype_file/method_$m/ -E -H 
		--pop1_country_id_ls $country1 --pop2_country_id_ls $country2 --pop1Header=$pop1 --pop2Header=$pop2
		-o workflow/popgen/Bootstrap_$n\Sampling_CompareAlleleFrequencyMethod$m\_sample$pop1ss\_$pop1\_vs_sample$pop2ss\_$pop2\_maxContigID$maxCID.xml
		-j hcondor -l hcondor -D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -z localhost -C 20 -x $maxCID -a 524 --pop1_sampleSize $pop1ss
		--pop2_sampleSize $pop2ss
		--plinkIBDCheckOutputFname PlinkIBDCheck/PlinkIBDCheck_Method38_W50Z20R0.5.2012.9.13T102232/ibdCheckIBDCheck/LDPrunedMerged_ibdCheck.tsv
		--maxIBDSharing 0.1 --minAbsDelta 0.2 --noOfSamplings $n;
	
Description:
	2012-10-03 workflow that investigates allele sharing between two different populations in VCF folder.
		
	#. each population is defined using tax-ID, site & country IDs
	#. split the VCF
	  #. extract individuals from individual population into separate VCF, each VCF with re-calculated AC/AF. ExtractSamplesFromVCF.py => run gatk's selectVariants to update AC/AF & DP, (run two this-kind jobs) 
	       #. arguments: tax-ID, site-ID, country-ID
	  #. JuxtaposeAlleleFrequencyFromMultiVCFInput.py VariantDiscovery.compareAlternativeAlleleFrequencyOfTwoCallSets(), output it in the format of Draw2DHistogramOfMatrix.py. header & 3 columns: AF_1, AF_2, count. 
	       #. arguments: two input file, 
	#. merge all the output and run Draw2DHistogramOfMatrix.py
	#. focus on SNPs whose AAF in either population is [0.1,0.9]
	     #. calculate how much around y=x, how much not
    
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from vervet.src.AbstractVervetWorkflow import AbstractVervetWorkflow
from CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder import CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder

class BootstrapCompareAlleleFrequencyOfTwoPopulation(CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder):
	__doc__ = __doc__
	option_default_dict = CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.option_default_dict.copy()
	option_default_dict.update({
			('noOfSamplings', 1, int): [5, '', 1, 'how many samplings to calculate the median AFS frequency correlation', ],\
			})
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.__init__(self, **keywords)
		
		"""
		listArgumentName_data_type_ls = [('pop1_site_id_ls', int), ("pop1_country_id_ls", int), \
								("pop1_tax_id_ls", int), \
								('pop2_site_id_ls', int), ("pop2_country_id_ls", int), \
								("pop2_tax_id_ls", int),]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[])
		"""
	def addAllJobs(self, workflow=None, inputVCFData=None, chr2IntervalDataLs=None, \
				GenomeAnalysisTKJar=None, samtools=None, \
				CreateSequenceDictionaryJava=None, CreateSequenceDictionaryJar=None, \
				BuildBamIndexFilesJava=None, BuildBamIndexFilesJar=None,\
				mv=None, \
				refFastaFList=None, \
				needFastaIndexJob=False, needFastaDictJob=False, \
				data_dir=None, no_of_gatk_threads = 1, \
				intervalSize=3000, intervalOverlapSize=0, \
				outputDirPrefix="", transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2012.10.15
			architect of the whole map-reduce framework
			call the parent's addAllJobs in a loop
		"""
		samplingReturnDataLs = []
		for i in xrange(self.noOfSamplings):
			oneSamplingReturnData = CompareAlleleFrequencyOfTwoPopulationFromOneVCFFolder.addAllJobs(self, \
					workflow=workflow, inputVCFData=inputVCFData, \
					chr2IntervalDataLs=chr2IntervalDataLs, samtools=samtools, \
				GenomeAnalysisTKJar=GenomeAnalysisTKJar, \
				CreateSequenceDictionaryJava=CreateSequenceDictionaryJava, CreateSequenceDictionaryJar=CreateSequenceDictionaryJar, \
				BuildBamIndexFilesJava=BuildBamIndexFilesJava, BuildBamIndexFilesJar=BuildBamIndexFilesJar,\
				mv=mv, \
				refFastaFList=refFastaFList,\
				needFastaIndexJob=needFastaIndexJob, needFastaDictJob=needFastaDictJob, \
				data_dir=data_dir, no_of_gatk_threads = 1, \
				intervalSize=intervalSize, intervalOverlapSize=intervalOverlapSize, \
				outputDirPrefix='%s_%s_'%(outputDirPrefix, i), transferOutput=transferOutput, job_max_memory=job_max_memory,\
				**keywords)
			samplingReturnDataLs.append(oneSamplingReturnData)
		
		topOutputDir = "%sFinalReduce"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		
		#a ReduceMatrixByAverageColumnsWithSameKey job
		outputFile = File(os.path.join(topOutputDir, 'medianAlleleSharingStatAcrossAllSampling.tsv'))
		medianReduceJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByAverageColumnsWithSameKey, \
						outputF=outputFile, extraArguments='-k 0 -v 1-8', parentJobLs=[topOutputDirJob], \
						extraDependentInputLs=None, transferOutput=True)
		#a MergeSameHeaderTablesIntoOne job
		outputFile = File(os.path.join(topOutputDir, 'alleleSharingStatAcrossAllSampling.tsv'))
		mergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.MergeSameHeaderTablesIntoOne, \
						outputF=outputFile, extraArguments=None, parentJobLs=[topOutputDirJob], \
						extraDependentInputLs=None, transferOutput=True)
		for oneSamplingReturnData in samplingReturnDataLs:
			self.addInputToStatMergeJob(workflow=workflow, statMergeJob=medianReduceJob, parentJobLs=[oneSamplingReturnData.estimateOutlierJob])
			self.addInputToStatMergeJob(workflow=workflow, statMergeJob=mergeJob, parentJobLs=[oneSamplingReturnData.estimateOutlierJob])
		
		outputFile = File( os.path.join(topOutputDirJob.output, 'outlierFraction_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[mergeJob.output], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="outlierFraction", whichColumnPlotLabel="outlierFraction", \
					logY=False, positiveLog=True, logCount=False, valueForNonPositiveYValue=-1,\
					minNoOfTotal=5,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[topOutputDirJob, mergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		outputFile = File( os.path.join(topOutputDirJob.output, 'AFS_cor_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[mergeJob.output], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="corr", whichColumnPlotLabel="AFSCorrelation", \
					logY=False, positiveLog=True, logCount=False, valueForNonPositiveYValue=-1,\
					minNoOfTotal=5,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[topOutputDirJob, mergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=True,  job_max_memory=2000)
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		
		
if __name__ == '__main__':
	main_class = BootstrapCompareAlleleFrequencyOfTwoPopulation
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()