#!/usr/bin/env python
"""
Examples:
	#2012.8.30 run on hoffman2 condor pool (hcondor), filtered sequences (-Q 1), alignment method 2 (-G 2)
	# site==VRC filtering (-S 447), -C 1 (no job clustering)
	# add "-I 690-720" to restrict the alignment input (in addition to the filter above) 
	%s -a 524 -j hcondor -l hcondor -u yh -z localhost
		-o workflow/VerifyBAM/VerifyBAM_Site447Filtered1Ref524AlnMethod2.xml
		-H -C 1 -S 447 -Q 1 -G2 
		-e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java  -f ~/NetworkData/vervet/db/genotype_file/method_10_26996_VCF_998chromosomes_998VCFIntoOne.vcf.gz
		#-I 690-720
	
Description:
	2012.8.30
		a pegasus workflow that estimates mixture fraction of BAM file given genotype (frequency) from input VCF.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], )	#, sys.argv[0], sys.argv[0]

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from InspectAlignmentPipeline import InspectAlignmentPipeline

class VerifyBAMWorkflow(InspectAlignmentPipeline):
	__doc__ = __doc__
	option_default_dict = InspectAlignmentPipeline.commonOptionDict.copy()
	option_default_dict.update({
						("inputVCFFname", 1, ): [None, 'f', 1, 'VCF file which contains genotype for verifyBamID'],\
						("verifyBamIDPath", 0, ): ['%s/bin/verifyBamID', '', 1, 'path to the verifyBamID binary'],\
						("needPerContigJob", 0, int): [0, 'P', 0, 'toggle to add DepthOfCoverage and VariousReadCount jobs for each contig.'],\
						})

	def __init__(self, **keywords):
		"""
		2011-11-4
		"""
		self.pathToInsertHomePathList.append('verifyBamIDPath')
		InspectAlignmentPipeline.__init__(self, **keywords)
	
	
	def addVerifyBamIDJob(self, workflow=None, executable=None, inputVCF=None, inputBAM=None, outputFnamePrefix=None,\
				doSiteEstimation=False, doSelfEstimation=False, doBestEstimateion=False,\
				doFreeMix=False, doFreeFull=False, doFreeRefBias=False, doFreeNone=False, \
				doChipMix=None, doChipFull=None, doChipRefBias=None, doChipNone=None, \
				minAF=0.01, genoError=1e-03, minCallRate=0.50, \
				ignoreRG=False, ignoreOverlapPair=False, noEOF=False, precise=False, \
				minMapQ=10, maxDepth=20, minQ=13, maxQ=40, grid=None, \
				refRef=None, refHet=None, refAlt=None, \
				verbose=False, \
				parentJobLs=None, extraDependentInputLs=None, transferOutput=False, \
				extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.8.30
			
polyacti@n7204:~/NetworkData/vervet/vervetPipeline$ verifyBamID
verifyBamID 1.0.0 -- verify identity and purity of sequence data
(c) 2010 Hyun Min Kang, Goo Jun, and Goncalo Abecasis


The following parameters are available.  Ones with "[]" are in effect:

--minQ : minimum base quality to include
--maxQ : maximum base quality to cap
--site : If set, use only site information in the VCF and do not compare with the actual genotypes
--self : Only compare the ID-matching individuals between the VCF and BAM file
--best : Find the best matching individuals (.bestSM and .bestRG files will be produced). 
			This option is substantially longer than the default option

Available Options
                             Input Files : --vcf [], --bam [], --subset [],
                                           --smID []
                    VCF analysis options : --genoError [1.0e-03],
                                           --minAF [0.01],
                                           --minCallRate [0.50]
   Individuals to compare with chip data : --site, --self, --best
          Chip-free optimization options : --free-none, --free-mix [ON],
                                           --free-refBias, --free-full
          With-chip optimization options : --chip-none, --chip-mix [ON],
                                           --chip-refBias, --chip-full
                    BAM analysis options : --ignoreRG, --ignoreOverlapPair,
                                           --noEOF, --precise, --minMapQ [10],
                                           --maxDepth [20], --minQ [13],
                                           --maxQ [40], --grid [0.05]
                 Modeling Reference Bias : --refRef [1.00], --refHet [0.50],
                                           --refAlt [0.00]
                          Output options : --out [], --verbose

			
		"""
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		
		extraArgumentList = []
		extraOutputLs = []
		key2ObjectForJob = {}
		suffixAndNameTupleList = []	# a list of tuples , in each tuple, 1st element is the suffix. 2nd element is the proper name of the suffix.
			#job.$nameFile will be the way to access the file.
			#if 2nd element (name) is missing, suffix[1:].replace('.', '_') is the name (dot replaced by _) 
		
		if inputVCF:
			extraDependentInputLs.append(inputVCF)
			extraArgumentList.extend(["--vcf", inputVCF])
		if inputBAM:
			extraDependentInputLs.append(inputBAM)
			extraArgumentList.extend(["--bam", inputBAM])
		
		if outputFnamePrefix:
			extraArgumentList.extend(["--out", outputFnamePrefix])
			suffixAndNameTupleList.append(['.log', 'log'])
			suffixAndNameTupleList.extend([['.selfSM', ], ['.selfRG', ], ['.depthSM', ], ['.depthRG', ]])
			
		else:
			sys.stderr.write("outputFnamePrefix (=%s now) must be provided.\n"%(outputFnamePrefix))	
			sys.exit(3)
		
		if doSiteEstimation:
			extraArgumentList.append("--site")
		elif doSelfEstimation:
			extraArgumentList.append("--self")
		elif doBestEstimateion:
			extraArgumentList.append("--best")
			suffixAndNameTupleList.extend([['.bestSM', ], ['.bestRG', ]])
		
		
		if doFreeMix:
			extraArgumentList.append("--free-mix")
		elif doFreeFull:
			extraArgumentList.append("--free-full")
		elif doFreeRefBias:
			extraArgumentList.append("--free-refBias")
		elif doFreeNone:
			extraArgumentList.append("--free-none")
		
		if doChipMix:
			extraArgumentList.append("--chip-mix")
		elif doChipFull:
			extraArgumentList.append("--chip-full")
		elif doChipRefBias:
			extraArgumentList.append("--chip-refBias")
		elif doChipNone:
			extraArgumentList.append("--chip-none")
			
		if minAF is not None:
			extraArgumentList.append("--minAF %s"%(minAF))
		if genoError is not None:
			extraArgumentList.append("--genoError %s"%(genoError))
		if minCallRate is not None:
			extraArgumentList.append("--minCallRate %s"%(minCallRate))
		if ignoreRG:
			extraArgumentList.append("--ignoreRG")
		if ignoreOverlapPair:
			extraArgumentList.append("--ignoreOverlapPair")
		if noEOF:
			extraArgumentList.append("--noEOF")
		if precise:
			extraArgumentList.append("--precise")
		if minMapQ is not None:
			extraArgumentList.append("--minMapQ %s"%(minMapQ))
		if maxDepth is not None:
			extraArgumentList.append("--maxDepth %s"%(int(maxDepth)))	#must convert it to integer. verifyBamID won't take 3.0, (=3). 
		if minQ is not None:
			extraArgumentList.append("--minQ %s"%(minQ))
		if maxQ is not None:
			extraArgumentList.append("--maxQ %s"%(maxQ))
		if grid is not None:
			extraArgumentList.append("--grid %s"%(grid))
		if refRef is not None:
			extraArgumentList.append("--refRef %s"%(refRef))
		if refHet is not None:
			extraArgumentList.append("--refHet %s"%(refHet))
		if refAlt is not None:
			extraArgumentList.append("--refAlt %s"%(refAlt))
		if verbose:
			extraArgumentList.append("--verbose")
			
		if extraArguments:
			extraArgumentList.append(extraArguments)
		
		self.setupMoreOutputAccordingToSuffixAndNameTupleList(outputFnamePrefix=outputFnamePrefix, suffixAndNameTupleList=suffixAndNameTupleList, \
													extraOutputLs=extraOutputLs, key2ObjectForJob=key2ObjectForJob)
		
		job= self.addGenericJob(executable=executable, inputFile=None, outputFile=None, \
				parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
				extraOutputLs=extraOutputLs,\
				transferOutput=transferOutput, \
				extraArgumentList=extraArgumentList, key2ObjectForJob=key2ObjectForJob, job_max_memory=job_max_memory, **keywords)
		return job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-25
			split out of run()
		"""
		InspectAlignmentPipeline.registerCustomExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		verifyBamID = Executable(namespace=namespace, name="verifyBamID", version=version, os=operatingSystem,\
								arch=architecture, installed=True)
		verifyBamID.addPFN(PFN("file://" + self.verifyBamIDPath, site_handler))
		executableClusterSizeMultiplierList.append((verifyBamID, 1))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	def addJobs(self, workflow=None, alignmentDataLs=None, refName2size=None, inputVCF=None, verifyBamID=None, \
				dataDir=None, needPerContigJob=False, needSSHDBTunnel=0, outputDirPrefix="",\
				transferOutput=True):
		"""
		2012.8.30
			
		"""
		if workflow is None:
			workflow = self
		sys.stderr.write("Adding jobs for %s references and %s alignments..."%(len(refName2size), len(alignmentDataLs)))
		if len(alignmentDataLs)==0:
			sys.stderr.write("No alignment for verifyBamID. Exit now.\n")
			sys.exit(0)
		no_of_jobs = 0
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		topOutputDir = "%sverifyBAMOutput"%(outputDirPrefix)
		topOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		
		mergedOutputDir = "%smergedOutput"%(outputDirPrefix)
		mergedOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=mergedOutputDir)
		no_of_jobs += 1
		
		plotOutputDir = "%splot"%(outputDirPrefix)
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=plotOutputDir)
		no_of_jobs += 1
		
		selfSampleMixupMergeFile = File(os.path.join(mergedOutputDir, 'selfSMMerge.tsv'))
		selfSampleMixupMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=selfSampleMixupMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[selfSampleMixupMergeJob], \
												fileList=[selfSampleMixupMergeFile]))
		no_of_jobs += 1
		"""
		output of *.selfSM from verifyBamID
		
#SEQ_ID	RG	CHIP_ID	#SNPS	#READS	AVG_DP	FREEMIX	FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF	RDPHET  RDPALT
1968_3017_2005001_GA_vs_524	ALL	NA	24414	22222	0.91	0.00003 13245.89	13305.27	0.55489 0.05202 NA	NA	NA	NA	NA	NA	NA	NA
		
		"""
		outputFile = File( os.path.join(plotOutputDir, 'freeMix_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selfSampleMixupMergeFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFraction", \
					logWhichColumn=False, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 1
		
		outputFile = File( os.path.join(plotOutputDir, 'refProbGivenHet_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selfSampleMixupMergeFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREE_RH", whichColumnPlotLabel="refAlleleProbGivenHet", \
					logWhichColumn=False, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 1
		
		outputFile = File( os.path.join(plotOutputDir, 'refProbGivenAlternative_Hist.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addDrawHistogramJob(workflow=workflow, executable=workflow.DrawHistogram, inputFileList=[selfSampleMixupMergeFile], \
							outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREE_RA", whichColumnPlotLabel="refAlleleProbGivenAlternativeAllele", \
					logWhichColumn=False, positiveLog=True, logCount=True, valueForNonPositiveYValue=-1,\
					minNoOfTotal=10,\
					figureDPI=100, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs += 1

		outputFile = File(os.path.join(plotOutputDir, 'freeMix_vs_chipMix.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFractionByHet", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="CHIPMIX", xColumnPlotLabel="mixAgainstSelfGenotype", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		outputFile = File(os.path.join(plotOutputDir, 'freeMix_vs_AVG_DP.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFractionByHet", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="AVG_DP", xColumnPlotLabel="avgDepth", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		outputFile = File(os.path.join(plotOutputDir, 'chipMix_vs_AVG_DP.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="CHIPMIX", whichColumnPlotLabel="mixAgainstSelfGenotype", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="AVG_DP", xColumnPlotLabel="avgDepth", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		outputFile = File(os.path.join(plotOutputDir, 'freeMix_vs_refProbGivenHet.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFraction", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="FREE_RH", xColumnPlotLabel="refAlleleProbGivenHet", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		outputFile = File(os.path.join(plotOutputDir, 'freeMix_vs_refProbGivenAlternativeAllele.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFraction", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=50,\
					xColumnHeader="FREE_RA", xColumnPlotLabel="refAlleleProbGivenAlternativeAllele", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		outputFile = File(os.path.join(plotOutputDir, 'refAlleleProbGivenHet_vs_refProbGivenAlternativeAllele.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[selfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREE_RH", whichColumnPlotLabel="refAlleleProbGivenHet", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
					xColumnHeader="FREE_RA", xColumnPlotLabel="refAlleleProbGivenAlternativeAllele", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, selfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		#subtract the FREELIK1 - FREELIK0
		likSubstractedSelfSampleMixupMergeFile = File(os.path.join(mergedOutputDir, 'selfSMMerge_likDelta.tsv'))
		likSubstractedSelfSampleMixupMergeJob = self.addStatMergeJob(workflow, \
					statMergeProgram=workflow.ReduceMatrixBySumSameKeyColsAndThenDivide, \
					outputF=likSubstractedSelfSampleMixupMergeFile, transferOutput=False, \
					extraArguments="--operatorType 2 -k 0 -v 7,8,6")
					#column 0 is sample ID, column 6 is FREEMIX. 7 is FREELIK1, 8 is FRELIK0
		self.addInputToStatMergeJob(workflow, statMergeJob=likSubstractedSelfSampleMixupMergeJob, inputF=selfSampleMixupMergeJob.output,\
						parentJobLs=[selfSampleMixupMergeJob])
		returnData.jobDataLs.append(PassingData(jobLs=[likSubstractedSelfSampleMixupMergeJob], \
								fileList=[likSubstractedSelfSampleMixupMergeFile]))
		no_of_jobs += 1
		
		outputFile = File(os.path.join(plotOutputDir, 'freeMix_vs_deltaMinusLogLikelihood.png'))
		#no spaces or parenthesis or any other shell-vulnerable letters in the x or y axis labels (whichColumnPlotLabel, xColumnPlotLabel)
		self.addAbstractPlotJob(workflow=workflow, executable=workflow.AbstractPlot, \
					inputFileList=[likSubstractedSelfSampleMixupMergeJob.output], \
					outputFile=outputFile, \
					whichColumn=None, whichColumnHeader="FREEMIX", whichColumnPlotLabel="mixFraction", \
					logWhichColumn=False, positiveLog=True, valueForNonPositiveYValue=-1,\
					xColumnHeader="FREELK1_by_FREELK0", xColumnPlotLabel="deltaMinusLogLikelihood", \
					minNoOfTotal=5,\
					figureDPI=150, samplingRate=1,\
					parentJobLs=[plotOutputDirJob, likSubstractedSelfSampleMixupMergeJob], \
					extraDependentInputLs=None, \
					extraArguments=None, transferOutput=transferOutput,  job_max_memory=2000)
		no_of_jobs +=1
		
		
		selfRGMixupMergeFile = File(os.path.join(mergedOutputDir, 'selfRGMerge.tsv'))
		selfRGMixupMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=selfRGMixupMergeFile, transferOutput=False)
		returnData.jobDataLs.append(PassingData(jobLs=[selfRGMixupMergeJob], \
												fileList=[selfRGMixupMergeFile]))
		no_of_jobs += 1
		
		#alignmentId2RGJobDataLs = returnData.alignmentId2RGJobDataLs
		i =0
		for alignmentData in alignmentDataLs:
			alignment = alignmentData.alignment
			
			bamF= alignmentData.bamF
			baiF = alignmentData.baiF
			if i ==0:	#need at least one log file
				transferOutputForThisJob = True
			else:
				transferOutputForThisJob = False
			i += 1
			
			outputFnamePrefix = os.path.join(topOutputDir, alignment.getReadGroup())
			verifyBamIDJob = self.addVerifyBamIDJob(executable=self.verifyBamID, inputVCF=inputVCF, inputBAM=bamF, \
									outputFnamePrefix=outputFnamePrefix,\
				doFreeFull=True,\
				doChipMix=None, doChipFull=None, doChipRefBias=None, doChipNone=None, \
				minAF=0.01, genoError=1e-03, minCallRate=0.50, \
				minMapQ=20, maxDepth=int(3*alignment.median_depth), minQ=13, maxQ=40, \
				parentJobLs=[topOutputDirJob]+alignmentData.jobLs, extraDependentInputLs=[baiF], \
				transferOutput=transferOutputForThisJob, \
				extraArguments=None, job_max_memory=5000)
				# pass transferOutput to it so to keep log
			no_of_jobs += 1
			
			
			self.addInputToStatMergeJob(workflow, statMergeJob=selfSampleMixupMergeJob, inputF=verifyBamIDJob.selfSMFile,\
						parentJobLs=[verifyBamIDJob])
			self.addInputToStatMergeJob(workflow, statMergeJob=selfRGMixupMergeJob, inputF=verifyBamIDJob.selfRGFile,\
						parentJobLs=[verifyBamIDJob])
		
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		
		#2012.8.30 gzip the final output
		newReturnData = self.addGzipSubWorkflow(workflow=workflow, inputData=returnData, transferOutput=transferOutput,\
						outputDirPrefix="%smergedOutputGzip"%(outputDirPrefix))
		return newReturnData
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		if not self.dataDir:
			self.dataDir = db_vervet.data_dir
		
		if not self.localDataDir:
			self.localDataDir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
		
		if self.needPerContigJob:
			refName2size = self.getTopNumberOfContigs(self.contigMaxRankBySize)
			#refName2size = set(['Contig149'])	#temporary when testing Contig149
			#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		else:
			refName2size = {}
		
		alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, dataDir=self.localDataDir)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
												individual_site_id_set=set(self.site_id_ls))
		
		self.registerJars(workflow)
		
		self.registerExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, dataDir=self.dataDir)
		
		inputVCF = self.registerOneInputFile(workflow=workflow, inputFname=self.inputVCFFname, \
											folderName=self.pegasusFolderName)
		self.addJobs(alignmentDataLs=alignmentDataLs, refName2size=refName2size, inputVCF=inputVCF, \
				verifyBamID=self.verifyBamID, \
				dataDir=self.dataDir, needPerContigJob=self.needPerContigJob, needSSHDBTunnel=self.needSSHDBTunnel, \
				outputDirPrefix="", transferOutput=True)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = VerifyBAMWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
	
