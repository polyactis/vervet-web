#!/usr/bin/env python
"""
Examples:
	
	# 2011.12.14 run all VRC on hoffman2, top 7559 contigs. "--noOfCallingJobsPerNode 3" controls clustering of calling programs.
	# "--clusters_size 30" controls clustering for other programs., "-S 447" dictates monkeys from VRC
	%s -u yh -a 524 -s 2 -z localhost -o dags/AlignmentToCall/AlignmentToTrioCallPipeline_VRC_top7559Contigs.xml -j hcondor -l hcondor 
		-N7559 --noOfCallingJobsPerNode 3 --clusters_size 30 -S 447 -e /u/home/eeskin/polyacti/ --data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --needSSHDBTunnel
	
	# 2013.06.19
	%s -I /Network/Data/vervet/db/genotype_file/method_36 -u yh -a 524  -z localhost
		--plinkIBDCheckOutputFname 
		-o  dags/RefineCall/BeagleTrioCallerOnMethod36Contig96_100.xml -j hcondor -l hcondor
		--noOfCallingJobsPerNode 1 --clusters_size 1
		--data_dir /Network/Data/vervet/db/ --local_data_dir /Network/Data/vervet/db/
		--minContigID 96 --maxContigID 100
	
	# 2012.8.15 run TrioCaller on method 14 samtools calls, contig ID from 96 to 100 (--minContigID 96 --maxContigID 100)
	# sequenced filtered (--sequence_filtered 1), alignment by method 2 (--alignment_method_id 2), onlyKeepBiAllelicSNP (--onlyKeepBiAllelicSNP) 
	# calling job clusters size=1, others =1 (--noOfCallingJobsPerNode 1 --clusters_size 1)
	# add -Y (not guess #loci from 1st number in filename) if the input VCF is not db affiliated
	# 3000 (-Z 3000) sites per unit, 500 overlapping between units (-U 500)
	# add "--treatEveryOneIndependent" if you want to treat everyone independent (no mendelian constraints from TrioCaller)
	%s --run_type 2 -I ~/NetworkData/vervet/db/genotype_file/method_14/
		-u yh -a 524  -z localhost -o  dags/AlignmentToCall/TrioCallerOnMethod14Contig96_100.xml
		-j hcondor -l hcondor
		--noOfCallingJobsPerNode 1 --clusters_size 1 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--minContigID 96 --maxContigID 100 --sequence_filtered 1 --alignment_method_id 2  --onlyKeepBiAllelicSNP --needSSHDBTunnel
		# -Y -Z 3000 -U 500 --treatEveryOneIndependent
	
	#2013.1.4 run polymutt on method 36, no overlapping between intervals, --intervalOverlapSize 0
	#(because polymutt runs on loci one by one, no dependency)
	%s --run_type 3 -I ~/NetworkData/vervet/db/genotype_file/method_36/ -u yh -a 524
		-z localhost -o dags/AlignmentToCall/PolymuttOnMethod36Contig96_100.xml
		-j hcondor -l hcondor --noOfCallingJobsPerNode 1 --clusters_size 1 -e /u/home/eeskin/polyacti/
		--data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/ --local_data_dir /u/home/eeskin/polyacti/NetworkData/vervet/db/
		--minContigID 96 --maxContigID 100 --sequence_filtered 1 --alignment_method_id 2  --onlyKeepBiAllelicSNP --needSSHDBTunnel
		--intervalOverlapSize 500 --intervalSize 4000
	
Description:
	2013.06.14
		a program which generates a 3-stage pedigree genotype-calling pegasus workflow dag (xml file).
			#. Beagle phasing/imputing on high-coverage members
			#. Beagle phasing/imputing on all members with step 1's output as reference panel
			#. TrioCaller on step 2's output.
		
		Sample IDs in --plinkIBDCheckOutputFname are in original VCF form (alignment.read_group, no replicate tag).
"""
import sys, os
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from pymodule import VCFFile
from Pegasus.DAX3 import *
from vervet.src import VervetDB, AbstractVervetWorkflow
from vervet.src.genotyping.AlignmentToTrioCallPipeline import AlignmentToTrioCallPipeline

parentClass=AlignmentToTrioCallPipeline


class BeagleAndTrioCallerOnVCFWorkflow(AbstractVervetWorkflow, parentClass):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetWorkflow.option_default_dict)
	option_default_dict.update(parentClass.commonCallPipelineOptionDict)
	#option_default_dict.pop(('ind_aln_id_ls', 0, ))
	#option_default_dict.pop(('ind_seq_id_ls', 0, ))
	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update({
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						("minCoverageForRefPanel", 0, int): [8, '', 1, 'minimum coverage for a pedigree member to be in reference panel'],\
						("treatEveryOneIndependent", 0, int): [0, '', 0, 'toggle this to treat everyone in the pedigree independent (thus no replicates)'],\
						("run_type", 1, int): [1, '', 1, '1: not used right now'],\
						("pedigreeFileFormat", 1, int): [2, '', 1, 'passed to OutputVRCPedigreeInTFAMGivenOrderFromFile.py. 1: plink; 2: TrioCaller/Merlin/GATK; \n\
	3: for polymutt, replicate certain individuals to make pedigree loop-free'],\
						("minProbForValidCall", 1, float): [0.6, '', 1, 'minimum probability for a call to be regarded as valid'],\
						('plinkIBDCheckOutputFname', 1, ): [None, '', 1, 'file that contains IBD check result, PI_HAT=relatedness.\n\
	at least 3-columns with header: IID1, IID2, PI_HAT. IID1 and IID2 should match the whichColumn (whichColumnHeader) of inputFname.\n\
	The sampling will try to avoid sampling close pairs, PI_HAT(i,j)<=maxIBDSharing'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
	#2012.8.26 run_type 2 lower the default. In the addTrioCallerJobsONVCFFiles(), these two numbers refer to the number of sites, not base pairs. 
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 1000
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
	def __init__(self,  **keywords):
		"""
		2013.05
		"""
		parentClass.__init__(self, **keywords)
		self.extractRefPanelSampleIDJob = None	#this job only needs to be run once during mapEachInterval(). 
		#can't put it inside preReduce because it requires a VCF with replicates among them
		
		#self.trioCallerPath =  self.insertHomePath(self.trioCallerPath, self.home_path)
		self.needSplitChrIntervalData = False	#2013.06.21 turn this off before setup_run() to not construct chr2IntervalDataLs
	
	def preReduce(self, workflow=None, outputDirPrefix="", passingData=None, transferOutput=True, **keywords):
		"""
		2013.05.01
			1. a job that outputs the pedigree from db, with members from the VCF file. used by various filter programs and TrioCaller
			2. a job that extracts the high-coverage individuals from the VCF file
			3. figure out the existence of beagle unrelated cohort, trio cohort, pair/duo cohort for high-coverage group and all individuals
				need the pedigree graph, a VCF file => all sample IDs and only high-coverage individuals
				
		"""
		returnData = AbstractVervetWorkflow.preReduce(self, workflow=workflow, outputDirPrefix=outputDirPrefix, \
													passingData=passingData, transferOutput=transferOutput, **keywords)
		
		self.statDirJob = self.addMkDirJob(outputDir="%sStat"%(outputDirPrefix))
		self.phaseByGATKDirJob = self.addMkDirJob(outputDir="%sPhaseByGATK"%(outputDirPrefix))
		self.highCoveragePanelDirJob = self.addMkDirJob(outputDir="%sHighCoveragePanel"%(outputDirPrefix))
		self.auxDirJob = self.addMkDirJob(outputDir="%sAuxilliary"%(outputDirPrefix))
		
		# self.reduceOutputDirJob would contain non-replicate VCF files
		#this folder would store all the reduced VCF files with replicates among samles. 
		self.replicateVCFDirJob = self.addMkDirJob(outputDir="%sReplicateVCF"%(outputDirPrefix))
		
		self.plinkIBDCheckOutputFile = self.registerOneInputFile(inputFname=self.plinkIBDCheckOutputFname, \
										folderName='aux')
		
		# output pedigree to get pedigree file (for GATK, TrioCaller, own programs) and sampleID2FamilyCountF 
		#		(for ReplicateVCFGenotypeColumns job)
		# find a way to cache this job (used for same set of samples, but different chromosome intervals)
		pedFile = File(os.path.join(self.auxDirJob.output, 'pedigree_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'sampleID2FamilyCount_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		if self.pedigreeFileFormat==3:	#for polymutt
			polymuttDatFile = File(os.path.join(self.auxDirJob.output, 'datFile_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		else:
			polymuttDatFile = None
		outputPedigreeJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=self.firstVCFJobData.file, outputFile=pedFile, sampleID2FamilyCountF=sampleID2FamilyCountF,\
				polymuttDatFile = polymuttDatFile,\
				outputFileFormat=self.pedigreeFileFormat, replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=self.firstVCFJobData.jobLs + [self.auxDirJob], \
				extraDependentInputLs=[], transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		self.outputPedigreeJob = outputPedigreeJob
		
		"""
		#ExtractSamplesFromVCF non-replicated samples with coverage >=min_coverage
		# this is for outputPedigreeOfHghCoverageSamplesJob. the output of sample IDs do not have replicate tags.
		outputFile = File(os.path.join(self.auxDirJob.output, '%s_minCoverage%s.tsv'%(self.pedigreeFileFormat, self.minCoverageForRefPanel)))
		extractHighCoverageSampleIDJob = self.addExtractSampleIDJob(workflow=workflow, \
							inputFile=self.firstVCFJobData.file, \
							outputFile=outputFile,\
							min_coverage=self.minCoverageForRefPanel, returnData=returnData,\
							transferOutput=True, \
							parentJobLs=self.firstVCFJobData.jobLs + [self.auxDirJob])
		
		
		# GATK SelectVariants: select High-coverage individuals out into a new VCF
		#	selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#	samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		#		?can it deal with Platypus output, which does not have AC/AF/DP?
		# selectHighCoverageSampleJob is needed here because a VCF file of high-coverage members is need 
		# 	for outputPedigreeOfHghCoverageSamplesJob
		#
		inputFBasenamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(os.path.basename(self.firstVCFJobData.file.name))
		highCoverageSampleVCF = File(self.auxDirJob.output, '%s_highCoverageSample.vcf'%(inputFBasenamePrefix))
		selectHighCoverageSampleJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				inputF=self.firstVCFJobData.file, \
				outputF=highCoverageSampleVCF, \
				refFastaFList=self.registerReferenceData.refFastaFList, sampleIDKeepFile=extractHighCoverageSampleIDJob.output,\
				parentJobLs=[self.auxDirJob, extractHighCoverageSampleIDJob]+self.firstVCFJobData.jobLs, \
				extraDependentInputLs=[], transferOutput=transferOutput, \
				extraArguments=None, job_max_memory=2000)
		
		# output a plink pedigree that contains these HC members only
		# output pedigree to get pedigree file (for GATK, TrioCaller, own programs) and sampleID2FamilyCountF (for ReplicateVCFGenotypeColumns job)
		# find a way to cache this job (used for same set of samples, but different chromosome intervals)
		pedFile = File(os.path.join(self.auxDirJob.output, 'highCoveragePedigree_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		sampleID2FamilyCountF = File(os.path.join(self.auxDirJob.output, 'highCoveragePedigree_sampleID2FamilyCount_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		if self.pedigreeFileFormat==3:	#for polymutt
			polymuttDatFile = File(os.path.join(self.auxDirJob.output, 'highCoverage_datFile_outputFileFormat%s.txt'%(self.pedigreeFileFormat)))
		else:
			polymuttDatFile = None
		self.outputPedigreeOfHghCoverageSamplesJob = self.addOutputVRCPedigreeInTFAMGivenOrderFromFileJob(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, \
				inputFile=selectHighCoverageSampleJob.output, outputFile=pedFile, sampleID2FamilyCountF=sampleID2FamilyCountF,\
				polymuttDatFile = None,\
				outputFileFormat=self.pedigreeFileFormat, replicateIndividualTag=self.replicateIndividualTag,\
				treatEveryOneIndependent=self.treatEveryOneIndependent,\
				parentJobLs=[self.auxDirJob, selectHighCoverageSampleJob], \
				extraDependentInputLs=[], transferOutput=True, \
				extraArguments=None, job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		"""
		
		#analyze the pedigree graph to figure out singletons, trios, duos
		self.alignmentLs = self.db.getAlignmentsFromVCFFile(inputFname=yh_pegasus.getAbsPathOutOfFile(self.firstVCFJobData.file))
		"""
		#2013.06.14 approach below does not work because pedigree of extracting-high-coverage + replication is different from that of replication + extracting-high-coverage (=reality).
		# some replicates might end up as singletons in the latter, while not so in the former.
		#
		self.highCoverageAlignmentLs = self.db.filterAlignments(alignmentLs=self.alignmentLs, min_coverage=self.minCoverageForRefPanel, \
			max_coverage=None, individual_site_id=None, \
			sequence_filtered=None, individual_site_id_set=None, \
			mask_genotype_method_id=None, parent_individual_alignment_id=None,\
			country_id_set=None, tax_id_set=None, excludeContaminant=False, excludeTissueIDSet=None,\
			local_realigned=None, reduce_reads=None, report=False)
		
		self.highCoveragePedigreeGraph = self.db.constructPedgreeGraphOutOfAlignments(self.highCoverageAlignmentLs).DG
		self.highCoverageMemberPedigreeSplitStructure = self.db.getPedigreeSplitStructure(pedigreeGraph=self.highCoveragePedigreeGraph, \
															removeFamilyFromGraph=False)
		"""
		pedigreeDataStructure = self.db.constructPedgreeGraphOutOfAlignments(alignmentLs=self.alignmentLs)
		self.pedigreeGraph = pedigreeDataStructure.DG
		individual_id2alignmentLs = pedigreeDataStructure.individual_id2alignmentLs
		
		# work on the pedigree graph to figure out if singleton, trio, duo file will exist.
		#figure out the singletons, duos, trios, using the function in AlignmentToTrioCall
		self.allMemberPedigreeSplitStructure = self.db.getPedigreeSplitStructure(pedigreeGraph=self.pedigreeGraph, \
															removeFamilyFromGraph=False)
		sys.stderr.write("Constructing pedigree split structure for hihg-coverage members ...")
		self.highCoverageMemberPedigreeSplitStructure = PassingData(familySize2familyLs={})
		no_of_families = 0
		for familySize, familyLs in self.allMemberPedigreeSplitStructure.familySize2familyLs.iteritems():
			newFamilyLs = []
			for family in familyLs:
				newFamily = []
				for member in family:
					alignment = individual_id2alignmentLs.get(member)[0]
					if alignment.median_depth>=self.minCoverageForRefPanel:
						newFamily.append(member)
				if len(newFamily)>0:
					newFamilyLs.append(newFamily)
					no_of_families += 1
			self.highCoverageMemberPedigreeSplitStructure.familySize2familyLs[familySize] = newFamilyLs
		sys.stderr.write(" %s families with %s different sizes.\n"%(no_of_families, \
												len(self.highCoverageMemberPedigreeSplitStructure.familySize2familyLs)))
		#a stat merge job (keeping track of how many mendel error sites were filtered)
		filterByRemoveMendelErrorSiteStatMergeFile = File(os.path.join(self.statDirJob.folder, 'filterByRemoveMendelErrorSiteStatMerge.tsv'))
		self.filterByRemoveMendelErrorSiteStatMergeJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
								outputF=filterByRemoveMendelErrorSiteStatMergeFile, \
								transferOutput=True, parentJobLs=[self.statDirJob],\
								extraArguments="-k 1 -v 2-4")	#column 1 is the chromosome length, which are set to be all same.
								#column 2-4 are #sitesInInput1, #sitesInInput2, #overlapping
		
		#a job that outputs alignment coverage (alignment.id, median_depth)
		inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(self.firstVCFJobData.file.name)
		alignmentDepthFile = File(os.path.join(self.auxDirJob.folder, '%s_alignmentDepth.tsv'%(inputFBasenamePrefix)))
		self.outputAlignmentDepthJob = self.addOutputVCFAlignmentDepthRangeJob(executable=self.OutputVCFAlignmentDepthRange, \
						inputFile=self.firstVCFJobData.file, \
						ref_ind_seq_id=self.ref_ind_seq_id, depthFoldChange=None, minGQ=None,\
						outputFile=alignmentDepthFile, outputFileFormat=1,\
						extraArgumentList=None,\
						parentJobLs=[self.auxDirJob]+self.firstVCFJobData.jobLs, \
						extraDependentInputLs=None, transferOutput=True, \
						job_max_memory=2000, sshDBTunnel=self.needSSHDBTunnel)
		
		#concordance stat reduce jobs
		#reduce the replicate concordance results from before TrioCaller (after beagle phasing)
		#
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.allSites.tsv'))
		reduceBeaglePhaseReplicateConcordanceJob_AllSites = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 2,3', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.homo.tsv'))
		reduceBeaglePhaseReplicateConcordanceJob_HomoOnly = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 5,6', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'beaglePhaseReplicateConcordance.tsv'))
		concatenateTwoBeaglePhaseConcordanceResultJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByMergeColumnsWithSameKey, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 2,3,4', transferOutput=False)
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoBeaglePhaseConcordanceResultJob, \
							parentJobLs=[reduceBeaglePhaseReplicateConcordanceJob_AllSites])
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoBeaglePhaseConcordanceResultJob, \
							parentJobLs=[reduceBeaglePhaseReplicateConcordanceJob_HomoOnly])
		returnData.jobDataLs.append(PassingData(jobLs=[concatenateTwoBeaglePhaseConcordanceResultJob], \
											fileLs=[concatenateTwoBeaglePhaseConcordanceResultJob.output]))
		#pass to self, as they will be used in reduceEachVCF()
		self.reduceBeaglePhaseReplicateConcordanceJob_AllSites = reduceBeaglePhaseReplicateConcordanceJob_AllSites
		self.reduceBeaglePhaseReplicateConcordanceJob_HomoOnly = reduceBeaglePhaseReplicateConcordanceJob_HomoOnly
		
		#reduce replicate concordance results from after-TrioCaller VCFs 
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.allSites.tsv'))
		reduceTrioCallerReplicateConcordanceJob_AllSites = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 2,3', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.homo.tsv'))
		reduceTrioCallerReplicateConcordanceJob_HomoOnly = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixBySumSameKeyColsAndThenDivide, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 5,6', transferOutput=False)
		outputFile = File(os.path.join(self.statDirJob.folder, 'trioCallerReplicateConcordance.tsv'))
		concatenateTwoTrioCallerConcordanceResultJob = self.addStatMergeJob(statMergeProgram=self.ReduceMatrixByMergeColumnsWithSameKey, \
							outputF=outputFile, \
							extraArguments='-k 0,1 -v 2,3,4', transferOutput=False)
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoTrioCallerConcordanceResultJob, \
							parentJobLs=[reduceTrioCallerReplicateConcordanceJob_AllSites])
		self.addInputToStatMergeJob(statMergeJob=concatenateTwoTrioCallerConcordanceResultJob, \
							parentJobLs=[reduceTrioCallerReplicateConcordanceJob_HomoOnly])
		returnData.jobDataLs.append(PassingData(jobLs=[concatenateTwoTrioCallerConcordanceResultJob], \
											fileLs=[concatenateTwoTrioCallerConcordanceResultJob.output]))
		#pass to self, as they will be used in reduceEachVCF()
		self.reduceTrioCallerReplicateConcordanceJob_AllSites = reduceTrioCallerReplicateConcordanceJob_AllSites
		self.reduceTrioCallerReplicateConcordanceJob_HomoOnly = reduceTrioCallerReplicateConcordanceJob_HomoOnly
		
		return returnData
	
	def addBeagle3Job(self, workflow=None, executable=None, BeagleJar=None, \
					phasedBeagleInputFile=None,\
					likelihoodBeagleInputFile=None, triosBeagleInputFile=None, pairsBeagleInputFile=None,\
					unphasedBeagleInputFile=None,\
					markersBeagleInputFile=None,
					outputFnamePrefix=None, noOfIterations=None, noOfSamplingHaplotypesPerSample=None, \
					parentJobLs=None, transferOutput=True, job_max_memory=2000,\
					frontArgumentList=None, extraArguments=None, extraArgumentList=None, extraOutputLs=None, \
					extraDependentInputLs=None, no_of_cpus=None, walltime=120, **keywords):
		"""
		i.e.
			
		2013.06.13 a generic function to add Beagle version 3 jobs
		
		crocea@vervetNFS:~/bin/Beagle/example/imputation$ 
		prefix=method_36_Contig791_replicated_phasedByGATK_noMendelError_unphased_minProb0.65;
		java -Xmx15g -XX:PermSize=1000m -XX:MaxPermSize=6000m
			-jar ~/bin/Beagle/beagle.jar like=$prefix\_familySize1.bgl
			pairs=$prefix\_familySize2.bgl trios=$prefix\_familySize3.bgl
			markers=$prefix.markers out=$prefix\_beagled missing=?
			unphased=...bgl
			
		Note:
			"missing=? " can not be used if all input is likelihood data.
		
		Beagle help http://faculty.washington.edu/browning/beagle/beagle.html
		
		niterations=<number of iterations> where <number of iterations> is a positive even 
			integer giving the number of iterations of the phasing algorithm. If an odd integer is 
			specified, the next even integer is used. The niterations argument is optional. The default 
			value is niterations=10. The default value typically gives good accuracy
		
		nsamples=<number of samples> where <number of samples> is positive integer giving 
			the number of haplotype pairs to sample for each individual during each iteration of the 
			phasing algorithm. The nsamples argument is optional. The default value is nsamples=4.
			If you are phasing an extremely large sample (say > 4000 individuals), you may want to 
			use a smaller nsamples parameter (e.g. 1 or 2) to reduce computation time. If you are 
			phasing a small sample (say < 200 individuals), you may want to use a larger nsamples
			parameter (say 10 or 20) to increase accuracy.
		"""
		if workflow is None:
			workflow = self
		if executable is None:
			executable = self.java
		if BeagleJar is None:
			BeagleJar = self.BeagleJar
		if frontArgumentList is None:
			frontArgumentList = []
		if extraArgumentList is None:
			extraArgumentList = []
		if extraDependentInputLs is None:
			extraDependentInputLs = []
		if extraOutputLs is None:
			extraOutputLs = []
		
		#place holder for key output files, corresponding to these input like=, unphased=, trios=, pairs=
		key2ObjectForJob = {"likelihoodPhasedOutputFile":None,\
						"singletonPhasedOutputFile":None,\
						"pairsPhasedOutputFile":None,\
						"triosPhasedOutputFile":None,\
						}
		
		memRequirementObject = self.getJVMMemRequirment(job_max_memory=job_max_memory, minMemory=2000)
		job_max_memory = memRequirementObject.memRequirement
		javaMemRequirement = memRequirementObject.memRequirementInStr
		
		frontArgumentList.extend([javaMemRequirement, '-jar', BeagleJar, "out=%s"%(outputFnamePrefix)])
		extraDependentInputLs.append(BeagleJar)
		
		if noOfIterations is not None:
			frontArgumentList.append("niterations=%s"%(noOfIterations))
		if noOfSamplingHaplotypesPerSample is not None:
			frontArgumentList.append("nsamples=%s"%(noOfSamplingHaplotypesPerSample))
		
		logFile=File('%s.log'%(outputFnamePrefix))
		extraOutputLs.append(logFile)
		key2ObjectForJob['logFile'] = logFile
		
		phasedBeagleFileSuffix = '.bgl.phased.gz'
		onlyLikelihoodInput = True	#determines whether missing=? should be added
		if likelihoodBeagleInputFile:
			frontArgumentList.extend(["like=%s"%(likelihoodBeagleInputFile.name)])
			extraDependentInputLs.append(likelihoodBeagleInputFile)
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(likelihoodBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['likelihoodPhasedOutputFile'] = phasedFile
			#additional output for the likelihood input
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if triosBeagleInputFile:
			frontArgumentList.extend(["trios=%s"%(triosBeagleInputFile.name)])
			extraDependentInputLs.append(triosBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(triosBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['triosPhasedOutputFile'] = phasedFile
		
		if pairsBeagleInputFile:
			frontArgumentList.extend(["pairs=%s"%(pairsBeagleInputFile.name)])
			extraDependentInputLs.append(pairsBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(pairsBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['pairsPhasedOutputFile'] = phasedFile
		
		if unphasedBeagleInputFile:
			frontArgumentList.extend(["unphased=%s"%(unphasedBeagleInputFile.name)])
			extraDependentInputLs.append(unphasedBeagleInputFile)
			onlyLikelihoodInput = False
			inputFBasenamePrefix = utils.getFileBasenamePrefixFromPath(unphasedBeagleInputFile.name)
			
			phasedFile = File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, phasedBeagleFileSuffix))
			extraOutputLs.append(phasedFile)
			key2ObjectForJob['singletonPhasedOutputFile'] = phasedFile
			#additional likelihood-related output
			for suffix in ['.bgl.dose.gz', '.bgl.gprobs.gz', '.bgl.r2']:
				_outputFile=File('%s.%s%s'%(outputFnamePrefix, inputFBasenamePrefix, suffix))
				extraOutputLs.append(_outputFile)
		
		if phasedBeagleInputFile:	#reference panel
			frontArgumentList.extend(["phased=%s"%(phasedBeagleInputFile.name)])
			extraDependentInputLs.append(phasedBeagleInputFile)
			
		
		if markersBeagleInputFile:
			frontArgumentList.extend(["markers=%s"%(markersBeagleInputFile.name)])
			extraDependentInputLs.append(markersBeagleInputFile)
		if not onlyLikelihoodInput:
			frontArgumentList.append("missing=?")
		
		job = self.addGenericJob(workflow=workflow, executable=executable, inputFile=None, \
					inputArgumentOption=None,  inputFileList=None,\
					argumentForEachFileInInputFileList=None,\
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
					extraOutputLs=extraOutputLs, \
					transferOutput=transferOutput, \
					frontArgumentList=frontArgumentList, extraArguments=extraArguments, \
					extraArgumentList=extraArgumentList, job_max_memory=job_max_memory,  sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, no_of_cpus=no_of_cpus, walltime=walltime, **keywords)
		
		return job
	
	def convertVCF2Beagle(self, workflow=None, VCFJobData=None, outputDirJob=None, \
						outputFileBasenamePrefix=None,\
						outputPedigreeJob=None, pedigreeSplitStructure=None,\
						transferOutput=False, \
						job_max_memory=None, walltime=None):
		"""
		2013.05.01
		"""
		
		#replicate the individuals involved in more than one trio/duo
		#2012.4.2 replicate individuals who appear in more than 1 families
		#note: remove individuals in pedigree file but not in VCF file
		
		#GATK ProduceBeagleInput
		beagleLikelihoodFile = File(os.path.join(outputDirJob.folder, "%s.bgl"%(outputFileBasenamePrefix)))
		produceBeagleInputJob = self.addGATKJob(executable=self.ProduceBeagleInputJava, \
					GATKAnalysisType="ProduceBeagleInput", \
					inputFile=VCFJobData.file, inputArgumentOption="--variant:VCF", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					outputFile=beagleLikelihoodFile, \
					parentJobLs=VCFJobData.jobLs + [outputDirJob], \
					transferOutput=transferOutput, \
					extraArguments=None, extraArgumentList=None, \
					extraOutputLs=None, extraDependentInputLs=None, \
					job_max_memory=job_max_memory, \
					no_of_cpus=None, walltime=walltime)
		
		#a SplitPedigreeVCFIntoBeagleTriosDuosFiles job
		#need pedigreeFile (from outputPedigreeJob)
		outputFnamePrefix = os.path.join(outputDirJob.folder, outputFileBasenamePrefix)
		
		key2File= {'size1File': None, 'size2File': None, 'size3File':None}	#to store output files
			#size1File: output for singletons (likelihood format)
			#size2File: output for duos (Beagle genotype format)
			#size3File: output for trios (Beagle genotype format)
		markerFile = File("%s.markers"%(outputFnamePrefix))
		extraOutputLs = []	#first is the markers file
		extraOutputLs.append(markerFile)
		key2File['markerFile'] = markerFile
		for familySize, familyLs in pedigreeSplitStructure.familySize2familyLs.iteritems():
			if familyLs:	#non-empty
				outputFile = File('%s_familySize%s.bgl'%(outputFnamePrefix, familySize))
				key2File['size%sFile'%(familySize)] = outputFile
				extraOutputLs.append(outputFile)
		
		splitPedigreeVCFIntoBeagleTriosDuosFilesJob = self.addGenericJob(executable=self.SplitPedigreeVCFIntoBeagleTriosDuosFiles, \
						inputFile=VCFJobData.file, \
						outputFile=None, \
						parentJobLs=VCFJobData.jobLs + [produceBeagleInputJob, outputDirJob, outputPedigreeJob], \
						extraDependentInputLs=[produceBeagleInputJob.output, outputPedigreeJob.output], \
						extraOutputLs=extraOutputLs,  transferOutput=transferOutput, \
						extraArguments=None, \
						extraArgumentList=["--gatkPrintBeagleFname", produceBeagleInputJob.output, \
										"--plinkPedigreeFname", outputPedigreeJob.output, \
										"--minProbForValidCall %s"%(self.minProbForValidCall), \
										"--dummyIndividualNamePrefix dummy", \
										"--outputFnamePrefix", outputFnamePrefix], \
						job_max_memory=job_max_memory, \
						sshDBTunnel=None, key2ObjectForJob=key2File, objectWithDBArguments=None, \
						no_of_cpus=None, walltime=walltime)
		
		#
		return splitPedigreeVCFIntoBeagleTriosDuosFilesJob
	
	def mapEachInterval(self, workflow=None, intervalData=None, chromosome=None, \
					VCFJobData=None, passingData=None, 
					mapEachChromosomeData=None, transferOutput=False, \
					**keywords):
		"""
		2013.04.30
		"""
		if workflow is None:
			workflow = self
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		
		topOutputDirJob = passingData.topOutputDirJob
		
		intervalFileBasenamePrefix = passingData.intervalFileBasenamePrefix
		span = passingData.span
		noOfIndividuals= passingData.noOfIndividuals
		SNPVCFFile = VCFJobData.file
		SNPVCFJobLs = VCFJobData.jobLs
		"""
		### 2013.06.19 intervalData does not exsit for input that is entirely composed of VCF files (SplitVCFFile job does not return intervals)
		if intervalData.file:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.file
		else:
			mpileupInterval = intervalData.interval
			bcftoolsInterval = intervalData.interval
		intervalFileBasenameSignature = intervalData.intervalFileBasenameSignature
		overlapInterval = intervalData.overlapInterval
		overlapFileBasenameSignature = intervalData.overlapIntervalFnameSignature
		span = intervalData.span
		"""
		if chromosome is None:
			chromosome = getattr(passingData, 'chromosome', None)
		
		#noOfIndividuals
		realInputVolume = noOfIndividuals * span
		baseInputVolume = 200*2000000
		#base is 200 individual X 2Mb region => 120 minutes
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value
		#base is 4X, => 5000M
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=10000).value
		
		splitVCFJob = passingData.mapEachVCFData.splitVCFJob
		
		##### stage 0, replicate individuals, simple phase by GATK, remove mendel error sites
		# 2013.06.11 replicate individuals who appear in more than 1 families
		round1_IndividualsReplicatedVCF = File( os.path.join(self.phaseByGATKDirJob.folder, \
											'%s.replicate.vcf'%(intervalFileBasenamePrefix)))
		replicateVCFGenotypeColumnsJob = self.addReplicateVCFGenotypeColumnsJob(workflow, \
					executable=workflow.ReplicateVCFGenotypeColumns, inputF=VCFJobData.file, \
					sampleID2FamilyCountF=self.outputPedigreeJob.sampleID2FamilyCountF, \
					outputF=round1_IndividualsReplicatedVCF, \
					replicateIndividualTag=self.replicateIndividualTag,\
					parentJobLs=[self.outputPedigreeJob, self.phaseByGATKDirJob]+VCFJobData.jobLs, \
					extraDependentInputLs=[], \
					transferOutput=False, \
					extraArguments=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=10000).value, \
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value,\
					)
		
		if self.extractRefPanelSampleIDJob is None:
			#ExtractSamplesFromVCF samples with coverage >=min_coverage
			# the input VCF contains replicates.
			outputFile = File(os.path.join(self.auxDirJob.output, 'pedigree.format%s.minCoverage%s.tsv'%(self.pedigreeFileFormat, self.minCoverageForRefPanel)))
			extractRefPanelSampleIDJob = self.addExtractSampleIDJob(workflow=workflow, inputFile=replicateVCFGenotypeColumnsJob.output, \
								outputFile=outputFile,\
								min_coverage=self.minCoverageForRefPanel, returnData=returnData,\
								transferOutput=True, \
								parentJobLs=[replicateVCFGenotypeColumnsJob, self.auxDirJob])
			self.extractRefPanelSampleIDJob = extractRefPanelSampleIDJob
		
		
		#  phase by GATK, need a pedigree file (plink format)
		"""
		java -Xmx4g -jar ~/script/gatk2/GenomeAnalysisTK.jar
			-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
			-T PhaseByTransmission -V method_36_Contig791_replicated.vcf
			-ped trioCaller_VRC.merlin -o method_36_Contig791_replicated_phasedByGATK.vcf
			--MendelianViolationsFile method36_Contig791_replicated_MendelViolation.txt
		"""
		phasedByGATKVCFFile = File(os.path.join(self.phaseByGATKDirJob.folder, '%s.phasedByGATK.vcf'%(intervalFileBasenamePrefix)))
		mendelianViolationsTextFile = File(os.path.join(self.phaseByGATKDirJob.folder, '%s.mendelViolation.txt'%(intervalFileBasenamePrefix))) 
		phaseByGATKJob = self.addGATKJob(executable=self.PhaseByGATK, GATKAnalysisType='PhaseByTransmission', \
					inputFile=replicateVCFGenotypeColumnsJob.output, \
					inputArgumentOption="--variant:VCF", \
					refFastaFList=self.registerReferenceData.refFastaFList,\
					outputFile=phasedByGATKVCFFile, \
					parentJobLs=[self.phaseByGATKDirJob, self.outputPedigreeJob, replicateVCFGenotypeColumnsJob], \
					transferOutput=False, \
					extraArguments=None, \
					extraArgumentList=["-ped", self.outputPedigreeJob.output, "--MendelianViolationsFile", mendelianViolationsTextFile], \
					extraOutputLs=[mendelianViolationsTextFile], \
					extraDependentInputLs=[self.outputPedigreeJob.output], no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		#find mendel error sites:
		"""
		# input is phased by GATK
		java -Xmx4g -jar ~/script/gatk2/GenomeAnalysisTK.jar
			-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
			-T SelectVariants --variant:VCF method_36_Contig791_replicated_phasedByGATK.vcf
			-ped trioCaller_VRC.merlin -mvq 0 --mendelianViolation
			-o method_36_Contig791_replicated_phasedByGATK_mendelViolations.vcf
		"""
		mendelViolationVCFFile = File(os.path.join(self.phaseByGATKDirJob.folder, '%s.mendelViolation.vcf'%(intervalFileBasenamePrefix)))
		findMendelErrorJob = self.addGATKJob(executable=self.FindMendelErrorByGATK, GATKAnalysisType="SelectVariants",\
				GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, inputFile=phaseByGATKJob.output, \
				inputArgumentOption="--variant:VCF",\
				outputFile=mendelViolationVCFFile, \
				refFastaFList=self.registerReferenceData.refFastaFList, \
				parentJobLs=[self.phaseByGATKDirJob, phaseByGATKJob, self.outputPedigreeJob], \
				extraDependentInputLs=[self.outputPedigreeJob.output], transferOutput=False, \
				extraArguments=None, \
				extraArgumentList=["-mvq 0", "--mendelianViolation", "-ped", self.outputPedigreeJob.output], \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		# remove mendel error sites
		"""
		java -Xmx4g -jar ~/script/gatk2/GenomeAnalysisTK.jar
			-R /Network/Data/vervet/db/individual_sequence/524_superContigsMinSize2000.fasta
			-T SelectVariants --variant method_36_Contig791_replicated_phasedByGATK.vcf
			--discordance method_36_Contig791_replicated_phasedByGATK_mendelViolations.vcf
			-o method_36_Contig791_replicated_phasedByGATK_noMendelError.vcf
		"""
		noMendelErrorVCFFile = File(os.path.join(self.phaseByGATKDirJob.folder, '%s_noMendelError.vcf'%(intervalFileBasenamePrefix)))
		removeMendelErrorSiteByGATKJob = self.addGATKJob(executable=self.RemoveMendelErrorSiteByGATK, \
					GATKAnalysisType="SelectVariants", \
					inputFile=phaseByGATKJob.output, inputArgumentOption="--variant:VCF", \
					refFastaFList = self.registerReferenceData.refFastaFList, \
					outputFile=noMendelErrorVCFFile, \
					extraArguments=None, extraArgumentList=["--discordance", findMendelErrorJob.output], \
					extraDependentInputLs=[findMendelErrorJob.output],
					parentJobLs=[phaseByGATKJob, findMendelErrorJob, self.phaseByGATKDirJob], \
					transferOutput=False, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value,\
					no_of_cpus=None)
		
		#count how many sites before and after removal
		outputF = File(os.path.join(self.statDirJob.output, '%s_BeforeAfterRemovingMendelErrorSiteStat.tsv'%(intervalFileBasenamePrefix)))
		self.addVCFBeforeAfterFilterStatJob(chromosome=chromosome, outputF=outputF, \
							lastVCFJob=phaseByGATKJob, \
							currentVCFJob=removeMendelErrorSiteByGATKJob,\
							statMergeJob=self.filterByRemoveMendelErrorSiteStatMergeJob, \
							parentJobLs=[phaseByGATKJob, removeMendelErrorSiteByGATKJob, self.statDirJob])
		
		#### Part 1 generate high-quality reference panel through Beagle on high-coverage individuals
		# extractRefPanelSampleIDJob outputs sample IDs with replicate tags
		# select the high-coverage members
		outputVCF = File(os.path.join(self.highCoveragePanelDirJob.output, \
									'%s.minCoverageForRefPanel%s.vcf'%(intervalFileBasenamePrefix, self.minCoverageForRefPanel)))
		#selectVariants would re-generate AC, AF so that TrioCaller could read it.
		#samtools uses 'AC1' instead of AC, 'AF1' instead of AF.
		selectHighCoverageSampleJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
				GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, inputF=removeMendelErrorSiteByGATKJob.output, outputF=outputVCF, \
				refFastaFList=self.registerReferenceData.refFastaFList, sampleIDKeepFile=self.extractRefPanelSampleIDJob.output,\
				parentJobLs=[self.highCoveragePanelDirJob, splitVCFJob, self.extractRefPanelSampleIDJob, removeMendelErrorSiteByGATKJob], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		#convert VCF to Beagle input
		outputFileBasenamePrefix = "%s.unphased"%(intervalFileBasenamePrefix)
		highCoverageVCF2BeagleInputJob = self.convertVCF2Beagle(VCFJobData=PassingData(file=selectHighCoverageSampleJob.output, \
																			jobLs=[selectHighCoverageSampleJob]), \
									outputDirJob=self.highCoveragePanelDirJob, \
									outputFileBasenamePrefix=outputFileBasenamePrefix,
									outputPedigreeJob=self.outputPedigreeJob, \
									pedigreeSplitStructure=self.highCoverageMemberPedigreeSplitStructure,\
									transferOutput=False,\
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		# run Beagle 
		outputFnamePrefix = os.path.join(self.highCoveragePanelDirJob.folder, "%s.highCoverageBeagled"%(intervalFileBasenamePrefix))
		beagleOnHighCoverageJob = self.addBeagle3Job(executable=self.BeagleJava, \
								likelihoodBeagleInputFile=highCoverageVCF2BeagleInputJob.size1File, \
								triosBeagleInputFile=highCoverageVCF2BeagleInputJob.size3File, \
								pairsBeagleInputFile=highCoverageVCF2BeagleInputJob.size2File, \
								markersBeagleInputFile=highCoverageVCF2BeagleInputJob.markerFile, \
								outputFnamePrefix=outputFnamePrefix, noOfIterations=30,\
								noOfSamplingHaplotypesPerSample=4,\
								parentJobLs=[self.highCoveragePanelDirJob, highCoverageVCF2BeagleInputJob], \
								transferOutput=False, \
								extraArguments=None, extraArgumentList=["lowmem=true"],\
								extraOutputLs=None, extraDependentInputLs=None, \
								no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=9000, maxJobPropertyValue=18000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		#a SelectDistantMembersFromGenotypeFile.py job to generate a ref panel for 2nd-round beagle
		# need the pedigree file
		phasedRefPanelGenotypeFile = File(os.path.join(self.highCoveragePanelDirJob.folder, '%s_PhasedRefPanel.bgl'%(intervalFileBasenamePrefix)))
		selectDistantMembersFromGenotypeFileJob = self.addGenericJob(executable=self.SelectDistantMembersFromGenotypeFile, \
						outputFile=phasedRefPanelGenotypeFile, outputArgumentOption="-o", \
						inputFileList=[beagleOnHighCoverageJob.likelihoodPhasedOutputFile, beagleOnHighCoverageJob.singletonPhasedOutputFile, \
									beagleOnHighCoverageJob.pairsPhasedOutputFile, beagleOnHighCoverageJob.triosPhasedOutputFile], \
						argumentForEachFileInInputFileList="", \
						parentJobLs=[self.highCoveragePanelDirJob, beagleOnHighCoverageJob, self.outputPedigreeJob], \
						extraDependentInputLs=[self.plinkIBDCheckOutputFile], \
						extraOutputLs=None, transferOutput=False, frontArgumentList=None, \
						extraArguments=None, \
						extraArgumentList=["--maxPairwiseKinship 0.2", "--sampleSize 40", \
							"--plinkIBDCheckOutputFname", self.plinkIBDCheckOutputFile, \
							"--replicateIndividualTag", self.replicateIndividualTag,\
							"--individualAlignmentCoverageFname", self.outputAlignmentDepthJob.output, \
							"--pedigreeFname", self.outputPedigreeJob.output], \
						no_of_cpus=None, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		##### Part 2 run Beagle on everyone with reference panel
		
		# convert VCF to Beagle input
		AllMemberVCFJobData = PassingData(file=removeMendelErrorSiteByGATKJob.output, jobLs=[removeMendelErrorSiteByGATKJob])
		outputFileBasenamePrefix = "%s.everyone.unphased"%(intervalFileBasenamePrefix)
		vcf2BeagleInputJob = self.convertVCF2Beagle(VCFJobData=AllMemberVCFJobData, \
							outputDirJob=self.mapDirJob, outputFileBasenamePrefix=outputFileBasenamePrefix, \
							outputPedigreeJob=self.outputPedigreeJob, \
							pedigreeSplitStructure=self.allMemberPedigreeSplitStructure, \
							transferOutput=False, \
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		# run Beagle
		outputFnamePrefix = os.path.join(self.mapDirJob.folder, '%s.beagled'%(intervalFileBasenamePrefix))
		beagleJob = self.addBeagle3Job(executable=self.BeagleJava, \
						phasedBeagleInputFile=selectDistantMembersFromGenotypeFileJob.output,\
						likelihoodBeagleInputFile=vcf2BeagleInputJob.size1File, \
						triosBeagleInputFile=vcf2BeagleInputJob.size3File, \
						pairsBeagleInputFile=vcf2BeagleInputJob.size2File, \
						markersBeagleInputFile=vcf2BeagleInputJob.markerFile, \
						outputFnamePrefix=outputFnamePrefix, noOfIterations=30,\
						noOfSamplingHaplotypesPerSample=4,\
						extraArguments=None, extraArgumentList=["lowmem=true"],\
						parentJobLs=[self.mapDirJob, selectDistantMembersFromGenotypeFileJob, vcf2BeagleInputJob], \
						transferOutput=False, no_of_cpus=None, \
						job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=10000, maxJobPropertyValue=25000).value,\
						walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value,\
						)
		
		#a CombinePhasedBeagleOutputsIntoVCF job
		phasedVCFFile = File(os.path.join(self.mapDirJob.folder, "%s.phased.vcf"%(intervalFileBasenamePrefix)))
		combinePhasedBeagleOutputsIntoVCFJob = self.addGenericJob(executable=self.CombinePhasedBeagleOutputsIntoVCF, \
						inputFile=removeMendelErrorSiteByGATKJob.output, inputArgumentOption="--originalVCFFname", \
						outputFile=phasedVCFFile, outputArgumentOption="-o", \
						inputFileList=[beagleJob.likelihoodPhasedOutputFile, beagleJob.singletonPhasedOutputFile, \
									beagleJob.pairsPhasedOutputFile, beagleJob.triosPhasedOutputFile], \
						argumentForEachFileInInputFileList="", \
						parentJobLs=[beagleJob, removeMendelErrorSiteByGATKJob, self.mapDirJob], \
						extraDependentInputLs=None, \
						extraOutputLs=None, transferOutput=False, \
						extraArgumentList=["--replicateIndividualTag", self.replicateIndividualTag], \
						no_of_cpus=None, \
						job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
						walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		"""
		#### 2013.06.19 non-overlap interval is not clear, so skip it
		# select interval job to remove the overlapping region
		#select non-overlap interval job to remove the overlapping region
		#select the variants to get rid of overlap, so that trioCallerReplicateConcordanceJob could be properly reduced
		nonOverlapBeaglePhaseOutputF = File(os.path.join(self.mapDirJob.folder, '%s.beaglePhased.nonoverlap.vcf'%intervalFileBasenameSignature))
		selectBeaglePhaseOutputJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
						inputF=combinePhasedBeagleOutputsIntoVCFJob.output, \
						outputF=nonOverlapBeaglePhaseOutputF, \
						refFastaFList=self.registerReferenceData.refFastaFList, \
						interval=mpileupInterval,\
						parentJobLs=[combinePhasedBeagleOutputsIntoVCFJob, self.mapDirJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		"""
		# a CheckGenotypeConcordanceAmongReplicates.py job
		beaglePhasedReplicateConcordanceFile = File(os.path.join(self.statDirJob.folder, \
								'%s.phased.concordance.tsv'%(intervalFileBasenamePrefix)))
		#specify GenomeAnalysisTKJar due to gatk1.6 (not gatk2), 
		returnData.beaglePhasedReplicateConcordanceJob = self.addGATKJob(executable=self.CalculateConcordanceJava, \
					GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
					GATKAnalysisType="CalculateConcordanceAmongReplicates",\
					inputFile=combinePhasedBeagleOutputsIntoVCFJob.output, inputArgumentOption="--variant", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					interval=None, \
					outputFile=beaglePhasedReplicateConcordanceFile, outputArgumentOption="--concordanceStatFname",\
					frontArgumentList=None, extraArguments="--replicateIndividualTag %s"%(self.replicateIndividualTag), \
					extraArgumentList=None, extraOutputLs=None, \
					parentJobLs=[self.statDirJob, combinePhasedBeagleOutputsIntoVCFJob], \
					transferOutput=False, \
					no_of_cpus=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		#TrioCaller
		refineGenotypeOutputF = File(os.path.join(self.mapDirJob.folder, \
												'%s.trioCaller.vcf'%(intervalFileBasenamePrefix)))
		refineGenotypeJob = self.addTrioCallerJob(trioCallerWrapper=self.trioCallerWrapper, \
				trioCallerPath=self.trioCallerPath, \
				inputVCF=combinePhasedBeagleOutputsIntoVCFJob.output,\
				pedFile=self.outputPedigreeJob.output, outputVCF=refineGenotypeOutputF, \
				parentJobLs=[self.mapDirJob, combinePhasedBeagleOutputsIntoVCFJob, self.outputPedigreeJob], \
				extraDependentInputLs=[], transferOutput=False, \
				extraArguments="--inputPhased", \
				job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
				walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)	#1.2G memory for 12K loci
		
		returnData.refineGenotypeJob = refineGenotypeJob
		
		"""
		#### 2013.06.19 non-overlap interval is not clear, so skip it
		#select non-overlap interval job to remove the overlapping region
		#select the variants to get rid of overlap, so that trioCallerReplicateConcordanceJob could be properly reduced
		nonOverlapTrioCallerOutputF = File(os.path.join(self.mapDirJob.folder, '%s.trioCaller.nonoverlap.vcf'%intervalFileBasenameSignature))
		selectTrioCallerOutputJob = self.addSelectVariantsJob(SelectVariantsJava=self.SelectVariantsJava, \
						inputF=refineGenotypeJob.output, outputF=nonOverlapTrioCallerOutputF, \
						interval=mpileupInterval,\
						refFastaFList=self.registerReferenceData.refFastaFList, parentJobLs=[refineGenotypeJob, self.mapDirJob], \
						extraDependentInputLs=[], transferOutput=False, \
						extraArguments=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		"""
		# a CheckGenotypeConcordanceAmongReplicates.py job
		trioCallerReplicateConcordanceFile = File(os.path.join(self.statDirJob.folder, \
								'%s.trioCaller.concordance.tsv'%(intervalFileBasenamePrefix)))
		returnData.trioCallerReplicateConcordanceJob = self.addGATKJob(executable=self.CalculateConcordanceJava, \
					GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
					GATKAnalysisType="CalculateConcordanceAmongReplicates",\
					inputFile=refineGenotypeJob.output, inputArgumentOption="--variant", \
					refFastaFList=self.registerReferenceData.refFastaFList, \
					interval=None, \
					outputFile=trioCallerReplicateConcordanceFile, outputArgumentOption="--concordanceStatFname",\
					frontArgumentList=None, extraArguments="--replicateIndividualTag %s"%(self.replicateIndividualTag), \
					extraArgumentList=None, extraOutputLs=None, \
					parentJobLs=[self.statDirJob, refineGenotypeJob], \
					transferOutput=False, \
					no_of_cpus=None, \
					job_max_memory=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value, \
					walltime=self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		#2013.06.14
		#merge replicates to generate consensus call
		# (not haplotype-based, as different recombination points across replicate haplotypes make it non-trivial )
		mergeReplicateOutputF = File(os.path.join(self.mapDirJob.folder, \
									'%s.replicatesMerged.vcf'%(intervalFileBasenamePrefix)))
		returnData.mergeVCFReplicateColumnsJob = self.addMergeVCFReplicateGenotypeColumnsJob(\
							executable=self.MergeVCFReplicateHaplotypesJava,\
							GenomeAnalysisTKJar=self.GenomeAnalysisTKJar, \
							inputF=refineGenotypeJob.output, outputF=mergeReplicateOutputF, \
							replicateIndividualTag=self.replicateIndividualTag, \
							refFastaFList=self.registerReferenceData.refFastaFList, parentJobLs=[self.mapDirJob, refineGenotypeJob], \
							extraDependentInputLs=[], transferOutput=False, \
							extraArguments=None, \
							analysis_type='MergeVCFReplicateGenotypeColumns',\
					job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=4000, \
							minJobPropertyValue=4000, maxJobPropertyValue=9000).value,\
					walltime= self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=1200).value)
		
		
		return returnData
	
	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, mapEachIntervalDataLs=None,\
					transferOutput=True, **keywords):
		"""
		2013.05.01
			#. concatenate all the sub-VCFs into one
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		returnData.mapEachIntervalDataLs = mapEachIntervalDataLs
		
		refineGenotypeJobLs = [pdata.refineGenotypeJob for pdata in mapEachIntervalDataLs]
		mergeVCFReplicateColumnsJobLs = [pdata.mergeVCFReplicateColumnsJob for pdata in mapEachIntervalDataLs]
		
		
		realInputVolume = passingData.jobData.file.noOfIndividuals * passingData.jobData.file.noOfLoci
		baseInputVolume = 200*2000000
		#base is 4X coverage in 20Mb region => 120 minutes
		walltime = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=60, \
							minJobPropertyValue=60, maxJobPropertyValue=500).value
		#base is 4X, => 5000M
		job_max_memory = self.scaleJobWalltimeOrMemoryBasedOnInput(realInputVolume=realInputVolume, \
							baseInputVolume=baseInputVolume, baseJobPropertyValue=2000, \
							minJobPropertyValue=2000, maxJobPropertyValue=8000).value
		self.concatenateOverlapIntervalsIntoVCF(chromosome=chromosome, passingData=passingData, \
						intervalJobLs=refineGenotypeJobLs, outputDirJob=self.replicateVCFDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime, \
						**keywords)
		self.concatenateOverlapIntervalsIntoVCF(chromosome=chromosome, passingData=passingData, \
						intervalJobLs=mergeVCFReplicateColumnsJobLs, outputDirJob=self.reduceOutputDirJob, \
						transferOutput=True, job_max_memory=job_max_memory, walltime=walltime,\
						**keywords)
		
		for pdata in mapEachIntervalDataLs:
			#add this output to the union job
			self.addInputToStatMergeJob(statMergeJob=self.reduceBeaglePhaseReplicateConcordanceJob_AllSites, \
							parentJobLs=[pdata.beaglePhasedReplicateConcordanceJob])
			self.addInputToStatMergeJob(statMergeJob=self.reduceBeaglePhaseReplicateConcordanceJob_HomoOnly, \
							parentJobLs=[pdata.beaglePhasedReplicateConcordanceJob])
			self.addInputToStatMergeJob(statMergeJob=self.reduceTrioCallerReplicateConcordanceJob_AllSites, \
							parentJobLs=[pdata.trioCallerReplicateConcordanceJob])
			self.addInputToStatMergeJob(statMergeJob=self.reduceTrioCallerReplicateConcordanceJob_HomoOnly, \
							parentJobLs=[pdata.trioCallerReplicateConcordanceJob])
		
		return returnData
		
		
	def concatenateOverlapIntervalsIntoVCF(self, chromosome=None, passingData=None, intervalJobLs=None,\
					outputDirJob=None,
					transferOutput=True, job_max_memory=None, walltime=None, **keywords):
		"""
		2013.06.14
			#. concatenate overlapping (share some loci) VCFs into one, used in reduceEachVCF
		"""
		
		#ligate vcf job (different segments of a chromosome into one chromosome) for replicate VCFs
		concatTrioCallerOutputFname = os.path.join(outputDirJob.folder, '%s.vcf'%(passingData.fileBasenamePrefix))
		concatTrioCallerOutputF = File(concatTrioCallerOutputFname)
		concatJob = self.addLigateVcfJob(executable=self.ligateVcf, \
									ligateVcfPerlPath=self.ligateVcfPerlPath, \
									outputFile=concatTrioCallerOutputF, \
									parentJobLs=[outputDirJob], \
									extraDependentInputLs=None, transferOutput=False, \
									extraArguments=None, job_max_memory=job_max_memory, walltime=walltime/2)
		
		for intervalJob in intervalJobLs:
			#add this output to the union job
			# 2012.6.1 done it through addInputToStatMergeJob()
			self.addInputToStatMergeJob(statMergeJob=concatJob, inputF=intervalJob.output, \
							parentJobLs=[intervalJob], extraDependentInputLs=[])
		
		
		#convert to vcf4 so that other vcftools software could be used.
		vcf4Filename = os.path.join(outputDirJob.folder, '%s.v4.vcf'%passingData.fileBasenamePrefix)
		vcf4File = File(vcf4Filename)
		vcf_convert_TrioCallerOutputJob = self.addVCFFormatConvertJob(vcf_convert=self.vcf_convert, \
					parentJob=concatJob, inputF=concatJob.output, \
					outputF=vcf4File, transferOutput=False, job_max_memory=job_max_memory, walltime=walltime/2)
		
		
		#bgzip and tabix the trio caller output
		trioGzipOutputF = File("%s.gz"%vcf4Filename)
		bgzip_tabix_job = self.addBGZIP_tabix_Job(parentJob=vcf_convert_TrioCallerOutputJob, \
								inputF=vcf4File, outputF=trioGzipOutputF, transferOutput=True,\
								job_max_memory=job_max_memory/4, walltime=walltime/4)
		
		return bgzip_tabix_job
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2013.06
		"""
		parentClass.registerCustomExecutables(self, workflow)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='CalculateConcordanceJava', \
															clusterSizeMultipler=0.6)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='ProduceBeagleInputJava', \
															clusterSizeMultipler=0.6)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='PhaseByGATK',\
															clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='FindMendelErrorByGATK', \
															clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=self.javaPath, name='RemoveMendelErrorSiteByGATK', \
															clusterSizeMultipler=0.5)
		
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/extractor/SelectDistantMembersFromGenotypeFile.py"), \
										name='SelectDistantMembersFromGenotypeFile', clusterSizeMultipler=0.2)
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/splitter/SplitPedigreeVCFIntoBeagleTriosDuosFiles.py"), \
										name='SplitPedigreeVCFIntoBeagleTriosDuosFiles', clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(self.pymodulePath, "pegasus/mapper/converter/CombinePhasedBeagleOutputsIntoVCF.py"), \
										name='CombinePhasedBeagleOutputsIntoVCF', clusterSizeMultipler=1)		#2013.05.02 no clustering for OutputVRCPedigreeInTFAMGivenOrderFromFile
		self.setOrChangeExecutableClusterSize(executable=self.OutputVRCPedigreeInTFAMGivenOrderFromFile, clusterSizeMultipler=0)
	
if __name__ == '__main__':
	main_class = BeagleAndTrioCallerOnVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()