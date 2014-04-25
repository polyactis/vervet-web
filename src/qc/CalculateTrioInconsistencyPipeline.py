#!/usr/bin/env python
"""
Examples:
	# 2011-9-29 (site depth >=1, --minDepth 1)
	%s --ref_ind_seq_id 524 -I ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o TrioInconsistency92VRC_20110922T2216.xml -l condorpool -j condorpool  --db_user yh --hostname uclaOffice --minDepth 1
	
	# 2011-9-29 (including depth=0 sites, --minDepth 0)
	%s --ref_ind_seq_id 524 -I ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o TrioInconsistency92VRC_20110922T2216.xml -l condorpool -j condorpool  --db_user yh --hostname uclaOffice --minDepth 0
	
	# 2011.12.16 run on hoffman2 condor (site depth >=1, --minDepth 1)  (turn on checkEmptyVCFByReading, --checkEmptyVCFByReading)
	# only contig ID<=100 (-x)
	%s --ref_ind_seq_id 524 -I AlignmentToTrioCallPipeline_VRC_Aln559_600_Trio620_632_648_top2Contigs.2011.12.14T1432/trioCaller/
		-o dags/MendelInconsistency/TrioInconsistency_TrioCall_VRC_Aln559_600_Trio620_632_648_top2Contigs.xml
		-l hcondor -j hcondor
		--db_user yh --hostname localhost
		-t ~/NetworkData/vervet/db/
		-D ~/NetworkData/vervet/db/  --clusters_size 1 --checkEmptyVCFByReading --minDepth 1 -x 100
	
	
Description:
	2011-9-29
		a program which generates a pegasus workflow dag (xml file) to 
		
		1. distribution of inconsistent rate
		2. inconsistent rate v.s. AAF
		3. inconsistent rate v.s. position of contig
		4. distribution of AAF (todo) 
				
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

from Pegasus.DAX3 import Executable, File, Link, Job
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from pymodule import VCFFile
from vervet.src import VervetDB, AbstractVervetWorkflow

parentClass = AbstractVervetWorkflow
class CalculateTrioInconsistencyPipeline(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update({
						('addTrioContigSpecificPlotJobs', 0, ): [0, '', 0, 'whether to add mendelian-inconsistency plotting jobs for each trio and each contig', ],\
						('addTrioSpecificPlotJobs', 0, ): [0, '', 0, 'whether to add mendelian-inconsistency plotting jobs for each trio', ],\
						('inputDir', 0, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		parentClass.__init__(self, **keywords)
		self.inputDir = os.path.abspath(self.inputDir)
	
	
	def getDuoTrioThatExistInVCF(self, db_vervet, VCFFilename):
		"""
		2011.12.13
			1. get all alignment IDs (their individual ID) from VCFFilename
			2. construct a directed graph from db and the given alignments
			3. find all trios/duos
		"""
		sys.stderr.write("Getting alignments based on VCF file %s ... "%(VCFFilename))
		vcfFile = VCFFile(inputFname=VCFFilename)
		#each sample id is '%s_%s_%s_%s_vs_%s'%(aln.id, ind_seq_id, individual_sequence.individual.code, sequencer, ref_ind_seq_id)
		alignmentLs = []
		for sample_id in vcfFile.sample_id_ls:
			try:
				alignment_id = int(sample_id.split('_')[0])
				alignment = VervetDB.IndividualAlignment.get(alignment_id)
				alignmentLs.append(alignment)
			except:
				sys.stderr.write('Except at sample_id=%s: %s\n'%(sample_id, repr(sys.exc_info())))
				import traceback
				traceback.print_exc()
		sys.stderr.write(" %s alignments.\n"%(len(alignmentLs)))
		return self.getDuoTrioFromAlignmentLs(db_vervet, alignmentLs)
		
	def getDuoTrioFromAlignmentLs(self, db_vervet, alignmentLs):
		"""
		2012.1.25
			append alignment.id, individual.code to individual_sequence.id (the current trio member ID)
		2012.1.9
			split out of getDuoTrioThatExistInVCF()
		"""
		sys.stderr.write("Getting duos and trios from %s  alignments... "%(len(alignmentLs)))
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		DG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		
		#find trios first
		trioLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=3, removeFamilyFromGraph=False)
		#find duos
		duoLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=2, removeFamilyFromGraph=False)
		familyLs = trioLs + duoLs
		familyStrLs = []
		for family in familyLs:
			familySize = len(family)
			trio_id = None
			if familySize==2:	#one parent missing
				parent1ID, offspring_id = family[:2]
				#output the single parent first
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				parent1Sex = parent1Alignment.individual_sequence.individual.codeSexInNumber()
				#output a fake 2nd parent, with opposite sex
				parent2ID = 0
				parent2Sex = 1-(parent1Sex-1)+1	#1 is father. 2 is mother.
				#output the offspring
				childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
				child_id = childAlignment.getCompositeID()
				if parent1Sex==1:
					father_id = parent1Alignment.getCompositeID()
					mother_id = 0
				else:
					father_id = 0
					mother_id = parent1Alignment.getCompositeID()
				trio_id = '%s,%s,%s'%(father_id, mother_id, child_id)
				
			elif familySize==3:
				parent1ID, parent2ID, offspring_id = family[:3]
				parent1Alignment = individual_id2alignmentLs.get(parent1ID)[0]
				parent1Sex = parent1Alignment.individual_sequence.individual.codeSexInNumber()
				parent2Alignment = individual_id2alignmentLs.get(parent2ID)[0]
				childAlignment = individual_id2alignmentLs.get(offspring_id)[0]
				child_id = childAlignment.getCompositeID()
				
				if parent1Sex==1:
					father_id = parent1Alignment.getCompositeID()
					mother_id = parent2Alignment.getCompositeID()
				else:
					father_id = parent2Alignment.getCompositeID()
					mother_id = parent1Alignment.getCompositeID()
				trio_id = '%s,%s,%s'%(father_id, mother_id, child_id)
			
			if trio_id:
				familyStrLs.append(trio_id)
		sys.stderr.write(" %s trios fetched.\n"%(len(familyStrLs)))
		return familyStrLs
	
	def getAllTrios(self, db_vervet, aln_ref_ind_seq_id):
		"""
		2013.09.04 deprecated
		2011-9-28
			find all trios, a list of "father,mother,child", all of which are identified by IndividualSequence.id.
			
			for each individual_alignment, -> ind_seq -> individual
				get its parents from ind2ind
				check if both parents have also been sequenced and aligned.
		"""
		sys.stderr.write("Finding trios, sequenced and aligned against IndividualSequence %s ..."%(aln_ref_ind_seq_id))
		#find out all available alignments.
		alignmentQuery = VervetDB.IndividualAlignment.query.filter_by(ref_ind_seq_id=aln_ref_ind_seq_id)
		individual_with_alignment_id2alignment = {}
		
		for individual_alignment in alignmentQuery:
			ind = individual_alignment.individual_sequence.individual
			individual_with_alignment_id2alignment[ind.id] = individual_alignment
		
		trioLs = []
		for ind_id, individual_alignment in individual_with_alignment_id2alignment.iteritems():
			#find the father first
			ind2ind = VervetDB.Ind2Ind.query.filter_by(individual2_id=ind_id).filter_by(relationship_type_id=1).first()
			if ind2ind:
				fatherIndividualID = ind2ind.individual1_id
			else:
				fatherIndividualID = None
			#find mother
			ind2ind = VervetDB.Ind2Ind.query.filter_by(individual2_id=ind_id).filter_by(relationship_type_id=2).first()
			if ind2ind:
				motherIndividualID = ind2ind.individual1_id
			else:
				motherIndividualID = None
			motherAlignment = individual_with_alignment_id2alignment.get(motherIndividualID)
			fatherAlignment = individual_with_alignment_id2alignment.get(fatherIndividualID)
			if fatherAlignment and motherAlignment:
				trio_id = '%s,%s,%s'%(fatherAlignment.ind_seq_id, motherAlignment.ind_seq_id, individual_alignment.ind_seq_id)
				trioLs.append(trio_id)
		sys.stderr.write(" %s trios fetched.\n"%(len(trioLs)))
		return trioLs
	
	def addTrioInconsistencyCalculationJob(self, workflow=None, executable=None, \
						inputFile=None, outputFnamePrefix=None, \
						trio_id=None, windowSize=200000, minDepth=1,\
						parentJobLs=None, extraArgumentList=None,\
						transferOutput=False, **keywords):
		"""
		2013.09.04 modernized
		2011-9-29
		"""
		if extraArgumentList is None:
			extraArgumentList = []
		key2ObjectForJob = {}
		extraOutputLs = []
		if minDepth is not None:
			extraArgumentList.append("--minDepth %s"%(minDepth))
		if windowSize is not None:
			extraArgumentList.append("--windowSize %s"%(windowSize))
		if trio_id:
			extraArgumentList.append("--trio_id %s"%(trio_id))
		extraArgumentList.append("--outputFnamePrefix %s"%outputFnamePrefix)
		
		
		summaryOutputF = File("%s.summary.tsv"%outputFnamePrefix)
		extraOutputLs.append(summaryOutputF)
		key2ObjectForJob["summaryOutputF"] = summaryOutputF
		
		windowOutputF = File("%s.window.%s.tsv"%(outputFnamePrefix, windowSize))
		extraOutputLs.append(windowOutputF)
		key2ObjectForJob["windowOutputF"] = windowOutputF
		
		
		frequencyOutputF = File("%s.frequency.tsv"%(outputFnamePrefix))
		extraOutputLs.append(frequencyOutputF)
		key2ObjectForJob["frequencyOutputF"] = frequencyOutputF
		
		
		depthOutputF = File("%s.vs.depth.tsv"%(outputFnamePrefix))
		extraOutputLs.append(depthOutputF)
		key2ObjectForJob["depthOutputF"] = depthOutputF
		
		job = self.addGenericJob(executable=executable, inputFile=inputFile, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=None, \
					extraOutputLs=extraOutputLs, \
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=extraArgumentList, \
					transferOutput=transferOutput, sshDBTunnel=None, \
					key2ObjectForJob=key2ObjectForJob, objectWithDBArguments=None, objectWithDBGenomeArguments=None,\
					no_of_cpus=None, job_max_memory=3000, walltime=180, \
					max_walltime=None)
		return job
	
	def addPlotJob(self, workflow=None, PlotExecutable=None, outputFnamePrefix=None, outputFnameToRegisterLs=None,\
				title=None, parentJob=None, parentJobLs=None, transferOutput=True):
		"""
		2013.09.04 modernized
		2011-10-20
			outputFname is the argument passed to plotJob.
			outputFnameToRegisterLs contains a list of all output from this job that need to be registered and transferred.
		2011-9-29
		"""
		extraOutputLs=[]
		extraArgumentList = ["-o", outputFnamePrefix, "-t", title]
		for outputFname in outputFnameToRegisterLs:
			outputF = File(outputFname)
			extraOutputLs.append(outputF)
		
		job = self.addGenericJob(executable=PlotExecutable, inputFile=None, inputArgumentOption="-i", \
					outputFile=None, outputArgumentOption="-o", \
					parentJob=None, parentJobLs=parentJobLs, extraDependentInputLs=None, \
					extraOutputLs=extraOutputLs, \
					frontArgumentList=None, extraArguments=None, \
					extraArgumentList=extraArgumentList, \
					transferOutput=transferOutput, sshDBTunnel=None, \
					no_of_cpus=None, job_max_memory=3000, walltime=180, \
					max_walltime=None)
		return job
	
	def addParentToPlotJob(self, workflow=None, parentJob=None, parentOutputF=None, plotJob=None):
		"""
		2011-9-29
		"""
		if workflow is None:
			workflow = self
		plotJob.addArguments(parentOutputF)
		plotJob.uses(parentOutputF, transfer=True, register=True, link=Link.INPUT)
		if parentJob:
			workflow.depends(parent=parentJob, child=plotJob)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		if workflow is None:
			workflow = self
		
		vervetSrcPath = self.vervetSrcPath
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "mapper/CalculateTrioInconsistency.py"),\
												name="CalculateTrioInconsistency", clusterSizeMultipler=1)
		
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencySummaryHist.py"),\
												name="PlotTrioInconsistencySummaryHist", clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyOverPosition.py"),\
												name="PlotTrioInconsistencyOverPosition", clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyOverFrequency.py"),\
												name="PlotTrioInconsistencyOverFrequency", clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyVsDepth.py"),\
												name="PlotTrioInconsistencyVsDepth", clusterSizeMultipler=1)
		
		self.addOneExecutableFromPathAndAssignProperClusterSize(path=os.path.join(vervetSrcPath, "reducer/ReduceTrioInconsistencyByPosition.py"),\
												name="ReduceTrioInconsistencyByPosition", clusterSizeMultipler=1)
		
	
	def addCalculateTrioInconsistencyJobs(self, workflow=None, inputVCFData=None, trioLs=None, \
				addTrioSpecificPlotJobs=None, addTrioContigSpecificPlotJobs=None,\
				outputDirPrefix="", **keywords):
		"""
		addTrioContigSpecificPlotJobs: whether to add mendelian-inconsistency plotting jobs for each trio and each contig
		addTrioSpecificPlotJobs: whether to add mendelian-inconsistency plotting jobs for each trio

		2012.10.19 rename it to addAllJobs()
		2011.1.8
			add outputDirPrefix to differentiate one run from another if multiple trio call workflows are run simultaneously
			outputDirPrefix could contain "/" to denote sub-folders.
		"""
		if workflow is None:
			workflow =self
		if addTrioSpecificPlotJobs is None:	#use class's own option to figure out what to do
			addTrioSpecificPlotJobs = getattr(self, 'addTrioSpecificPlotJobs', False)
		if addTrioContigSpecificPlotJobs is None:
			addTrioContigSpecificPlotJobs = getattr(self, 'addTrioContigSpecificPlotJobs', False)
		
		if not trioLs:
			#2012.1.9 try to guess it from the first VCF file
			#trioLs = self.getAllTrios(self.db_vervet, aln_ref_ind_seq_id=self.aln_ref_ind_seq_id)
			if inputVCFData and inputVCFData.jobDataLs:
				firstInputF = inputVCFData.jobDataLs[0].file
				if hasattr(firstInputF, 'abspath') and os.path.isfile(firstInputF.abspath):
					trioLs = self.getDuoTrioThatExistInVCF(self.db_vervet, firstInputF.abspath)
		if not trioLs:
			#still empty. stop here
			sys.stderr.write("No trios specified/found. Couldn't add trio inconsistency calculation jobs.\n")
			return None
		sys.stderr.write("Adding trio-inconsistency calculation jobs for %s trios ..."%(len(trioLs)))
		
		returnJobData = PassingData()
		
		contig_id2trioInconsistencyJobLs = {}
		
		
		topOutputDir = "%sTrioInconsistency"%(outputDirPrefix)
		topOutputDirJob = self.addMkDirJob(outputDir=topOutputDir)
		
		#each contig in each trio gets a summary.
		trioInconsistencyByContigSummaryFile = File(os.path.join(topOutputDir, 'trio_inconsistency_by_contig_homo_het.tsv'))
		trioInconsistencyByContigSummaryMergeJob = self.addStatMergeJob(statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
							outputF=trioInconsistencyByContigSummaryFile, transferOutput=True, parentJobLs=[topOutputDirJob])
		
		#reduce the trio consistency by trio (each trio gets a summary)
		trioInconsistencySummaryFile = File(os.path.join(topOutputDir, 'trio_inconsistency.tsv'))
		trioInconsistencySummaryJob = self.addStatMergeJob(statMergeProgram=workflow.ReduceMatrixBySumSameKeyColsAndThenDivide, \
						outputF=trioInconsistencySummaryFile, transferOutput=True, extraArguments='-k 0 -v 4,5', parentJobLs=[topOutputDirJob])
		self.addInputToStatMergeJob(statMergeJob=trioInconsistencySummaryJob, \
								inputF=trioInconsistencyByContigSummaryMergeJob.output, \
								parentJobLs=[trioInconsistencyByContigSummaryMergeJob])
		returnJobData.trioInconsistencySummaryJob = trioInconsistencySummaryJob
		
		windowSize = 200000	#trio inconsistency calculation windows
		for trio_id in trioLs:
			# Add a mkdir job for the call directory.
			#letting numerou genotype call jobs detect&create this directory runs into race condition.
			trioReprInFilename = trio_id.replace(",", "_")
			trioDir = "%sTrio_%s"%(outputDirPrefix, trioReprInFilename)
			trioDirJob = self.addMkDirJob(outputDir=trioDir)
			
			#homoOnlyDir = os.path.join(trioDir, "homoOnly")
			#homoOnlyDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=homoOnlyDir)
			#workflow.depends(parent=trioDirJob, child=homoOnlyDirJob)
			
			allSitesDir = os.path.join(trioDir, "HomoHet")
			allSitesDirJob = self.addMkDirJob(outputDir=allSitesDir)
			workflow.depends(parent=trioDirJob, child=allSitesDirJob)
			
			#depthOutputDir = os.path.join(trioDir, "depth")
			#depthOutputDirJob = self.addMkDirJob(outputDir=depthOutputDir)
			
			if addTrioContigSpecificPlotJobs:
				positionOutputDir = os.path.join(trioDir, "Position")
				positionOutputDirJob = self.addMkDirJob(outputDir=positionOutputDir)
			
			# no space in a single argument
			title = "trio-%s-%s-files"%(trio_id, len(inputVCFData.jobDataLs))
			if addTrioSpecificPlotJobs:
				"""
				#2012.8.1 no more homo-site only jobs
				
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_summary_hist_homo_only.png'%(trioReprInFilename))
				summaryHomoOnlyPlotJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencySummaryHist, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
							parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_summary_hist_homo_het.png'%(trioReprInFilename))
				summaryAllSitesPlotJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencySummaryHist, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
							parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_over_position_homo_only.png'%(trioReprInFilename))
				inconsistencyOverPositionHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverPosition, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
							parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_over_position_homo_het.png'%(trioReprInFilename))
				inconsistencyOverPositionAllSitesJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverPosition, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_over_frequency_homo_only.png'%(trioReprInFilename))
				inconsistencyOverFrequencyHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverFrequency, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFname = os.path.join(topOutputDir, '%s_inconsistency_over_frequency_homo_het.png'%(trioReprInFilename))
				inconsistencyOverFrequencyAllSitesJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverFrequency, \
							outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
							parentJobLs=[topOutputDirJob],\
							title=title)
				"""
				outputFnamePrefix = '%s_homo_only_inconsistency_over'%(trioReprInFilename)
				fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
				fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
				fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
				fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
				mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
				mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
	
				inconsistencyOverDepthHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyVsDepth, \
							outputFnamePrefix=outputFnamePrefix, outputFnameToRegisterLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
									fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
							namespace=namespace, version=version, parentJob=trioDirJob,\
							title=title)
				
				outputFnamePrefix = '%s_homo_het_inconsistency_over'%(trioReprInFilename)
				fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
				fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
				fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
				fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
				mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
				mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
				inconsistencyOverDepthAllSitesJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyVsDepth, \
							outputFnamePrefix=outputFnamePrefix, outputFnameToRegisterLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
										fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
							namespace=namespace, version=version, parentJob=trioDirJob,\
							title=title)
				"""
			for jobData in inputVCFData.jobDataLs:
				inputFile = jobData.file
				contig_id = self.getContigIDFromFname(inputFile.name)
				if self.maxContigID:
					try:
						contig_id = int(contig_id)
						if contig_id>self.maxContigID:	#skip the small contigs
							continue
					except:
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
				if contig_id not in contig_id2trioInconsistencyJobLs:
					contig_id2trioInconsistencyJobLs[contig_id] = []
				
				"""
				#2012.8.1 no more homo-site only jobs
				
				outputFnamePrefix  = os.path.join(homoOnlyDir, '%s.inconsistency.homoOnly'%(os.path.basename(inputFile.name)))
				trioInconsistencyCaculationHomoOnlyJob = self.addTrioInconsistencyCalculationJob(workflow, \
							executable=workflow.CalculateTrioInconsistency, \
							inputFile=inputFile, outputFnamePrefix=outputFnamePrefix, \
							trio_id=trio_id,\
							parentJobLs=[homoOnlyDirJob] + jobData.jobLs, \
							windowSize=windowSize, minDepth=self.minDepth)	#homoOnly
				"""
				outputFnamePrefix  = os.path.join(allSitesDir, '%s.inconsistency.homoHet'%(os.path.basename(inputFile.name)))
				trioInconsistencyCaculationAllSitesJob = self.addTrioInconsistencyCalculationJob(\
							executable=workflow.CalculateTrioInconsistency, \
							inputFile=inputFile, outputFnamePrefix=outputFnamePrefix, \
							trio_id=trio_id,\
							windowSize=windowSize, minDepth=self.minDepth,\
							parentJobLs=[allSitesDirJob] + jobData.jobLs)	#all sites
				contig_id2trioInconsistencyJobLs[contig_id].append(trioInconsistencyCaculationAllSitesJob)
				#add output to some reduce job
				self.addInputToStatMergeJob(statMergeJob=trioInconsistencyByContigSummaryMergeJob, \
								inputF=trioInconsistencyCaculationAllSitesJob.summaryOutputF, \
								parentJobLs=[trioInconsistencyCaculationAllSitesJob])
				
				if addTrioSpecificPlotJobs:
					"""
					#2012.8.1 no more homo-site only jobs
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.summaryOutputF, \
									plotJob=summaryHomoOnlyPlotJob)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.windowOutputF, \
									plotJob=inconsistencyOverPositionHomoOnlyJob)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.frequencyOutputF, \
									plotJob=inconsistencyOverFrequencyHomoOnlyJob)
					"""
					"""
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.depthOutputF, \
									plotJob=inconsistencyOverDepthHomoOnlyJob)
					"""
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.summaryOutputF, \
									plotJob=summaryAllSitesPlotJob)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.windowOutputF, \
									plotJob=inconsistencyOverPositionAllSitesJob)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.frequencyOutputF, \
									plotJob=inconsistencyOverFrequencyAllSitesJob)
					"""
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.depthOutputF, \
									plotJob=inconsistencyOverDepthAllSitesJob)
					"""
				
				if addTrioContigSpecificPlotJobs:
					newTitle = "trio-%s-chr-%s"%(trio_id, contig_id)
					"""
					#2012.8.1 no more homo-site only jobs
					outputFname = os.path.join(positionOutputDir, '%s_chromosome_%s_inconsistency_over_position_homo_only.png'%(trioDir, contig_id))
					contigInconsistencyOverPositionHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverPosition, \
								outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
								parentJob=positionOutputDirJob, title=newTitle)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.windowOutputF, \
									plotJob=contigInconsistencyOverPositionHomoOnlyJob)
					"""
					outputFname = os.path.join(positionOutputDir, '%s_chromosome_%s_inconsistency_over_position_homo_het.png'%(trioDir, contig_id))
					contigInconsistencyOverPositionAllSitesJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyOverPosition, \
								outputFnamePrefix=outputFname, outputFnameToRegisterLs=[outputFname], \
								parentJob=positionOutputDirJob, title=newTitle)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.windowOutputF, \
									plotJob=contigInconsistencyOverPositionAllSitesJob)
					"""
					outputFnamePrefix = os.path.join(depthOutputDir, '%s_chromosome_%s_homo_only_inconsistency_over'%(trioDir, contig_id))
					fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
					fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
					mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
					fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
					fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
					mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
					contigInconsistencyOverDepthHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyVsDepth, \
								outputFnamePrefix=outputFnamePrefix, outputFnameToRegisterLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
									fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
								namespace=namespace, version=version, parentJob=depthOutputDirJob,\
								title=newTitle)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
									parentOutputF=trioInconsistencyCaculationHomoOnlyJob.depthOutputF, \
									plotJob=contigInconsistencyOverDepthHomoOnlyJob)
					
					outputFnamePrefix = os.path.join(depthOutputDir, '%s_chromosome_%s_homo_het_inconsistency_over'%(trioDir, contig_id))
					fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
					fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
					mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
					fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
					fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
					mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
					contigInconsistencyOverDepthAllSitesJob = self.addPlotJob(workflow, PlotExecutable=workflow.PlotTrioInconsistencyVsDepth, \
								outputFnamePrefix=outputFnamePrefix, outputFnameToRegisterLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,
										fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
								namespace=namespace, version=version, parentJob=depthOutputDirJob,\
								title=newTitle)
					self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
									parentOutputF=trioInconsistencyCaculationAllSitesJob.depthOutputF, \
									plotJob=contigInconsistencyOverDepthAllSitesJob)
					"""
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		
		sys.stderr.write(" \t Reducing %s chromosome trio inconsistency calculation jobs  ..."%(len(contig_id2trioInconsistencyJobLs)))
		
		if len(contig_id2trioInconsistencyJobLs)>0:
			#2011.12.16 reduce job to calculate average duo/trio inconsistency per position 
			avgTrioInconsistencyByPositionFile = File(os.path.join(topOutputDir, 'avgTrioInconsistencyByPosition.tsv'))
			avgTrioInconsistencyByPositionMergeJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
								outputF=avgTrioInconsistencyByPositionFile, transferOutput=False, parentJobLs=[topOutputDirJob])
			
			contig_id_ls = contig_id2trioInconsistencyJobLs.keys()
			contig_id_ls.sort()
			for contig_id in contig_id_ls:
				trioInconsistencyJobLs = contig_id2trioInconsistencyJobLs.get(contig_id)
				
			#reduce job to calculate average trio inconsistency for each chromosome.
				oneChrAvgTrioInconsistencyByPositionFile = File(os.path.join(topOutputDir, 'chr_%s_avgTrioInconsistencyByPosition.tsv'%(contig_id)))
				oneChrAvgTrioInconsistencyByPositionJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceTrioInconsistencyByPosition, \
									outputF=oneChrAvgTrioInconsistencyByPositionFile, transferOutput=False, parentJobLs=[topOutputDirJob])
				for job in trioInconsistencyJobLs:
					self.addInputToStatMergeJob(workflow, statMergeJob=oneChrAvgTrioInconsistencyByPositionJob, inputF=job.depthOutputF, \
								parentJobLs=[job])
				
				self.addInputToStatMergeJob(workflow, statMergeJob=avgTrioInconsistencyByPositionMergeJob, \
										inputF=oneChrAvgTrioInconsistencyByPositionJob.output, \
										parentJobLs=[oneChrAvgTrioInconsistencyByPositionJob])
			#bgzip and index the avgTrioInconsistencyByPositionFile
			avgTrioInconsistencyByPosGzipFile = File("%s.gz"%avgTrioInconsistencyByPositionFile.name)
			avgTrioInconsistencyByPosBGZipTabixJob = self.addBGZIP_tabix_Job(workflow, bgzip_tabix=workflow.bgzip_tabix, \
								parentJob=avgTrioInconsistencyByPositionMergeJob, inputF=avgTrioInconsistencyByPositionFile, \
								outputF=avgTrioInconsistencyByPosGzipFile, parentJobLs=[topOutputDirJob],\
								transferOutput=True, tabixArguments="-s 1 -b 2 -e 2")
			
			returnJobData.avgTrioInconsistencyByPosBGZipTabixJob = avgTrioInconsistencyByPosBGZipTabixJob
		
		sys.stderr.write("%s jobs.\n"%(self.no_of_jobs))
		return returnJobData

	addAllJobs = addCalculateTrioInconsistencyJobs
if __name__ == '__main__':
	main_class = CalculateTrioInconsistencyPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
