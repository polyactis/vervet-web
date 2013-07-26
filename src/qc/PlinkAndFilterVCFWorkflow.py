#!/usr/bin/env python
"""
Examples:
	
	#2011.12.19 run on hoffman2's condorpool
	%s ....
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		
	#2012.5.1 filter trioCaller output with total depth (no minGQ filter anymore), minMAC=10 (-n 10), maxSNPMismatchRate=1 (-L 1.0)
	# minMAF=0.05 (-f 0.05), no depth-band filter (-A 0)
	%s -I AlignmentToTrioCallPipeline_VRC_top7559Contigs.2011.12.15T0057/trioCaller/ -l condorpool -j condorpool
		-z uclaOffice -u yh -q ./alnStatForFilter.2012.5.1T1430.tsv  -n 10 -L 1.0 -f 0.05 -A 0 -a 524 -C 50 -E -H
		-o FilterVCF_trioCallerTop7559Contigs.xml
	
	#2012.8.1 FilterGenotypeMethod5_ByMethod7Sites (-S ...) NoDepthFilter (-A 0) MaxSNPMissing0.5 (-L 0.5)
	%s -I ~/NetworkData/vervet/db/genotype_file/method_5/ -q ./aux/alnStatForFilter.2012.7.30T1542.tsv
		-L 0.5 -a 524  -E -H -o workflow/FilterGenotypeMethod5_ByMethod7Sites_NoDepthFilter_MaxSNPMissing0.5.xml  -l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh -C 5 -S ./method7_sites.tsv -A 0
	
	#2012.8.1 FilterGenotypeMethod6_ByMaskingZeroDPSite (-Z 1) 2FoldDepthFilter (-A 2) MaxSNPMissing1.0 (-L 1.0)
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	# "-g 5" to the minimum distance between neighboring SNPs
	# "-H" to require db ssh tunnel for db interacting jobs
	# add "-y2" to do pruning, then filter based on those pruned sites
	# for "-y2", also add "-W 50 --LDPruneWindowShiftSize 20 -R $R" for LD pruning parameters
	# W=500; Z=200; R=0.5 ;
	%s -I ~/NetworkData/vervet/db/genotype_file/method_6/ -q ./aux/alnStatForFilter.2012.8.1T1805.tsv
		-L 1.0 -a 524  -E -H -o workflow/PlinkLDPrune/LDPrune_Method38_2VCF_AndLD_contig96_100_W$W\Z$Z\R$R\.xml  
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -C 5 -Z 1 -A 2 -g 5
		#-V 90 -x 100 --locusSamplingRate 1
		#-y1 --maxMendelError 6
		#-y2 --mergeListFname ./aux/Method38_LDPrune_merge_list.2012.9.13T1339 -W $W --LDPruneWindowShiftSize $Z -R $R
	
Description:
	2012.9.13 pipeline that does vcf2plink, plink mendel or LD pruning, select sites (<=maxMendelError), filter VCF
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import Executable, PFN, File
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from vervet.src.qc.FilterVCFPipeline import FilterVCFPipeline
from vervet.src.PlinkOnVCFWorkflow import PlinkOnVCFWorkflow

parentClass = FilterVCFPipeline
class PlinkAndFilterVCFWorkflow(parentClass, PlinkOnVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('vcf2Dir', 0, ))
	option_default_dict.update({
				('maxMendelError', 1, int): [6, '', 1, 'sites with #mendel errors beyond this number are discarded.', ],\
				})
	option_default_dict.update(
						{
						('LDPruneMinR2', 0, float): [0.4, 'R', 1, 'minimum r2 for LD pruning', ],\
						('locusSamplingRate', 0, float): [0.01, 'c', 1, 'how many loci to sample', ],\
						('mergeListFname', 0, ): [None, '', 1, 'the file to contain the merge-list for plink, required for run_type>1', ],\
						('run_type', 1, int): [1, 'y', 1, 'which run_type to run. \n\
		1: plink mendel, select sites <=maxMendelError, filter VCF\n\
		2: LD pruning, select sites, filter VCF \n\
		3: plink mendel, select sites <=maxMendelError, filter VCF, LD pruning, select sites, filter VCF\n\
		any other run_type is simply filter VCF.\n'
		],\
						('LDPruneWindowSize', 1, int): [50, 'W', 1, ' window size (in the number of SNPs, not bp) for plink LD pruning'],\
						('LDPruneWindowShiftSize', 1, int): [20, '', 1, 'adjacent window shift (in the number of SNPs), not bp '],\
					})
	#2013.7.17 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 10000
	option_default_dict[('max_walltime', 1, int)][0] = 1320	#under 23 hours

	def __init__(self,  **keywords):
		"""
		"""
		parentClass.__init__(self, **keywords)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		PlinkOnVCFWorkflow.registerCustomExecutables(self, workflow=workflow)
		#CalculateVCFStatPipeline.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		
		OutputSitesBelowMaxMendelError = Executable(namespace=namespace, name="OutputSitesBelowMaxMendelError", version=version, \
							os=operatingSystem, arch=architecture, installed=True)
		OutputSitesBelowMaxMendelError.addPFN(PFN("file://" + os.path.join(self.pymodulePath, "pegasus/mapper/OutputSitesBelowMaxMendelError.py"), \
							site_handler))
		executableClusterSizeMultiplierList.append((OutputSitesBelowMaxMendelError, 0))
		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
	def run(self):
		"""
		"""
		
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, db_user=self.db_user,
							db_passwd=self.db_passwd, hostname=self.hostname, dbname=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		"""
		2012.9.13 chr2size is only needed for LD calculation jobs later (commented out right now)
		#chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
		#									sequence_type_id=9)
		"""
		if self.run_type==2 or self.run_type==3:	#run_type that involve LD-pruning at one step or another
			#chrOrder=2 means chromosomes are not ordered alphabetically but by their sizes (descendingly)
			oneGenomeData = db_genome.getOneGenomeData(tax_id=self.ref_genome_tax_id, chr_gap=0, chrOrder=2, \
													sequence_type_id=self.ref_genome_sequence_type_id,\
													maxPseudoChrSize=1000000000)	#plink could not handle a chromosome of >2^31 bp length
			chr_id2cumu_chr_start = oneGenomeData.chr_id2cumu_chr_start
			if self.run_type==4:
				#fake selected contigs as X chromosomes and ignore others.
				sexContigList = ['Contig83', 'Contig149', 'Contig193']
				for contig in sexContigList:
					if contig in chr_id2cumu_chr_start:
						newChr, cumuStart = chr_id2cumu_chr_start.get(contig)
						chr_id2cumu_chr_start[contig] = ['X', cumuStart]	#assign some contigs to X
				#2012.9.12 contigs beyond 193 is a waste of time.
				self.maxContigID=195
				sys.stderr.write("Warning: maxContigID is manually set to 195 since contigs with ID bigger than 193 are useless.\n")
			if not self.mergeListFname:
				sys.stderr.write("Error: mergeListFname %s is nothing. required for this run_type %s.\n"%(self.mergeListFname,\
																										self.run_type))
				sys.exit(3)
			
			#for LD pruning, need to get rid files with too few loci, used in self.setup_run()
			#files with too few loci (like 1) cause problem by becoming empty after LD-prune (why? at least one SNP).
			self.minNoOfLociInVCF = 2
		else:
			chr_id2cumu_chr_start = None
			self.minNoOfLociInVCF = 0
		
		#below is for the 1st round of vcf2plink 
		if self.run_type==2:	#all run_type that involve LD-pruning
			ModifyTPEDRunType = 3	#plink LD-prune will skip non-recognizable chromosomes (1-24, human), so fake the chromosome
				#using chr_id2cumu_chr_start
			#plink mendel job needs full VRC pedigree, but LD pruning does not
			treatEveryOneIndependent = True
		else:
			ModifyTPEDRunType = 1	#plink mendel doesn't skip non-human chromosomes
			#plink mendel job needs full VRC pedigree
			treatEveryOneIndependent = False
		
		if self.run_type in [1,3]:
			addUngenotypedDuoParents = True
		else:
			addUngenotypedDuoParents = False
		
		
		self.needToKnowNoOfLoci = True
		self.setup_run()
		
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		"""
		

		"""
		if self.depthFoldChange and self.depthFoldChange>0:
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
			alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
		else:
			alnStatForFilterF = None
		"""
		alnStatForFilterF = None
		
		if not self.vcf1Dir:
			sys.stderr.write("ERRROR vcf1Dir %s is not available.\n"%(self.vcf1Dir))
			sys.exit(3)
		# 2012.5.1 filter only on the 1st vcf folder
		vcf2PlinkJobData1 = self.addVCF2PlinkJobs(inputData=self.inputData, db_vervet=self.db_vervet, minMAC=None, minMAF=None,\
						maxSNPMissingRate=None, transferOutput=False,\
						maxContigID=self.maxContigID, outputDirPrefix="vcf2plinkRuntype%s_s1"%(self.run_type), \
						outputPedigreeAsTFAM=True,\
						treatEveryOneIndependent=treatEveryOneIndependent,\
						returnMode=2, ModifyTPEDRunType=ModifyTPEDRunType, chr_id2cumu_chr_start=chr_id2cumu_chr_start,\
						addUngenotypedDuoParents=addUngenotypedDuoParents)
		
		if self.run_type ==1 or self.run_type==3:
			mendelWorkflowData = self.addPlinkMendelErrorJobs(inputData=vcf2PlinkJobData1, transferOutput=True,\
						maxContigID=self.maxContigID, outputDirPrefix="mendelRuntype%s_s2"%(self.run_type), \
						locusSamplingRate=self.locusSamplingRate, 
						returnMode=2)
		
			lmendelMergeJob = mendelWorkflowData.lmendelMergeJob	#last job from PlinkMendelWorkflow is the locus-mendel merge job.
			keepSNPPosF = File('sitesWithMax%sMendelError.s2.tsv'%(self.maxMendelError))
			outputSitesBelowMaxMendelJob = self.addGenericJob(executable=self.OutputSitesBelowMaxMendelError, \
				inputFile=lmendelMergeJob.output, \
				inputArgumentOption="-i", \
				outputFile=keepSNPPosF, outputArgumentOption="-o", \
				parentJobLs=[lmendelMergeJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
				extraArguments="-m %s"%(self.maxMendelError), extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
				key2ObjectForJob=None)
			keepSNPPosParentJobLs = [outputSitesBelowMaxMendelJob]
		elif self.run_type==2:
			mergeListFile = self.registerOneInputFile(inputFname=self.mergeListFname, folderName=self.pegasusFolderName,\
													checkFileExistence=False)
			LDPruneWorkflowData = self.addPlinkLDPruneJobs(inputData=vcf2PlinkJobData1, transferOutput=True,\
					maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
					LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
					outputDirPrefix="ldPruneRuntype%s_s2"%(self.run_type), returnMode=1, \
					mergeListFile=mergeListFile)
			#plinkMergeJobData = LDPruneWorkflowData.jobDataLs[-1]	#last job from LD-prune workflow is the plink merge job.
			#plinkMergeJob = plinkMergeJobData.jobLs[0]
			plinkMergeJob = LDPruneWorkflowData.plinkMergeJob
			keepSNPPosF = File('LDPrunedSitesW%sZ%sR%s.s2.tsv'%(self.LDPruneWindowSize, self.LDPruneWindowShiftSize, self.LDPruneMinR2))
			convertJob = self.addGenericJob(executable=self.ConvertPlinkBIM, inputFile=plinkMergeJob.bimFile, \
				inputArgumentOption="-i", \
				outputFile=keepSNPPosF, outputArgumentOption="-o", \
				parentJobLs=[plinkMergeJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
				extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
				key2ObjectForJob=None)
			keepSNPPosParentJobLs = [convertJob]
		
		#select sites and filter
		filterReturnData = self.addJobsToFilterOneVCFDir(inputData=self.inputData, \
								registerReferenceData=self.registerReferenceData, \
								alnStatForFilterF=alnStatForFilterF, keepSNPPosF=keepSNPPosF, \
								onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
								minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate,\
								minDepthPerGenotype=self.minDepthPerGenotype, outputDirPrefix="filterRunType%s_s3"%(self.run_type),\
								minNeighborDistance=self.minNeighborDistance, transferOutput=True,\
								keepSNPPosParentJobLs=keepSNPPosParentJobLs)
		if self.run_type==3:	#add LD pruning
			#VCF 2 plink for filtered VCFs
			# outputPedigreeAsTFAM=True because addUngenotypedDuoParents=True needs it
			# treatEveryOneIndependent=True because LD pruning use only founders genotype, if False, will use 2 or 3 founders (plink definition) only
			# ModifyTPEDRunType=3
			# addUngenotypedDuoParents = True (because Mark mendel error genotype missing needs it)
			vcf2PlinkJobData2 = self.addVCF2PlinkJobs(inputData=filterReturnData, db_vervet=self.db_vervet, minMAC=None, minMAF=None,\
							maxSNPMissingRate=None, transferOutput=False,\
							maxContigID=self.maxContigID, outputDirPrefix="vcf2plinkRuntype%s_s4"%(self.run_type),\
							outputPedigreeAsTFAM=True,\
							treatEveryOneIndependent=True,\
							returnMode=2, ModifyTPEDRunType=3, chr_id2cumu_chr_start=chr_id2cumu_chr_start,\
							addUngenotypedDuoParents=True)
			
			#copy the family-file jobs because their files (famFile, tfamFile) are modified inside markMendelErrorMissingJobData
			LDPruneFamilyFileJob = copy.deepcopy(vcf2PlinkJobData2.famJob)
			LDPruneTFamilyFileJob = copy.deepcopy(vcf2PlinkJobData2.tfamJob)
			
			#2013.07.24 use vcf2PlinkJobData1's family job because it has a real pedigree. Everyone is independent in vcf2PlinkJobData2 's family file
			# for LD pruning. 
			markMendelErrorMissingJobData = self.markMendelErrorLociMissingSubWorkflow(inputData=vcf2PlinkJobData2, \
												transferOutput=False, \
												maxContigID=self.maxContigID, \
												famJob=vcf2PlinkJobData1.famJob, tfamJob=vcf2PlinkJobData1.tfamJob, \
												outputDirPrefix='markMendelErrorGenotypeMissing_s5', returnMode=2)
			
			#LD pruning
			mergeListFile = self.registerOneInputFile(inputFname="%s.merge.s6"%self.mergeListFname, \
													folderName=self.pegasusFolderName,\
													checkFileExistence=False)
			#2013.07.24 use vcf2PlinkJobData2's family job because everyone is independent in it
			LDPruneWorkflowData = self.addPlinkLDPruneJobs(inputData=markMendelErrorMissingJobData, transferOutput=True,\
							famJob=LDPruneFamilyFileJob, tfamJob=LDPruneTFamilyFileJob, \
							maxContigID=self.maxContigID, LDPruneMinR2=self.LDPruneMinR2, \
							LDPruneWindowSize=self.LDPruneWindowSize, LDPruneWindowShiftSize=self.LDPruneWindowShiftSize, \
							outputDirPrefix="ldPruneRunType%s_s6"%(self.run_type), returnMode=1, \
							mergeListFile=mergeListFile)
			#plinkMergeJobData = LDPruneWorkflowData.jobDataLs[-1]	#last job from LD-prune workflow is the plink merge job.
			plinkMergeJob = LDPruneWorkflowData.plinkMergeJob
			keepSNPPosF = File('LDPrunedSitesW%sZ%sR%s.s6.tsv'%(self.LDPruneWindowSize, self.LDPruneWindowShiftSize, self.LDPruneMinR2))
			convertJob = self.addGenericJob(executable=self.ConvertPlinkBIM, inputFile=plinkMergeJob.bimFile, \
				inputArgumentOption="-i", \
				outputFile=keepSNPPosF, outputArgumentOption="-o", \
				parentJobLs=[plinkMergeJob], extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
				extraArguments=None, extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
				key2ObjectForJob=None)
			keepSNPPosParentJobLs = [convertJob]
			
			#pick the LD-pruned sites
			self.addJobsToFilterOneVCFDir(inputData=filterReturnData, \
						registerReferenceData=self.registerReferenceData, \
						keepSNPPosF=keepSNPPosF, \
						onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
						minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate,\
						minDepthPerGenotype=self.minDepthPerGenotype, outputDirPrefix="filterRunType%s_s7"%(self.run_type),\
						minNeighborDistance=self.minNeighborDistance, transferOutput=True,\
						keepSNPPosParentJobLs=keepSNPPosParentJobLs)
		
		"""
		import random
		statCalculationInputData = PassingData()
		statCalculationInputData.jobDataLs = []
		#randomly choose 6 contigs
		chosenIndexLs = random.sample(range(0,len(filterReturnData.jobDataLs)), 6)
		for i in chosenIndexLs:
			statCalculationInputData.jobDataLs.append(filterReturnData.jobDataLs[i])
		#statCalculationInputData = filterReturnData
		#caclualte LD on some data
		returnData = self.addStatCalculationJobs(workflow=workflow, inputData=statCalculationInputData, \
							registerReferenceData=self.registerReferenceData, \
								chr2size=chr2size, windowSize=100000, minChrLengthForPlot=500000, \
								minChrSize=100000, LDWindowSize=100000, outputDirPrefix="VCFStat",\
								transferOutput=True,\
								samplingRate=self.locusSamplingRate, minSiteGap=30000)
		"""
		self.end_run()
		


	
if __name__ == '__main__':
	main_class = PlinkAndFilterVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()