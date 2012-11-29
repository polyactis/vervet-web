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
	%s -I ~/NetworkData/vervet/db/genotype_file/method_6/ -q ./aux/alnStatForFilter.2012.8.1T1805.tsv
		-L 1.0 -a 524  -E -H -o workflow/FilterGenotypeMethod6_ByMaskingZeroDPSite_2FoldDepthFilter_MaxSNPMissing1.0.xml 
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -C 5 -Z 1 -A 2 -g 5
		#-V 90 -x 100 --maxMendelError 6
	
Description:
	2012.9.13 pipeline that does vcf2plink, plink mendel, select sites <=maxMendelError, filter VCF
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from vervet.src import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from vervet.src.AbstractVervetWorkflow import AbstractVervetWorkflow
#from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
from FilterVCFPipeline import FilterVCFPipeline
from vervet.src.PlinkOnVCFWorkflow import PlinkOnVCFWorkflow

class FilterVCFByMendelErrorWorkflow(FilterVCFPipeline, PlinkOnVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = FilterVCFPipeline.option_default_dict.copy()
	option_default_dict.pop(('vcf2Dir', 0, ))
	option_default_dict.update({
				('maxMendelError', 1, int): [6, '', 1, 'sites with mendel errors below this number are discarded.', ],\
				})

	def __init__(self,  **keywords):
		"""
		"""
		FilterVCFPipeline.__init__(self, **keywords)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		FilterVCFPipeline.registerCustomExecutables(self, workflow=workflow)
		PlinkOnVCFWorkflow.registerCustomExecutables(self, workflow=workflow)
		
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
		if self.depthFoldChange and self.depthFoldChange>0:
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=self.depthFoldChange, minGQ=self.minGQ)
			alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname))
		else:
			alnStatForFilterF = None
		
		
		if self.vcf1Dir:
			# 2012.5.1 filter only on the 1st vcf folder
			#a relative-path name for vcf1Dir
			vcf1Name = self.findProperVCFDirIdentifier(self.vcf1Dir, defaultName='vcf1')
			inputData = self.registerAllInputFiles(workflow, self.vcf1Dir, input_site_handler=self.input_site_handler, \
												checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
												pegasusFolderName="%s_%s"%(self.pegasusFolderName, vcf1Name), \
												maxContigID=self.maxContigID, \
												minContigID=self.minContigID)
			
			vcf2PlinkJobData = self.addVCF2PlinkJobs(workflow, inputData=inputData, db_vervet=db_vervet, minMAC=None, minMAF=None,\
							maxSNPMissingRate=None, transferOutput=False,\
							maxContigID=self.maxContigID, outputDirPrefix="vcf2plink", outputPedigreeAsTFAM=True,\
							treatEveryOneIndependent=False,\
							returnMode=2, ModifyTPEDRunType=1, chr_id2cumu_chr_start=None)
			mendelJobData = self.addPlinkMendelErrorJobs(inputData=vcf2PlinkJobData, transferOutput=True,\
							maxContigID=self.maxContigID, outputDirPrefix="mendel", locusSamplingRate=0.1, 
							returnMode=2)
			
			locusMendelJobData = mendelJobData.jobDataLs[-1]	#last job from PlinkMendelWorkflow is the locus-mendel merge job.
			keepSNPPosF = File('sitesWithMax%sMendelError.tsv'%(self.maxMendelError))
			outputSitesBelowMaxMendelJob = self.addGenericJob(executable=self.OutputSitesBelowMaxMendelError, inputFile=locusMendelJobData.file, \
					inputArgumentOption="-i", \
					outputFile=keepSNPPosF, outputArgumentOption="-o", \
					parentJobLs=locusMendelJobData.jobLs, extraDependentInputLs=None, extraOutputLs=None, transferOutput=True, \
					extraArguments="-m %s"%(self.maxMendelError), extraArgumentList=None, job_max_memory=2000,  sshDBTunnel=None, \
					key2ObjectForJob=None)
			
			self.addJobsToFilterOneVCFDir(workflow, inputData=inputData, refFastaFList=refFastaFList, \
									alnStatForFilterF=alnStatForFilterF, keepSNPPosF=keepSNPPosF, \
									onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
									minMAC=self.minMAC, minMAF=self.minMAF, maxSNPMissingRate=self.maxSNPMissingRate,\
									minDepthPerGenotype=self.minDepthPerGenotype, outputDirPrefix="filter",\
									minNeighborDistance=self.minNeighborDistance, keepSNPPosParentJobLs=[outputSitesBelowMaxMendelJob])
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = FilterVCFByMendelErrorWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()