#!/usr/bin/env python
"""
Examples:
	
	#2011.12.19 run on hoffman2's condorpool
	%s ....
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti/
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-H
		
	#2012.5.1 filter trioCaller output with total depth (no minGQ filter anymore), minMAC=10 (-B 10), maxSNPMismatchRate=1 (-R 1.0)
	# minMAF=0.05 (-M 0.05), no depth-band filter (-A 0)
	%s -i AlignmentToTrioCallPipeline_VRC_top7559Contigs.2011.12.15T0057/trioCaller/ -l condorpool -j condorpool
		-z uclaOffice -u yh -q ./alnStatForFilter.2012.5.1T1430.tsv  -B 10 -R 1.0 -M 0.05 -A 0 -a 524 -C 50 -E
		-o FilterVCF_trioCallerTop7559Contigs.xml
	
	#2012.8.1 FilterGenotypeMethod5_ByMethod7Sites (-S ...) NoDepthFilter (-A 0) MaxSNPMissing0.5 (-R 0.5)
	%s -i ~/NetworkData/vervet/db/genotype_file/method_5/ -q ./aux/alnStatForFilter.2012.7.30T1542.tsv
		-R 0.5 -a 524  -E -o workflow/FilterGenotypeMethod5_ByMethod7Sites_NoDepthFilter_MaxSNPMissing0.5.xml  -l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/ -u yh -C 5 -S ./method7_sites.tsv -A 0 -H
	
	#2012.8.1 FilterGenotypeMethod6_ByMaskingZeroDPSite (-Z 1) 2FoldDepthFilter (-A 2) MaxSNPMissing1.0 (-R 1.0)
	%s -i ~/NetworkData/vervet/db/genotype_file/method_6/ -q ./aux/alnStatForFilter.2012.8.1T1805.tsv
		-R 1.0 -a 524  -E -o workflow/FilterGenotypeMethod6_ByMaskingZeroDPSite_2FoldDepthFilter_MaxSNPMissing1.0.xml 
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-u yh -C 5 -Z 1 -A 2 -H
	
	#2012.8.2 try different combinations of MAC (-B) and maxSNPMissingRateLs (-R) on method 9
	%s -i ~/NetworkData/vervet/db/genotype_file/method_9/
		-R 0.4,0.8,1.0 -a 524 -o workflow/BootstrapFilterCalculateVCF_Method9_R.4.8.1.0_B_2_8_12.xml
		-l hcondor -j hcondor
		-e /u/home/eeskin/polyacti/  -t /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-z localhost -u yh -C 60  -B 2,8,12 -H
	
	# 2012.8.2 try different combinations of MAC (-B) and maxSNPMissingRateLs (-R)
	# different minDepthPerGenotype (-Z 0,3), depthFoldChange=3 (-A) 
	# on method 6
	# remove non-biallelic SNPs (-K)
	# "-V 90 -x 100" are used to restrict contig IDs between 90 and 100.
	# LD within 400k (-L ...), all other stats (snp density, pi, ti/tv) are calculated in 200kb windows (-w 200000).
	# max minChrLengthForPlot=2million (-P), minChrSiz for computing =500k (-n)
	# need ssh tunnel (-H) 
	%s -i ~/NetworkData/vervet/db/genotype_file/method_6 -q ./aux/alnStatForFilter.2012.8.2T2250.tsv
		-R 0,0.8 -a 524 -o workflow/BootstrapFilterCalculateVCFMethod6_R.0.8._B_2_8_A3_Z0_3.xml
		-l hcondor -j hcondor -e /u/home/eeskin/polyacti/ 
		-t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-z localhost -u yh -C 60  -B 2,8 -A 3 -Z 0,3 -K
		-w 200000 -L 400000 -P 2000000 -n 500000 -H
		#-V 90 -x 100 
	
Description:
	2012.8.14 a program that combines FilterVCFPipeline.py and CalculateVCFStatPipeline.py into one workflow.
	2012.7.30 set depthFoldChange=0 to skip the filter by depth.
	
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], \
				sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from FilterVCFPipeline import FilterVCFPipeline
from CalculateVCFStatPipeline import CalculateVCFStatPipeline

class BootstrapFilterCalculateVCFWorkflow(FilterVCFPipeline, CalculateVCFStatPipeline):
	__doc__ = __doc__
	option_default_dict = CalculateVCFStatPipeline.option_default_dict.copy()
	option_default_dict.update(FilterVCFPipeline.common_filter_option_dict)
	option_default_dict.update({
						('minGQ', 1, int): [50, 'G', 1, 'minimum GQ/GenotypeQuality for one genotype. 2012.5.1 no longer enforced in FilterVCFByDepth.java', ],\
						('depthFoldChangeLs', 0, ): [None, 'A', 1, 'coma/dash-separated list of a variant is retained if its depth within this fold change of meanDepth,\
				set this to 0 or below to eliminate this step of filtering.', ],\
						("minDepthPerGenotypeLs", 0, ): [None, 'Z', 1, 'coma/dash-separated list of mask genotype with below this depth as ./. (other fields retained), \
	esp. necessary for SAMtools, which output homozygous reference if no read for one sample.'],\
						("minMACLs", 0, ): [None, 'B', 1, 'coma/dash-separated list of minimum MinorAlleleCount (by chromosome)'],\
						("minMAFLs", 0, ): [None, 'f', 1, 'coma/dash-separated list of minimum MinorAlleleFrequency (by chromosome)'],\
						("maxSNPMissingRateLs", 0, ): [None, 'R', 1, 'coma/dash-separated list of maximum SNP missing rate in one vcf (denominator is #chromosomes)'],\
						('inputDir', 1, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						})
	"""
	# 2012.8.2 
		('windowSize', 1, int): [200000, 'w', 1, 'window size for TiTv, pi, snpDensity calculation by vcftools', ],\
		('LDWindowSize', 0, int): [0, 'L', 1, 'window size for LD calculation by vcftools, set it to 0 to skip LD jobs.', ],\
		('minChrLengthForPlot', 1, int): [4000000, 'P', 1, 'minimum chromosome size for a chromosome to be included in plot', ],\
		('minChrSize', 1, int): [1000000, 'n', 1, 'minimum chromosome size for any computing ', ],\
		('includeIndelVCF', 1, int): [0, '', 1, 'toggle this to include indel VCF, filename with "indel"'],\
	"""
	def __init__(self,  **keywords):
		"""
		"""
		FilterVCFPipeline.__init__(self, **keywords)
		self.inputDir = os.path.abspath(self.inputDir)
		
		
		listArgumentName_data_type_ls = [('depthFoldChangeLs', int), ("minDepthPerGenotypeLs", int), \
								("minMACLs", int), ("minMAFLs", float), ("maxSNPMissingRateLs", float)]
		listArgumentName2hasContent = self.processListArguments(listArgumentName_data_type_ls, emptyContent=[None])
		if listArgumentName2hasContent['depthFoldChangeLs'] and not self.alnStatForFilterFname:
			sys.stderr.write("Error: alnStatForFilterFname (%s) is nothing while depthFoldChangeLs=%s.\n"%\
							(self.alnStatForFilterFname, repr(self.depthFoldChangeLs)))
			sys.exit(3)
		
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		FilterVCFPipeline.registerCustomExecutables(self, workflow=workflow)
		CalculateVCFStatPipeline.registerCustomExecutables(self, workflow=workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#have to be in front of the db_vervet connection code. Otherwise schema "genome" wont' be default path and its visible will not be visible.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
											sequence_type_id=9)
		
		
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
		
		depthFoldChange = self.depthFoldChangeLs[0]
		if depthFoldChange>0:
			self.outputAlignmentDepthAndOthersForFilter(db_vervet=db_vervet, outputFname=self.alnStatForFilterFname, \
												ref_ind_seq_id=self.ref_ind_seq_id, \
												foldChange=depthFoldChange, minGQ=self.minGQ)
			alnStatForFilterF = self.registerOneInputFile(inputFname=os.path.abspath(self.alnStatForFilterFname),\
														folderName=self.pegasusFolderName)
		else:
			alnStatForFilterF = None
		if self.keepSNPPosFname:
			keepSNPPosF = self.registerOneInputFile(inputFname=os.path.abspath(self.keepSNPPosFname),\
														folderName=self.pegasusFolderName)
		else:
			keepSNPPosF = None
		# 2012.5.1 filter only on the 1st vcf folder
		inputData = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler, \
									checkEmptyVCFByReading=self.checkEmptyVCFByReading,\
									pegasusFolderName="%s"%(self.pegasusFolderName), maxContigID=self.maxContigID, \
									minContigID=self.minContigID)
		counter =0
		for minMAC in self.minMACLs:
			for minMAF in self.minMAFLs:
				for maxSNPMissingRate in self.maxSNPMissingRateLs:
					for minDepthPerGenotype in self.minDepthPerGenotypeLs:
						folderSignature='minMAC%s_minMAF%s_maxSNPMissingRate%s_minDepthPerGenotype%s_depthFoldChange%s_keepSNPFile%s_onlyBiAllelic%s'%\
							(minMAC, minMAF, maxSNPMissingRate, minDepthPerGenotype, depthFoldChange, getattr(keepSNPPosF, 'name', None),\
							self.onlyKeepBiAllelicSNP)
						filterJobData = self.addJobsToFilterOneVCFDir(workflow, inputData=inputData, refFastaFList=refFastaFList, \
												alnStatForFilterF=alnStatForFilterF, keepSNPPosF=keepSNPPosF, \
												onlyKeepBiAllelicSNP=self.onlyKeepBiAllelicSNP,\
												minMAC=minMAC, minMAF=minMAF, maxSNPMissingRate=maxSNPMissingRate,\
												minDepthPerGenotype=minDepthPerGenotype, transferOutput=False, \
												outputDirPrefix="run%s_%s"%(counter, folderSignature))
						returnData = self.addStatCalculationJobs(workflow=workflow, inputData=filterJobData, refFastaFList=refFastaFList, \
												chr2size=chr2size, windowSize=self.windowSize, minChrLengthForPlot=self.minChrLengthForPlot, \
												minChrSize=self.minChrSize, LDWindowSize=self.LDWindowSize, transferOutput=True, \
												outputDirPrefix="run%s_%s"%(counter, folderSignature),\
												samplingRate=self.samplingRate, minSiteGap=30000)
						counter += 1
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)

if __name__ == '__main__':
	main_class = BootstrapFilterCalculateVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
