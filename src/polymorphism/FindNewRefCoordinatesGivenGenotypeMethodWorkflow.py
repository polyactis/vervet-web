#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	%s -a 524 -I Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5/
		-C 1 -H -m1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/ -u yh -z localhost
		-o workflow/PairwiseDistance/PairwiseDistance_Combine_FilterVCF_VRC_SK_Nevis_FilteredSeq_top1000Contigs_and_VariantsOf36RNA-SeqMonkeysFromNam_minDepth5.xml
		#-U 0 -Z 3000

Description:
	2012-10-03 unfinished child of pymodule/polymorphism/FindNewRefCoordinatesGivenVCFFolderWorkflow()
    #. mapEachVCF:
    		split VCF into several small ones, N-SNP each
      #. mapEachInterval
         #. extract flanking sequences from the input VCF ():
         	ExtractFlankingSequenceForVCFLoci.py -i vcf -a reference fasta -o output.fasta -f flankingLength (=24)
         #. blast them
         #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
             #. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
         #. output a mapping file between old SNP and new SNP coordinates.
         		#. reduce this thing by combining everything
         #. make a new VCF file based on the input split VCF file
            (replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)
	(#. merge same contig small VCF files into one big file (ligateVCF job))
    #. reduce()
         #. merge all small VCF files into one big one
         #. sort it by chromosome position
         #. split the big one into individual chromosomes
         #. add the contig-only file into db, with new reference and new genotype method ID (parent=genotype-method-id) (AddVCF2DB.py)
         
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from vervet.src import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from vervet.src.AbstractVervetWorkflow import AbstractVervetWorkflow
from pymodule.polymorphism.FindNewRefCoordinatesGivenVCFFolderWorkflow import FindNewRefCoordinatesGivenVCFFolderWorkflow

class FindNewRefCoordinatesGivenGenotypeMethodWorkflow(FindNewRefCoordinatesGivenVCFFolderWorkflow):
	__doc__ = __doc__
	option_default_dict = FindNewRefCoordinatesGivenVCFFolderWorkflow.option_default_dict.copy()
	option_default_dict.update(AbstractVervetWorkflow.db_option_dict.copy())
	option_default_dict.update({
						})
	
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		FindNewRefCoordinatesGivenVCFFolderWorkflow.__init__(self, **keywords)

	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		if workflow==None:
			workflow=self
		
		FindNewRefCoordinatesGivenVCFFolderWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)		
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	def mapEachInterval(self, workflow=None, \
					VCFJobData=None, passingData=None, transferOutput=False, **keywords):
		"""
		2012.10.3
			#. extract flanking sequences from the input VCF (ref sequence file => contig ref sequence)
			#. blast them
			#. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
				#. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
			#. output a mapping file between old SNP and new SNP coordinates.
				#. reduce this thing by combining everything
			#. make a new VCF file based on the input split VCF file
				(replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)

		"""
		# a flanking sequence extraction job
		
		# a blast job
		# a FindSNPPositionOnNewRefFromFlankingBlastOutput job
		
	def reduceEachChromosome(self, workflow=None, chromosome=None, passingData=None, transferOutput=True, \
						**keywords):
		"""
		2012.10.3
			#. merge all small VCF files
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData

	def reduceEachVCF(self, workflow=None, chromosome=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.10.3
			#. merge all small VCF files
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData
	
	def reduce(self, workflow=None, passingData=None, transferOutput=True, **keywords):
		"""
		2012.10.3
			#. merge all small VCF files into one big one
			#. sort it by chromosome position
			#. split the big one into individual chromosomes
			#. add the contig-only file into db, with new reference and new genotype method ID (parent=genotype-method-id) (AddVCF2DB.py)
		
		"""
		returnData = PassingData(no_of_jobs = 0)
		returnData.jobDataLs = []
		return returnData

if __name__ == '__main__':
	main_class = FindNewRefCoordinatesGivenGenotypeMethodWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()