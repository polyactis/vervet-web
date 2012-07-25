#!/usr/bin/env python
"""
Examples:
	%s
	
	dirPrefix=AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top4Contigs_single_sample_condor_20111105T0143/556_
	%s -i $dirPrefix\gatk/ -I $dirPrefix\samtools/ -l condorpool -j condorpool
		-o 4HighCovVRC_isq_15_18_vs_524_top804Contigs_gatk_vs_samtools_overlap_stat.xml -z uclaOffice -u yh -k genome
		-C 100
	
Description:
	2011-11-7 pegasus workflow that compares overlap between two vcf files (mapper/CheckTwoVCFOverlap.py),
			calculate mismatch rate, pi statistics based on the intersection
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, GenomeDB, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow

class CheckTwoVCFOverlapPipeline(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('vcf1Dir', 1, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						('vcf2Dir', 1, ): ['', 'I', 1, 'input folder that contains vcf or vcf.gz files', ],\
						})

	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
	
	def registerCustomExecutables(self, workflow):
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
		
		checkTwoVCFOverlap = Executable(namespace=namespace, name="CheckTwoVCFOverlap", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		checkTwoVCFOverlap.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "mapper/CheckTwoVCFOverlap.py"), site_handler))
		checkTwoVCFOverlap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(checkTwoVCFOverlap)
		workflow.checkTwoVCFOverlap = checkTwoVCFOverlap
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		"""
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, sequence_type_id=9)
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		self.registerJars(workflow)
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		statOutputDir = "statDir"
		statOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=statOutputDir)
		
		import re
		chr_pattern = re.compile(r'(\w+\d+).*')
		input_site_handler = self.input_site_handler
		
		overlapStatF = File('overlapStat.tsv')
		overlapStatMergeJob=self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
					outputF=overlapStatF, parentJobLs=[], \
					extraDependentInputLs=[], transferOutput=True, extraArguments=None)
		overlapPosFile = File("overlapPos.tsv")
		overlapPosMergeJob=self.addStatMergeJob(workflow, statMergeProgram=workflow.mergeSameHeaderTablesIntoOne, \
					outputF=overlapPosFile, parentJobLs=[], \
					extraDependentInputLs=[], transferOutput=True, extraArguments=None)
		
		overlapStatSumF = File('overlapStat.sum.tsv')
		overlapStatSumJob = self.addStatMergeJob(workflow, statMergeProgram=workflow.ReduceMatrixByChosenColumn, \
						outputF=overlapStatSumF, extraDependentInputLs=[], transferOutput=True, \
						extraArguments='-k 1000000 -v 1-1000')	#1000 is a random big upper limit. 100 monkeys => 330 columns
		self.addInputToStatMergeJob(workflow, statMergeJob=overlapStatSumJob, inputF=overlapStatF, \
							parentJobLs=[overlapStatMergeJob])
		
		counter = 0
		vcfFileID2path_1 = self.getVCFFileID2path(self.vcf1Dir)
		vcfFileID2path_2 = self.getVCFFileID2path(self.vcf2Dir)
		sharedVCFFileIDSet = set(vcfFileID2path_1.keys())&set(vcfFileID2path_2.keys())
		sys.stderr.write("%s shared vcf files.\n"%(len(sharedVCFFileIDSet)))
		
		for vcfFileID in sharedVCFFileIDSet:
			gatkVCFAbsPath = vcfFileID2path_1.get(vcfFileID)
			samtoolsVCFAbsPath = vcfFileID2path_2.get(vcfFileID)
			if not NextGenSeq.isVCFFileEmpty(gatkVCFAbsPath) and not NextGenSeq.isVCFFileEmpty(samtoolsVCFAbsPath, \
									checkContent=self.checkEmptyVCFByReading):	#make sure the samtools vcf is not empty
				gatkVCFFileBaseName = os.path.basename(gatkVCFAbsPath)
				chr = vcfFileID
				chr_size = chr2size.get(chr)
				if chr_size is None:
					sys.stderr.write("size for chr %s is unknown. set it to 1000.\n"%(chr))
					chr_size = 1000
				
				vcf1 = File(os.path.join('vcf1', gatkVCFFileBaseName))	#relative path
				vcf1.addPFN(PFN("file://" + gatkVCFAbsPath, input_site_handler))
				workflow.addFile(vcf1)
				
				vcf2 = File(os.path.join('vcf2', os.path.basename(samtoolsVCFAbsPath)))	#relative path
				vcf2.addPFN(PFN("file://" + samtoolsVCFAbsPath, input_site_handler))
				workflow.addFile(vcf2)
				
				
				outputFnamePrefix = os.path.join(statOutputDir, os.path.splitext(gatkVCFFileBaseName)[0])
				checkTwoVCFOverlapJob = self.addCheckTwoVCFOverlapJob(workflow, executable=workflow.checkTwoVCFOverlap, \
						vcf1=vcf1, vcf2=vcf2, chromosome=chr, chrLength=chr_size, \
						outputFnamePrefix=outputFnamePrefix, parentJobLs=[statOutputDirJob], \
						extraDependentInputLs=[], transferOutput=False, extraArguments=" -m %s "%(self.minDepth), job_max_memory=1000)
				
				self.addInputToStatMergeJob(workflow, statMergeJob=overlapStatMergeJob, inputF=checkTwoVCFOverlapJob.output, \
							parentJobLs=[checkTwoVCFOverlapJob], extraDependentInputLs=[])
				self.addInputToStatMergeJob(workflow, statMergeJob=overlapPosMergeJob, inputF=checkTwoVCFOverlapJob.overlapSitePosF, \
							parentJobLs=[checkTwoVCFOverlapJob], extraDependentInputLs=[])
				counter += 1
		sys.stderr.write("%s jobs.\n"%(counter+1))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = CheckTwoVCFOverlapPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
