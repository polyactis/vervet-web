#!/usr/bin/env python
"""
Examples:
	# 8 genomes versus top 156 contigs
	%s -o workflow_8GenomeVsTop156Contigs.xml -u yh  -a 128 -i 1-8 -N 156 -y2 -s2
	
	# 2011-7-21 use GATK + coverage filter, top 2 contigs
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8 -N 2 -y1  -s2
	
	# 2011-7-21 use GATK + coverage filter on hoffman2 and site_handler, top 5 contigs
	%s -o workflow_8GenomeVsTop2Contig_GATK.xml -u yh  -a 120 -i 1-8
		-N 5 -y1 -l hoffman2 -e /u/home/eeskin/polyacti -t /u/home/eeskin/polyacti/NetworkData/vervet/db  -s2
	
	#2011-8-31 work 10 VWP and one VRC ref monkeys, variants only
	%s -a 9 -I 495,498-507 -u yh  
		-l condorpool -y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC.xml -s2 -q /tmp/all_isq_coverage.tsv
	
	%s -a 120 -I 34,38 -u yh -l hoffman2
	-y1 -o AlignmentToCallPipeline_10VWP_VRC_ref_vs_1Mb_BAC_hoffman2.xml  -s2 -e /u/home/eeskin/polyacti
	-t /u/home/eeskin/polyacti/NetworkData/vervet/db -N 4 -q /tmp/all_isq_coverage.tsv
	
	#2011-9-14 top 25 contigs, variants only, run on uschpc cluster
	%s -I 559-656 -j uschpc -l uschpc -u yh -a 524 -s 2 -e /home/cmb-03/mn/yuhuang -t /home/cmb-03/mn/yuhuang/NetworkData/vervet/db/
		-z 10.8.0.10 -o ./AlignmentToCallPipeline_559_656_vs_524_top_25Contigs_uschpc.xml
		-D ~/mnt/hpc-cmb_home/NetworkData/vervet/db/ -N25
	
	# 2011-11-4 run GATK/samtools on single-sample at a time, for 4 high-coverage VRC monkeys, top 804 contigs
	%s -a 524 -i 15-18 -u yh -l condorpool -j condorpool -s 2 -N 804
		-o AlignmentToCallPipeline_4HighCovVRC_isq_15_18_vs_524_top804Contigs_single_sample_condor.xml -z uclaOffice -n 2
	
	# 2011-11-22 prepare a workflow to run on condorpool
	%s -a 524 -u yh -j condorpool -l condorpool -N 3 -z uclaOffice -o AlignmentToCallLowPass_top3Contigs.xml
	
	#2011-11-22 prepare a low-coverage (coverage<=15) workflow to run on hoffman2's condorpool
	%s -a 524 -u yh -j hcondor -l hcondor -N 7559 -z localhost -o AlignmentToCallLowPass_top7559Contigs.xml
		-e /u/home/eeskin/polyacti/ -t /u/home/eeskin/polyacti/NetworkData/vervet/db/ -D /u/home/eeskin/polyacti/NetworkData/vervet/db/
		-J ~/bin/jdk/bin/java -w 15 -S 447 -C 30 -O 10
	
Description:
	2011-11-22
		a derivative of AlignmentToCallPipeline.py for monkeys given a max-coverage and site-id filter.
"""
import sys, os, math
from AlignmentToCallPipeline import AlignmentToCallPipeline
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0],\
				sys.argv[0], sys.argv[0], sys.argv[0])

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *

class AlignmentToCallPipelineForLowPass(AlignmentToCallPipeline):
	__doc__ = __doc__
	option_default_dict = AlignmentToCallPipeline.option_default_dict.copy()
	option_default_dict.update({
			("max_coverage", 1, int): [15, 'w', 1, 'maximum coverage (individual_sequence.coverage) for an sequence run to be considered as low-pass'],\
			("site_id", 1, int): [447, 'S', 1, 'individuals must come from this site.'],\
			})
	def __init__(self,  **keywords):
		"""
		2011-11-22
		"""
		AlignmentToCallPipeline.__init__(self, **keywords)
	
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
		
		refName2size = self.getTopNumberOfContigs(self.topNumberOfContigs, contigMinRankBySize=self.contigMinRankBySize)
		#refName2size = set(['Contig149'])	#temporary when testing Contig149
		#refName2size = set(['1MbBAC'])	#temporary when testing the 1Mb-BAC (formerly vervet_path2)
		refNameLs = refName2size.keys()
		
		alignmentLs = self.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										aln_method_id=2, dataDir=self.localDataDir)
		
		#site id 447 is the VRC site
		alignmentLs = self.filterAlignments(alignmentLs, max_coverage=self.max_coverage, individual_site_id=self.site_id)
		
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = self.initiateWorkflow(workflowName)
		
		refSequence = VervetDB.IndividualSequence.get(self.ref_ind_seq_id)
		refFastaFname = os.path.join(self.dataDir, refSequence.path)
		refFastaFList = yh_pegasus.registerRefFastaFile(workflow, refFastaFname, registerAffiliateFiles=True, \
							input_site_handler=self.input_site_handler,\
							checkAffiliateFileExistence=True)
		refFastaF = refFastaFList[0]
		
		self.registerJars(workflow)
		
		self.registerCommonExecutables(workflow)
		self.registerCustomExecutables(workflow)
		
		
		if self.run_type==1:	#multi-sample calling
			dirPrefix = ""
			# Add a mkdir job for the call directory.
			callOutputDir = "%scall"%(dirPrefix)
			callOutputDirJob = self.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=callOutputDir, \
											namespace=workflow.namespace, version=workflow.version)
			gatkDir = "%sgatk"%(dirPrefix)
			gatkDirJob = self.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=gatkDir, \
										namespace=workflow.namespace, version=workflow.version)
			samtoolsDir = "%ssamtools"%(dirPrefix)
			samtoolsDirJob = self.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=samtoolsDir, \
											namespace=workflow.namespace, version=workflow.version)
			unionDir = "%sgatk_samtools_union"%(dirPrefix)
			unionDirJob = self.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=unionDir, \
									namespace=workflow.namespace, version=workflow.version)
			intersectionDir = "%sgatk_samtools_intersection"%(dirPrefix)
			intersectionDirJob = self.addMkDirJob(workflow, mkdir=workflow.mkdirWrap, outputDir=intersectionDir, \
												namespace=workflow.namespace, version=workflow.version)
			
			
			self.addGenotypeCallJobs(workflow, alignmentLs, refName2size, samtools=workflow.samtools, \
					genotyperJava=workflow.genotyperJava, genomeAnalysisTKJar=workflow.genomeAnalysisTKJar, \
					addOrReplaceReadGroupsJava=workflow.addOrReplaceReadGroupsJava, addOrReplaceReadGroupsJar=workflow.addOrReplaceReadGroupsJar, \
					createSequenceDictionaryJava=workflow.createSequenceDictionaryJava, createSequenceDictionaryJar=workflow.createSequenceDictionaryJar, \
					mergeSamFilesJar=workflow.mergeSamFilesJar, \
					BuildBamIndexFilesJava=workflow.BuildBamIndexFilesJava, BuildBamIndexFilesJar=workflow.BuildBamIndexFilesJar, \
					mv=workflow.mv, CallVariantBySamtools=workflow.CallVariantBySamtools,\
					bgzip_tabix=workflow.bgzip_tabix, vcf_convert=workflow.vcf_convert, vcf_isec=workflow.vcf_isec, vcf_concat=workflow.vcf_concat, \
					concatGATK=workflow.concatGATK, concatSamtools=workflow.concatSamtools,\
					genotypeCallByCoverage=workflow.genotypeCallByCoverage, refFastaFList=refFastaFList, bamListF=None, \
					callOutputDirJob =callOutputDirJob, gatkDirJob=gatkDirJob, samtoolsDirJob=samtoolsDirJob, unionDirJob=unionDirJob, intersectionDirJob=intersectionDirJob,\
					namespace=workflow.namespace, version=workflow.version, site_handler=self.site_handler, input_site_handler=self.input_site_handler,\
					seqCoverageF=None, \
					needFastaIndexJob=self.needFastaIndexJob, needFastaDictJob=self.needFastaDictJob, \
					chunkSize=2000000, site_type=self.site_type, dataDir=self.dataDir)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
	
if __name__ == '__main__':
	main_class = AlignmentToCallPipelineForLowPass
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
