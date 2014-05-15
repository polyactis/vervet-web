#!/usr/bin/env python
"""
Examples:
	# scaffold (120) is used as reference in alignment. get genome sequence id from 1 to 8.
	%s -o workflow.xml -a 120 --ind_seq_id_ls 1-8 -y2 -s2
	
	# 1Mb-BAC (9) is used as reference.
	%s -o workflow.xml -a 9 --ind_seq_id_ls 1-4 -y2 -s2
	
	# 8 genomes versus top 156 contigs
	%s -o workflow_8GenomeVsTop156Contigs.xml -u yh  -a 128 --ind_seq_id_ls 1-8 -N 156 -y2 -s2
	
	
Description:
	2013.2.26
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from Pegasus.DAX3 import *
from pymodule import ProcessOptions, PassingData, AbstractWorkflow, utils
from pymodule.pegasus import yh_pegasus
from pymodule.popgen.PopGenSimulationWorkflow import PopGenSimulationWorkflow

from vervet.src import AbstractVervetWorkflow, VervetDB, AbstractVervetAlignmentAndVCFWorkflow

parentClass = PopGenSimulationWorkflow
class SimulatePedigreeSequenceWorkflow(PopGenSimulationWorkflow):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('inputDir', 0, ))
	
	
	commonCallPipelineOptionDict= {
						#('seqCoverageFname', 0, ): ['', 'q', 1, 'The sequence coverage file. tab/comma-delimited: individual_sequence.id coverage'],\
						("pop_gen_simulation_id", 1, int): [None, '', 1, 'id from popGenSimulation table, use as founders'],\
						
						
						("sampleSize", 1, int ): [50, '', 1, 'how many individuals from the pop-gen simulation to constitute the sample (final output)'],\
						("noOfLociToSimulate", 1, int): [1, '', 1, 'how many loci to simulate'],\
						("simulateLocusLengthList", 1, ): ['1000000', '', 1, 'a coma-separated list of locus length to simulate'],\
						("sfs_code_path", 1, ): ['%s/script/repository/sfscode/bin/sfs_code', '', 1, 'path to the sfs_code pop-gen forward simulator'],\
						("slim_path", 1, ): ['%s/script/repository/slim/slim', '', 1, 'path to the slim pop-gen forward simulator'],\
						
						("sequencingErrorRate", 1, float ): [0.01, '', 1, 'sequencing error rate for each raw read'],\
						("pairedEndInsertSizeMedian", 1, int): [400, '', 1, 'median insert sizes of paired end reads'],\
						("pedigreeSequencingCoverageStrategy", 0, int): [0, '', 1, 'what coverage strategy to use to sequence pedigree members. \
							0: existing pedigree coverage'],\
						
						
						
						("simulateIndels", 0, int): [0, '', 1, 'toggle to simulate insertion/deletion polymorphism as well in pop-gen simulation'],\
						("noOfReplicates", 1, int): [5, '', 1, 'how many replicates for this simulation setup'],\
						
						}
	
	commonCallPipelineOptionDict.update(parentClass.partitionWorkflowOptionDict.copy())
	commonCallPipelineOptionDict.update(parentClass.commonAlignmentWorkflowOptionDict.copy())
	option_default_dict.update(commonCallPipelineOptionDict.copy())
	option_default_dict.update({
						("run_type", 1, int): [1, 'n', 1, '1: multi-sample calling, 2: single-sample one by one'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.9.17 call parentClass.__init__() directly
		2011-7-11
		"""
		self.pathToInsertHomePathList.extend(['sfs_code_path'])
		
		parentClass.__init__(self, **keywords)
		
		if hasattr(self, 'simulateLocusLengthList', None):
			self.simulateLocusLengthList = utils.getListOutOfStr(self.simulateLocusLengthList, data_type=int)
		else:
			self.simulateLocusLengthList = []
	
	def addAllJobs(self, workflow=None, db_250k=None, association_result_ls=None, \
				data_dir=None, min_MAF=None, \
				neighbor_distance=None, max_neighbor_distance=None, \
				min_score_ls=None, min_overlap_ratio_ls=None, ground_score=None,\
				peakPadding=None, tax_id=None, \
				outputDirPrefix="", transferOutput=True, job_max_memory=2000, **keywords):
		"""
		2013.2.27
			run ms
			estimate parameters from ms
			forward simulator with estimated ms-parameters or take the output of ms as input
			
			
		"""
		if workflow is None:
			workflow = self
		
		sys.stderr.write("Adding jobs for pop-gen & pedigree sequence simulation #jobs=%s... \n"%\
							(self.no_of_jobs))
		
		returnData = PassingData()
		returnData.jobDataLs = []
		
		passingData = PassingData(fileBasenamePrefix=None, \
					outputDirPrefix=outputDirPrefix, \
					jobData=None,\
					preReduceReturnData=None,\
					association_group_key2orderIndex = {},\
					association_group_key2resultList = {},\
					association_group_key2reduceAssociationPeakJobMatrix = {},\
					association_group_key2countAssociationLocusJobList = {},\
					resultID2defineLandscapeJobData = {},
					)
		
		preReduceReturnData = self.preReduce(workflow=workflow, outputDirPrefix=outputDirPrefix, \
									passingData=passingData, transferOutput=False,\
									**keywords)
		
		mapDirJob = preReduceReturnData.mapDirJob
		plotOutputDirJob = preReduceReturnData.plotOutputDirJob
		countAssociationLocusOutputDirJob = preReduceReturnData.countAssociationLocusOutputDirJob
		reduceOutputDirJob = preReduceReturnData.reduceOutputDirJob
		
		passingData.preReduceReturnData = preReduceReturnData
		
		#add output pedigree job
		
		for i in xrange(self.noOfReplicates):
			popGenSimulationFolderJob = self.addMkDirJob(outputDir=os.path.join(mapDirJob.output, 'popGenSim%s'%(i)), \
														parentJobLs=[mapDirJob])
			popSimulationJob = self.addPopGenSimulationJob()
			
			#. add pop-gen output 2 polymorphism-table file
			
			#. add pedigree haplotype simulator
			
			#. add polymorphismTableFile2fasta convertor
			
			#. add short-read simulator
			
		
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2013.2.26
		"""
		parentClass.registerCustomExecutables(self, workflow=workflow)
		
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2013.3.8
		self.addOneExecutableFromPathAndAssignProperClusterSize(\
											path=os.path.join(self.vervetSrcPath, 'pedigree/SimulatePedigreeHaplotype.py'), \
												name="SimulatePedigreeHaplotype", clusterSizeMultipler=0.2)
		
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = self.db_vervet
		
		if not self.data_dir:
			self.data_dir = db_vervet.data_dir
		
		if not self.local_data_dir:
			self.local_data_dir = db_vervet.data_dir
		
		workflow = self.initiateWorkflow()
		
		alignmentLs = db_vervet.getAlignments(self.ref_ind_seq_id, ind_seq_id_ls=self.ind_seq_id_ls, ind_aln_id_ls=self.ind_aln_id_ls,\
										alignment_method_id=self.alignment_method_id, data_dir=self.local_data_dir,\
										individual_sequence_file_raw_id_type=self.individual_sequence_file_raw_id_type)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=self.sequence_filtered, \
									individual_site_id_set=set(self.site_id_ls),\
									mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=set(self.country_id_ls), tax_id_set=set(self.tax_id_ls),\
									excludeContaminant=self.excludeContaminant)
		
		cumulativeMedianDepth = db_vervet.getCumulativeAlignmentMedianDepth(alignmentLs=alignmentLs, \
										defaultSampleAlignmentDepth=self.defaultSampleAlignmentDepth)
		
		refFastaFList = self.getReferenceSequence(workflow=self)
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		alignmentDataLs = self.registerAlignmentAndItsIndexFile(workflow, alignmentLs, data_dir=self.data_dir)
		
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		



if __name__ == '__main__':
	main_class = SimulatePedigreeSequenceWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
