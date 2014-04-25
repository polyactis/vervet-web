#!/usr/bin/env python
"""
Examples:
	
Description:
	2013.1.25 abstract workflow that takes alignment as input.
"""
import sys, os, math
__doc__ = __doc__%()


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, AbstractAlignmentWorkflow
from vervet.src import VervetDB
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow


class AbstractVervetAlignmentWorkflow(AbstractAlignmentWorkflow, AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractAlignmentWorkflow.option_default_dict)
	#change the short option to nothing to avoid conflict
	#option_default_dict[('inputDir', 0, )] = ['', 'L', 1, 'input folder that contains vcf or vcf.gz files', ]

	#option_default_dict.pop(('inputDir', 0, ))
	option_default_dict.update(AbstractAlignmentWorkflow.commonAlignmentWorkflowOptionDict.copy())
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		AbstractAlignmentWorkflow.__init__(self, **keywords)
		self.db = self.db_vervet	#2013.1.25 main db
		
		#2013.2.4 this will screw up the arguments like self.ind_seq_id_ls, that require extra processing in AbstractAlignmentWorkflow.__init__
		#AbstractVervetWorkflow.__init__(self, **keywords)
	
	getReferenceSequence=AbstractVervetWorkflow.getReferenceSequence
	
	connectDB =	AbstractVervetWorkflow.connectDB
	
		
	def registerCustomExecutables(self, workflow=None):
		
		"""
		"""
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow=workflow)
		AbstractAlignmentWorkflow.registerCustomExecutables(self, workflow=workflow)
		
		if workflow is None:
			workflow = self
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		executableClusterSizeMultiplierList = []
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
	
	
	#run=AbstractAlignmentWorkflow.run

if __name__ == '__main__':
	main_class = AbstractVervetAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
