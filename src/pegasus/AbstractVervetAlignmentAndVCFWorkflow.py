#!/usr/bin/env python
"""
Examples:
	
Description:
	2012.9.17 abstract workflow that takes both alignment and VCF as input. re-factored out of AlignmentToCallPipeline.py 
"""
import sys, os, math
__doc__ = __doc__%()


#bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
import subprocess, cStringIO
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus,\
	AbstractAlignmentAndVCFWorkflow
from vervet.src import VervetDB
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow
from vervet.src.pegasus.AbstractVervetAlignmentWorkflow import AbstractVervetAlignmentWorkflow


class AbstractVervetAlignmentAndVCFWorkflow(AbstractVervetWorkflow, AbstractAlignmentAndVCFWorkflow, AbstractVervetAlignmentWorkflow):
	__doc__ = __doc__
	option_default_dict = copy.deepcopy(AbstractVervetWorkflow.option_default_dict)
	option_default_dict.update(copy.deepcopy(AbstractAlignmentAndVCFWorkflow.option_default_dict))
	#change the short option to nothing to avoid conflict
	#option_default_dict[('inputDir', 0, )] = ['', 'L', 1, 'input folder that contains vcf or vcf.gz files', ]
	
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		#if self.option_default_dict.get('ligateVcfPerlPath'):
		#	self.pathToInsertHomePathList.append('ligateVcfPerlPath')
		AbstractVervetWorkflow.__init__(self, **keywords)
		
	
	def registerCustomExecutables(self, workflow=None):
		
		"""
		2011-11-28
		"""
		AbstractVervetWorkflow.registerCustomExecutables(self, workflow=workflow)
		
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
	
	def run(self):
		"""
		2013.2.14 overwrite the default : AbstractVervetWorkflow.run()
		"""
		AbstractAlignmentAndVCFWorkflow.run(self)

if __name__ == '__main__':
	main_class = AbstractVervetAlignmentAndVCFWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()