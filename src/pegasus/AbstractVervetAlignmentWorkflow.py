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
from Pegasus.DAX3 import *
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, AbstractAlignmentWorkflow
from vervet.src import VervetDB
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow


class AbstractVervetAlignmentWorkflow(AbstractVervetWorkflow, AbstractAlignmentWorkflow):
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

if __name__ == '__main__':
	main_class = AbstractVervetAlignmentWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
