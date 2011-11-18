#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011-9-29
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib
from PlotTrioInconsistencyOverFrequency import PlotTrioInconsistencyOverFrequency

class PlotTrioInconsistencyOverPosition(PlotTrioInconsistencyOverFrequency):
	__doc__ = __doc__
	option_default_dict = PlotTrioInconsistencyOverFrequency.option_default_dict

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		PlotTrioInconsistencyOverFrequency.__init__(self, inputFnameLs, **keywords)
		#Super(PlotTrioInconsistencyOverPosition, self).__init__(inputFnameLs, **keywords)
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
		"""
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		import pylab
		pylab.clf()
		
		if len(self.inputFnameLs)>1:	#multiple input file (multiple contigs). use lines.
			self.formatString = '-'
		else:	#only one input file (one contig). use dot.
			self.formatString = '.'
		for inputFname in self.inputFnameLs:
			self.trioInconsistentRateFileWalker(inputFname, processFunc=self.plotXY, minNoOfTotal=self.minNoOfTotal,\
											run_type=2)
		
		if self.title is None:
			title = " %s refs"%(len(self.inputFnameLs))
		else:
			title = self.title
		
		pylab.title(title)
		pylab.xlabel("position")
		pylab.ylabel("inconsistent rate")
		
		pylab.savefig(self.outputFname, dpi=self.figureDPI)
		sys.stderr.write("Done.\n")
		
		

if __name__ == '__main__':
	main_class = PlotTrioInconsistencyOverPosition
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
