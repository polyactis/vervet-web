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
import pylab


class PlotTrioInconsistencyOverFrequency(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
						('minNoOfTotal', 1, int): [100, '', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
						('title', 1, ): [None, 't', 1, 'title for the figure.'],\
						('figureDPI', 1, int): [200, '', 1, 'dpi for the output figures (png)'],\
						('formatString', 1, ): ['-', '', 1, 'formatString passed to matplotlib plot'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	@classmethod
	def trioInconsistentRateFileWalker(cls, inputFname, processFunc=None, minNoOfTotal=100, run_type=1):
		"""
		2011-9-30
		"""
		try:
			reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
			inconsistent_rate_index = col_name2index.get("inconsistency")
			if run_type==1:
				index_of_x_data = col_name2index.get("stopFrequency")
			elif run_type==2:
				index_of_x_data = col_name2index.get("stop")
			else:
				sys.stderr.write("Unsupported run_type %s in trioInconsistentRateFileWalker().\n"%(run_type))
				sys.exit(3)
			index_of_no_of_total = col_name2index.get("no_of_total")
			inconsistent_rate_ls = []
			x_ls = []
			for row in reader:
				no_of_total = int(float(row[index_of_no_of_total]))
				if no_of_total<=minNoOfTotal:
					continue
				inconsistency = float(row[inconsistent_rate_index])
				inconsistent_rate_ls.append(inconsistency)
				x_data = float(row[index_of_x_data])
				x_ls.append(x_data)
			processFunc(x_ls, inconsistent_rate_ls)
			del reader
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
	
	def plotXY(self, x_ls, y_ls, ):
		"""
		2011-9-30
		"""
		pylab.plot(x_ls, y_ls, self.formatString)
		
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		pylab.clf()
		
		for inputFname in self.inputFnameLs:
			if os.path.isfile(inputFname):
				self.trioInconsistentRateFileWalker(inputFname, processFunc=self.plotXY, minNoOfTotal=self.minNoOfTotal,\
											run_type=1)
		
		if self.title is None:
			title = " %s refs"%(len(self.inputFnameLs))
		else:
			title = self.title
		
		pylab.title(title)
		pylab.xlabel("frequency")
		pylab.ylabel("inconsistent rate")
		
		pylab.savefig(self.outputFname, dpi=self.figureDPI)
		sys.stderr.write("Done.\n")


if __name__ == '__main__':
	main_class = PlotTrioInconsistencyOverFrequency
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
