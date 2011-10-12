#!/usr/bin/env python
"""
Examples:
	%s -o /tmp/test.png -t anything /tmp/Contig119.trio.75_17_86.het.homo.inconsistency.summary.tsv 
		/tmp/Contig119.trio.75_17_86.het.homo.inconsistency.summary.tsv ...
	
	%s 
	

Description:
	2011-9-29
		a list of input files are appended in the end
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
import csv, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib


class PlotTrioInconsistencySummaryHist(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
						('title', 1, ): [None, 't', 1, 'title for the figure.'],\
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
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		inconsistent_rate_ls = []
		for inputFname in self.inputFnameLs:
			reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
			inconsistent_rate_index = col_name2index.get("inconsistency")
			for row in reader:
				inconsistency = float(row[inconsistent_rate_index])
				inconsistent_rate_ls.append(inconsistency)
			del reader
		
		
		if self.title is None:
			title = "histogram of inconsistent rate from %s refs"%(len(inconsistent_rate_ls))
		else:
			title = self.title
		if len(inconsistent_rate_ls)>10:
			medianInconsistentRate = numpy.median(inconsistent_rate_ls)
			title += " median %.4f"%(medianInconsistentRate)
		yh_matplotlib.drawHist(inconsistent_rate_ls, title=title, \
									xlabel_1D="Inconsistent Rate", xticks=None, outputFname=self.outputFname, min_no_of_data_points=20, needLog=False, \
									dpi=200)
		

if __name__ == '__main__':
	main_class = PlotTrioInconsistencySummaryHist
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
