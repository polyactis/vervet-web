#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011-11-4
		input should be the *.sample_interval_summary output file of GATK's DepthOfCoverageWalker.
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


class ReduceDepthOfCoverage(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-11-4
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		writer.writerow(['#sampleID', 'chromosome', 'meanDepth', 'medianDepth'])
		for inputFname in self.inputFnameLs:
			inputFile = utils.openGzipFile(inputFname)
			delimiter = figureOutDelimiter(inputFile)
			reader = csv.reader(inputFile, delimiter=delimiter)
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header)
			
			intervalIDIndex = col_name2index.get("Target")
			#only the first read group among the output (so don't run the DepthOfCoverageWalker over multi-read-group bam files
			avgCoverageIndex = 4
			sampleID = header[avgCoverageIndex][:-9]	#this column header is like $sampleID_mean_cvg. so get rid of _mean_cvg
			medianCoverageIndex = 6
			
			for row in reader:
				intervalID = row[intervalIDIndex]
				writer.writerow([sampleID, intervalID, row[avgCoverageIndex], row[medianCoverageIndex]])
		del writer
		sys.stderr.write("Done.\n")


if __name__ == '__main__':
	main_class = ReduceDepthOfCoverage
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
