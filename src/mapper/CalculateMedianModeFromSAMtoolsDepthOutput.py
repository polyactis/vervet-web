#!/usr/bin/env python
"""
Examples:
	%s -i /tmp/depth.tsv.gz -a 100 -n 0 -o /tmp/depth.mode.tsv.gz -f 0.3
	
	%s -i /tmp/depth.tsv.gz -a 100 -n 0 -o /tmp/depth.mode.tsv -f 0.3

Description:
	2012.5.6
		Given a tsv-format matrix, this program calculates the mean/median/mode of one column's statistics.
		input & output file could be either gzipped or not.
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, getColName2IndexFromHeader 
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import random, numpy, scipy.stats


class CalculateMedianModeFromSAMtoolsDepthOutput(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('alignmentID', 1, ): [None, 'a', 1, 'ID of this alignment from which all the stats are extracted.'],\
						('fractionToSample', 1, float): [0.001, 'f', 1, 'fraction of rows to be included in calculating mean/median/mode'],\
						('whichColumn', 1, int ): [2, 'w', 1, 'which column of inputFname is the target stat.'],\
						('noOfLinesInHeader', 0, int): [1, 'n', 1, 'how many lines in the header'],\
						('maxNumberOfSamplings', 0, float): [1E7, 'm', 1, 'max number of samples to take into memory for median/mode/mean calculation to avoid memory blowup'],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		"""
		if self.whichColumnLs:
			self.whichColumnLs = getListOutOfStr(self.whichColumnLs)
		else:
			self.whichColumnLs = self.whichColumnLs
		"""
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		inf = utils.openGzipFile(self.inputFname, openMode='r')
		
		reader = csv.reader(inf, delimiter=figureOutDelimiter(inf))
		header = None
		for i in xrange(self.noOfLinesInHeader):
			if i==0:
				header = reader.next()	#first line is taken as header
			else:
				reader.next()
		if header is not None:
			colName2Index = getColName2IndexFromHeader(header)
		
		newHeader = ['alignmentID', 'total_base_count', 'sampled_base_count', 'meanDepth', 'medianDepth', 'modeDepth']
		inputStatLs = []
		
		writer = csv.writer(utils.openGzipFile(self.outputFname, openMode='w'), delimiter='\t')
		writer.writerow(newHeader)
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			if real_counter <= self.maxNumberOfSamplings:
				r = random.random()
				if r<=self.fractionToSample and real_counter<=self.maxNumberOfSamplings:
					inputStatLs.append(float(row[self.whichColumn]))
					real_counter += 1
		
		meanDepth = numpy.mean(inputStatLs)
		medianDepth = numpy.median(inputStatLs)
		modeDepth = scipy.stats.mode(inputStatLs)[0][0]
		outputRow = [self.alignmentID, counter, real_counter, meanDepth, medianDepth, modeDepth]
		writer.writerow(outputRow)
		del writer

if __name__ == '__main__':
	main_class = CalculateMedianModeFromSAMtoolsDepthOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()