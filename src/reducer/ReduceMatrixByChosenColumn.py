#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s 
	
Description:
	2011-11-12
		This program sums individual columns from all input files based on keys from the keyColumn.
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

from pymodule import ProcessOptions, figureOutDelimiter, utils
import csv

class ReduceMatrixByChosenColumn(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						("keyColumn", 1, int): [0, 'k', 1, 'data is keyed by this column'],\
						('valueColumnLs', 1, ):["1", 'v', 1, 'comma/tab-separated list, specifying columns from which to aggregate total value by key'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
		if self.valueColumnLs:
			self.valueColumnLs = utils.getListOutOfStr(self.valueColumnLs, data_type=int)
		else:
			self.valueColumnLs = []
		
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
			
		newHeader = None
		key2dataLs = {}	#key is the keyColumn, dataLs corresponds to the sum of each column from valueColumnLs 
		delimiter = None
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				continue
			try:
				delimiter = figureOutDelimiter(inputFname)
				if not delimiter:
					delimiter='\t'
				reader = csv.reader(open(inputFname), delimiter=delimiter)
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			
			try:
				header = reader.next()
				originalHeaderLength = len(header)
				if newHeader is None:
					newHeader = []
					newHeader.append(header[self.keyColumn])
					for columnIndex in self.valueColumnLs:
						if columnIndex<len(header):
							newHeader.append(header[columnIndex])
			except:	#in case something wrong (i.e. file is empty)
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			
			for row in reader:
				try:
					key = row[self.keyColumn]
					if key not in key2dataLs:
						key2dataLs[key] = [0]*(len(newHeader)-1)	#should not be longer than the newHeader
					for i in xrange(len(self.valueColumnLs)):
						columnIndex = self.valueColumnLs[i]
						if columnIndex<len(row):
							value = float(row[columnIndex])
							key2dataLs[key][i] = key2dataLs[key][i] + value
				except:	#in case something wrong (i.e. file is empty)
					sys.stderr.write('Ignore this row: %s.\n'%repr(row))
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			del reader
			
		
		if key2dataLs and delimiter and newHeader:
			writer = csv.writer(open(self.outputFname, 'w'), delimiter=delimiter)
			writer.writerow(newHeader)
			keyLs = key2dataLs.keys()
			keyLs.sort()
			for key in keyLs:
				dataLs = key2dataLs.get(key)
				writer.writerow([key] + dataLs)
			del writer

if __name__ == '__main__':
	main_class = ReduceMatrixByChosenColumn
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
