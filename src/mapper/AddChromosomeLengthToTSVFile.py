#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s 
	
Description:
	2011-11-12

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

from pymodule import ProcessOptions, figureOutDelimiter
import csv

class AddChromosomeLengthToTSVFile(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'input tsv file with header'],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						("chromosome", 0, ): [None, 'c', 1, 'chromosome name for this reference (not important.'],\
						("chrLength", 1, int): [0, 'l', 1, 'length of the reference for this inputFname.'],\
						('divideByLength', 0, int):[0, 'v', 0, 'toggle to divide a bunch stats by chrLength'],\
						('divideStartingColumn', 1, int):[2, 's', 1, 'which column to start divide its value by chrLength'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		try:
			delimiter = figureOutDelimiter(self.inputFname)
			if not delimiter:
				delimiter='\t'
			reader = csv.reader(open(self.inputFname), delimiter=delimiter)
			writer = csv.writer(open(self.outputFname, 'w'), delimiter=delimiter)
			lengthName = 'chrLength'
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			print sys.exc_info()
			sys.exit(0)
		try:
			header = reader.next()
			originalHeaderLength = len(header)
			header.append(lengthName)
			if self.divideByLength:
				i = self.divideStartingColumn
				while (i<originalHeaderLength):
					statColumnHeader = header[i]
					header.append("%s_div_by_%s"%(statColumnHeader, lengthName))
					i += 1;
			writer.writerow(header)
			for row in reader:
				new_data_row = row[:]
				new_data_row.append(self.chrLength)
				if self.divideByLength:
					i = self.divideStartingColumn
					while (i<originalHeaderLength):
						if self.chrLength>0:
							try:
								statData = float(row[i])
								rate = statData/self.chrLength
							except:
								rate = -1
						else:
							rate = -1
						new_data_row.append(rate)
						i += 1;
				writer.writerow(new_data_row)
			del reader
			del writer
		except:	#in case something wrong (i.e. file is empty)
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			print sys.exc_info()
			sys.exit(0)
		

if __name__ == '__main__':
	main_class = AddChromosomeLengthToTSVFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
