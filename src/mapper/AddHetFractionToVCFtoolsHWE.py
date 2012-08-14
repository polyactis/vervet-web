#!/usr/bin/env python
"""
Examples:
	%s  -o /tmp/ccc.tsv.gz -i /tmp/call_1.tsv -c Contig0 -l 14324143 -v -s 3
	
	%s 
	
Description:
	2012.8.6 it does same work as AddChromosomeLengthToTSVFile.py + calculate/add a column of site heterozygosity.

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	#sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, figureOutDelimiter
import csv
from AddChromosomeLengthToTSVFile import AddChromosomeLengthToTSVFile
class AddHetFractionToVCFtoolsHWE(AddChromosomeLengthToTSVFile):
	__doc__ = __doc__
	option_default_dict = AddChromosomeLengthToTSVFile.option_default_dict
	option_default_dict.update({
							("homoHetVectorHeader", 0, ): ["OBS(HOM1/HET/HOM2)", '', 1, 'the header for the homo/het vector column from input, i.e. 13/3/0'],\
							("hetFractionHeader", 0, ): ["hetFraction", '', 1, 'the header for the new heterozygote fraction column,'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		2011-7-12
		"""
		AddChromosomeLengthToTSVFile.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def processHeader(self, reader=None, extendHeader=None, chrLengthHeader = 'chrLength'):
		"""
		2012.8.7
		"""
		header = AddChromosomeLengthToTSVFile.processHeader(self, reader=reader, extendHeader=extendHeader, chrLengthHeader=chrLengthHeader)
		header.append(self.hetFractionHeader)
		return header
	
	def processRow(self, row=None):
		"""
		2012.8.7
			OBS(HOM1/HET/HOM2) (i.e. 13/3/0 )
		"""
		new_data_row = AddChromosomeLengthToTSVFile.processRow(self, row=row)
		col_index = self.col_name2index.get(self.homoHetVectorHeader, None)
		yValue = new_data_row[col_index]
		vector = yValue.split('/')
		vector = map(int, vector)
		noOfHomo1, noOfHet, noOfHomo2 = vector
		noOfTotal = sum(vector)
		if noOfTotal>0:
			yValue = float(noOfHet/float(noOfTotal))
		else:
			yValue = -1
		new_data_row.append(yValue)
		return new_data_row
	
	
if __name__ == '__main__':
	main_class = AddHetFractionToVCFtoolsHWE
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()