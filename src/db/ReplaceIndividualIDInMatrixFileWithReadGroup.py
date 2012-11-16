#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/Contig315_StKitts_vs_Nevis.tsv --replaceColumnHeaderLs=IID1,IID2
		-s 1.0 -o /tmp/Contig315_StKitts_vs_Nevis.2D.tsv
	

Description:
	2012.10.15
		This program replaces Individual.id in some columns with the read group IDs supplied in self.readGroupFname.
			Samples that are not in self.readGroupFname will be discarded.
	If "-i ..." is given, it is regarded as one of the input files (plus the ones in trailing arguments). 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, SNP
import numpy, random
from pymodule.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from pymodule import MatrixFile
from SelectRowsWithinCoverageRange import SelectRowsWithinCoverageRange

class ReplaceIndividualIDInMatrixFileWithReadGroup(SelectRowsWithinCoverageRange):
	__doc__ = __doc__
	option_default_dict = AbstractMatrixFileWalker.option_default_dict.copy()
	option_default_dict.update(AbstractMatrixFileWalker.db_option_dict)
	option_default_dict[('minNoOfTotal', 1, int)][0] = 0
	option_default_dict.update({
						('readGroupFname', 1, ): [None, '', 1, 'This file contains a column with read groups in it. has a header.'],\
						('readGroupHeader', 0, ): [None, '', 1, 'The header of the read group in the readGroupFname. \n\
	If not given, column 0 is the read group header.'],\
						('replaceColumnHeaderLs', 1, ): [None, '', 1, 'a coma/dash-separated list of headers for the columns\n\
	whose individual IDs are to be replaced',],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		SelectRowsWithinCoverageRange.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		#self.connectDB() called within its __init__()
		self.replaceColumnHeaderLs = getListOutOfStr(self.replaceColumnHeaderLs, data_type=str)
	
	def setup(self, **keywords):
		"""
		"""
		AbstractMatrixFileWalker.setup(self, **keywords)
		
		#construct a individualCode2readGroup from readGroupFname
		self.invariantPData.individualCode2readGroup = {}
		reader = MatrixFile(inputFname=self.readGroupFname)
		reader.constructColName2IndexFromHeader()
		if self.readGroupHeader:
			readGroupIndex = reader.getColIndexGivenColHeader(self.readGroupHeader)
		else:
			readGroupIndex = 0
		for row in reader:
			readGroup = row[readGroupIndex]
			individualAlignment = self.db_vervet.parseAlignmentReadGroup(readGroup).individualAlignment
			if individualAlignment:
				individual_code = individualAlignment.individual_sequence.individual.code
				self.invariantPData.individualCode2readGroup[individual_code] = readGroup
		del reader
		return 1
	
	def processRow(self, row=None, pdata=None):
		"""
		2012.10.7
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		if col_name2index and self.replaceColumnHeaderLs:
			includeThisRow = True
			for columnHeader in self.replaceColumnHeaderLs:
				columnIndex = col_name2index.get(columnHeader, None)
				if columnIndex is not None:
					columnValue = row[columnIndex]
					readGroup = self.invariantPData.individualCode2readGroup.get(columnValue)
					if readGroup is None:	#some columns won't be translated. so skip the whole row.
						includeThisRow = False
						break
					row[columnIndex] = readGroup
			
			if includeThisRow:
				self.invariantPData.writer.writerow(row)
				returnValue = 1
		return returnValue

if __name__ == '__main__':
	main_class = ReplaceIndividualIDInMatrixFileWithReadGroup
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()