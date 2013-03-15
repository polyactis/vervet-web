#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i 

Description:
	2012.3.14
		count the number of reads/bases of fastq or fasta input file.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class CountFastqReadBaseCount(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.update({
							('isq_id', 0, int): [0, 'q', 1, 'IndividualSequence.id'],\
							('isqf_id', 0, int): [0, 'f', 1, 'IndividualSequenceFile.id if applicable'],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		self.inputFnameLs = inputFnameLs	#2012.3.19 not used
	
	@classmethod
	def getReadBaseCount(cls, inputFname, ignore_set = set(['>', '+', '@']), onlyForEmptyCheck=False):
		"""
		2012.3.19
			inputFname could be fastq or fasta
		"""
		inf = utils.openGzipFile(inputFname, openMode='r')
		read_count = 0
		base_count = 0
		
		for line in inf:
			if line[0] in ignore_set:
				if line[0]=='+':	#skip the quality-score line right after this "+" line
					inf.next()
				continue
			read_count += 1
			base_count += len(line.strip())
			if onlyForEmptyCheck:	#2012.3.19 one read is enough.
				break
		
		del inf
		return PassingData(read_count=read_count, base_count=base_count)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		baseCountData = self.getReadBaseCount(self.inputFname)
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['isq_id', 'isqf_id', 'read_count', 'base_count']
		writer.writerow(header)
		data_row = [self.isq_id, self.isqf_id, baseCountData.read_count, baseCountData.base_count]
		writer.writerow(data_row)
		del writer
	
if __name__ == '__main__':
	main_class = CountFastqReadBaseCount
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()