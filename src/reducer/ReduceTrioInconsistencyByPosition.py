#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011.12.16
		input is the depthOutputWriter (one output) of CalculateTrioInconsistency.py
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


class ReduceTrioInconsistencyByPosition(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-12-16
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#['trio_set', 'chromosome', 'pos', 'depthOfFather','depthOfMother', 'depthOfChild', 'isInconsistent']
		
		chr_pos2inconsistentData = {}	#key is (chr,pos),
		#value is (noOfInconsistencyInTrio, noOfTotalInTrio, noOfInconsistencyInDuo, noOfTotalInDuo)
		sys.stderr.write("Reading from %s files ...\n"%(len(self.inputFnameLs)))
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				continue
			reader = None
			trioSetStrIndex = None
			chromosomeIndex = None
			posIndex = None
			isInconsistentIndex = None
			try:
				reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
				header = reader.next()
				col_name2index = getColName2IndexFromHeader(header)
				
				trioSetStrIndex = col_name2index.get("#trio_set")
				chromosomeIndex = col_name2index.get("chromosome")
				posIndex = col_name2index.get("pos")
				isInconsistentIndex = col_name2index.get("isInconsistent")
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
			if reader is not None and isInconsistentIndex is not None:
				for row in reader:
					trio_set_str = row[trioSetStrIndex]
					chromosome = row[chromosomeIndex]
					pos = int(row[posIndex])
					isInconsistent = int(row[isInconsistentIndex])
					chr_pos = (chromosome, pos)
					if chr_pos not in chr_pos2inconsistentData:
						chr_pos2inconsistentData[chr_pos] = [0, 0, 0, 0]
					#trio_set_ls = trio_set_str.split(',')
					if trio_set_str.find("0")==0 or trio_set_str.find(",0")!=-1:	#it's a duo. one parent is missing.
						chr_pos2inconsistentData[chr_pos][2] += isInconsistent
						chr_pos2inconsistentData[chr_pos][3] += 1
					else:	#it's a trio
						chr_pos2inconsistentData[chr_pos][0] += isInconsistent
						chr_pos2inconsistentData[chr_pos][1] += 1
						
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Outputting ...")
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		writer.writerow(['#chromosome', 'pos', 'noOfInconsistencyInTrio', 'noOfTotalInTrio', 'inconsistencyRateInTrio',\
						'noOfInconsistencyInDuo', 'noOfTotalInDuo', 'inconsistencyRateInDuo'])
		chr_pos_ls = chr_pos2inconsistentData.keys()
		chr_pos_ls.sort()
		for chr_pos in chr_pos_ls:
			chromosome, pos = chr_pos
			noOfInconsistencyInTrio, noOfTotalInTrio, noOfInconsistencyInDuo, noOfTotalInDuo = chr_pos2inconsistentData.get(chr_pos)
			if noOfTotalInTrio>0:
				inconsistencyRateInTrio = noOfInconsistencyInTrio/float(noOfTotalInTrio)
			else:
				inconsistencyRateInTrio = -1
			if noOfTotalInDuo>0:
				inconsistencyRateInDuo = noOfInconsistencyInDuo/float(noOfTotalInDuo)
			else:
				inconsistencyRateInDuo = -1
			writer.writerow([chromosome, pos, noOfInconsistencyInTrio, noOfTotalInTrio, inconsistencyRateInTrio,\
							noOfInconsistencyInDuo, noOfTotalInDuo, inconsistencyRateInDuo])
		
		del writer
		sys.stderr.write("Done.\n")


if __name__ == '__main__':
	main_class = ReduceTrioInconsistencyByPosition
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
