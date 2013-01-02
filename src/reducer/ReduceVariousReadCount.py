#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011-11-4
		input should be the output file of GATK's VariousReadCountWalker (custom written).
		The input name contains the chromosome ID. like xxx_chr2.tsv
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


class ReduceVariousReadCount(object):
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
		writer.writerow(['#sampleID', 'chromosome', 'length', 'noOfReadsAlignedByLength', 'noOfSingletonsByLength', \
						'noOfPairsOnSameContigByLength',\
						'meanInferInsertSize', 'noOfPairsOnDifferentContigsByLength'])
		for inputFname in self.inputFnameLs:
			inputFile = utils.openGzipFile(inputFname)
			delimiter = figureOutDelimiter(inputFile)
			reader = csv.reader(inputFile, delimiter=delimiter)
			header = reader.next()
			col_name2index = getColName2IndexFromHeader(header)
			
			sampleIDIndex = col_name2index.get("readGroup")
			chromosomeIndex = col_name2index.get("firstReferenceName")
			chromosomeLengthIndex = col_name2index.get("firstReferenceLength")
			
			numberOfReadsIndex = col_name2index.get("numberOfReads")
			numberOfReadsAlignedIndex = col_name2index.get("numberOfReadsAligned")
			numberOfSingletonsMappedIndex = col_name2index.get("numberOfSingletonsMapped")
			numberOfPairsOnSameContigIndex = col_name2index.get("numberOfPairsOnSameContig")
			numberOfPairsOnDifferentContigsIndex = col_name2index.get("numberOfPairsOnDifferentContigs")
			meanInsertSizeIndex = col_name2index.get("meanInsertSize")
			
			for row in reader:
				sampleID = row[sampleIDIndex]
				chromosome = row[chromosomeIndex]
				chromosomeLength = int(row[chromosomeLengthIndex])
				
				numberOfReads = float(row[numberOfReadsIndex])
				numberOfReadsAligned = float(row[numberOfReadsAlignedIndex])
				numberOfSingletonsMapped = float(row[numberOfSingletonsMappedIndex])
				numberOfPairsOnSameContig = float(row[numberOfPairsOnSameContigIndex])
				numberOfPairsOnDifferentContigs = float(row[numberOfPairsOnDifferentContigsIndex])
				meanInsertSize = row[meanInsertSizeIndex]
				
				writer.writerow([sampleID, chromosome, chromosomeLength, numberOfReadsAligned/chromosomeLength, \
								numberOfSingletonsMapped/chromosomeLength,\
								numberOfPairsOnSameContig/chromosomeLength, meanInsertSize, \
								numberOfPairsOnDifferentContigs/chromosomeLength])
		del writer
		sys.stderr.write("Done.\n")


if __name__ == '__main__':
	main_class = ReduceVariousReadCount
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
