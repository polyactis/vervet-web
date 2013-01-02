#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i 180samples.topoff.included.readCount.xls  -o 180samples.topoff.included.readCount.normalized.tsv -m 180

Description:
	2012.5.8
		Input is Gene X Sample expression matrix. The 1st column is gene ID. The 1st row is the header (geneID + list of sample ID).
		For normalization:
			1. divide each sample's gene expression by its total sum.
			2. for each gene, derive its mean/stdev across all samples
			3. new value = (exprValue-mean)/stdev 
		Output is same format. tab-delimited.
	

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, getColName2IndexFromHeader 
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import random, numpy

class NormalizeGeneExprForEskinGroup(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('minExprSumPerGene', 1, float): [None, 'm', 1, 'mininum sum of gene expression values across all samples. a filter before any normalization is applied..'],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
	
	def readDataMatrix(self, inputFname, minExprSumPerGene=180):
		"""
		2012.5.8
		"""
		sys.stderr.write("Reading the gene expression matrix from %s ..."%(inputFname))
		
		suffix = os.path.splitext(inputFname)[1]
		if suffix=='.gz':
			import gzip
			inf = gzip.open(inputFname, 'r')
		else:
			inf = open(inputFname, 'r')
		
		reader = csv.reader(inf, delimiter=figureOutDelimiter(inf))
		header = reader.next()	#first line is taken as header
		colName2Index = getColName2IndexFromHeader(header)
		data_matrix = []
		row_id_ls = []
		counter = 0
		real_counter = 0
		for row in reader:
			data_row = row[1:]
			data_row = map(float, data_row)
			exprSumPerGene = sum(data_row)
			counter += 1
			if exprSumPerGene>=minExprSumPerGene:
				real_counter += 1
				row_id_ls.append(row[0])
				data_matrix.append(data_row)
		data_matrix = numpy.array(data_matrix)
		sys.stderr.write("%s rows out of %s selected. %s rows , %s columns.\n"%(real_counter, counter, \
																	len(row_id_ls), len(header)-1))
		return PassingData(row_id_ls=row_id_ls, header=header, data_matrix=data_matrix)
	
	def outputMatrixData(self, matrixData, outputFname=None):
		"""
		2012.5.8
		"""
		sys.stderr.write("Outputting matrix data %s ..."%(repr(matrixData.data_matrix.shape)))
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		writer.writerow(matrixData.header)
		for i in xrange(len(matrixData.row_id_ls)):
			writer.writerow([matrixData.row_id_ls[i]]+ list(matrixData.data_matrix[i,:]))
		del writer
		sys.stderr.write("Done.\n")
	
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		matrixData = self.readDataMatrix(self.inputFname, minExprSumPerGene=self.minExprSumPerGene)
		
		sys.stderr.write("Normalizing ...")
		colSumArray = numpy.sum(matrixData.data_matrix, axis=0)
		matrixData.data_matrix = matrixData.data_matrix/colSumArray
		
		rowMeanArray = numpy.mean(matrixData.data_matrix, axis=1)
		rowStdArray = numpy.std(matrixData.data_matrix, axis=1, dtype=numpy.float, ddof=1)
		
		def subtractMeanDivideByStd(dataRow, meanArray, stdArray):
			return (dataRow-meanArray)/stdArray
		
		subtractMeanDivideByStdVec = numpy.vectorize(subtractMeanDivideByStd)
		
		#iteration of a 2D matrix goes row by row (dimension: no of columns). However the rowMeanArray's dimension is number of rows.
		#Applying the transpose() solves the problem. 
		matrixData.data_matrix = subtractMeanDivideByStdVec(matrixData.data_matrix.transpose(), rowMeanArray, rowStdArray).transpose()
		sys.stderr.write("Done.\n")
		
		self.outputMatrixData(matrixData, outputFname=self.outputFname)
	
if __name__ == '__main__':
	main_class = NormalizeGeneExprForEskinGroup
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()