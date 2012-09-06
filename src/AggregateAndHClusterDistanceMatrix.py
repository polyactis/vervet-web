#!/usr/bin/env python
"""
Examples:
	#
	%s -f dendrogram -o /tmp/454_illu_6_sub_vs_1MbBAC.tsv pairwiseDistMatrix/*.tsv
	
	
	%s 
	
Description:
	Aggregate pairwise distance matrices, output by CalculatePairwiseDistanceOutOfSNPXStrainMatrix.py
		into one distance matrix. Then run PCA and hclustering on this matrix.
	Output:
		1. .png .svg: hcluster dendrogram output figures
		2. _PCA.tsv: matrix with top 6 PCs from running PCA. It could be uploaded to custom motionchart web app.
"""
import os, sys, numpy
__doc__ = __doc__%(sys.argv[0], sys.argv[0])
#2007-03-05 common codes to initiate database connection
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib


class AggregateAndHClusterDistanceMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 0, ): ['', 'o', 1, 'output the combined distance matrix if given', ],\
						('figureFnamePrefix', 1, ): ['', 'f', 1, 'dendrogram output filename prefix. .png and .svg will be outputted.', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-10-19
		This class is the entry to all others.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.inputFnameLs = inputFnameLs
	
	def massageSampleId(self, col_id_ls):
		"""
		2012.9.5
			added takeLastSplit further massage the sample ID
		"""
		
		take1stSplit = lambda x: x.split("_GA_vs_")[0]
		takeLastSplit = lambda x: x.split("_")[-1]
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			newXLabel = take1stSplit(col_id)
			if newXLabel[:7] == 'VRC_ref':
				newXLabel = col_id.split('_vs_')[0]	#a different split for VRC_ref to retain the sequencer
				newXLabel = newXLabel[4:]
			else:
				newXLabel = takeLastSplit(newXLabel)
			
			col_id_ls[i] = newXLabel
		return col_id_ls
	
	def outputMismatchData(self, outputFname, samplePair2data, distanceMatrix, sampleId2index, sampleIdLs=None):
		"""
		2011-10-19
		"""
		sys.stderr.write("Outputting aggregated distance data to %s ..."%(outputFname))
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		samplePairLs = samplePair2data.keys()
		samplePairLs.sort()
		
		for samplePair in samplePairLs:
			data = samplePair2data.get(samplePair)
			no_of_mismatches, no_of_total_non_NA = data[:2]
			distance = no_of_mismatches/no_of_total_non_NA
			writer.writerow([samplePair[0], samplePair[1], distance, no_of_mismatches, no_of_total_non_NA])
		
		header = [''] + sampleIdLs
		writer.writerow(header)
		no_of_rows = len(sampleId2index)
		for i in range(no_of_rows):
			row_id = sampleIdLs[i]
			data_row = [row_id] + list(distanceMatrix[i])
			writer.writerow(data_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def readInput(self, inputFnameLs, ):
		sys.stderr.write("Reading distance data from %s files ..."%(len(inputFnameLs)))
		sampleId2index = {}
		samplePair2data = {}	#value is [no_of_mismatches, no_of_total_non_NA]
		for inputFname in self.inputFnameLs:
			reader = csv.reader(open(inputFname, ), delimiter=figureOutDelimiter(inputFname))
			matrixStart = False
			for row in reader:
				if row[0]=='':
					matrixStart = True
					break
				
				sample1Id = row[0]
				if sample1Id not in sampleId2index:
					sampleId2index[sample1Id] = len(sampleId2index)
				sample2Id = row[1]
				if sample2Id not in sampleId2index:
					sampleId2index[sample2Id] = len(sampleId2index)
				no_of_mismatches = float(row[-2])
				no_of_total_non_NA = float(row[-1])
				samplePair = (sample1Id, sample2Id)
				if samplePair not in samplePair2data:
					samplePair2data[samplePair] = [0, 0]
				samplePair2data[samplePair][0] += no_of_mismatches
				samplePair2data[samplePair][1] += no_of_total_non_NA
			
			del reader
		sys.stderr.write("Done.\n")
		return sampleId2index, samplePair2data
	
	def runPCAOnDistanceMatrix(self, distanceMatrix, col_id_ls=[], outputFname=None):
		"""
		2011-10-21
		"""
		sys.stderr.write("Running PCA on distance matrix and output ...")
		#import pca_module
		from pymodule.PCA import PCA
		PCA_writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		#T, P, explained_var = pca_module.PCA_svd(phenData_trans.data_matrix, standardize=True)
		T, P, explained_var = PCA.eig(distanceMatrix, normalize=True)	#normalize=True causes missing value in the covariance matrix
		header = ['MonkeyID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		PCA_writer.writerow(header)
		for i in range(len(distanceMatrix)):
			col_id = col_id_ls[i]
			T[i,0] = -T[i,0]	#reverse the sign of PC1 to match geography better
			data_row = [col_id] + list(T[i,0:6])
			PCA_writer.writerow(data_row)
		del PCA_writer
		sys.stderr.write("Done.\n")
	
	def run(self,):
		if self.debug:	# 2010-4-18 enter debug mode "~/.../variation/misc.py -b"
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		sampleId2index, samplePair2data = self.readInput(self.inputFnameLs)
		
		sys.stderr.write("Calculating distance matrix for aggregated data ...")
		distanceMatrix = numpy.zeros([len(sampleId2index), len(sampleId2index)])
		for samplePair, data in samplePair2data.iteritems():
			no_of_mismatches, no_of_total_non_NA = data[:2]
			distance = no_of_mismatches/no_of_total_non_NA
			sample1Id, sample2Id = samplePair[:2]
			sample1Index = sampleId2index[sample1Id]
			sample2Index = sampleId2index[sample2Id]
			distanceMatrix[sample1Index, sample2Index] = distance
			distanceMatrix[sample2Index, sample1Index] = distance
		sys.stderr.write("Done.\n")
		
		sampleIdLs = sampleId2index.keys()
		for sampleId, list_index in sampleId2index.iteritems():
			sampleIdLs[list_index] = sampleId
		
		if self.outputFname:
			self.outputMismatchData(self.outputFname, samplePair2data, distanceMatrix, sampleId2index, sampleIdLs)
		
		massagedSampleIDLs = self.massageSampleId(sampleIdLs)
		
		#2012.9-6 stop massaging sample IDs for PCA output. mapper/AppendInfo2SmartPCAOutput.py could be applied to this.
		self.runPCAOnDistanceMatrix(distanceMatrix, col_id_ls=sampleIdLs, outputFname='%s_PCA.tsv'%(self.figureFnamePrefix))
		
		import pylab
		from hcluster import pdist, linkage, dendrogram
		pylab.clf()
		Z=linkage(distanceMatrix, 'single')
		yh_matplotlib.setFontAndLabelSize(base_size=3)
		dendrogram(Z, color_threshold=0.001, labels=massagedSampleIDLs, orientation='right', leaf_font_size=None)	#leaf_font_size=1 or 5 has no effect
		pylab.savefig('%s.svg'%self.figureFnamePrefix, dpi=200)
		pylab.savefig('%s.png'%self.figureFnamePrefix, dpi=300)
		sys.exit(0)
	
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = AggregateAndHClusterDistanceMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
