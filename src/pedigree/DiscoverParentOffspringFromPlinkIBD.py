#!/usr/bin/env python
"""
Examples:
	# 2012.8.21
	%s -i ~/NetworkData/vervet/Kinx2Apr2012.txt -u yh
		-l PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/LDPrunedMerged.genome
		-n PlinkSexCheck/PlinkSexCheck_Method18_W100Z10R0.8_maxContigID195.2012.8.21T1207/sexCheckSexCheck/LDPrunedMerged.sexcheck
		-o PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/crossLabelCheck.tsv
	
	#2012.8.22 greedily removes the worst contaminant
	%s -i ~/NetworkData/vervet/Kinx2Apr2012.txt
		-l PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/LDPrunedMerged.genome
		-n PlinkSexCheck/PlinkSexCheck_Method18_W100Z10R0.8_maxContigID195.2012.8.21T1207/sexCheckSexCheck/LDPrunedMerged.sexcheck
		-o PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/wrongLabelChiSq_greedy.tsv
		-u yh  -A
		

Description:
	2012.8.21
		Program that plots kinship vs. PI_hat from plinkIBD result.
		inputFname is the kinship file from Sue.
		Besides the png scatterplot (_scatter.png), a tsv file (_table.tsv) will be created to hold pedigree kinship vs. IBD PI_HAT. 
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

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, SNPData
from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB, utils
from pymodule.MatrixFile import MatrixFile
from pymodule import SNP
import numpy, random, pylab
import numpy as np
from vervet.src import VervetDB
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
#from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
import heapq
from scipy import stats
import networkx as nx

class DiscoverParentOffspringFromPlinkIBD(AbstractMapper):
	__doc__ = __doc__
#						
	option_default_dict = AbstractMapper.option_default_dict
	option_default_dict.update({
						('maxDistanceToPOVector', 1, float): [0.04, 'm', 1, 'the IBD0,IBD1,IBD2 sharing probability vector for PO is 0,1,0. \
			This cutoff controls how far the observed data should be from that vector.'], \
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def constructPedigreeGraphFromPlinkIBD(self, inputFname=None, maxDistanceToPOVector=0.04, drawDistribution=False, outputFnamePrefix=None):
		"""
		2012.8.14
		"""
		sys.stderr.write("Constructing pedigree-graph out of plink-ibd %s ..."%(inputFname))
		DG=nx.DiGraph()
		childNodeSet = set()
		reader = MatrixFile(inputFname)
		reader.constructColName2IndexFromHeader()
		
		monkey1IDIndex = reader.getColIndexGivenColHeader("IID1")
		monkey2IDIndex = reader.getColIndexGivenColHeader("IID2")
		Z0Index = reader.getColIndexGivenColHeader("Z0")
		Z1Index = reader.getColIndexGivenColHeader("Z1")
		Z2Index = reader.getColIndexGivenColHeader("Z2")
		
		poVector = numpy.array([0,1,0.0])
		counter = 0
		real_counter = 0
		
		data_ls = []
		for row in reader:
			monkey1ID = int(row[monkey1IDIndex])	#turn it into integer so could compare age
			monkey2ID = int(row[monkey2IDIndex])
			Z0 = float(row[Z0Index])
			Z1 = float(row[Z1Index])
			Z2 = float(row[Z2Index])
			ZVector = numpy.array([Z0, Z1, Z2])
			dist = numpy.linalg.norm(poVector-ZVector)
			if drawDistribution and outputFnamePrefix:
				data_ls.append(dist)
			if dist<=maxDistanceToPOVector:
				if monkey1ID>monkey2ID:
					childID = monkey1ID
					parentID = monkey2ID
				else:
					childID = monkey2ID
					parentID = monkey1ID
				DG.add_edge(parentID, childID, weight=dist)
				childNodeSet.add(childID)
				real_counter += 1
			counter += 1
		del reader
		sys.stderr.write("%s out of %s lines become PO pairs. %s children, %s nodes. %s edges. %s connected components.\n"%(\
							real_counter, counter, len(childNodeSet), DG.number_of_nodes(), DG.number_of_edges(), \
							nx.number_connected_components(DG.to_undirected())))
		if drawDistribution and outputFnamePrefix:
			outputFname = '%s_IBDVector2POVectorDist_hist.png'%(outputFnamePrefix)
			yh_matplotlib.drawHist(data_ls, title='', \
								xlabel_1D="dist(ZVector,POVector)", xticks=None, \
								outputFname=outputFname, min_no_of_data_points=10, \
								needLog=True, \
								dpi=200, min_no_of_bins=25)
		return PassingData(DG=DG, childNodeSet=childNodeSet)
	
	def outputPedigreeGraph(self, DG=None, outputFname=None):
		"""
		2012.8.23
		"""
		sys.stderr.write("Outputting graph %s nodes. %s edges. %s connected components..." %\
						(DG.number_of_nodes(), DG.number_of_edges(), \
							nx.number_connected_components(DG.to_undirected())))
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['parentID', 'childID', 'distToPOVector']
		writer.writerow(header)
		counter = 0
		for e in DG.edges_iter(data=True):
			parentID = e[0]
			childID = e[1]
			data = e[2]
			weight = data['weight']
			data_row = [parentID, childID, weight]
			writer.writerow(data_row)
			counter +=1
		sys.stderr.write(" %s edges outputted. \n"%(counter))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		pdata = self.constructPedigreeGraphFromPlinkIBD(inputFname=self.inputFname, \
								maxDistanceToPOVector=self.maxDistanceToPOVector, drawDistribution=False, \
								outputFnamePrefix=self.outputFnamePrefix)
		self.outputPedigreeGraph(DG=pdata.DG, outputFname=self.outputFname)
		

if __name__ == '__main__':
	main_class = DiscoverParentOffspringFromPlinkIBD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
