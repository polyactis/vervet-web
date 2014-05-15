#!/usr/bin/env python
"""
Examples:
	%s -i input.vcf -o -o selectedStKitts.vcf --country_id_ls 144 --tax_id_ls 60711
	
	%s -i input.vcf -o -o selectedNevis.vcf --country_id_ls 148 --tax_id_ls 60711

	%s -i input.vcf.gz -o selected.vcf --country_id_ls 135,136,144,148,151 --tax_id_ls 60711
		

Description:
	2012.10.5 program that extracts samples from a VCF and form a new VCF.
		need to re-calculate the AC/AF values of each variant.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])


sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, utils
from vervet.src import VervetDB
#from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from pymodule.yhio.AbstractMatrixFileWalker import AbstractMatrixFileWalker
from pymodule.yhio.PlinkPedigreeFile import PlinkPedigreeFile
from pymodule.yhio.MatrixFile import MatrixFile
from AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey import AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey
parentClass = AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey

class AnalyzeTrioInconsistencyInRelationToPerturbedMonkey(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict
	option_default_dict.update({
				('originalTrioInconsistencyFname', 1, ): ['', '', 1, 'file contains original (before downsampling) trio inconsistency result'],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def setup(self, **keywords):
		"""
		2013.09.30
		"""
		parentClass.setup(self, **keywords)
		
		#. read in alignment coverage data
		trioInconsistencyFile = MatrixFile(inputFname=self.originalTrioInconsistencyFname)
		trioInconsistencyFile.constructColName2IndexFromHeader()
		self.trioID2dataList = trioInconsistencyFile.constructDictionary(keyColumnIndexList=[0], valueColumnIndexList=[1,2,3],\
														keyUniqueInInputFile=True)
		trioInconsistencyFile.close()
		
		self.alignmentID2alignmentReadGroup = {}
		for alignmentReadGroup in self.alignmentReadGroup2individualID:
			alignmentID = alignmentReadGroup.split('_')[0]
			self.alignmentID2alignmentReadGroup[alignmentID] = alignmentReadGroup
		
	def getAlignmentReadGroupFromCompositeID(self, alignmentCompositeID=None, alignmentID2alignmentReadGroup=None):
		"""
		2013.10.04
		"""
		alignmentID = alignmentCompositeID.split('_')[0]
		return alignmentID2alignmentReadGroup.get(alignmentID)
	
	def processHeader(self, header=None, pdata=None, rowDefinition=None):
		"""
		2013.09.30
		"""
		header = header + ["monkeyDownsampled", "coverage", "parentalCoverage", "coAncestryPathLength", "coAncestryPath", "trioInconsistencyDelta"]
		self._writeHeader(header=header, pdata=pdata, rowDefinition=rowDefinition)
	
	def processRow(self, row=None, pdata=None):
		"""
		2013.09.30
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index is not None:
			trioSetIndex = col_name2index.get("#trio_set")
			trioInconsistencyIndex = col_name2index.get("no_of_inconsistent_by_no_of_total")
			trioSet = row[trioSetIndex]	#in alignment read-group
			trioCoverage = []
			shortestCoAncestryPath = None
			for alignmentCompositeID in trioSet.split(','):
				alignmentReadGroup = self.getAlignmentReadGroupFromCompositeID(alignmentCompositeID, \
																			alignmentID2alignmentReadGroup=self.alignmentID2alignmentReadGroup)
				
				coverage = self.alignmentReadGroup2coverage.get(alignmentReadGroup, 0)
				trioCoverage.append(coverage)
				if alignmentCompositeID=='0':
					continue
				
				targetIndividualID = self.alignmentReadGroup2individualID.get(alignmentReadGroup)
				for perturbedMonkey in self.perturbedMonkeyList:
					sourceIndividualID = self.alignmentReadGroup2individualID.get(perturbedMonkey)
					coAncestryPath = self.pedigreeGraph.getShortestPathInUndirectedVersion(sourceIndividualID, targetIndividualID)
					if shortestCoAncestryPath is None or len(coAncestryPath)<len(shortestCoAncestryPath):
						shortestCoAncestryPath = coAncestryPath
			
			parentalCoverage = sum(trioCoverage[:2])
			trioCoverage = map(repr, trioCoverage)
			shortestCoAncestryPathEdgeDirectionLs = []
			#edge[2] for edge in shortestCoAncestryPath ]
			for edge in shortestCoAncestryPath:
				if edge[2]==1:
					shortestCoAncestryPathEdgeDirectionLs.append('+')
				else:
					shortestCoAncestryPathEdgeDirectionLs.append("-")
			trioInconsistency = float(row[trioInconsistencyIndex])
			originalTrioInconsistency = float(self.trioID2dataList[trioSet][2])
			data_row = row + [','.join(self.perturbedMonkeyList), ','.join(trioCoverage),  parentalCoverage, len(shortestCoAncestryPath), ''.join(shortestCoAncestryPathEdgeDirectionLs),
						trioInconsistency-originalTrioInconsistency]
			if self.invariantPData.writer:
				self.invariantPData.writer.writerow(data_row)
				returnValue = 1
		return returnValue
	

if __name__ == '__main__':
	main_class = AnalyzeTrioInconsistencyInRelationToPerturbedMonkey
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()