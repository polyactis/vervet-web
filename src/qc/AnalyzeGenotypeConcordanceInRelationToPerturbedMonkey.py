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

parentClass = AbstractMatrixFileWalker

class AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict
	option_default_dict.update({
				('perturbedMonkeyList', 1, ): ['', '', 1, 'a coma-separated list of individual-alignment.read_group for all perturbed monkeys'],\
				('individualAlignmentCoverageFname', 1, ): ['', '', 1, 'file contains two columns, individual-alignment.read_group, coverage, individual.id.'],\
				('pedigreeFname', 1, ): ['', '', 1, 'pedigree (trios/duos/singletons) file in plink format with all members included. ID is Individual.id\n\
	This is used to figure out the shortest path between two members.'],\
				})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.perturbedMonkeyList = utils.getListOutOfStr(list_in_str=self.perturbedMonkeyList, data_type=str, separator1=',', separator2=None)
	
	def setup(self, **keywords):
		"""
		2013.09.30
		"""
		parentClass.setup(self, **keywords)
		#. read in pedigree graph
		self.plinkPedigreeFile = PlinkPedigreeFile(inputFname=self.pedigreeFname)
		self.pedigreeGraph = self.plinkPedigreeFile.pedigreeGraph
		
		#. read in alignment coverage data
		alignmentCoverageFile = MatrixFile(inputFname=self.individualAlignmentCoverageFname)
		alignmentCoverageFile.constructColName2IndexFromHeader()
		self.alignmentReadGroup2coverage = alignmentCoverageFile.constructDictionary(keyColumnIndexList=[0], valueColumnIndexList=[1],\
														keyUniqueInInputFile=True, valueDataType=float)
		
		#. also get map from read_group to individual.id (used in graph)
		alignmentCoverageFile._resetInput()
		self.alignmentReadGroup2individualID = alignmentCoverageFile.constructDictionary(keyColumnIndexList=[0], valueColumnIndexList=[2], \
														keyUniqueInInputFile=True)
		#. also get map from individual.id to read_group (used in graph)
		alignmentCoverageFile._resetInput()
		self.individualID2alignmentReadGroup = alignmentCoverageFile.constructDictionary(keyColumnIndexList=[2], valueColumnIndexList=[0], \
																		keyUniqueInInputFile=True)
		alignmentCoverageFile.close()
	
	def processHeader(self, header=None, pdata=None, rowDefinition=None):
		"""
		2013.09.30
		"""
		header = header + ["monkeyDownsampled", "coverage", "parentalCoverage", "coAncestryPathLength", "coAncestryPath"]
		self._writeHeader(header=header, pdata=pdata, rowDefinition=rowDefinition)
	
	def getParentalCoverage(self, alignmentReadGroup=None, alignmentReadGroup2coverage=None, \
						pedigreeGraph=None, alignmentReadGroup2individualID=None, individualID2alignmentReadGroup=None):
		"""
		2013.10.04
		"""
		individualID = alignmentReadGroup2individualID.get(alignmentReadGroup)
		parents = pedigreeGraph.predecessors(individualID)
		parentCoverageLs = []
		for parentID in parents:
			parentReadGroup = individualID2alignmentReadGroup.get(parentID)
			coverage = alignmentReadGroup2coverage.get(parentReadGroup, 0)
			parentCoverageLs.append(coverage)
		return sum(parentCoverageLs)
	
	def processRow(self, row=None, pdata=None):
		"""
		2013.09.30
		"""
		returnValue = 0
		col_name2index = getattr(pdata, 'col_name2index', None)
		y_ls = getattr(pdata, 'y_ls', None)
		if col_name2index is not None:
			monkeyIDIndex = col_name2index.get("sample_id", None)
			alignmentReadGroup = row[monkeyIDIndex]	#alignment read-group
			coverage = self.alignmentReadGroup2coverage.get(alignmentReadGroup)
			shortestCoAncestryPath = None
			targetIndividualID = self.alignmentReadGroup2individualID.get(alignmentReadGroup)
			parentalCoverage = self.getParentalCoverage(alignmentReadGroup, self.alignmentReadGroup2coverage, \
									self.pedigreeGraph, self.alignmentReadGroup2individualID, self.individualID2alignmentReadGroup)
			for perturbedMonkey in self.perturbedMonkeyList:
				sourceIndividualID = self.alignmentReadGroup2individualID.get(perturbedMonkey)
				coAncestryPath = self.pedigreeGraph.getShortestPathInUndirectedVersion(sourceIndividualID, targetIndividualID)
				if shortestCoAncestryPath is None or len(coAncestryPath)<len(shortestCoAncestryPath):
					shortestCoAncestryPath = coAncestryPath
			shortestCoAncestryPathEdgeDirectionLs = []
			#edge[2] for edge in shortestCoAncestryPath ]
			for edge in shortestCoAncestryPath:
				if edge[2]==1:
					shortestCoAncestryPathEdgeDirectionLs.append('+')
				else:
					shortestCoAncestryPathEdgeDirectionLs.append("-")
			data_row = row + [','.join(self.perturbedMonkeyList), coverage, parentalCoverage, len(shortestCoAncestryPath), ''.join(shortestCoAncestryPathEdgeDirectionLs)]
			if self.invariantPData.writer:
				self.invariantPData.writer.writerow(data_row)
				returnValue = 1
		return returnValue
	

if __name__ == '__main__':
	main_class = AnalyzeGenotypeConcordanceInRelationToPerturbedMonkey
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()