#!/usr/bin/env python
"""
Examples:
	%s  -o /tmp/inconsistentParentsOfSharedKids.tsv -y 2
		~/../sservice/Vervet/SeqData2012/Linda.fam 
		PlinkIBDCheck/PlinkIBDCheck_Method14_W100Z10R0.4.2012.8.13T1413/vcf2plinkVCF2PlinkMerged/pedigree.2.tfam
		/u/home/eeskin/sservice/Vervet/Expression/SolarPed.txt

	%s -o /tmp/inconsistentParentsOfSharedIndividuals.tsv
		~/../sservice/Vervet/SeqData2012/Linda.fam 
		PlinkIBDCheck/PlinkIBDCheck_Method14_W100Z10R0.4.2012.8.13T1413/vcf2plinkVCF2PlinkMerged/pedigree.2.tfam
		/u/home/eeskin/sservice/Vervet/Expression/SolarPed.txt

Description:
	2012.8.14 program that compares pedigrees from multiple input files.
		Input files should have the child/father/mother on the same columns.
		Doesn't mater which delimiter is used.
		If first line's first cell contains characters [a-df-zA-DF-Z\-] in the end, 
			#no 'e' or 'E', as could be used in scientific number,
			first line is regarded as header and skipped.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter
from pymodule import SNP
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
import networkx as nx

class ComparePedigreeFromMultipleInput(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						('childColumnIndex', 1, int):[1, '', 1, 'index of the child ID column, 0-based. same across input.'],\
						('fatherColumnIndex', 1, int):[2, '', 1, 'index of the father ID column, 0-based'],\
						('motherColumnIndex', 1, int):[3, '', 1, 'index of the mother ID column, 0-based'],\
						('run_type', 1, int):[1, 'y', 1, 'run type \
		1: check parents of all individuals shared among all input pedigrees,\
		2: check parents of only the children that are shared among all input pedigrees'],\
						
						})
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def constructPedigreeGraphFromOneFile(self, inputFname=None):
		"""
		2012.8.14
		"""
		sys.stderr.write("Constructing pedigree-graph out of %s ..."%(inputFname))
		DG=nx.DiGraph()
		reader = None
		childNodeSet = set()
		try:
			inputFile = utils.openGzipFile(inputFname)
			delimiter = figureOutDelimiter(inputFname)
			isCSVReader = True
			if delimiter=='\t' or delimiter==',':
				reader = csv.reader(inputFile, delimiter=delimiter)
			else:
				reader = inputFile
				isCSVReader = False
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.exit(3)
		
		counter = 0
		if reader is not None:
			for row in reader:
				if not isCSVReader:
					row = row.strip().split()
				if counter ==0 and self.p_char.search(row[0]):	#character in 1st cell of 1st line, it's header skip.
					continue
				childID = row[self.childColumnIndex]
				DG.add_node(childID)	#in case this guy has no parents, then won't be added via add_edge()
				childNodeSet.add(childID)
				fatherID = row[self.fatherColumnIndex]
				if fatherID!='0':
					DG.add_edge(fatherID, childID)
				motherID = row[self.motherColumnIndex]
				if motherID!='0':
					DG.add_edge(motherID, childID)
				counter += 1
			del reader
		sys.stderr.write("%s children, %s nodes. %s edges. %s connected components.\n"%(\
										len(childNodeSet), DG.number_of_nodes(), DG.number_of_edges(), \
										nx.number_connected_components(DG.to_undirected())))
		return PassingData(DG=DG, childNodeSet=childNodeSet)
	
	def traverse(self):
		"""
		2012.1.9
		"""
		newHeader = []
		inputFname2DG = {}
		
		commonChildNodeSet = None
		
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				sys.stderr.write("Warning: file %s doesn't exist.\n"%(inputFname))
				continue
			graphData = self.constructPedigreeGraphFromOneFile(inputFname=inputFname)
			inputFname2DG[inputFname] = graphData.DG
			if commonChildNodeSet is None:
				commonChildNodeSet = graphData.childNodeSet
			else:
				commonChildNodeSet = commonChildNodeSet.intersection(graphData.childNodeSet)
		
		returnData = PassingData(inputFname2DG=inputFname2DG, commonChildNodeSet=commonChildNodeSet)
		return returnData
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		returnData = self.traverse()
		inputFname2DG = returnData.inputFname2DG
		if self.run_type==1:
			commonNodeSet = None
			for inputFname, DG in inputFname2DG.iteritems():
				if commonNodeSet is None:
					commonNodeSet = set(DG.nodes())
				else:
					commonNodeSet = commonNodeSet.intersection(set(DG.nodes()))
		elif self.run_type==2:
			#only use the child node set
			commonNodeSet = returnData.commonChildNodeSet
		else:
			sys.stderr.write("run_type %s  not supported.\n"%(self.run_type))
			sys.exit(4)
		sys.stderr.write("%s common nodes across %s files.\n"%(len(commonNodeSet), len(inputFname2DG)))
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ["childID", ]
		for inputFname in self.inputFnameLs:
			header.append("parentsIn%s"%(os.path.basename(inputFname)))
		
		writer.writerow(header)
		counter = 0
		real_counter = 0
		for node_id in commonNodeSet:
			parentTuple2Count = {}
			inputFname2parentTuple = {}
			for inputFname, DG in inputFname2DG.iteritems():
				if node_id in DG:
					parents = DG.predecessors(node_id)
					if parents is None:
						parentTuple = ""
					else:
						parents.sort()
						parentTuple = tuple(parents)
				else:
					parentTuple = "notInGraph"
				if parentTuple not in parentTuple2Count:
					parentTuple2Count[parentTuple] = 0
				parentTuple2Count[parentTuple] += 1
				inputFname2parentTuple[inputFname] = parentTuple
			counter += 1
			if len(parentTuple2Count)>1:
				real_counter += 1
				#sys.stderr.write("inconsistent parents for node %s. %s.\n"%(node_id, repr(inputFname2parentTuple)))
				data_row = [node_id]
				for inputFname in self.inputFnameLs:
					data_row.append(repr(inputFname2parentTuple[inputFname]))
				writer.writerow(data_row)
		del writer
		sys.stderr.write("%s out of %s nodes have inconsistent parents.\n"%(real_counter, counter))


if __name__ == '__main__':
	main_class = ComparePedigreeFromMultipleInput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()