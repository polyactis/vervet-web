#!/usr/bin/env python
"""
Examples:
	# 2012.8.21
	%s -i
	
	#2012.8.22 greedily removes the worst contaminant
	%s -i wrongLabelChiSq_greedy_4_minAbsDelta0.2.tsv -O ~/script/vervet/doc/figures/labelContaminationChecking/pedigreeWithOutlier/method22Outlier
		-u yh -z uclaOffice
		

Description:
	2012.9.4 program that plots whole pedigree with the addition of outlier edges, colored by IBD-kinsip. 
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
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
import networkx as nx
import matplotlib as mpl


class PlotPedigree(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
						('kinshipIBDOutlierFname', 1, ): [None, 'i', 1, 'tsv file containing the outlier pairs,\
	output of DetectWrongLabelByCompKinshipVsIBD.py.'], \
						('minAbsDeltaForOutlier', 1, float): [0.2, 'm', 1, 'if not 0, this will be used as minAbsDelta in outlier frequency counting'], \
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getMonkeyID2Coverage(self, inputFname):
		"""
		2012.9.4
			copied from vervet/src/misc.py
		2012.2.10
			inputFname is output of SequencingStrategy.assignVRCSequencePriorityBasedOnPedigree() + manual change of top ones
		"""
		sys.stderr.write("Reading the list of ranked monkeys from %s ..."%(inputFname))
		reader = MatrixFile(inputFname)
		reader.constructColName2IndexFromHeader()
		
		monkey_id_index = reader.getColIndexGivenColHeader("UCLAID")
		pre_set_coverage_index = reader.getColIndexGivenColHeader("pre-set-coverage")
		future_coverage_index = reader.getColIndexGivenColHeader("future coverage")
		to_sequence_monkey_id2coverage = {}
		for row in reader:
			monkey_id = row[monkey_id_index]
			pre_set_coverage = row[pre_set_coverage_index]
			if pre_set_coverage:
				pre_set_coverage = float(pre_set_coverage)
			else:
				pre_set_coverage = 0
			future_coverage = 0
			if len(row)>=future_coverage_index+1:
				future_coverage = float(row[future_coverage_index])
			to_sequence_monkey_id2coverage[monkey_id] = max(future_coverage, pre_set_coverage)
		del reader
		sys.stderr.write(" %s monkeys are to-be-sequenced.\n"%(len(to_sequence_monkey_id2coverage)))
		return to_sequence_monkey_id2coverage
		
		"""
		#below is to read from a more simple file 
		from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
		import VervetDB
		import csv
		
		sys.stderr.write("Reading from %s ... "%(inputFname))
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey_id_index = col_name2index.get("UCLAID")
		coverage_index = col_name2index.get("coverage")
		if coverage_index is None:
			coverage_index = col_name2index.get("targetCoverage")
		monkey_id2preDeterminedCoverage = {}
		for row in reader:
			monkey_id = row[monkey_id_index]
			coverage = float(row[coverage_index])
			monkey_id2preDeterminedCoverage[monkey_id] = coverage
		del reader
		sys.stderr.write("%s monkeys.\n"%(len(monkey_id2preDeterminedCoverage)))
		return monkey_id2preDeterminedCoverage
		"""
	
	def getOutlierPair(self, inputFname=None, minAbsDeltaForOutlier=0.2):
		"""
		2012.9.4
		"""
		sys.stderr.write("Getting outlier pairs from %s with minAbsDeltaForOutlier=%s ..."%(inputFname, minAbsDeltaForOutlier))
		reader = MatrixFile(inputFname)
		reader.constructColName2IndexFromHeader()
		monkey1IDIndex = reader.getColIndexGivenColHeader("monkey1ID")
		monkey2IDIndex = reader.getColIndexGivenColHeader("monkey2ID")
		pedigreeKinshipIndex = reader.getColIndexGivenColHeader("pedigreeKinship")
		IBDIndex = reader.getColIndexGivenColHeader("IBD_PI_HAT")
		monkeyPair2data  = {}
		counter = 0
		real_counter = 0
		minDelta = None
		maxDelta = None
		monkeyID2monkeyPair2data = {}
		for row in reader:
			monkey1ID = row[monkey1IDIndex]
			monkey2ID = row[monkey2IDIndex]
			kinship = float(row[pedigreeKinshipIndex])
			IBD = float(row[IBDIndex])
			absDelta = abs(kinship-IBD)
			if absDelta>=minAbsDeltaForOutlier:
				delta = IBD-kinship
				data = PassingData(kinship=kinship, IBD=IBD, delta=delta, absDelta=absDelta)
				key = (monkey1ID, monkey2ID)
				if key in monkeyPair2data:
					sys.stderr.write("WARNING: key %s has value %s in monkeyPair2data already. value ovewritten with %s.\n"%\
									(repr(key), monkeyPair2data.get(key), data))
				monkeyPair2data[key] = data
				if minDelta is None or delta<minDelta:
					minDelta = delta
				if maxDelta is None or delta>maxDelta:
					maxDelta = delta
				real_counter += 1
				if monkey1ID not in monkeyID2monkeyPair2data:
					monkeyID2monkeyPair2data[monkey1ID] ={}
				monkeyID2monkeyPair2data[monkey1ID][key] = data
				
				if monkey2ID not in monkeyID2monkeyPair2data:
					monkeyID2monkeyPair2data[monkey2ID] ={}
				monkeyID2monkeyPair2data[monkey2ID][key] = data
				
			counter += 1
		sys.stderr.write(" %s/%s outlier pairs, minDelta=%s, maxDelta=%s.\n"%(real_counter, counter, minDelta, maxDelta))
		return PassingData(monkeyPair2data=monkeyPair2data, minDelta=minDelta, maxDelta=maxDelta,\
						monkeyID2monkeyPair2data=monkeyID2monkeyPair2data)
	
	sex2node_property_list = {}		# in node_property_list, each entry is (node, size, color)
		#size depends on whether it's deep-sequenced (30X), 4 for 30X, 8 for the REF, 2 for all others. 
		#color depends on it's sequenced or not
	def drawPedigree(self, db_vervet=None, outputFnamePrefix=None, baseNodeSize=4, monkey_id2coverage=None, \
					monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None):
		"""
		2012.9.4
			copied from vervet/src/misc.py
		2012.2.10
			add argument monkeyCoverageFname to override coverage data from db
		2011-5-6
			
		"""
		DG = db_vervet.constructPedgree()
		
		"""
		sys.stderr.write("Plotting out degree histogram ....")
		out_degree_ls = []
		#for n in DG.nodes_iter():
		#	out_degree_ls.append(DG.out_degree(n))
		#outputFname = '%s_outDegreeHist.png'%(outputFnamePrefix)
		#yh_matplotlib.drawHist(out_degree_ls, title="Histogram of no. of children per monkey", xlabel_1D="no. of children", \
		#					outputFname=outputFname, min_no_of_data_points=50, needLog=True)
		"""
		sys.stderr.write("Assigning each node with different size/color ...")
		
		sex2node_property_list = self.sex2node_property_list
		if not self.sex2node_property_list:
			for v in DG:
				individual = VervetDB.Individual.get(v)
				sex = individual.sex
				
				if individual.sex not in sex2node_property_list:
					sex2node_property_list[sex] = []
				individual_sequence = VervetDB.IndividualSequence.query.filter_by(sequencer='GA').\
					filter_by(individual_id=individual.id).first()
				node_color = 0.8	#color for the sequenced
				
				if monkey_id2coverage:
					coverage = monkey_id2coverage.get(individual.code)
				else:
					coverage = None
				if coverage is None and individual_sequence is not None:
					coverage = individual_sequence.coverage
				
				if coverage is None:
					node_size = baseNodeSize
					node_color = 0	#not sequenced in a different color
				else:
					node_size = baseNodeSize
					if individual.id==1:	#the reference
						node_size = baseNodeSize*108
						node_size = baseNodeSize*4
					"""
					elif coverage<2:
						node_size = baseNodeSize*3
					elif coverage<8:
						node_size = baseNodeSize*12
					elif coverage>=20:
						node_size = baseNodeSize*36
					else:
						node_size = baseNodeSize
					"""
				
				if individual.vrc_founder:
					node_color = 0.25
					if individual_sequence is None:
						node_size = baseNodeSize*12
				
				sex2node_property_list[sex].append((v, node_size, node_color))
			sys.stderr.write("Done.\n")
		
		#nx.draw_circular(DG,with_labels=False, alpha=0.5)
		pylab.clf()
		pylab.axis("off")
		pylab.figure(figsize=(100, 60))
		layout = 'dot'
		pos = nx.graphviz_layout(DG, prog=layout)
		nx.draw_networkx_edges(DG, pos, edgelist=DG.edges(), \
							alpha=0.2, width=2.0,\
							arrows=False)
		G = DG.to_undirected()
		
		def drawOutlierEdge(DG, pos=None, monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None):
			phenotype_cmap = mpl.cm.jet
			max_phenotype = maxEdgeColor
			min_phenotype = minEdgeColor
			phenotype_gap = max_phenotype - min_phenotype
			phenotype_jitter = phenotype_gap/10.
			phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
			axe_x_offset4 = 0
			axe_y_offset1 = 0.9
			axe_width4 = 0.1
			axe_height4 = 0.2
			
			sys.stderr.write("Drawing outlier edges  ...")
			#2012.9.4 draw the outlier as additional edges
			edgelist = []
			edge_color = []
			for monkeyPair, data in monkeyPair2data.iteritems():
				dbID1 = db_vervet.getIndividualDBEntry(ucla_id=monkeyPair[0]).id
				dbID2 = db_vervet.getIndividualDBEntry(ucla_id=monkeyPair[1]).id
				edge_color.append(data.delta)
				edgelist.append([dbID1, dbID2])
			nx.draw_networkx_edges(DG, pos, edgelist=edgelist, edge_color=edge_color, edge_cmap=phenotype_cmap, \
								edge_vmin=min_phenotype, edge_vmax=max_phenotype, \
								alpha=1.0, width=0.5, style='solid',\
								arrows=False)
			sys.stderr.write(".\n")
		
		def drawEdgeColorLegend(DG, pos=None, monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None):
			phenotype_cmap = mpl.cm.jet
			max_phenotype = maxEdgeColor
			min_phenotype = minEdgeColor
			phenotype_gap = max_phenotype - min_phenotype
			phenotype_jitter = phenotype_gap/10.
			phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
			axe_x_offset4 = 0
			axe_y_offset1 = 0.9
			axe_width4 = 0.1
			axe_height4 = 0.2
			#2012.9.4 draw the legend
			sys.stderr.write("Drawing edge color legend ...")
			axe_map_phenotype_legend = pylab.axes([axe_x_offset4+0.02, axe_y_offset1, axe_width4-0.02, axe_height4], frameon=False)
			cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
										norm=phenotype_norm,
										orientation='horizontal')
			cb.set_label('Legend for the edge color, IBD-kinship')
			sys.stderr.write(".\n")
		
		def drawGraphNodes(G, pos, sex2node_property_list):
			import matplotlib as mpl
			for sex, node_property_ls in sex2node_property_list.iteritems():
				if sex=='M':
					node_shape = 's'
				else:
					node_shape = 'o'
				node_list = []
				node_color_list = []
				node_size_list = []
				for node_property in node_property_ls:
					node, node_size, node_color = node_property[:3]
					node_list.append(node)
					node_size_list.append(node_size)
					node_color_list.append(node_color)
				node_size_ar = numpy.array(node_size_list)
				if sex=='M':
					node_size_ar = node_size_ar*1.5	#by default, the square is smaller than a circle icon.
				
				nx.draw_networkx_nodes(DG, pos, nodelist=node_list, node_color=node_color_list, node_size=node_size_ar, \
									node_shape=node_shape, alpha=0.5, width=0, linewidths=0, cmap =mpl.cm.jet, vmin=0, vmax=1.0)
		
		drawGraphNodes(DG, pos, sex2node_property_list)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			drawOutlierEdge(DG, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			#draw outlier edges in the end because the new axes() would change the current figure in pylab
			drawEdgeColorLegend(DG, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=100)
		
		"""
		###spectral layout output
		pylab.clf()
		pylab.axis("off")
		pylab.figure(figsize=(60, 60))
		pos = nx.spectral_layout(G)
		nx.draw_networkx_edges(G, pos, alpha=0.9, width=0.5)
		
		#nx.draw_networkx_nodes(G, pos, node_color='r', node_size=baseNodeSize, alpha=0.8, width=0, linewidths=0)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		drawGraphNodes(G, pos, sex2node_property_list)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			drawOutlierEdge(G, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			#draw outlier edges in the end because the new axes() would change the current figure in pylab
			drawEdgeColorLegend(G, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		pylab.savefig('%s_spectral_layout_graph.png'%(outputFnamePrefix), dpi=150)
		"""
		"""
		layout = 'twopi'
		pylab.clf()
		nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=300)
		
		layout = 'circo'
		pylab.clf()
		nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=300)
		"""
	
	"""
		
		#2011-5-6
		outputFnamePrefix = os.path.expanduser('~/script/vervet/data/pedigree')
		DBVervet.drawPedigree(db_vervet, outputFnamePrefix=outputFnamePrefix)
		sys.exit(0)
		
		#2012.2.10
		outputFnamePrefix = os.path.expanduser('~/script/vervet/data/pedigree_723Highlight')
		monkeyCoverageFname = "/tmp/723MonkeysForWGS_to_be_sequenced.tsv"
		DBVervet.drawPedigree(db_vervet, outputFnamePrefix=outputFnamePrefix, monkeyCoverageFname=monkeyCoverageFname)
		sys.exit(0)
		
	"""
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#if self.monkeyCoverageFname and os.path.isfile(self.monkeyCoverageFname):
		#	monkey_id2coverage = cls.getMonkeyID2Coverage(self.monkeyCoverageFname)
		#else:
		monkey_id2coverage = {}
		
		outlierData = self.getOutlierPair(inputFname=self.kinshipIBDOutlierFname, minAbsDeltaForOutlier=self.minAbsDeltaForOutlier)
		self.drawPedigree(db_vervet=self.db_vervet, outputFnamePrefix=self.outputFnamePrefix, baseNodeSize=40, \
						monkey_id2coverage=monkey_id2coverage, \
						monkeyPair2data=outlierData.monkeyPair2data, minEdgeColor=outlierData.minDelta, \
						maxEdgeColor=outlierData.maxDelta)
		
		for monkeyID, monkeyPair2data in outlierData.monkeyID2monkeyPair2data.iteritems():
			noOfOutliers = len(monkeyPair2data)
			if noOfOutliers>=10:
				outputFnamePrefix = '%s_%soutlier_%s'%(self.outputFnamePrefix, noOfOutliers, monkeyID)
				self.drawPedigree(db_vervet=self.db_vervet, outputFnamePrefix=outputFnamePrefix, baseNodeSize=40, \
						monkey_id2coverage=monkey_id2coverage, \
						monkeyPair2data=monkeyPair2data, minEdgeColor=outlierData.minDelta, \
						maxEdgeColor=outlierData.maxDelta)
		

if __name__ == '__main__':
	main_class = PlotPedigree
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
