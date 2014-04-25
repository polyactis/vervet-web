#!/usr/bin/env python
"""
Examples:
	#2012.8.22 draw pedigree with focus on IBD-kinship outlier edges. 
	%s -i ~/NetworkData/vervet/Kinx2Aug2012.txt.gz -l LDPrunedMerged_ibdCheck.tsv
		-O ~/script/vervet/doc/figures/labelContaminationChecking/pedigreeWithOutlier/method22Outlier
		-u yh -z uclaOffice
	
	#2012.11.30 draw connected components among sequenced monkeys
	%s -i ~/NetworkData/vervet/Kinx2Aug2012.txt.gz -l LDPrunedMerged_ibdCheck.tsv
		-O ~/script/vervet/doc/figures/VRC_CC/Pedigree -z uclaOffice -u yh --db_passwd secret
		#these settings are for drawing the largest component with 662 nodes.
		--defaultFontLabelSize 30 --defaultNodeSize 2000 --defaultEdgeWidth 10 --defaultFigureWidth 400

Description:
	2012.9.4 program that draws co-ancestry pedigree graph of two individuals with extra edges describing (IBD-kinship) difference.
		outlier edges will be colored according to IBD-kinsip.
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
#set in the very beginning , otherwise, no effect.
from pymodule import yh_matplotlib

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, SNPData
from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB, utils
from pymodule import MatrixFile
from pymodule import SNP
import numpy, random, pylab
import numpy as np
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
import networkx as nx
import matplotlib as mpl
from DetectWrongLabelByCompKinshipVsIBD import DetectWrongLabelByCompKinshipVsIBD

class PlotPedigree(AbstractVervetMapper, DetectWrongLabelByCompKinshipVsIBD):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.update({
						('plinkIBDCheckOutputFname', 1, ): ["", 'l', 1, 'file that contains IBD check result, output of plink IBD'], \
						('kinshipFname', 1, ): [None, 'i', 1, 'the header-less 3-column tab/coma-delimited kinship filename.'], \
						('minAbsDeltaForOutlier', 1, float): [0.2, 'm', 1, 'if not 0, used to identify abs(kinship-IBD) outliers'], \
						('monkey1ID', 0, ): [None, 'M', 1, 'if both monkey1ID and monkey2ID are given, only plot ancestry graph for these two'], \
						('monkey2ID', 0, ): [None, 'N', 1, 'if both monkey1ID and monkey2ID are given, only plot ancestry graph for these two'], \
						('defaultFontLabelSize', 1, int): [70, '', 1, 'default font & label size on the plot'], \
						('defaultNodeSize', 1, int): [3500, '', 1, 'default node size in the pedigree plot'], \
						('defaultEdgeWidth', 1, float): [25, '', 1, 'default edge width in the pedigree plot'], \
						('defaultFigureWidth', 1, float): [120, '', 1, 'default figure width in the pedigree plot'], \
						('defaultFigureHeight', 1, float): [60, '', 1, 'default figure height in the pedigree plot'], \
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		if self.monkey1ID and self.monkey2ID:
			self.monkeyPairDesignated = [self.monkey1ID, self.monkey2ID]
			self.monkeyPairDesignated.sort()
			self.monkeyPairDesignated = tuple(self.monkeyPairDesignated)
		else:
			self.monkeyPairDesignated = None
		yh_matplotlib.setFontAndLabelSize(self.defaultFontLabelSize)
		yh_matplotlib.setDefaultFigureSize((self.defaultFigureWidth, self.defaultFigureHeight))
		
		#DetectWrongLabelByCompKinshipVsIBD.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
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
	
	def getOutlierPair(self, kinshipData=None, minAbsDeltaForOutlier=0.2, monkeyPair2IBDVector=None,\
					monkeyIDInIBDDataSet=None, monkeyPairDesignated=None):
		"""
		2012.9.10
			add argument monkeyPairDesignated, and read data from kinship file, rather than outlier output from DetectWrongLabelByCompKinshipVsIBD.py
			
			
		2012.9.4
		"""
		
		
		sys.stderr.write("Getting outlier pairs from kinshipData (%s X %s) and %s pairs of IBD vectors (%s monkeys) minAbsDeltaForOutlier=%s ..."%\
						(kinshipData.data_matrix.shape[0], kinshipData.data_matrix.shape[1], \
						len(monkeyPair2IBDVector), len(monkeyIDInIBDDataSet), minAbsDeltaForOutlier))
		
		sharedMonkeyIDSet = set(kinshipData.row_id_ls)&set(monkeyIDInIBDDataSet)
		monkey_id_ls = list(sharedMonkeyIDSet)
		monkey_id_ls.sort()
		
		monkeyPair2data  = {}
		counter = 0
		real_counter = 0
		minDelta = None
		maxDelta = None
		monkeyID2monkeyPair2data = {}
		row_id2row_index = kinshipData.row_id2row_index
		no_of_monkeys = len(monkey_id_ls)
		for i in xrange(no_of_monkeys):
			monkey1ID = monkey_id_ls[i]
			for j in xrange(i+1, no_of_monkeys):
				counter += 1
				monkey2ID = monkey_id_ls[j]
				kinship = kinshipData.getCellDataGivenRowColID(monkey1ID, monkey2ID)
				monkey_id_pair = [monkey1ID, monkey2ID]
				monkey_id_pair.sort()
				key = tuple(monkey_id_pair)
				IBD = monkeyPair2IBDVector.get(key).IBD
				if kinship is not None and not hasattr(kinship, 'mask') and IBD is not None:
					delta = IBD-kinship
					absDelta = abs(delta)
					if absDelta>=minAbsDeltaForOutlier or key==self.monkeyPairDesignated:
						if monkeyPair2IBDVector and key in monkeyPair2IBDVector:
							IBDVector = monkeyPair2IBDVector.get(key).IBDVector
							IBDVectorStr = monkeyPair2IBDVector.get(key).IBDVectorStr
						else:
							IBDVector = None
							IBDVectorStr = None
					
						data = PassingData(kinship=kinship, IBD=IBD, delta=delta, absDelta=absDelta, IBDVector=IBDVector, IBDVectorStr=IBDVectorStr)
						if key in monkeyPair2data:
							sys.stderr.write("WARNING: key %s has value %s in monkeyPair2data already. value overwritten with %s.\n"%\
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
				
		
		sys.stderr.write(" %s/%s outlier pairs, minDelta=%s, maxDelta=%s.\n"%(real_counter, counter, minDelta, maxDelta))
		return PassingData(monkeyPair2data=monkeyPair2data, minDelta=minDelta, maxDelta=maxDelta,\
						monkeyID2monkeyPair2data=monkeyID2monkeyPair2data)
	
	def getMonkeyPair2IBDVector(self, inputFname=None):
		"""
		2012.9.10
			return monkeyIDSet as well
		2012.9.6
		"""
		sys.stderr.write("Getting monkey pair 2 IBD vector from %s  ..."%(inputFname))
		reader = MatrixFile(inputFname)
		reader.constructColName2IndexFromHeader()
		monkey1IDIndex = reader.getColIndexGivenColHeader("IID1")
		monkey2IDIndex = reader.getColIndexGivenColHeader("IID2")
		IBDIndex = reader.getColIndexGivenColHeader("PI_HAT")
		Z0Index = reader.getColIndexGivenColHeader("Z0")
		Z1Index = reader.getColIndexGivenColHeader("Z1")
		Z2Index = reader.getColIndexGivenColHeader("Z2")
		formatFunc = lambda x: '%.2f'%(x)
		monkeyPair2IBDVector = {}
		counter = 0
		monkeyIDSet = set()
		for row in reader:
			monkey1ID = row[monkey1IDIndex]
			monkey2ID = row[monkey2IDIndex]
			monkey_id_pair = [monkey1ID, monkey2ID]
			monkey_id_pair.sort()
			key = tuple(monkey_id_pair)
			Z0 = float(row[Z0Index])
			Z1 = float(row[Z1Index])
			Z2 = float(row[Z2Index])
			IBD = float(row[IBDIndex])
			IBDVector = [Z0, Z1, Z2]
			IBDVector = map(formatFunc, IBDVector)
			IBDVectorStr = ','.join(IBDVector)
			data = PassingData(IBD=IBD, IBDVector=IBDVector, IBDVectorStr=IBDVectorStr)
			if key in monkeyPair2IBDVector:
				sys.stderr.write("WARNING: key %s has value %s in monkeyPair2IBDVector already. value overwritten with %s.\n"%\
									(repr(key), monkeyPair2IBDVector.get(key), data))
			monkeyPair2IBDVector[key] = data
			monkeyIDSet.add(monkey1ID)
			monkeyIDSet.add(monkey2ID)
			counter += 1
		sys.stderr.write(" %s pairs of IBD vectors for %s unique monkeys.\n"%(len(monkeyPair2IBDVector), len(monkeyIDSet)))
		return PassingData(monkeyPair2IBDVector=monkeyPair2IBDVector, monkeyIDSet=monkeyIDSet)
	
	def getSex2NodePropertyList(self, DG=None,db_vervet=None, monkey_id2coverage=None, baseNodeSize=1, \
							nodeSizeInProportionToCoverage=True):
		"""
		2012.9.6
			
		"""
		sys.stderr.write("Assigning each node with different size/color ...")
		sex2NodePropertyList = {}		# in node_property_list, each entry is (node, size, color)
		#size depends on whether it's deep-sequenced (30X), 4 for 30X, 8 for the REF, 2 for all others. 
		#color depends on it's sequenced or not
		for v in DG:
			individual = db_vervet.getIndividualDBEntryViaDBID(v)
			sex = individual.sex
			
			if individual.sex not in sex2NodePropertyList:
				sex2NodePropertyList[sex] = []
			individual_sequence = VervetDB.IndividualSequence.query.filter_by(sequencer_id=2).\
				filter_by(individual_id=individual.id).first()
			node_color = 0.8	#color for the sequenced
			
			if monkey_id2coverage:
				coverage = monkey_id2coverage.get(individual.code)
			else:
				coverage = None
			if coverage is None and individual_sequence is not None:
				coverage = individual_sequence.coverage
			
			if coverage is None:	#un-sequenced monkeys
				node_size = baseNodeSize
				node_color = 0	#not sequenced in a different color, 0 is blue. [0-1] depends on mpl.cm.jet.
				#2014.03.04 temporary
				node_color = 0.2
			else:
				node_size = baseNodeSize
				if individual.id==1:	#the reference
					node_size = baseNodeSize*8
					#node_size = baseNodeSize*2
				elif nodeSizeInProportionToCoverage:
					if coverage<2:
						node_size = baseNodeSize*1
					elif coverage<8:
						node_size = baseNodeSize*2
					elif coverage>=20:
						node_size = baseNodeSize*4
					else:
						node_size = baseNodeSize
				
			
			if individual.vrc_founder:
				node_color = 0.2	#[0-1] depends on mpl.cm.jet.
				if individual_sequence is None:
					node_size = baseNodeSize*2
					#2014.03.04 temporary
					node_size = baseNodeSize
					
			
			sex2NodePropertyList[sex].append((v, node_size, node_color))
		sys.stderr.write("Done.\n")
		return sex2NodePropertyList
	
			
	def drawOutlierEdge(self, DG=None, db_vervet=None, pos=None, monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None, alpha=0.5, \
					edgeWidth=1.0):
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
			if DG.has_node(dbID1) and DG.has_node(dbID2):	#must be in the graph already
				edge_color.append(data.delta)
				edgelist.append([dbID1, dbID2])
		nx.draw_networkx_edges(DG, pos, edgelist=edgelist, edge_color=edge_color, edge_cmap=phenotype_cmap, \
							edge_vmin=min_phenotype, edge_vmax=max_phenotype, \
							alpha=alpha, width=edgeWidth, style='solid',\
							arrows=False)
		sys.stderr.write(".\n")
	
	def drawEdgeColorLegend(self, DG=None, pos=None,  minEdgeColor=None, maxEdgeColor=None):
		phenotype_cmap = mpl.cm.jet
		max_phenotype = maxEdgeColor
		min_phenotype = minEdgeColor
		phenotype_gap = max_phenotype - min_phenotype
		phenotype_jitter = phenotype_gap/10.
		phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
		axe_x_offset4 = 0
		axe_y_offset1 = 0.97
		axe_width4 = 0.6
		axe_height4 = 0.03
		#2012.9.4 draw the legend
		sys.stderr.write("Drawing edge color legend ...")
		axe_map_phenotype_legend = pylab.axes([axe_x_offset4+0.02, axe_y_offset1, axe_width4-0.02, axe_height4], frameon=False)
		cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
									norm=phenotype_norm,
									orientation='horizontal')
		cb.set_label('Legend for the edge color, IBD-kinship')
		sys.stderr.write(".\n")
	
	def drawGraphNodes(self, G=None, pos=None, sex2NodePropertyList=None, alpha=0.9):
		"""
		2014.03.04, temporary, change alpha to 0.9
		
		default alpha=0.3
		"""
		import matplotlib as mpl
		for sex, node_property_ls in sex2NodePropertyList.iteritems():
			if sex=='M':
				node_shape = 's'
			else:
				node_shape = 'o'
			node_list = []
			node_color_list = []
			node_size_list = []
			for node_property in node_property_ls:
				node, node_size, node_color = node_property[:3]
				if G.has_node(node):
					node_list.append(node)
					node_size_list.append(node_size)
					node_color_list.append(node_color)
			node_size_ar = numpy.array(node_size_list)
			if sex=='M':
				node_size_ar = node_size_ar*1.5	#by default, the square is smaller than a circle icon.
			
			nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=node_color_list, node_size=node_size_ar, \
								node_shape=node_shape, alpha=alpha, width=0, linewidths=0, cmap =mpl.cm.jet, vmin=0, vmax=1.0)
	
	
	
	def drawPedigree(self, DG=None, db_vervet=None, outputFnamePrefix=None, monkey_id2coverage=None, \
					monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None, sex2NodePropertyList=None):
		"""
		2012.9.4
			copied from vervet/src/misc.py
		2012.2.10
			add argument monkeyCoverageFname to override coverage data from db
		2011-5-6
			
		"""
		
		"""
		sys.stderr.write("Plotting out degree histogram ....")
		out_degree_ls = []
		#for n in DG.nodes_iter():
		#	out_degree_ls.append(DG.out_degree(n))
		#outputFname = '%s_outDegreeHist.png'%(outputFnamePrefix)
		#yh_matplotlib.drawHist(out_degree_ls, title="Histogram of no. of children per monkey", xlabel_1D="no. of children", \
		#					outputFname=outputFname, min_no_of_data_points=50, needLog=True)
		"""
		
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
		self.drawGraphNodes(DG, pos, sex2NodePropertyList)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			self.drawOutlierEdge(DG, db_vervet=db_vervet, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
			#if monkeyPair2data and minEdgeColor and maxEdgeColor:
			#draw outlier edges in the end because the new axes() would change the current figure in pylab
			self.drawEdgeColorLegend(DG, pos=pos, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=25)
		
		"""
		###spectral layout output
		pylab.clf()
		pylab.axis("off")
		pylab.figure(figsize=(60, 60))
		pos = nx.spectral_layout(G)
		nx.draw_networkx_edges(G, pos, alpha=0.9, width=0.5)
		
		#nx.draw_networkx_nodes(G, pos, node_color='r', node_size=baseNodeSize, alpha=0.8, width=0, linewidths=0)
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		self.drawGraphNodes(G, pos, sex2NodePropertyList)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			self.drawOutlierEdge(G, db_vervet=db_vervet, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			#draw outlier edges in the end because the new axes() would change the current figure in pylab
			self.drawEdgeColorLegend(G, pos=pos, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
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
	
	def drawListOfMonkeysWithAncestorsAndDescendants(self, db_vervet=None, DG=None, reverseDG=None, \
				outputFnamePrefix=None, node_list=None, \
				sex2NodePropertyList=None, defaultEdgeWidth=25, addParentsForDescendants=False, \
				**keywords):
		"""
		2012.11.29
		"""
		sys.stderr.write("Drawing list of %s monkeys with ancestors and descendants ...  "%(len(node_list)))
		
		pylab.clf()
		pylab.axis("off")
		axe_pvalue = pylab.axes([-0, -0, 1, 0.93], frameon=False)	#left gap, bottom gap, width, height.
		pylab.figure(axe_pvalue.figure.number)#, figsize=(100, 60))	#figure size was set in the beginning of the program.
		
		#fig = matplotlib.pyplot.gcf()
		#fig.set_size_inches(185, 60)
		#pylab.figure(figsize=(100, 60))
		
		nodeSet = set(node_list)
		#add all ancestors
		for dbID in node_list:
			for node in db_vervet.getDescendantSet(dbID):
				nodeSet.add(node)
				if addParentsForDescendants:
					for parentNode in DG.predecessors(node):
						nodeSet.add(parentNode)
			for node in db_vervet.getAncestorSet(dbID):
				nodeSet.add(node)
			"""
			for node, bridgeList in nx.predecessor(reverseDG, source=dbID, target=None).iteritems():	#all ancestors of dbID.
				#predecessor in the path from dbID to all nodes in reverseDG
				#bridgeList only contains the the immediate predecessor to node, not the whole path from dbID to node.
				nodeSet.add(node)
			for node, bridgeList in nx.predecessor(DG, source=dbID, target=None).iteritems():	#all descendants of dbID.
				#predecessor in the path from dbID to all nodes in DG
				#bridgeList only contains the the immediate predecessor to node, not the whole path from dbID to node.
				nodeSet.add(node)
			"""
		nodeList = list(nodeSet)
		subDG = DG.subgraph(nodeList)
		sys.stderr.write(" %s nodes %s edges ...  "%(subDG.number_of_nodes(), subDG.number_of_edges()))
		if len(node_list)==1:
			individual = VervetDB.Individual.get(node_list[0])
			title = 'Monkey %s with ancestors&descendants (%s total)'%(individual.ucla_id, len(nodeList))
		else:
			title = '%s nodes with ancestors&descendants (%s total)'%(len(node_list), len(nodeList))
		
		pylab.title(title, fontsize=80)
		layout = 'dot'
		pos = nx.graphviz_layout(subDG, prog=layout)
		
		nx.draw_networkx_edges(subDG, pos, edgelist=subDG.edges(), \
							alpha=0.2, width=defaultEdgeWidth, style='dashed',\
							arrows=False)
		self.drawGraphNodes(subDG, pos, sex2NodePropertyList)
		
		labels = {}	#dicitonary to pass into draw_networkx_labels
		for n in subDG.nodes():
			individual = db_vervet.getIndividualDBEntryViaDBID(n)
			labels[n] = individual.code
		nx.draw_networkx_labels(subDG, pos, labels=labels, font_size=50, font_color='k', font_family='sans-serif', \
							font_weight='normal', alpha=0.5, ax=None)
		
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=30)
		sys.stderr.write(".\n")
	
	def drawListOfMonkeysAndTheirParents(self, db_vervet=None, DG=None, reverseDG=None, \
				outputFnamePrefix=None, node_list=None, \
				sex2NodePropertyList=None, defaultEdgeWidth=25,\
				**keywords):
		"""
		2012.11.29	
		"""
		sys.stderr.write("Drawing list of %s monkeys along with their parents ...  "%(len(node_list)))
		
		pylab.clf()
		pylab.axis("off")
		axe_pvalue = pylab.axes([-0, -0, 1, 0.93], frameon=False)	#left gap, bottom gap, width, height.
		pylab.figure(axe_pvalue.figure.number)#, figsize=(100, 60))	#figure size was set in the beginning of the program.
		
		#fig = matplotlib.pyplot.gcf()
		#fig.set_size_inches(185, 60)
		#pylab.figure(figsize=(100, 60))
		
		nodeSet = set(node_list)
		#add all ancestors
		for dbID in node_list:
			for parentNode in DG.predecessors(dbID):	#parents of each node
				nodeSet.add(parentNode)
		
		nodeList = list(nodeSet)
		subDG = DG.subgraph(nodeList)
		sys.stderr.write(" %s nodes %s edges ...  "%(subDG.number_of_nodes(), subDG.number_of_edges()))
		
		title = '%s nodes with parents (%s total)'%(len(node_list), len(nodeList))
		
		pylab.title(title, fontsize=80)
		layout = 'dot'
		pos = nx.graphviz_layout(subDG, prog=layout)
		
		nx.draw_networkx_edges(subDG, pos, edgelist=subDG.edges(), \
							alpha=0.2, width=defaultEdgeWidth, style='dashed',\
							arrows=False)
		self.drawGraphNodes(subDG, pos, sex2NodePropertyList)
		
		labels = {}	#dicitonary to pass into draw_networkx_labels
		for n in subDG.nodes():
			individual = db_vervet.getIndividualDBEntryViaDBID(n)
			labels[n] = individual.code
		nx.draw_networkx_labels(subDG, pos, labels=labels, font_size=50, font_color='k', font_family='sans-serif', \
							font_weight='normal', alpha=0.5, ax=None)
		
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=30)
		sys.stderr.write(".\n")
	
	def drawSubPedigree(self, db_vervet=None, DG=None, reverseDG=None, outputFnamePrefix=None, monkeyPair=None, monkeyPairData=None, \
					monkeyPair2data=None, minEdgeColor=None, maxEdgeColor=None, sex2NodePropertyList=None,  defaultEdgeWidth=25,\
					**keywords):
		"""
		2012.9.4
			copied from vervet/src/misc.py
		2012.2.10
			add argument monkeyCoverageFname to override coverage data from db
		2011-5-6
			
		"""
		sys.stderr.write("Drawing sub-pedigree ...  ")
		#nx.draw_circular(DG,with_labels=False, alpha=0.5)
		pylab.clf()
		pylab.axis("off")
		axe_pvalue = pylab.axes([-0, -0, 1, 0.93], frameon=False)	#left gap, bottom gap, width, height.
		pylab.figure(axe_pvalue.figure.number)#, figsize=(100, 60))	#figure size was set in the beginning of the program.
		
		#fig = matplotlib.pyplot.gcf()
		
		#fig.set_size_inches(185, 60)
		
		#pylab.figure(figsize=(100, 60))
		monkey1ID = monkeyPair[0]
		monkey2ID = monkeyPair[1]
		dbID1 = db_vervet.getIndividualDBEntry(ucla_id=monkey1ID).id
		dbID2 = db_vervet.getIndividualDBEntry(ucla_id=monkey2ID).id
		
		nodeSet = set([dbID1, dbID2])
		#add all ancestors
		for dbID in [dbID1, dbID2]:
			for node, bridgeList in nx.predecessor(reverseDG, source=dbID, target=None).iteritems():	#all ancestors of dbID.
				#predecessor in the path from dbID to all nodes in reverseDG
				#bridgeList only contains the the immediate predecessor to node, not the whole path from dbID to node.
				nodeSet.add(node)
			for node in DG.successors(dbID):	#children of dbID
				nodeSet.add(node)
				for parentNode in DG.predecessors(node):	#other parents of children of dbID
					nodeSet.add(parentNode)
		
		nodeList = list(nodeSet)
		subDG = DG.subgraph(nodeList)
		sys.stderr.write(" %s nodes %s edges ...  "%(subDG.number_of_nodes(), subDG.number_of_edges()))
		
		title = 'outlier: %s-%s, kinship=%.3f, IBD=%.3f, IBDVector=%s'%(monkey1ID, monkey2ID, monkeyPairData.kinship,\
														monkeyPairData.IBD, monkeyPairData.IBDVectorStr)
		pylab.title(title, fontsize=80)
		layout = 'dot'
		pos = nx.graphviz_layout(subDG, prog=layout)
		nx.draw_networkx_edges(subDG, pos, edgelist=subDG.edges(), \
							alpha=0.2, width=self.defaultEdgeWidth, style='dashed',\
							arrows=False)
		self.drawGraphNodes(subDG, pos, sex2NodePropertyList)
		
		labels = {}	#dicitonary to pass into draw_networkx_labels
		for n in subDG.nodes():
			individual = db_vervet.getIndividualDBEntryViaDBID(n)
			labels[n] = individual.code
		nx.draw_networkx_labels(subDG, pos, labels=labels, font_size=50, font_color='k', font_family='sans-serif', \
							font_weight='normal', alpha=0.5, ax=None)
		
		#nx.draw_graphviz(DG, prog=layout,with_labels=False, alpha=0.5)
		#draw outlier edges in the end because the new axes() would change the current figure in pylab
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			self.drawOutlierEdge(subDG, db_vervet=db_vervet, pos=pos, monkeyPair2data=monkeyPair2data, minEdgeColor=minEdgeColor, \
								maxEdgeColor=maxEdgeColor, alpha=0.6, edgeWidth=4.0)
		if monkeyPair2data and minEdgeColor and maxEdgeColor:
			self.drawEdgeColorLegend(subDG, pos=pos, minEdgeColor=minEdgeColor, maxEdgeColor=maxEdgeColor)
		
		pylab.savefig('%s_graphviz_%s_graph.png'%(outputFnamePrefix, layout), dpi=30)
		sys.stderr.write(".\n")
	
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#if self.monkeyCoverageFname and os.path.isfile(self.monkeyCoverageFname):
		#	monkey_id2coverage = cls.getMonkeyID2Coverage(self.monkeyCoverageFname)
		#else:
		db_vervet = self.db_vervet
		monkey_id2coverage = {}
		
		DG = self.db_vervet.constructPedgree()
		reverseDG = self.db_vervet.constructPedgree(directionType=2)
		sex2NodePropertyList = self.getSex2NodePropertyList(DG=DG, db_vervet=self.db_vervet, monkey_id2coverage=None, \
											baseNodeSize=self.defaultNodeSize,\
											nodeSizeInProportionToCoverage=True)
		#2014.03.04 temporarily draw teh whole pedigree
		self.drawPedigree(DG=DG, db_vervet=self.db_vervet, outputFnamePrefix=self.outputFnamePrefix, \
						monkey_id2coverage=monkey_id2coverage, \
						sex2NodePropertyList=sex2NodePropertyList)
		#2014.03.14 halts after whole pedigree is drawn
		sys.exit(0)
		
		###########  2012.11.30 draw particular monkey with its ancestors & descendants
		for monkeyID in ['1986085', '2005038', '1982001', '1984010', '1986007', '1986002', '1983087', '1978007', '1987011']:
			outputFnamePrefix = '%s_monkey_%s_with_ancestors_descendants'%(self.outputFnamePrefix,monkeyID)
			self.drawListOfMonkeysWithAncestorsAndDescendants(db_vervet=db_vervet, DG=DG, reverseDG=reverseDG, \
				outputFnamePrefix=outputFnamePrefix, node_list=[db_vervet.getIndividualDBEntry(ucla_id=monkeyID).id], \
				sex2NodePropertyList=sex2NodePropertyList, defaultEdgeWidth=self.defaultEdgeWidth, \
				addParentsForDescendants=False)
		sys.exit(0)
		
		######## 2012.11.29 draw connected components in sequenced VRC
		alignmentLs = db_vervet.getAlignments(ref_ind_seq_id=524, ind_seq_id_ls=None, ind_aln_id_ls=None,\
								alignment_method_id=2, individual_sequence_file_raw_id_type=0, excludeAlignmentWithoutLocalFile=False)
		alignmentLs = db_vervet.filterAlignments(alignmentLs, sequence_filtered=None, \
												individual_site_id_set=set([447]),\
												mask_genotype_method_id=None, parent_individual_alignment_id=None,\
									country_id_set=None, tax_id_set=None,\
									excludeContaminant=True)
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		sequencedDG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		i = 0
		for cc in nx.connected_components(sequencedDG.to_undirected()):
			i += 1
			if len(cc)>0:
				outputFnamePrefix = '%sCC%size%s'%(self.outputFnamePrefix, i, len(cc))
				self.drawListOfMonkeysAndTheirParents(db_vervet=db_vervet, DG=DG, reverseDG=reverseDG, \
									outputFnamePrefix=outputFnamePrefix, node_list=cc, \
									sex2NodePropertyList=sex2NodePropertyList, defaultEdgeWidth=self.defaultEdgeWidth)
		sys.exit(0)
		
		
		
		########  2012.11.29 below for outlier-edge drawing
		monkeyPair2IBDVectorData = self.getMonkeyPair2IBDVector(inputFname=self.plinkIBDCheckOutputFname)
		monkeyPair2IBDVector = monkeyPair2IBDVectorData.monkeyPair2IBDVector
		monkeyIDInIBDDataSet = monkeyPair2IBDVectorData.monkeyIDSet
		
		kinshipData = self.getMonkeyKinshipData(self.kinshipFname)
		outlierData = self.getOutlierPair(kinshipData=kinshipData, minAbsDeltaForOutlier=self.minAbsDeltaForOutlier,\
									monkeyPair2IBDVector=monkeyPair2IBDVector, monkeyIDInIBDDataSet=monkeyIDInIBDDataSet,\
									monkeyPairDesignated=self.monkeyPairDesignated)
		"""
		self.drawPedigree(DG=DG, db_vervet=self.db_vervet, outputFnamePrefix=self.outputFnamePrefix, \
						monkey_id2coverage=monkey_id2coverage, \
						monkeyPair2data=outlierData.monkeyPair2data, minEdgeColor=outlierData.minDelta, \
						maxEdgeColor=outlierData.maxDelta, sex2NodePropertyList=sex2NodePropertyList)
		"""
		
		#for monkeyID, monkeyPair2data in outlierData.monkeyID2monkeyPair2data.iteritems():
		for monkeyPair, monkeyPairData in outlierData.monkeyPair2data.iteritems():
			if self.monkeyPairDesignated and monkeyPair!=self.monkeyPairDesignated:
				continue
			monkeyPairList = list(monkeyPair)
			monkeyPairList.sort()
			monkey1ID = monkeyPairList[0]
			monkey2ID = monkeyPairList[1]
			dbID1 = self.db_vervet.getIndividualDBEntry(ucla_id=monkey1ID).id
			dbID2 = self.db_vervet.getIndividualDBEntry(ucla_id=monkey2ID).id
			if DG.has_node(dbID1) and DG.has_node(dbID2):
				outputFnamePrefix = '%s_outlier_%s_vs_%s_delta_%.3f'%(self.outputFnamePrefix, monkeyPairList[0], monkeyPairList[1], \
													monkeyPairData.delta)
				self.drawSubPedigree(db_vervet=self.db_vervet, DG=DG, reverseDG=reverseDG, outputFnamePrefix=outputFnamePrefix, \
					monkeyPair = monkeyPair, monkeyPairData=monkeyPairData, \
					monkeyPair2data=outlierData.monkeyPair2data, minEdgeColor=outlierData.minDelta, \
					maxEdgeColor=outlierData.maxDelta, sex2NodePropertyList=sex2NodePropertyList, defaultEdgeWidth=self.defaultEdgeWidth)
			"""
			noOfOutliers = len(monkeyPair2data)
			if noOfOutliers>=10:
				outputFnamePrefix = '%s_%soutlier_%s'%(self.outputFnamePrefix, noOfOutliers, monkeyID)
				self.drawPedigree(DG=DG, db_vervet=self.db_vervet, outputFnamePrefix=outputFnamePrefix, \
						monkey_id2coverage=monkey_id2coverage, \
						monkeyPair2data=monkeyPair2data, minEdgeColor=outlierData.minDelta, \
						maxEdgeColor=outlierData.maxDelta, sex2NodePropertyList=sex2NodePropertyList)
			"""

if __name__ == '__main__':
	main_class = PlotPedigree
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
