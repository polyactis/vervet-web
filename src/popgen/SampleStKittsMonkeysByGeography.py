#!/usr/bin/env python
"""
Examples:
	%s -i
	
	#2012.11.27
	%s  -i /tmp/StKittsMonkeys.txt -O /tmp/StKittsMonkeys_sampling  -u yh --db_passwd secret --maxDist 1
		 -z uclaOffice   --noOfMonkeysToChoose 60

Description:
	2012.11.26 The goal is to maximize the unrelatedness among the final set.
		The way to do that without a priori genetic information is to base the sampling strategy on the migratory distance.
		For example, since there is a mountain range in the middle of St Kitts which is not cross-able for the monkeys (to our knowledge),
			the migratory distance should be measured around the coast of St Kitts, reflected in the sampling (check attached figure).
		So first, I'm going to construct a graph of sites to recapitulate the migratory path.
		Within this graph, each node is a site and two sites are connected if their direct geographic distance is within 2km 
			(a quarter of the width of the island). Sites across the island won't be connected in this graph.
		The migratory distance between any two sites is the cumulative geographic distance along their shortest path in this graph.
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
#yh_matplotlib.setDefaultFigureSize((100,60))
#yh_matplotlib.setFontAndLabelSize(70)

import csv
import numpy, random, pylab
import networkx as nx
import matplotlib as mpl


from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, SNPData
from pymodule import yh_matplotlib, GenomeDB, utils
from pymodule import MatrixFile
from pymodule import RBDict
from pymodule import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper

class SampleStKittsMonkeysByGeography(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.update({
						('minShortestDistance', 1, float): [0.4, '', 1, 'minimum shortest distance allowed for new monkeys to be sampled'], \
						('maxDist', 1, float): [1, '', 1, 'max distance for two sites to be connected in graph. unit in kilometers.'], \
						('maxLongitude', 1, float): [-62.672, '', 1, 'longitude of newly sampled monkeys be <=maxLongitude.\
			used to exclude monkeys on the SE peninsula of StKitts, which disrupts the graph representation .'], \
						('sequencedMonkeyCountryIDList', 1, ): ["144,148", '', 1, 'both St. Kitts and Nevis'], \
						('newSampleMonkeyCountryIDList', 1, ): ["144", '', 1, 'only from st Kitts'], \
						('noOfMonkeysToChoose', 1, int): [60, '', 1, 'total number of monkeys to choose (+chosen ones, + sequenced ones)'], \
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.sequencedMonkeyCountryIDList = utils.getListOutOfStr(self.sequencedMonkeyCountryIDList)
		self.newSampleMonkeyCountryIDList = utils.getListOutOfStr(self.newSampleMonkeyCountryIDList)
	
	def readInSequencedMonkeys(self, db_vervet=None, countryIDList=None):
		"""
		2012.11.26
		"""
		sys.stderr.write("Fetching sequenced monkeys of country %s from db ..."%\
							(utils.getStrOutOfList(countryIDList)))
		
		query = db_vervet.metadata.bind.execute("select * from view_individual_sequence where country_id in (%s)"%\
								(utils.getStrOutOfList(countryIDList)))
		monkey_UCLAID_set = set()
		for row in query:
			monkey_UCLAID_set.add(row.ucla_id)
		
		sys.stderr.write("%s monkeys.\n"%(len(monkey_UCLAID_set)))
		return monkey_UCLAID_set
	
	def readInChosenOnes(self, inputFname=None):
		"""
		2012.11.26
		"""
		sys.stderr.write("Reading the chosen monkeys from %s ... "%(inputFname))
		reader = MatrixFile(inputFname)
		chosenMonkeyIDSet = set()
		for row in reader:
			chosenMonkeyIDSet.add(row[0])
		del reader
		sys.stderr.write("%s monkeys.\n"%(len(chosenMonkeyIDSet)))
		return chosenMonkeyIDSet
	
	def readInMonkeysFromDB(self, db_vervet=None, countryIDList=None, maxLongitude=None):
		"""
		2012.12.26
		"""
		sys.stderr.write("Fetching all monkeys of country %s from db , maxLongitude = %s ..."%\
						(utils.getStrOutOfList(countryIDList), maxLongitude))
		
		monkeyID2Info = {}
		
		query_string = "select * from view_individual where country_id in (%s) and latitude is not null and longitude is not null "%\
				(utils.getStrOutOfList(countryIDList))
		
		if maxLongitude is not None:
			query_string += " and longitude<=%s "%(maxLongitude)
		
		query = db_vervet.metadata.bind.execute(query_string)
		for row in query:
			monkeyID2Info[row.ucla_id] = PassingData(latitude=row.latitude, longitude=row.longitude, sex=row.sex, db_id=row.id,\
												country=row.country, site_name=row.site_name, \
												alignment_id=None, alignment_depth=None)
			properAlignment = db_vervet.getProperAlignmentGivenIndividualID(ucla_id=row.ucla_id)
			if properAlignment:
				monkeyID2Info[row.ucla_id].alignment_id = properAlignment.alignment_id
				monkeyID2Info[row.ucla_id].alignment_depth = properAlignment.median_depth
		sys.stderr.write("%s monkeys.\n"%(len(monkeyID2Info)))
		return monkeyID2Info
	
	def calculateShortestGeographicDistanceBetweenTwoComponents(self, cc1=None, cc2=None, monkeyID2Info=None, \
															monkeyPair2Distance=None):
		"""
		2012.11.27
		"""
		shortest_pair = None
		shortest_dist = None
		for i in xrange(len(cc1)):
			monkey1ID = cc1[i]
			monkey1Info = monkeyID2Info.get(monkey1ID)
			for j in xrange(len(cc2)):
				monkey2ID = cc2[j]
				monkey2Info = monkeyID2Info.get(monkey2ID)
				monkeyIDPairLs = [monkey1ID, monkey2ID]
				monkeyIDPairLs.sort()
				monkeyIDPair = tuple(monkeyIDPairLs)
				distance = monkeyPair2Distance.get(monkeyIDPair)
				#distance = utils.calGreatCircleDistance(lat1=monkey1Info.latitude, lon1=monkey1Info.longitude, \
				#									lat2=monkey2Info.latitude, lon2=monkey2Info.longitude, earth_radius=6372.795)
				if shortest_dist is None or distance<shortest_dist:
					shortest_dist = distance
					shortest_pair = (monkey1ID, monkey2ID)
		return PassingData(shortest_pair=shortest_pair, shortest_dist=shortest_dist)
	
	def constructNeighborGraph(self, monkeyID2Info=None, maxDist=2):
		"""
		2012.11.26
		
		"""
		sys.stderr.write("Constructing geography-neighbor graph for %s monkeys..."%(len(monkeyID2Info)))
		graph = nx.Graph()
		monkeyIDList = monkeyID2Info.keys()
		monkeyPair2Distance = {}	#record keeping
		for i in xrange(len(monkeyIDList)):
			monkey1ID = monkeyIDList[i]
			monkey1Info = monkeyID2Info.get(monkey1ID)
			graph.add_node(monkey1ID)
			
			for j in xrange(i+1, len(monkeyIDList)):
				monkey2ID = monkeyIDList[j]
				monkey2Info = monkeyID2Info.get(monkey2ID)
				distance = utils.calGreatCircleDistance(lat1=monkey1Info.latitude, lon1=monkey1Info.longitude, \
													lat2=monkey2Info.latitude, lon2=monkey2Info.longitude, earth_radius=6372.795)
				monkeyIDPairLs = [monkey1ID, monkey2ID]
				monkeyIDPairLs.sort()
				monkeyIDPair = tuple(monkeyIDPairLs)
				if monkeyIDPair not in monkeyPair2Distance:
					monkeyPair2Distance[monkeyIDPair] = distance
				
				if distance<=maxDist:
					graph.add_edge(monkey1ID, monkey2ID, weight=distance)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(graph.number_of_nodes(), graph.size(), \
															nx.number_connected_components(graph)))
		
		sys.stderr.write("Adding extra edges to make the graph one connected component ...")
		#for each component, find the two components that are closest, and connect the components
		#by connecting, it's connecting the two nodes from two components that are closest.
		# supposedly, the orientation of new edges on the geographic map should help to decide
		# whether some components should only have one neighbor or two (i.e. one CC could be from a narrow isthmus).
		cc_ls = nx.connected_components(graph)
		cc_graph = nx.Graph()	#a graph of CCs. each node is index to one CC, edge's weight is shortest distance in graph.
		# edge's other attribute shortest_pair is the two nodes that are shortest.
		#calculate pairwise distance among all CCs 
		ccIndex2ShortestCCList = {}
		for i in xrange(len(cc_ls)):
			cc1 = cc_ls[i]
			cc_graph.add_node(i)
			if i not in ccIndex2ShortestCCList:
				ccIndex2ShortestCCList[i] = []
			for j in xrange(i+1, len(cc_ls)):
				cc2 = cc_ls[j]
				cc_graph.add_node(j)
				
				shortestData = self.calculateShortestGeographicDistanceBetweenTwoComponents(cc1=cc1, cc2=cc2, \
													monkeyID2Info=monkeyID2Info, monkeyPair2Distance=monkeyPair2Distance)
				cc_graph.add_edge(i, j, weight=shortestData.shortest_dist, shortest_pair = shortestData.shortest_pair)
				ccIndex2ShortestCCList[i].append((shortestData.shortest_dist, j))
				if j not in ccIndex2ShortestCCList:
					ccIndex2ShortestCCList[j] = []
				ccIndex2ShortestCCList[j].append((shortestData.shortest_dist, i))
		
		counter = 0
		#for each CC, find its two closest CCs and connect them in graph
		for ccIndex in ccIndex2ShortestCCList:
			ccIndex2ShortestCCList[ccIndex].sort()	#sort all other CCs by shortest distance
			for i in xrange(2):
				shortest_dist, cc2Index = ccIndex2ShortestCCList[ccIndex][i]
				shortest_pair = cc_graph[ccIndex][cc2Index]['shortest_pair']
				graph.add_edge(shortest_pair[0], shortest_pair[1], weight=shortest_dist)
				counter += 1
		sys.stderr.write("%s edges added.\n"%(counter))
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(graph.number_of_nodes(), graph.size(), \
															nx.number_connected_components(graph)))
		
		#trim the graph, get rid of extra long edges of nodes which already have >=4 edges
		
		#..
		return graph
	
	def getPairwiseDistanceWithinGraphOfChosenMonkey(self, graph=None, chosenMonkeyIDDict=None):
		"""
		2012.11.27
		"""
		sys.stderr.write("Getting list of pairwise distance of all pairs in graph for %s chosen monkeys ..."%\
						(len(chosenMonkeyIDDict)))
		distance_ls = []
		chosenMonkeyIDList = list(chosenMonkeyIDDict)
		for i in xrange(len(chosenMonkeyIDList)):
			monkey1ID = chosenMonkeyIDList[i]
			for j in xrange(i+1, len(chosenMonkeyIDList)):
				monkey2ID = chosenMonkeyIDList[j]
				try:
					distance = nx.astar_path_length(graph, source=monkey1ID, target=monkey2ID, weight='weight')
				except:
					distance = -5
				distance_ls.append(distance)
				#distance_ls.append(lengthStructure[monkey1ID][monkey2ID])
		sys.stderr.write("%s pairs .\n"%(len(distance_ls)))
		return distance_ls
	
	def constructNewMonkeyToChosenSetDistanceVector(self, graph=None, preChosenMonkeyIDSet=None, minShortestDistance=0.4):
		"""
		2012.11.26
			each monkey is assigned a probability mass based on its geographic distance to the closest monkey
				in the preChosenMonkeyIDSet.
		"""
		sys.stderr.write("Constructing distance vector from new monkey to %s chosen monkeys (%s total monkeys) ..."%\
						(len(preChosenMonkeyIDSet), len(graph)))
		counter = 0
		real_counter = 0
		
		spanStartPos = 0
		unChosenMonkeyIDSet = set(graph.nodes()) - preChosenMonkeyIDSet
		preChosenAndInGraphMonkeyIDSet = set(graph.nodes()) - unChosenMonkeyIDSet
		
		returnData = PassingData(shortestDistanceToChosenSet_monkeyID_ls = [])
		for monkey1ID in unChosenMonkeyIDSet:
			counter += 1
			shortestDistance = None
			shortestDistanceToThisChosenMonkey = None
			for monkey2ID in preChosenAndInGraphMonkeyIDSet:
				#get the shortest path
				#distance = nx.shortest_path_length(graph, source=monkey1ID, target=monkey2ID, weight='weight')
				distance = nx.astar_path_length(graph, source=monkey1ID, target=monkey2ID, weight='weight')
				#short path in graph, sum all the edge weights
				if shortestDistance is None or distance<shortestDistance:
					shortestDistance = distance
					shortestDistanceToThisChosenMonkey = monkey2ID
			if shortestDistance is not None and shortestDistance>=minShortestDistance:	#ignore monkeys from same sites
				returnData.shortestDistanceToChosenSet_monkeyID_ls.append((shortestDistance, monkey1ID, shortestDistanceToThisChosenMonkey))
				real_counter += 1
				spanStartPos += shortestDistance	#increase the distance.
		returnData.totalDistance = spanStartPos
		sys.stderr.write("%s (out of %s candidate) monkeys added into candidate vector, total distance to the chosen set is %s.\n"%\
						(real_counter, counter, spanStartPos))
		return returnData
	
	def updateNewMonkeyToChosenSetDistanceVector(self, graph=None, oldShortestDistanceVectorData=None, \
												newlyChosenMonkeyID=None, minShortestDistance=0.4):
		"""
		2012.11.26
			supplement to constructNewMonkeyToChosenSetDistanceVector()
		"""
		oldShortestDistanceToChosenSet_monkeyID_ls = oldShortestDistanceVectorData.shortestDistanceToChosenSet_monkeyID_ls
		
		sys.stderr.write("Updating the shortest-distance to chosen set vector (%s elements) because monkeys %s has been added into the chosen set ..."%\
						(len(oldShortestDistanceToChosenSet_monkeyID_ls), newlyChosenMonkeyID))
		
		returnData = PassingData(shortestDistanceToChosenSet_monkeyID_ls = [])
		
		counter = 0
		real_counter = 0
		spanStartPos = 0
		
		monkey2ID = newlyChosenMonkeyID
		for element in oldShortestDistanceToChosenSet_monkeyID_ls:
			shortestDistance, monkey1ID, shortestDistanceToThisChosenMonkey = element[:3]
			if monkey1ID!=monkey2ID:	#skip the chosen monkey
				counter += 1
				#get the shortest path
				#short path in graph, sum all the edge weights, A* star algorithm seems to be much faster (>100%) than default shortest_path.
				#distance = nx.shortest_path_length(graph, source=monkey1ID, target=monkey2ID, weight='weight')
				distance = nx.astar_path_length(graph, source=monkey1ID, target=monkey2ID, weight='weight')
				if distance<shortestDistance:
					shortestDistance = distance
					shortestDistanceToThisChosenMonkey = monkey2ID
					element = (shortestDistance, monkey1ID, shortestDistanceToThisChosenMonkey)
					real_counter += 1
				if shortestDistance>=minShortestDistance:
					returnData.shortestDistanceToChosenSet_monkeyID_ls.append(element)
					spanStartPos += shortestDistance	#increase the distance.
		returnData.totalDistance = spanStartPos
		sys.stderr.write("%s (out of %s) monkeys changed their shortest distance to the chosen set (%s elements). \n\
	total distance to the chosen set is %s.\n"%\
						(real_counter, counter, len(returnData.shortestDistanceToChosenSet_monkeyID_ls), spanStartPos))
		
		return returnData
	
	
	def chooseExtraSamples(self, graph=None, preChosenMonkeyIDSet=None, noOfMonkeysToChoose=60, shortestDistanceVectorData=None,\
						minShortestDistance=0.4):
		"""
		2012.11.26
			construct an RBTree with span=probability span (within 0-1).
			probability span is proportional to the shortest distance from one new monkey to any one in the preChosenMonkeyIDSet.
			
		"""
		#draw a random uniform variable
		#compareIns = CNVCompare(min_reciprocal_overlap=0.1)
		
		finalChosenMonkeyIDDict = {}
		for monkeyID in preChosenMonkeyIDSet:
			finalChosenMonkeyIDDict[monkeyID] = [-1,-1,-1]	#the value is a tuple of (shortestDistance, monkey1ID, shortestDistanceToThisChosenMonkey)
		
		
		while len(finalChosenMonkeyIDDict)<noOfMonkeysToChoose and len(shortestDistanceVectorData.shortestDistanceToChosenSet_monkeyID_ls)>0:
			#segmentKey = CNVSegmentBinarySearchTreeKey(chromosome="1", \
			#					span_ls=[u*probabilitySpanRBDict.totalDistance, u*probabilitySpanRBDict.totalDistance],\
			#					min_reciprocal_overlap=1)
			shortestDistanceToChosenSet_monkeyID_ls = shortestDistanceVectorData.shortestDistanceToChosenSet_monkeyID_ls
			shortestDistanceToChosenSet_monkeyID_ls.sort()
			element = shortestDistanceToChosenSet_monkeyID_ls[-1]
			
			newlyChosenMonkeyID = element[1]
			finalChosenMonkeyIDDict[newlyChosenMonkeyID] = element
			shortestDistanceVectorData = self.updateNewMonkeyToChosenSetDistanceVector(graph=graph, \
						oldShortestDistanceVectorData=shortestDistanceVectorData, newlyChosenMonkeyID=newlyChosenMonkeyID, \
						minShortestDistance=minShortestDistance)
		return finalChosenMonkeyIDDict
	
	def outputChosenMonkeys(self, monkeyID2Info=None, chosenMonkeyIDDict=None, outputFname=None):
		"""
		2012.11.26
			output include:
				monkeyID, dbID, alignmentDepth, country
		"""
		sys.stderr.write("Outputting %s chosen monkeys to %s ...  "%(len(chosenMonkeyIDDict), outputFname))
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		header = ['monkey_id', 'sex', 'latitude', 'longitude', 'site_name','country', 'alignment_id', 'alignment_depth', \
				'shortestDistanceToChosenSet', 'shortestDistanceToThisChosenMonkey']
		writer.writerow(header)
		for ucla_id, element  in chosenMonkeyIDDict.iteritems():
			monkeyInfo = monkeyID2Info.get(ucla_id)
			row = [ucla_id, monkeyInfo.sex, monkeyInfo.latitude, monkeyInfo.longitude, monkeyInfo.site_name, monkeyInfo.country, monkeyInfo.alignment_id,\
				monkeyInfo.alignment_depth, element[0], element[2]]
			writer.writerow(row)
		del writer
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
		
		sequencedMonkeyIDSet = self.readInSequencedMonkeys(db_vervet, countryIDList=self.sequencedMonkeyCountryIDList)
		
		if self.inputFname:
			preChosenMonkeyIDSet = self.readInChosenOnes(self.inputFname)
		else:
			preChosenMonkeyIDSet = set()
		
		preChosenMonkeyIDSet |= sequencedMonkeyIDSet
		
		#monkey with latitudes/longitudes
		monkeyID2Info = self.readInMonkeysFromDB(db_vervet, countryIDList=self.newSampleMonkeyCountryIDList,\
								maxLongitude=self.maxLongitude)
		
		#allMonkeyID2Info is for output purpose
		allMonkeyID2Info = self.readInMonkeysFromDB(db_vervet, countryIDList=self.sequencedMonkeyCountryIDList,\
								maxLongitude=None)
		
		graph = self.constructNeighborGraph(monkeyID2Info, maxDist=self.maxDist)
		
		allMonkeyGraph = self.constructNeighborGraph(allMonkeyID2Info, maxDist=self.maxDist)
		
		"""
		#draw it and check how many monkeys have degree=1
		pos=nx.graphviz_layout(graph, prog="neato")
		#nx.draw_shell(graph)
		nx.draw(graph, pos, with_labels=True
			)
		#node_size=40,
		pylab.savefig('%s_graphNeatoLayout.png'%(self.outputFnamePrefix), dpi=150)
		"""
		
		shortestDistanceVectorData = self.constructNewMonkeyToChosenSetDistanceVector(graph=graph, preChosenMonkeyIDSet=preChosenMonkeyIDSet, \
												minShortestDistance=self.minShortestDistance)
		#probabilitySpanRBDict = self.constructNewMonkeyToChosenSetDistanceVector(graph=graph, preChosenMonkeyIDSet=preChosenMonkeyIDSet)
		#sampling for 10 times
		for i in xrange(1):
			finalChosenMonkeyIDDict = self.chooseExtraSamples(graph, preChosenMonkeyIDSet=preChosenMonkeyIDSet, \
												noOfMonkeysToChoose=self.noOfMonkeysToChoose, \
												shortestDistanceVectorData=shortestDistanceVectorData, \
												minShortestDistance=self.minShortestDistance)
			self.outputChosenMonkeys(monkeyID2Info=allMonkeyID2Info, chosenMonkeyIDDict=finalChosenMonkeyIDDict, \
									outputFname='%s_sample%s_%sMonkeys.tsv'%(self.outputFnamePrefix, i, len(finalChosenMonkeyIDDict)))
			distance_ls = self.getPairwiseDistanceWithinGraphOfChosenMonkey(graph=allMonkeyGraph, \
																		chosenMonkeyIDDict=finalChosenMonkeyIDDict)
			
			yh_matplotlib.drawHist(data_ls=distance_ls, title=None, \
					xlabel_1D="pairwise distance within graph", \
					xticks=None, outputFname='%s_sample%s_%sMonkeys_pairwise_distance_hist.png'%\
						(self.outputFnamePrefix, i, len(finalChosenMonkeyIDDict)), \
					min_no_of_data_points=10, needLog=True, \
					dpi=200, max_no_of_bins=40)
			# draw Histogram of pairwise shortest-path distance among the chosen ones  
		

if __name__ == '__main__':
	main_class = SampleStKittsMonkeysByGeography
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
