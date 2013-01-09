#!/usr/bin/env python
"""
Examples:
	#for plink
	%s -i AlignmentToCall_AllVRC_vs_524_AllContig195Sites.2012.8.5T1549/samtools/Contig195.vcf.gz
		-o /tmp/723VRC.tfam -u yh

	#2013.1.3 for TrioCaller
	%s -i /Network/Data/vervet/db/genotype_file/method_25/25999_VCF_24115_VCF_Contig520.filterByMaxSNPMissingRate.recode.vcf.gz
		--outputFileFormat 2 -u yh -o /tmp/trioCaller_VRC.merlin --sampleID2FamilyCountFname=/tmp/sampleID2familyCount_y2.tsv
	
	#2013.1.3 output for polymutt
	%s -i /Network/Data/vervet/db/genotype_file/method_25/25999_VCF_24115_VCF_Contig520.filterByMaxSNPMissingRate.recode.vcf.gz
		--outputFileFormat 3 -u yh -o /tmp/polymutt_VRC.merlin --polymuttDatFname /tmp/datfile
		--sampleID2FamilyCountFname=/tmp/sampleID2familyCount_y3.tsv --db_passwd secret
	
Description:
	2012.8.9
		outputs VRC pedigrees (only the monkeys in inputFname/vcf format) in plink tfam format.
		parse the order of monkeys from VCF file (inputFname) and preserve it in TFAM file.
	six columns:
	 Family ID
	 Individual ID
	 Paternal ID
	 Maternal ID
	 Sex (1=male; 2=female; other=unknown)
	 Phenotype
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import heapq, copy
import networkx as nx
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, NextGenSeq, Genome
from pymodule import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB

class OutputVRCPedigreeInTFAMGivenOrderFromFile(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						('outputFileFormat', 1, int ):[1, '', 1, '1: for plink; 2: for TrioCaller; \n\
	3: for polymutt, replicate certain individuals to make pedigree loop-free'],\
						('sampleID2FamilyCountFname', 0, ): ['', '', 1, 'a tab-delimited file that records how many families in which each individual occurs'],\
						("treatEveryOneIndependent", 0, int): [0, '', 0, 'toggle this to treat everyone in the pedigree independent (parents=0)'],\
						('replicateIndividualTag', 1, ): ['copy', 'T', 1, 'the tag that separates the true ID and its replicate count'],\
						('polymuttDatFname', 0, ): ['', '', 1, 'if present, add "T\tGLF_Index" into the file. required for polymutt'],\
						})
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getIndividual(self, db_vervet=None, individual_id=None, individual_id2individual=None):
		"""
		2012.8.14
		"""
		if isinstance(individual_id, tuple):
			individual_id = individual_id[0]
		else:
			individual_id = individual_id
		if individual_id2individual:
			individual = individual_id2individual.get(individual_id)
			if individual:
				return individual
			else:
				individual = VervetDB.Individual.get(individual_id)
				individual_id2individual[individual_id] = individual
				return individual
		else:
			individual = VervetDB.Individual.get(individual_id)
			return individual
	
	def outputPedigreeForPlink(self, DG=None, db_vervet=None, inputFname=None, outputFname=None, treatEveryOneIndependent=None):
		"""
		2013.1.2
			copied from run()
		"""
		alignmentLs = db_vervet.getAlignmentsFromVCFFile(inputFname =inputFname)
		"""
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		DG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		"""
		individual_id2individual = {}
		
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		counter = 0
		family_id= 1	#all in one family
		for alignment in alignmentLs:
			nodeID = alignment.individual_sequence.individual_id
			individual = self.getIndividual(db_vervet=db_vervet, individual_id=nodeID, \
										individual_id2individual=individual_id2individual)
			if treatEveryOneIndependent:
				father_id = 0
				mother_id = 0
			elif nodeID in DG:
				parents = DG.predecessors(nodeID)
				if len(parents)==2:
					parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
												individual_id2individual=individual_id2individual)
					parent2 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[1], \
												individual_id2individual=individual_id2individual)
					parent1Sex = parent1.codeSexInNumber()
					if parent1Sex==2:
						#swap the father and mother row
						tmp = parent1
						parent1 = parent2
						parent2 = tmp
					
					father_id = parent1.ucla_id
					mother_id = parent2.ucla_id
				elif len(parents)==1:
					parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
										individual_id2individual=individual_id2individual)
					parent1Sex = parent1.codeSexInNumber()
					if parent1Sex==2:
						father_id = 0
						mother_id = parent1.ucla_id
					else:
						father_id = parent1.ucla_id
						mother_id = 0
				elif len(parents)==0:
					father_id = 0
					mother_id = 0
				else:
					sys.stderr.write("Error: number of parents (%s) for %s is %s.\n"%(repr(parents), nodeID, len(parents)))
					sys.exit(3)
			else:
				father_id = 0
				mother_id = 0
			data_row = [family_id, individual.ucla_id, father_id, mother_id, \
					individual.codeSexInNumber(), 1]
			writer.writerow(data_row)
			counter += 1
		sys.stderr.write("%s alignments outputted.\n"%(counter))
		del writer
		
	
	def writeOneLineToPedigreeFile(self, writer=None, family_id=None, individual_alignment=None, father_id=0, mother_id=0,
								sex=None, replicateIndividualTag='copy', sampleID2FamilyCount=None):
		"""
		2012.3.29
		2011.12.5
			used by outputPedigreeForTrioCaller()
		"""
		if type(individual_alignment)==str:
			sample_id = individual_alignment
			sexByGuess = 'x'
		else:
			#sample_id = individual_alignment.getReadGroup()
			sample_id = self.getSampleIDWithReplicateCount(individual_alignment, replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
			sexByGuess = individual_alignment.individual_sequence.individual.codeSexInNumber()
		if sex is None:
			sex = sexByGuess
		
		data_row = [family_id, sample_id, father_id, mother_id, sex]
		writer.writerow(data_row)
	
	def getSampleIDWithReplicateCount(self, alignment, replicateIndividualTag=None, sampleID2FamilyCount=None):
		"""
		2012.3.29
			use alignment.getReadGroup() + someReplicateIndicator to identify each VCF column
		"""
		hashedID = '%s%s%s'%(alignment.getReadGroup(), replicateIndividualTag, \
							sampleID2FamilyCount.get(alignment.getReadGroup()))
		return hashedID

	def adjustSampleID2FamilyCountFromAlignmentDBEntry(self, alignment, sampleID2FamilyCount=None):
		"""
		2012.3.29
		"""
		sample_id = alignment.getReadGroup()
		if sample_id not in sampleID2FamilyCount:
			sampleID2FamilyCount[sample_id] = 0
		sampleID2FamilyCount[sample_id] += 1
		
	def outputPedigreeForTrioCaller(self, db_vervet=None, inputFname=None, pedigreeOutputFname=None, \
								replicateIndividualTag='copy',\
								treatEveryOneIndependent=False, sampleID2FamilyCountFname=None):
		"""
		2013.1.2 moved from AlignmentToTrioCallPipeline.py. merlin format.
		2012.9.26 output is wrong if some members of duo or trio have >1 alignments
			need to fix this bug 
		2012.9.26 deal with the case that one individual appear in multiple alignments
		2012.8.28
			add argument treatEveryOneIndependent, which removes the family structure and also no replicates
		2012.3.29
			every individual gets replicated by the number of families it appears in.
		2011-11-27
			
		2011-11-3
			find all trios, a list of "father,mother,child", all of which are identified by IndividualSequence.id.
			
			for each individual_alignment, -> ind_seq -> individual
				get its parents from ind2ind
				check if both parents have also been sequenced and aligned.
		"""
		alignmentLs = db_vervet.getAlignmentsFromVCFFile(inputFname =inputFname)
		
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		DG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		
		#2012.8.29
		if treatEveryOneIndependent:
			removeFamilyFromGraph = True
		else:
			removeFamilyFromGraph = False
		#find trios first
		trioLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=3, removeFamilyFromGraph=removeFamilyFromGraph)
		#find duos
		duoLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=2, removeFamilyFromGraph=removeFamilyFromGraph)
		#find singletons (familySize=1 => noOfIncomingEdges=0, noOfIncomingEdges=0 => will not be parents of others)
		singletonLs = db_vervet.findFamilyFromPedigreeGivenSize(DG, familySize=1, removeFamilyFromGraph=removeFamilyFromGraph, noOfOutgoingEdges=0)
		#[[node] for node in DG.nodes()]. this will be good only when removeFamilyFromGraph=True for all sentences above.
		
		familyLs = trioLs + duoLs + singletonLs
		sys.stderr.write("Outputting %s families in units of trios (Merlin format)..."%(len(familyLs)))
		#time to output
		sampleID2FamilyCount = {}	#2012.3.29
		writer = csv.writer(open(pedigreeOutputFname, 'w'), delimiter=' ')
		noOfFamilies = 0
		for family in familyLs:
			familySize = len(family)
			for memberID in family:	#2012.3.29 make sure to record the occurrence of each member before writing them out
				# as the occurrence is used in output 
				for alignment in individual_id2alignmentLs.get(memberID):
					self.adjustSampleID2FamilyCountFromAlignmentDBEntry(alignment, sampleID2FamilyCount)
				
			if familySize==1:
				individual_id = family[0]
				for alignment in individual_id2alignmentLs.get(individual_id):
					sampleID = self.getSampleIDWithReplicateCount(alignment, replicateIndividualTag=replicateIndividualTag, \
												sampleID2FamilyCount=sampleID2FamilyCount)
					family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[alignment.getReadGroup()])
					self.writeOneLineToPedigreeFile(writer, family_id, alignment, father_id=0, mother_id=0,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
					noOfFamilies += 1
			elif familySize==2:	#2012.9.26 output below is wrong if some members of duo or trio have >1 alignments 
				parent1ID, offspring_id = family[:2]
				#output the single parent first
				for parent1Alignment in individual_id2alignmentLs.get(parent1ID):
					parent1SampleID = self.getSampleIDWithReplicateCount(parent1Alignment, replicateIndividualTag=replicateIndividualTag, \
												sampleID2FamilyCount=sampleID2FamilyCount)
					family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[parent1Alignment.getReadGroup()])	#2012.3.29 add count to family ID
					self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
					#output the offspring
					for childAlignment in individual_id2alignmentLs.get(offspring_id):
						family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[parent1Alignment.getReadGroup()])	#2012.3.29 add count to family ID
						if treatEveryOneIndependent:	#2012.8.29
							father_id = 0
							mother_id = 0
						else:
							#output a fake 2nd parent, with opposite sex
							sndParentID = 'dummy'
							sndParentSex = 1-(parent1Alignment.individual_sequence.individual.codeSexInNumber()-1)+1
							self.writeOneLineToPedigreeFile(writer, family_id, sndParentID, father_id=0, mother_id=0, sex=sndParentSex,\
													replicateIndividualTag=replicateIndividualTag, \
													sampleID2FamilyCount=sampleID2FamilyCount)
							if sndParentSex==1:
								father_id = sndParentID
								mother_id = parent1SampleID
							else:
								father_id = parent1SampleID
								mother_id = sndParentID
						self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
													replicateIndividualTag=replicateIndividualTag, \
													sampleID2FamilyCount=sampleID2FamilyCount)
						noOfFamilies += 1
			elif familySize==3:
				parent1ID, parent2ID, offspring_id = family[:3]
				#output one parent
				for parent1Alignment in individual_id2alignmentLs.get(parent1ID):
					parent1SampleID = self.getSampleIDWithReplicateCount(parent1Alignment, replicateIndividualTag=replicateIndividualTag, \
												sampleID2FamilyCount=sampleID2FamilyCount)
					family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[parent1Alignment.getReadGroup()])
					self.writeOneLineToPedigreeFile(writer, family_id, parent1Alignment, father_id=0, mother_id=0,\
											replicateIndividualTag=replicateIndividualTag, \
											sampleID2FamilyCount=sampleID2FamilyCount)
					#output 2nd parent
					for parent2Alignment in individual_id2alignmentLs.get(parent2ID):
						parent2SampleID = self.getSampleIDWithReplicateCount(parent2Alignment, replicateIndividualTag=replicateIndividualTag, \
													sampleID2FamilyCount=sampleID2FamilyCount)
						family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[parent1Alignment.getReadGroup()])
						self.writeOneLineToPedigreeFile(writer, family_id, parent2Alignment, father_id=0, mother_id=0,\
												replicateIndividualTag=replicateIndividualTag, \
												sampleID2FamilyCount=sampleID2FamilyCount)
						#output offspring
						for childAlignment in individual_id2alignmentLs.get(offspring_id):
							if treatEveryOneIndependent:	#2012.8.29
								#output the offspring
								father_id = 0
								mother_id = 0
							else:
								parent2Sex = parent2Alignment.individual_sequence.individual.codeSexInNumber()
								if parent2Sex==1:
									father_id = parent2SampleID
									mother_id = parent1SampleID
								else:
									father_id = parent1SampleID
									mother_id = parent2SampleID
							family_id = "F%s_%s"%(noOfFamilies, sampleID2FamilyCount[parent1Alignment.getReadGroup()])
							self.writeOneLineToPedigreeFile(writer, family_id, childAlignment, father_id=father_id, mother_id=mother_id,\
														replicateIndividualTag=replicateIndividualTag, \
														sampleID2FamilyCount=sampleID2FamilyCount)
							noOfFamilies += 1
		del writer
		sys.stderr.write("%s families.\n"%(noOfFamilies))
		
		self.outputSampleID2FamilyCount(sampleID2FamilyCount, outputFname=sampleID2FamilyCountFname)
		return sampleID2FamilyCount
		
	def outputSampleID2FamilyCount(self, sampleID2FamilyCount=None, outputFname=None):
		"""
		2013.1.6 output only when outputFname is something.
		2012.3.29
		"""
		if outputFname:
			sys.stderr.write("Writing sampleID2FamilyCount (%s entries) to %s ..."%(len(sampleID2FamilyCount), outputFname))
			header= ['individualID', 'familyCount']
			writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
			writer.writerow(header)
			count = 0
			for individualID, familyCount in sampleID2FamilyCount.iteritems():
				data_row = [individualID, familyCount]
				writer.writerow(data_row)
				count +=1 
			del writer
			sys.stderr.write("%s entries outputted.\n"%(len(sampleID2FamilyCount)))
		
	
	def rankNodesByOutDegree(self, graph=None, node2HierarchyLevel=None):
		"""
		2013.1.2
			use a heap to store, the data needs to include the children of each ranked node + itself.
			rank by out-degree and hierarchy-level (distance to leaf-nodes) 
		"""
		queue = []
		for n in graph.nodes():
			level = node2HierarchyLevel.get(n)
			out_degree = graph.out_degree(n)
			nodeData = PassingData(node=n, children=graph.successors(n))
			heapq.heappush(queue, [(-out_degree, -level), nodeData])
		return queue
		
	
	def breakLoopByReplicatingOneNode(self, loopGraph=None, fullGraph=None, node=None, nodeChildrenInLoopGraph=None, replicateIndividualTag=None,\
								nodeID2ReplicateCount=None):
		"""
		2013.1.2
			Leave other children of the node that are not in the loop-graph untouched.
			
			nodeChildrenInLoopGraph is the children of node in the loop-graph (a subset of its children in the fullGraph).
		"""
		if node not in nodeID2ReplicateCount:
			nodeID2ReplicateCount[node] = 1	#start from 1
		for child in nodeChildrenInLoopGraph[1:]:	#keep relationship between node and one child intact.
			#the child with fewst number of ancestors should be that child. right now, random.
			fullGraph.remove_edge(node, child)
			loopGraph.remove_edge(node, child)	#remove edge for the loopGraph
			
			nodeID2ReplicateCount[node] += 1
			replicateNodeID = (node, nodeID2ReplicateCount[node])	#the 2nd value is replicateIndex (start from 1)
			fullGraph.add_node(replicateNodeID)
			fullGraph.add_edge(replicateNodeID, child)
		
	
	def identifyLoopBreakers(self, loopGraph=None, fullGraph=None, replicateIndividualTag=None,\
							nodeID2ReplicateCount=None, node2HierarchyLevel=None):
		"""
		2013.1.2
		"""
		loopGraph.recursiveRemoveUniDegreeNodes()
		if loopGraph.number_of_edges()>0:
			nodeHeap = self.rankNodesByOutDegree(graph=loopGraph, node2HierarchyLevel=node2HierarchyLevel)
			topNodeData = heapq.heappop(nodeHeap)[1]
			nodeChildrenInLoopGraph = topNodeData.children
			topNode = topNodeData.node
			self.breakLoopByReplicatingOneNode(loopGraph=loopGraph, fullGraph=fullGraph, node=topNode, nodeChildrenInLoopGraph=nodeChildrenInLoopGraph,\
											replicateIndividualTag=replicateIndividualTag,\
											nodeID2ReplicateCount=nodeID2ReplicateCount)
			self.identifyLoopBreakers(loopGraph=loopGraph, fullGraph=fullGraph, replicateIndividualTag=replicateIndividualTag,\
							nodeID2ReplicateCount=nodeID2ReplicateCount, node2HierarchyLevel=node2HierarchyLevel)
	
	def convertReplicateNodeID2VCFSampleID(self, nodeID=None, replicateIndividualTag=None, nodeID2sampleID=None):
		"""
		2013.1.3
		"""
		if isinstance(nodeID, tuple):
			trueNodeID, replicateIndex = nodeID[:2]
		else:
			trueNodeID = nodeID
			replicateIndex = 1	#starting from 1
		if trueNodeID==0:	#missing individuals (or parents of founders)
			vcfSampleID = 0
		else:
			sampleID = nodeID2sampleID.get(trueNodeID, trueNodeID)
			vcfSampleID = '%s%s%s'%(sampleID, replicateIndividualTag, replicateIndex)
		return vcfSampleID
	
	def writeOutPolymuttPedigree(self, DG=None, db_vervet=None, outputFname=None, replicateIndividualTag=None,\
								nodeID2sampleID=None, treatEveryOneIndependent=None):
		"""
		2012.1.4
		"""
		sys.stderr.write("Writing pedigree to %s ..."%(outputFname))
		writer = csv.writer(open(outputFname, 'w'), delimiter=' ')
		counter = 0
		family_id= 0
		individual_id2individual = {}
		cc_nodes_list = nx.connected_component_subgraphs(DG.to_undirected())
		for cc_nodes in cc_nodes_list:	#each component is one family
			cc_subgraph = nx.subgraph(DG, cc_nodes)
			family_id += 1
			for nodeID in cc_subgraph.nodes():	#this nodeID could be (individual_id, replicate_index)
				individual = self.getIndividual(db_vervet=db_vervet, individual_id=nodeID, \
											individual_id2individual=individual_id2individual)
				if treatEveryOneIndependent:
					father_id = 0
					mother_id = 0
				elif nodeID in cc_subgraph:
					parents = cc_subgraph.predecessors(nodeID)
					if len(parents)==2:
						parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
													individual_id2individual=individual_id2individual)
						parent2 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[1], \
													individual_id2individual=individual_id2individual)
						parent1Sex = parent1.codeSexInNumber()
						if parent1Sex==2:
							#swap the father and mother row
							father_id = parents[1]
							mother_id = parents[0]
						else:
							father_id = parents[0]
							mother_id = parents[1]
					elif len(parents)==1:
						parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
											individual_id2individual=individual_id2individual)
						parent1Sex = parent1.codeSexInNumber()
						if parent1Sex==2:
							father_id = 0
							mother_id = parents[0]
						else:
							father_id = parents[0]
							mother_id = 0
					elif len(parents)==0:
						father_id = 0
						mother_id = 0
					else:
						sys.stderr.write("Error: number of parents (%s) for %s is %s.\n"%(repr(parents), nodeID, len(parents)))
						sys.exit(3)
				else:
					father_id = 0
					mother_id = 0
				thisNodeVCFSampleID = self.convertReplicateNodeID2VCFSampleID(nodeID=nodeID, replicateIndividualTag=replicateIndividualTag, \
															nodeID2sampleID=nodeID2sampleID)
				fatherVCFSampleID = self.convertReplicateNodeID2VCFSampleID(nodeID=father_id, replicateIndividualTag=replicateIndividualTag, \
															nodeID2sampleID=nodeID2sampleID)
				motherVCFSampleID = self.convertReplicateNodeID2VCFSampleID(nodeID=mother_id, replicateIndividualTag=replicateIndividualTag, \
															nodeID2sampleID=nodeID2sampleID)
				data_row = ['fam%s'%family_id, thisNodeVCFSampleID,\
						fatherVCFSampleID, motherVCFSampleID, individual.codeSexInNumber(), 1]
					#the sixth column is glf index if input is glf.
					#otherwise should be useless
				writer.writerow(data_row)
				counter += 1
		del writer
		sys.stderr.write("%s nodes outputed.\n"%(counter))
	
	def outputPedigreeForPolymutt(self, DG=None, db_vervet=None, inputFname=None, outputFname=None, replicateIndividualTag=None,\
								sampleID2FamilyCountFname=None, treatEveryOneIndependent=False):
		"""
		2013.1.2 for polymutt http://genome.sph.umich.edu/wiki/Polymutt
		"""
		
		alignmentLs = db_vervet.getAlignmentsFromVCFFile(inputFname =inputFname)
		sequencedIndividualIDSet = set()
		nodeID2sampleID = {}	#sampleID is the name for each sample in VCF file.
		for alignment in alignmentLs:
			nodeID = alignment.individual_sequence.individual.id
			sequencedIndividualIDSet.add(nodeID)
			nodeID2sampleID[nodeID] = alignment.sampleID
		
		DG.recursiveTrimOutOfSetLeafNodes(nodeIdSet=sequencedIndividualIDSet)
		node2HierarchyLevel = DG.calculateNodeHierarchyLevel()
		graphCopy = copy.deepcopy(DG)
		graphCopy2 = copy.deepcopy(DG)
		nodeID2ReplicateCount = {}
		self.identifyLoopBreakers(loopGraph=graphCopy,\
								fullGraph=DG, replicateIndividualTag=replicateIndividualTag,\
								nodeID2ReplicateCount=nodeID2ReplicateCount, node2HierarchyLevel=node2HierarchyLevel)
		
		if sampleID2FamilyCountFname:
			sampleID2FamilyCount = {}
			for nodeID in graphCopy2.nodes():
				count = nodeID2ReplicateCount.get(nodeID, 1)
				if nodeID in nodeID2sampleID:
					sampleID = nodeID2sampleID.get(nodeID)
				else:
					sampleID = nodeID
				sampleID2FamilyCount[sampleID] = count
			self.outputSampleID2FamilyCount(sampleID2FamilyCount, outputFname=sampleID2FamilyCountFname)
		
		#output the modified DG
		#if one node has sampleID from VCF file, use that instead of node ID.
		self.writeOutPolymuttPedigree(DG=DG, db_vervet=db_vervet, outputFname=outputFname, replicateIndividualTag=replicateIndividualTag,\
									nodeID2sampleID=nodeID2sampleID, treatEveryOneIndependent=treatEveryOneIndependent)
	
	def outputPolymuttDatFile(self, outputFname=None):
		"""
		2013.1.4 the datfile for polymutt
		"""
		if outputFname:
			sys.stderr.write("Writing to polymutt datfil %s ..."%(outputFname))
			outf = open(outputFname, 'w')
			outf.write("T\tGLF_Index\n")
			outf.close()
			sys.stderr.write("\n")
	
	
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		db_vervet = self.db_vervet
		session.begin()
		
		DG = db_vervet.constructPedgree()
		
		if self.outputFileFormat==1:
			self.outputPedigreeForPlink(DG=DG, db_vervet=db_vervet, inputFname=self.inputFname, outputFname=self.outputFname, \
									treatEveryOneIndependent=self.treatEveryOneIndependent)
		elif self.outputFileFormat==3:
			self.outputPedigreeForPolymutt(DG=DG, db_vervet=db_vervet, inputFname=self.inputFname, \
								outputFname=self.outputFname, replicateIndividualTag=self.replicateIndividualTag,\
								sampleID2FamilyCountFname=self.sampleID2FamilyCountFname, \
								treatEveryOneIndependent=self.treatEveryOneIndependent)
		elif self.outputFileFormat==2:
			self.outputPedigreeForTrioCaller(db_vervet=db_vervet, inputFname=self.inputFname, \
												pedigreeOutputFname=self.outputFname, replicateIndividualTag=self.replicateIndividualTag,\
												treatEveryOneIndependent=self.treatEveryOneIndependent,\
												sampleID2FamilyCountFname=self.sampleID2FamilyCountFname)
		else:
			sys.stderr.write("Error: un-supported output file format %s.\n"%(self.outputFileFormat))
			sys.exit(3)
		
		if self.polymuttDatFname:	#2013.1.4 the datfile for polymutt
			self.outputPolymuttDatFile(self.polymuttDatFname)
	
if __name__ == '__main__':
	main_class = OutputVRCPedigreeInTFAMGivenOrderFromFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()