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
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src.plot.PlotPedigreeKinshipVsGeneticIBD import PlotPedigreeKinshipVsGeneticIBD
import heapq
from scipy import stats
import rpy

class DetectWrongLabelByCompKinshipVsIBD(PlotPedigreeKinshipVsGeneticIBD):
	__doc__ = __doc__
#						
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.update({
						('plinkIBDCheckOutputFname', 1, ): ["", 'l', 1, 'file that contains IBD check result'], \
						('plinkSexCheckOutputFname', 0, ): ["", 'n', 1, 'file that contains plink sex check result'], \
						('kinshipMonkeyIDSetFname', 0, ): ["", 's', 1, 'restrict j in the kinship(i,j)-ibd(i,j) comparison. \
			j must be in this file. temporary 2012.8.24 four columns:monkeyID	noOfMismatches	noOfNonMissing	mismatchFraction.\
			and only when noOfMismatches=0 and noOfNonMissing>=30, this monkey is included.'], \
						('iterativeAlgorithm', 0, ): ["", 'A', 0, 'toggle to iteratively update chiSqStat. keep removing the top one'], \
						('minAbsDeltaForOutlier', 0, float): [0, 'm', 1, 'if not 0, this will be used as minAbsDelta in outlier frequency counting'], \
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		PlotPedigreeKinshipVsGeneticIBD.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def createDeltaMatrix(self, kinshipData=None, ibdData=None, takeAbs=False, onlySharedMonkeys=False, defaultValue=None):
		"""
		2012.8.22
		"""
		sys.stderr.write("Creating kinship-delta matrix ...")
		if onlySharedMonkeys:
			sharedMonkeyIDSet = set(kinshipData.row_id_ls)&set(ibdData.row_id_ls)
			monkey_id_ls = list(sharedMonkeyIDSet)
			monkey_id_ls.sort()
			row_id2row_index = {}
			for i in xrange(len(monkey_id_ls)):
				monkey_id = monkey_id_ls[i]
				row_id2row_index[monkey_id] = i
		else:
			monkey_id_ls = kinshipData.row_id_ls
			row_id2row_index = kinshipData.row_id2row_index
		no_of_monkeys = len(monkey_id_ls)
		kinshipIBDDeltaMatrix = numpy.zeros([no_of_monkeys, no_of_monkeys], dtype=numpy.float32)
		if defaultValue is None:
			defaultValue = numpy.nan
		kinshipIBDDeltaMatrix[:] = defaultValue
		for i in xrange(no_of_monkeys):
			monkey1_id = monkey_id_ls[i]
			monkey1_index = row_id2row_index.get(monkey1_id)
			for j in xrange(i+1, no_of_monkeys):
				monkey2_id = monkey_id_ls[j]
				monkey2_index = row_id2row_index.get(monkey2_id)
				kinship = kinshipData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				ibd = ibdData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				if kinship is not None and ibd is not None:
					if not hasattr(kinship, 'mask') and not hasattr(ibd, 'mask'):
						#if it has attribute mask, it's not a float value, and not defined in the array
						delta = kinship-ibd
						if takeAbs:
							delta = abs(delta)
						kinshipIBDDeltaMatrix[monkey1_index][monkey2_index] = delta
						kinshipIBDDeltaMatrix[monkey2_index][monkey1_index] = delta
		
		maskedMatrix = numpy.ma.array(kinshipIBDDeltaMatrix, mask=numpy.isnan(kinshipIBDDeltaMatrix))
		kinshipIBDDeltaData = SNPData(row_id_ls=monkey_id_ls, col_id_ls=monkey_id_ls, \
									data_matrix=maskedMatrix)
		sys.stderr.write(" %s rows&columns .\n"%(len(kinshipIBDDeltaData.row_id_ls)))
		return kinshipIBDDeltaData
	
	def cutOffKinshipIBDDeltaAndOutput(self, db_vervet=None, kinshipData=None, ibdData=None, \
					outputFnamePrefix=None, minAbsDelta=0.1, kinshipMonkeyIDSet=None):
		"""
		2012.8.23
			choose the pairs that have abs(kinship-ibd)>=minAbsDelta and output them
			count how many times one monkey appears in all these "outlier" pairs and output that as well.
			
			if kinshipMonkeyIDSet is given, the columns would only be among them.
		"""
		sys.stderr.write("Choosing pairs whose abs(kinship-ibd) >=%s ..."%(minAbsDelta))
		
		tableOutputWriter = csv.writer(open("%s_minAbsDelta%s.tsv"%(outputFnamePrefix, minAbsDelta), 'w'), delimiter='\t')
		header = ['monkey1ID', 'monkey2ID', 'pedigreeKinship', 'IBD_PI_HAT', 'kinshipIBDDelta','ageDifference']
		tableOutputWriter.writerow(header)
		
		monkey_id_pair_ls = []
		
		monkey_id2index = {}	#2012.8.21
		monkey_id_pair2kinship_ibd_delta = {}
		monkeyID2AbsDeltaVector = {}
		shared_monkey_id_set = set(ibdData.row_id_ls)&set(kinshipData.row_id_ls)
		shared_monkey_id_ls = list(shared_monkey_id_set)
		no_of_monkeys = len(shared_monkey_id_ls)
		for i in xrange(no_of_monkeys):
			monkey1_id = shared_monkey_id_ls[i]
			if kinshipMonkeyIDSet is None:
				kinshipMonkeyIDSet = shared_monkey_id_ls[i+1:]
			for monkey2_id in kinshipMonkeyIDSet:
				monkey_id_pair = [monkey1_id, monkey2_id]
				monkey_id_pair.sort()
				monkey_id_pair  = tuple(monkey_id_pair)
				if monkey_id_pair in monkey_id_pair2kinship_ibd_delta:	#already in , skip it. only happens when kinshipMonkeyIDSet is given
					continue
				ibd = ibdData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				kinship = kinshipData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				#if (not hasattr(ibd, 'mask')) and (not numpy.isnan(ibd)) and (not hasattr(kinship, 'mask'))\
				#		and (not numpy.isnan(kinship)):
				if ibd is not None and kinship is not None and (not numpy.isnan(ibd)) and (not numpy.isnan(kinship)) \
						and (not hasattr(ibd, 'mask')) and (not hasattr(kinship, 'mask')):
					delta = kinship-ibd
					absDelta = abs(delta)
					if absDelta>=minAbsDelta:
						monkey1 = self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=monkey1_id)
						if monkey1_id not in monkey_id2index:
							monkey_id2index[monkey1_id] = len(monkey_id2index)
						monkey2 =  self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=monkey2_id)
						if monkey2_id not in monkey_id2index:
							monkey_id2index[monkey2_id] = len(monkey_id2index)
						
						if monkey1 and monkey2 and monkey1.getCurrentAge() and monkey2.getCurrentAge():
							ageDelta = abs(monkey1.getCurrentAge() - monkey2.getCurrentAge())
						else:
							ageDelta = ''
						data_row = [monkey_id_pair[0], monkey_id_pair[1], kinship, ibd, delta, ageDelta]
						tableOutputWriter.writerow(data_row)
						monkey_id_pair2kinship_ibd_delta[monkey_id_pair] = kinship-ibd
						for j in xrange(len(monkey_id_pair)):
							monkeyID = monkey_id_pair[j]
							if monkeyID not in monkeyID2AbsDeltaVector:
								monkeyID2AbsDeltaVector[monkeyID] = []
							for k in xrange(len(monkey_id_pair)):	#2012.9.1 there is only one other monkey. but add them all.
								monkey2ID = monkey_id_pair[k]
								if monkeyID!=monkey2ID:
									monkeyID2AbsDeltaVector[monkeyID].append([monkey2ID, absDelta])
		
		del tableOutputWriter
		sys.stderr.write("  %s outlier pairs, covering %s monkeys.\n"%(len(monkey_id_pair2kinship_ibd_delta),\
															len(monkeyID2AbsDeltaVector)))
		
		monkeyIDList = monkeyID2AbsDeltaVector.keys()
		monkeyIDList.sort()
		queue = []	#a heap queue by the number of outlier pairs one monkey appears in 
		for monkeyID in monkeyIDList:
			data_ls = monkeyID2AbsDeltaVector.get(monkeyID)
			noOfOutliers = len(data_ls)
			#create a one-row X multi-column SNPData-like matrix structure
			row_id_ls = [monkeyID]
			col_id_ls = []	#the IDs of monkeys in those outlier pairs
			data_row = []	#store the absDelta value
			for data in data_ls:
				col_id_ls.append(data[0])
				data_row.append(data[1])
			
			deltaData = SNP.SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=numpy.array([data_row]))
			heapq.heappush(queue, [-noOfOutliers, monkeyID, deltaData])	#watch the -
		#return queue
		self.outputOutlierFrequencyQueueData(outputFnamePrefix=outputFnamePrefix, minAbsDelta=minAbsDelta, \
											outlierMonkeyQueue=queue)
	
	def updateOutlierMonkeyQueue(self, queue=None, monkeyToRemove=None):
		"""
		2012.9.1
		
		"""
		no_of_monkeys = len(queue)
		sys.stderr.write("Updating a outlierMonkeyQueue for %s monkeys ..."%\
						(no_of_monkeys))
		newQueue = []
		
		while len(queue)>0:
			minusNoOfOutliers, monkeyID, deltaData = heapq.heappop(queue)[:3]
			#if monkeyToRemove in deltaData.col_id2col_index:
			monkeyToRemoveIndex = deltaData.col_id2col_index.get(monkeyToRemove)
			if monkeyToRemoveIndex is not None:
				keptColIDLs = deltaData.col_id_ls[:monkeyToRemoveIndex] + deltaData.col_id_ls[monkeyToRemoveIndex+1:]
				deltaData = deltaData.keepColsByColID(deltaData, col_id_ls=keptColIDLs, dataType=numpy.float)
				minusNoOfOutliers += 1
			heapq.heappush(newQueue, [minusNoOfOutliers, monkeyID, deltaData])
		sys.stderr.write("  %s monkeys in the queue.\n"%(len(newQueue)))
		return PassingData(queue=newQueue)
		
	
	def outputOutlierFrequencyQueueData(self, outputFnamePrefix=None, minAbsDelta=None, outlierMonkeyQueue=None):
		"""
		2012.9.1
		"""
		sys.stderr.write("Outputting outlierMonkeyQueue ... \n")
		monkeyFrequencyInOutlierWriter = csv.writer(open('%s_monkeyFrequencyInMinAbsDelta%sPairs.tsv'%(outputFnamePrefix, minAbsDelta), 'w'),\
												delimiter='\t')
		header = ['monkeyID', 'outlierFrequency', 'medianAbsDeltaAmongOutliers']
		monkeyFrequencyInOutlierWriter.writerow(header)
		i = 0
		while len(outlierMonkeyQueue)>0:
			minusNoOfOutliers, monkeyID, deltaData = heapq.heappop(outlierMonkeyQueue)[:3]
			noOfOutliers = -minusNoOfOutliers
			
			medianAbsDeltaAmongOutliers = numpy.median(deltaData.data_matrix[0,:])
			data_row = [monkeyID, noOfOutliers, medianAbsDeltaAmongOutliers]
			monkeyFrequencyInOutlierWriter.writerow(data_row)
			pdata = self.updateOutlierMonkeyQueue(outlierMonkeyQueue, monkeyToRemove=monkeyID)
			outlierMonkeyQueue = pdata.queue
			i+= 1
		del monkeyFrequencyInOutlierWriter
		sys.stderr.write(" %s monkeys outputted.\n"%(i))
	
	
	def createKinshipIBDDeltaQueue(self, kinshipIBDDeltaData=None):
		"""
		2012.8.22
			a queue sorted by the -median(abs(list of kinship-ibd) for each monkey
			put - in front of it because heapq is a min queue.
		"""
		no_of_monkeys = len(kinshipIBDDeltaData.row_id2row_index)
		sys.stderr.write("Creating a  kinshipIBDDeltaQueue from a %sX%s kinshipIBDDeltaData matrix ..."%\
						(no_of_monkeys, no_of_monkeys))
		kinshipIBDDeltaQueue = []
		monkey_id2medianAbsDelta = {}
		maskedMatrix = kinshipIBDDeltaData.data_matrix
		for monkeyID , index in kinshipIBDDeltaData.row_id2row_index.iteritems():
			
			noOfNonMissing = len(maskedMatrix[index,:]) - numpy.ma.sum(maskedMatrix[index,:].mask)
			if noOfNonMissing>0:
				medianAbsDelta = numpy.ma.median(numpy.ma.abs(maskedMatrix[index,:]))
				heapq.heappush(kinshipIBDDeltaQueue, [-medianAbsDelta, monkeyID, noOfNonMissing])	#watch the -
				monkey_id2medianAbsDelta[monkeyID] = PassingData(medianAbsDelta=medianAbsDelta, noOfNonMissing=noOfNonMissing)
		sys.stderr.write("  %s monkeys in the kinshipIBDDeltaQueue.\n"%(len(kinshipIBDDeltaQueue)))
		return PassingData(kinshipIBDDeltaQueue=kinshipIBDDeltaQueue, monkey_id2medianAbsDelta=monkey_id2medianAbsDelta)
	
	def calculateMedianAbsDelta(self, kinshipData=None, kinshipDataMonkeyID=None, ibdData=None, ibdDataMonkeyID=None):
		"""
		2012.8.22
		"""
		medianAbsDelta = None
		noOfNonMissing = None
		kinshipVector = kinshipData.getRowVectorGivenRowID(row_id= kinshipDataMonkeyID)
		ibdVector = ibdData.getRowVectorGivenRowID(row_id=ibdDataMonkeyID)
		if kinshipVector is not None and ibdVector is not None:
			no_of_monkeys = len(kinshipData.col_id2col_index)
			kinshipIBDDeltaVector = numpy.zeros([no_of_monkeys], dtype=numpy.float32)
			kinshipIBDDeltaVector[:] = numpy.nan
			for monkeyID, kinship_index in kinshipData.col_id2col_index.iteritems():
				kinship = kinshipVector[kinship_index]
				if kinship is not None and not hasattr(kinship, 'mask'):
					ibdIndex = ibdData.col_id2col_index.get(monkeyID)
					if ibdIndex is not None:
						try:
							ibd = ibdVector[ibdIndex]
							if ibd is not None and not hasattr(ibd, 'mask') and not numpy.isnan(ibd):
								kinshipIBDDeltaVector[kinship_index] = kinship-ibd
						except:
							sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
							import traceback
							traceback.print_exc()
							import pdb
							pdb.set_trace()
			
			maskedVector = numpy.ma.array(kinshipIBDDeltaVector, mask=numpy.isnan(kinshipIBDDeltaVector))
			noOfNonMissing = len(maskedVector) - numpy.ma.sum(maskedVector.mask)
			if noOfNonMissing>0:
				medianAbsDelta = numpy.ma.median(numpy.ma.abs(maskedVector))
		return PassingData(medianAbsDelta=medianAbsDelta, noOfNonMissing=noOfNonMissing)
	
	def createKinshipIBDDeltaChiSqStatQueue(self, kinshipData=None, ibdData=None, \
								mean=None, std=None, given_row_id_ls=None):
		"""
		2012.8.22
			a queue sorted by the -\sum 2*log(Prob(abs(list of kinship-ibd))) for each monkey
			put - in front of it because heapq is a min queue.
		"""
		no_of_monkeys = len(ibdData.row_id2row_index)
		sys.stderr.write("Creating a  kinshipIBDDeltaChiSqStatQueue for %s monkeys ..."%\
						(no_of_monkeys))
		queue = []
		monkey_id2queueData = {}
		maskedMatrix = ibdData.data_matrix
		if given_row_id_ls is None:
			given_row_id_ls = set(kinshipData.row_id_ls) & set(ibdData.row_id_ls)
		for monkeyID in given_row_id_ls:
			chiSqStatData = self.calculateChiSqStatOfDeltaVector(kinshipData=kinshipData, kinshipDataMonkeyID=monkeyID, \
						ibdData=ibdData, ibdDataMonkeyID=monkeyID,\
						mean=mean, std=std)
			noOfNonMissing = chiSqStatData.noOfNonMissing
			chiSqStat = chiSqStatData.chiSqStat
			chiSqPvalue = chiSqStatData.chiSqPvalue
			if noOfNonMissing>0:
				heapq.heappush(queue, [-chiSqStat, monkeyID, noOfNonMissing, chiSqPvalue])	#watch the -
				monkey_id2queueData[monkeyID] = PassingData(chiSqStat=chiSqStat, noOfNonMissing=noOfNonMissing,
													chiSqPvalue=chiSqPvalue)
		sys.stderr.write("  %s monkeys in the queue.\n"%(len(queue)))
		return PassingData(queue=queue, monkey_id2queueData=monkey_id2queueData)
	
	
	def calculateChiSqStatOfDeltaVector(self, kinshipData=None, kinshipDataMonkeyID=None, ibdData=None, ibdDataMonkeyID=None,\
								mean=None, std=None):
		"""
		2012.8.22
			fisher's method in combining pvalue of z-score (normalized abs(kinship-ibd)).
		"""
		medianAbsDelta = None
		noOfNonMissing = None
		kinshipVector = kinshipData.getRowVectorGivenRowID(row_id= kinshipDataMonkeyID)
		ibdVector = ibdData.getRowVectorGivenRowID(row_id=ibdDataMonkeyID)
		chiSqStat = 0
		noOfNonMissing = 0
		chiSqMinusLogPvalue = None
		if kinshipVector is not None and ibdVector is not None:
			noOfNonMissingKinship = len(kinshipVector)- numpy.ma.sum(kinshipVector.mask)
			noOfNonMissingIBD = len(ibdVector)- numpy.ma.sum(ibdVector.mask)
			
			if noOfNonMissingKinship>0 and noOfNonMissingIBD>0:
				sharedMonkeyIDSet = set(kinshipData.col_id_ls) & set(ibdData.col_id_ls)
				for monkeyID in sharedMonkeyIDSet:
					kinship_index = kinshipData.col_id2col_index.get(monkeyID)
					#for monkeyID, kinship_index in kinshipData.col_id2col_index.iteritems():
					kinship = kinshipVector[kinship_index]
					if kinship is not None and not hasattr(kinship, 'mask'):
						ibdIndex = ibdData.col_id2col_index.get(monkeyID)
						if ibdIndex is not None:
							try:
								ibd = ibdVector[ibdIndex]
								if ibd is not None and not hasattr(ibd, 'mask') and not numpy.isnan(ibd):
									delta = kinship-ibd	#2012.8.23 no abs
									Z = abs((delta-mean)/std)	#2012.8.23 take absolute value, since it's always P(X>a), negative z-score gets wrong portion.
									logPvalue = rpy.r.pnorm(Z, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)	#the latter log is natural log.
									#should use two-tail, rather than one-tail
									logPvalue += math.log(2)
									
									#pvalue = utils.getZScorePvalue(Z)	#numerical underflow, pvalue=0
									chiSqStat += -2*logPvalue	#should use natural log
									noOfNonMissing += 1
							except:
								sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
								import traceback
								traceback.print_exc()
								import pdb
								pdb.set_trace()
				if noOfNonMissing>0:
					chiSqMinusLogPvalue = -rpy.r.pchisq(chiSqStat, 2*noOfNonMissing, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)
					#chiSqPvalue = stats.chisqprob(chiSqStat, df=2*noOfNonMissing)
					#equivalent to = 1- stats.chi2.cdf(chiSqStat, df) = Prob (x<chiSqStat)
		return PassingData(chiSqStat=chiSqStat, noOfNonMissing=noOfNonMissing, chiSqPvalue=chiSqMinusLogPvalue)
	
	def updateKinshipIBDDeltaChiSqStatQueue(self, queue=None, kinshipData=None, ibdData=None, \
								mean=None, std=None, dropMonkeyID=None):
		"""
		2012.8.23 updates the queue after one monkey is dropped.
		
			a queue sorted by the -\sum 2*log(Prob(abs(list of kinship-ibd))) for each monkey
			put - in front of it because heapq is a min queue.
		"""
		no_of_monkeys = len(queue)
		sys.stderr.write("Updating a kinshipIBDDeltaChiSqStatQueue for %s monkeys ..."%\
						(no_of_monkeys))
		newQueue = []
		
		monkey_id2queueData = {}
		maskedMatrix = ibdData.data_matrix
		while len(queue)>0:
			minusChiSqStat, monkeyID, noOfNonMissing, chiSqPvalue = heapq.heappop(queue)[:4]
			chiSqStat = -minusChiSqStat
			oldChiSqStatData = PassingData(chiSqStat=chiSqStat, noOfNonMissing=noOfNonMissing, chiSqPvalue=chiSqPvalue)
			
			chiSqStatData = self.updateChiSqStatOfDeltaVector(chiSqStatData=oldChiSqStatData, kinshipData=kinshipData, \
															kinshipDataMonkeyID=monkeyID, \
						ibdData=ibdData, ibdDataMonkeyID=monkeyID,\
						mean=mean, std=std, dropMonkeyID=dropMonkeyID)
			noOfNonMissing = chiSqStatData.noOfNonMissing
			chiSqStat = chiSqStatData.chiSqStat
			chiSqPvalue = chiSqStatData.chiSqPvalue
			if noOfNonMissing>0:
				heapq.heappush(newQueue, [-chiSqStat, monkeyID, noOfNonMissing, chiSqPvalue])	#watch the -
				monkey_id2queueData[monkeyID] = PassingData(chiSqStat=chiSqStat, noOfNonMissing=noOfNonMissing,
													chiSqPvalue=chiSqPvalue)
		sys.stderr.write("  %s monkeys in the queue.\n"%(len(newQueue)))
		return PassingData(queue=newQueue, monkey_id2queueData=monkey_id2queueData)
	
	def updateChiSqStatOfDeltaVector(self, chiSqStatData=None, kinshipData=None, kinshipDataMonkeyID=None, ibdData=None, ibdDataMonkeyID=None,\
								mean=None, std=None, dropMonkeyID=None):
		"""
		2012.8.22
			fisher's method in combining pvalue of z-score (normalized abs(kinship-ibd)).
		"""
		kinshipVector = kinshipData.getRowVectorGivenRowID(row_id= kinshipDataMonkeyID)
		ibdVector = ibdData.getRowVectorGivenRowID(row_id=ibdDataMonkeyID)
		oldNoOfNonMissing = chiSqStatData.noOfNonMissing
		if kinshipVector is not None and ibdVector is not None and chiSqStatData.noOfNonMissing>0:
			kinship_index = kinshipData.col_id2col_index.get(dropMonkeyID)
			#for monkeyID, kinship_index in kinshipData.col_id2col_index.iteritems():
			kinship = kinshipVector[kinship_index]
			if kinship is not None and not hasattr(kinship, 'mask'):
				ibdIndex = ibdData.col_id2col_index.get(dropMonkeyID)
				if ibdIndex is not None:
					try:
						ibd = ibdVector[ibdIndex]
						if not hasattr(ibd, 'mask') and not numpy.isnan(ibd):
							delta = kinship-ibd	#2012.8.23 no abs
							Z = abs((delta-mean)/std)	#2012.8.23 take absolute value, since it's always P(X>a), negative z-score gets wrong portion.
							logPvalue = rpy.r.pnorm(Z, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)	#the latter log is natural log.
							#should use two-tail, rather than one-tail
							logPvalue += math.log(2)
							#pvalue = utils.getZScorePvalue(Z)	#numerical underflow, pvalue=0
							chiSqStatData.chiSqStat -= -2*logPvalue	#subtract this guy out
							chiSqStatData.noOfNonMissing -= 1	#decrease the noOfNonMissing
					except:
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
						import traceback
						traceback.print_exc()
						import pdb
						pdb.set_trace()
			if oldNoOfNonMissing!=chiSqStatData.noOfNonMissing and chiSqStatData.noOfNonMissing>0:
				chiSqStatData.chiSqPvalue = -rpy.r.pchisq(chiSqStatData.chiSqStat, \
											2*chiSqStatData.noOfNonMissing, lower_tail = rpy.r.FALSE, log=rpy.r.TRUE)
				#chiSqPvalue = stats.chisqprob(chiSqStat, df=2*noOfNonMissing)
				#equivalent to = 1- stats.chi2.cdf(chiSqStat, df) = Prob (x<chiSqStat)
		return chiSqStatData
	
	
	def drawKinshipIBDDeltaVectorHistogram(self, kinshipIBDDeltaData=None, row_id=None, outputFnamePrefix=None):
		"""
		2012.8.22
		"""
		vector = kinshipIBDDeltaData.getRowVectorGivenRowID(row_id=row_id)
		if vector is not None:
			data_ls = []
			for i in xrange(len(vector)):
				if vector.mask[i]==False:
					data_ls.append(vector[i])
			if len(data_ls)>10:
				outputFname = '%s_monkey_%s_kinship_ibd_hist.png'%(outputFnamePrefix, row_id)
				yh_matplotlib.drawHist(data_ls, title='', \
								xlabel_1D="%s kinship-ibd"%(row_id), xticks=None, \
								outputFname=outputFname, min_no_of_data_points=10, \
								needLog=True, \
								dpi=200, min_no_of_bins=25)
	
	def estimateAbsDeltaMeanStd(self, kinshipIBDDeltaData=None, excludeTopFraction=0.2):
		"""
		2012.8.22
		"""
		sys.stderr.write("Estimating mean&std of abs(kinship-delta) using the middle %.1f of data ..."%((1-excludeTopFraction)*100))
		
		noOfRows = len(kinshipIBDDeltaData.row_id2row_index)
		data_ls = []
		for i in xrange(noOfRows):
			for j in xrange(i+1, noOfRows):
				if kinshipIBDDeltaData.data_matrix.mask[i][j]==False:
					data_ls.append(kinshipIBDDeltaData.data_matrix[i][j])
		
		# 2012.8.22 draw some histogram to check what data looks like
#		if len(data_ls)>10:
#			outputFname = '%s_kinship_ibd_hist.png'%(self.outputFnamePrefix)
#			yh_matplotlib.drawHist(data_ls, title='', \
#							xlabel_1D="kinship-ibd", xticks=None, \
#							outputFname=outputFname, min_no_of_data_points=10, \
#							needLog=True, \
#							dpi=200, min_no_of_bins=25)
		#data_ls = map(abs, data_ls)	#2012.8.23 no abs
		data_ls.sort()
		startIndex = min(0, int(len(data_ls)*(excludeTopFraction/2))-1)
		stopIndex = int(len(data_ls)*(1-excludeTopFraction/2))
		data_ls = data_ls[startIndex:stopIndex]
		
		data_mean = numpy.mean(data_ls)
		data_std = numpy.std(data_ls)
		
		sys.stderr.write(" mean=%.3f, std=%.3f.\n"%(data_mean, data_std))
		return PassingData(mean=data_mean, std=data_std)
	
	def PCAOnAbsKinshipIBDDeltaMatrix(self, kinshipData=None, ibdData=None, outputFnamePrefix=None):
		"""
		2012.8.24
			temporarily run PCA on abs(kinship-ibd) data matrix and then check clustering.
		"""
		#2012.8.24 temporarily doing something different
		kinshipIBDDeltaData = self.createDeltaMatrix(kinshipData=kinshipData, ibdData=ibdData, takeAbs=True, \
											onlySharedMonkeys=True, defaultValue=0)	#default is 0, not numpy.nan
		
		outputFname = '%s_PCAOnAbsKinshipIBDDelta.tsv'%outputFnamePrefix
		pcaData = kinshipIBDDeltaData.PCAOnDataMatrix(outputFname=outputFname)
		
		PC1_row_id_ls = []
		for i in xrange(len(kinshipIBDDeltaData.row_id_ls)):
			row_id = kinshipIBDDeltaData.row_id_ls[i]
			PC1 = pcaData.T[i,0]
			PC1_row_id_ls.append([PC1, row_id])
		PC1_row_id_ls.sort()
		row_id_sortedByPCLs = [row[1] for row in PC1_row_id_ls]
		#reorder the columns & rows of kinshipIBDDeltaData according to PC1 and then output
		kinshipIBDDeltaData = kinshipIBDDeltaData.reOrderRowsGivenNewRowIDList(row_id_ls=row_id_sortedByPCLs)
		kinshipIBDDeltaData = kinshipIBDDeltaData.reOrderColsGivenNewCOlIDList(col_id_ls = row_id_sortedByPCLs)
		
		outputFname = '%s_orderByPC1.tsv'%(outputFnamePrefix)
		kinshipIBDDeltaData.tofile(outputFname)
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname, schema=self.schema, port=self.port)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		
		kinshipData = self.getMonkeyKinshipData(inputFname=self.inputFname)
		#set kinshipData diagonal to 1
		ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.plinkIBDCheckOutputFname, rowIDHeader="IID1", colIDHeader="IID2", \
										rowIDIndex=None, colIDIndex=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		
		if self.minAbsDeltaForOutlier>0:
			#2012.8.23 cut data off for Sue
			if self.kinshipMonkeyIDSetFname:
				monkeyID2dataTuple = SNP.getKey2ValueFromMatrixLikeFile(inputFname=self.kinshipMonkeyIDSetFname, keyHeaderLs=['monkeyID'], \
									valueHeaderLs=['noOfMismatches', 'noOfNonMissing'], keyIndexLs=None, valueIndexLs=None, \
									hasHeader=True, valueDataType=float)
				kinshipMonkeyIDSet = set()
				for monkeyID, dataTuple in monkeyID2dataTuple.iteritems():
					if dataTuple[0]==0 and dataTuple[1]>30:
						kinshipMonkeyIDSet.add(monkeyID)
				sys.stderr.write("%s monkeys in kinshipMonkeyIDSet.\n"%(len(kinshipMonkeyIDSet)))
			else:
				kinshipMonkeyIDSet = None
			if self.outputFnamePrefix:
				self.cutOffKinshipIBDDeltaAndOutput(db_vervet=db_vervet, kinshipData=kinshipData, ibdData=ibdData, \
						outputFnamePrefix=self.outputFnamePrefix, minAbsDelta=self.minAbsDeltaForOutlier, kinshipMonkeyIDSet=kinshipMonkeyIDSet)
		
		#2012.8.24 output the delta matrix in PC1 order
		self.PCAOnAbsKinshipIBDDeltaMatrix(kinshipData=kinshipData,  ibdData=ibdData, outputFnamePrefix=self.outputFnamePrefix)
		
		if self.plinkSexCheckOutputFname:
			monkey_id2plinkSex = SNP.getKey2ValueFromMatrixLikeFile(inputFname=self.plinkSexCheckOutputFname, \
								keyHeaderLs=['IID'], valueHeaderLs=['SNPSEX'], keyIndexLs=None, valueIndexLs=None, \
								hasHeader=True, valueDataType=int)
		else:
			monkey_id2plinkSex = {}
		
		kinshipIBDDeltaData = self.createDeltaMatrix(kinshipData=kinshipData, ibdData=ibdData, takeAbs=False)
		
		meanStdData = self.estimateAbsDeltaMeanStd(kinshipIBDDeltaData=kinshipIBDDeltaData, excludeTopFraction=0.2)
		
		queueData = self.createKinshipIBDDeltaChiSqStatQueue(kinshipData=kinshipData, ibdData=ibdData, \
													mean=meanStdData.mean, std=meanStdData.std)
		
		queue = queueData.queue
		monkey_id2queueData = queueData.monkey_id2queueData
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['rank', 'monkeyID', 'chiSqStat', 'noOfNonMissing', 'chiSqPvalue', 'monkeySex','monkeyPlinkSex']
		#if self.iterativeAlgorithm:
		#	header.extend(['chiSqStatIter', 'noOfNonMissingIter', 'chiSqPvalueIter'])
		writer.writerow(header)
		
		i=0
		while i<5000 and len(queue)>0:
			
			minusChiSqStat, sourceMonkeyID, noOfNonMissing, chiSqPvalue = heapq.heappop(queue)[:4]
			chiSqStat = -minusChiSqStat
			sourceMonkeyDBEntry = self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=sourceMonkeyID)
			if sourceMonkeyDBEntry:
				sourceMonkeySex = sourceMonkeyDBEntry.codeSexInNumber()
			else:
				sourceMonkeySex = None
			sourceMonkeyPlinkSex = monkey_id2plinkSex.get(sourceMonkeyID)
			data_row = [i, sourceMonkeyID, chiSqStat, noOfNonMissing, chiSqPvalue, sourceMonkeySex, sourceMonkeyPlinkSex]
			
			
			if self.iterativeAlgorithm:
				"""
				if i>0:	#calculate the new chisq stat and p-value.
					chiSqStatData = self.calculateChiSqStatOfDeltaVector(kinshipData=kinshipData, kinshipDataMonkeyID=sourceMonkeyID, \
						ibdData=ibdData, ibdDataMonkeyID=sourceMonkeyID,\
						mean=meanStdData.mean, std=meanStdData.std)
					noOfNonMissing = chiSqStatData.noOfNonMissing
					chiSqStat = chiSqStatData.chiSqStat
					chiSqPvalue = chiSqStatData.chiSqPvalue
					data_row.extend([chiSqStat,noOfNonMissing,  chiSqPvalue])
				"""
				
				queueData = self.updateKinshipIBDDeltaChiSqStatQueue(queue=queue, kinshipData=kinshipData, ibdData=ibdData, \
								mean=meanStdData.mean, std=meanStdData.std, dropMonkeyID=sourceMonkeyID)
				#2012.8.23 old way not very efficient
				#remove itself.
#				ibdDataIndex = ibdData.row_id2row_index.get(sourceMonkeyID)
#				if ibdDataIndex:
#					ibdData.data_matrix[ibdDataIndex, :] = numpy.nan
#					ibdData.data_matrix[:, ibdDataIndex] = numpy.nan
#					ibdData.data_matrix.mask[ibdDataIndex, :] = True
#					ibdData.data_matrix.mask[:, ibdDataIndex] = True
				#queueData = self.createKinshipIBDDeltaChiSqStatQueue(kinshipData=kinshipData, ibdData=ibdData, \
				#									mean=meanStdData.mean, std=meanStdData.std,\
				#									given_row_id_ls=[row[1] for row in queue])
		
				queue = queueData.queue
				monkey_id2queueData = queueData.monkey_id2queueData
			writer.writerow(data_row)
			i+= 1
		del writer
	

if __name__ == '__main__':
	main_class = DetectWrongLabelByCompKinshipVsIBD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
