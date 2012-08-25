#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.21
	%s -i ~/NetworkData/vervet/Kinx2Apr2012.txt -u yh
		-l PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/LDPrunedMerged.genome
		-n PlinkSexCheck/PlinkSexCheck_Method18_W100Z10R0.8_maxContigID195.2012.8.21T1207/sexCheckSexCheck/LDPrunedMerged.sexcheck
		-o PlinkIBDCheck/PlinkIBDCheck_Method18_W100Z10R0.4_vsKinship.2012.8.21T1905/ibdCheckIBDCheck/crossLabelCheck.tsv
	

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
from pymodule import yh_matplotlib, GenomeDB
from pymodule.MatrixFile import MatrixFile
from pymodule import SNP
import numpy, random, pylab
import numpy as np
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src.plot.DetectWrongLabelByCompKinshipVsIBD import DetectWrongLabelByCompKinshipVsIBD
import heapq

class DetectLabelSwapByCompKinshipVsIBD(DetectWrongLabelByCompKinshipVsIBD):
	__doc__ = __doc__
#						
	option_default_dict = DetectWrongLabelByCompKinshipVsIBD.option_default_dict
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		DetectWrongLabelByCompKinshipVsIBD.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
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
		ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=self.plinkIBDCheckOutputFname, id1Header="IID1", id2Header="IID2", id1Index=None, id2Index=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		monkey_id2plinkSex = SNP.getKey2ValueFromMatrixLikeFile(inputFname=self.plinkSexCheckOutputFname, \
								keyHeaderLs=['IID'], valueHeaderLs=['SNPSEX'], keyIndexLs=None, valueIndexLs=None, \
								hasHeader=True, valueDataType=int)
		
		kinshipIBDDeltaData = self.createDeltaMatrix(kinshipData=kinshipData, ibdData=ibdData)
		queueData = self.createKinshipIBDDeltaQueue(kinshipIBDDeltaData)
		kinshipIBDDeltaQueue = queueData.kinshipIBDDeltaQueue
		monkey_id2medianAbsDelta = queueData.monkey_id2medianAbsDelta
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		header = ['sourceMonkeyID', 'medianAbsDelta', 'noOfNonMissing', 'sourceMonkeySex', 'sourceMonkeyPlinkSex', \
				'sourceMonkeyMedianAbsDeltaDropAfterSwap', \
			'targetMonkeyID', 'sourceMonkeyNoOfNonMissingAfterSwap', \
			'targetMonkeyMedianAbsDelta', 'targetMonkeyNoOfNonMissing', 'targetMonkeySex', 'targetMonkeyPlinkSex', \
			'targetMonkeyMedianAbsDeltaAfterSwap', 'targetMonkeyNoOfNonMissingAfterSwap']
		writer.writerow(header)
		
		i=0
		while i <50 and len(kinshipIBDDeltaQueue)>0:
			negativeMedianAbsDelta, sourceMonkeyID, noOfNonMissing = heapq.heappop(kinshipIBDDeltaQueue)[:3]
			medianAbsDelta = -negativeMedianAbsDelta
			sourceMonkeyDBEntry = self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=sourceMonkeyID)
			
			# 2012.8.22 draw some histogram to check what data looks like
			#self.drawKinshipIBDDeltaVectorHistogram(kinshipIBDDeltaData=kinshipIBDDeltaData, row_id=sourceMonkeyID, \
			#							outputFnamePrefix=self.outputFnamePrefix)
			
			medianAbsDeltaIncreaseQueue = []
			for targetMonkeyID in kinshipData.row_id_ls:
				if targetMonkeyID!=sourceMonkeyID:
					targetMonkeyDBEntry = self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=targetMonkeyID)
					#get the updated Median Delta for sourceMonkeyID
					pdata = self.calculateMedianAbsDelta(kinshipData=kinshipData, \
										kinshipDataMonkeyID=targetMonkeyID, ibdData=ibdData, ibdDataMonkeyID=sourceMonkeyID)
					sourceMonkeyMedianAbsDeltaAfterSwap = pdata.medianAbsDelta
					sourceMonkeyNoOfNonMissingAfterSwap  = pdata.noOfNonMissing
					
					#get the updated Median Delta for targetMonkeyID
					pdata = self.calculateMedianAbsDelta(kinshipData=kinshipData, \
										kinshipDataMonkeyID=sourceMonkeyID, ibdData=ibdData, ibdDataMonkeyID=targetMonkeyID)
					targetMonkeyMedianAbsDeltaAfterSwap = pdata.medianAbsDelta
					targetMonkeyNoOfNonMissingAfterSwap = pdata.noOfNonMissing
					
					if sourceMonkeyMedianAbsDeltaAfterSwap is not None:	#add to the queue
						#add the candidate monkey and how much median delta drops into the queue
						pdata = monkey_id2medianAbsDelta.get(targetMonkeyID)
						if pdata:
							targetMonkeyMedianAbsDelta = pdata.medianAbsDelta
							targetMonkeyNoOfNonMissing = pdata.noOfNonMissing
						else:
							targetMonkeyMedianAbsDelta = None
							targetMonkeyNoOfNonMissing = None
						item = [sourceMonkeyMedianAbsDeltaAfterSwap-medianAbsDelta, targetMonkeyID, sourceMonkeyNoOfNonMissingAfterSwap, \
							targetMonkeyMedianAbsDelta, targetMonkeyNoOfNonMissing, targetMonkeyMedianAbsDeltaAfterSwap, \
							targetMonkeyNoOfNonMissingAfterSwap]
						heapq.heappush(medianAbsDeltaIncreaseQueue, item)
				
			#the target monkey that increase the least (or drop the most) for the median delta is the prime candidate for label-swap 
			i+=1
			#output the top 5 candidates for each source monkey
			#output db sex for all monkeys and the plink sex check result
			j = 0
			while j<5 and len(medianAbsDeltaIncreaseQueue)>0:
				sourceMonkeyMedianAbsDeltaDropAfterSwap, targetMonkeyID, sourceMonkeyNoOfNonMissingAfterSwap, \
						targetMonkeyMedianAbsDelta, targetMonkeyNoOfNonMissing, targetMonkeyMedianAbsDeltaAfterSwap, targetMonkeyNoOfNonMissingAfterSwap =\
							heapq.heappop(medianAbsDeltaIncreaseQueue)[:7]
				sourceMonkeySex = sourceMonkeyDBEntry.codeSexInNumber()
				sourceMonkeyPlinkSex = monkey_id2plinkSex.get(sourceMonkeyID)
				
				targetMonkeySex = targetMonkeyDBEntry.codeSexInNumber()
				targetMonkeyPlinkSex = monkey_id2plinkSex.get(targetMonkeyID)
				
				data_row = [sourceMonkeyID, medianAbsDelta, noOfNonMissing, sourceMonkeySex, sourceMonkeyPlinkSex,\
						sourceMonkeyMedianAbsDeltaDropAfterSwap, targetMonkeyID, sourceMonkeyNoOfNonMissingAfterSwap,\
						targetMonkeyMedianAbsDelta, targetMonkeyNoOfNonMissing, targetMonkeySex, targetMonkeyPlinkSex, targetMonkeyMedianAbsDeltaAfterSwap,\
						targetMonkeyNoOfNonMissingAfterSwap]
				writer.writerow(data_row)
				j+= 1
		del writer
	

if __name__ == '__main__':
	main_class = DetectLabelSwapByCompKinshipVsIBD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
