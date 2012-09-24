#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.21
	%s -i ~/NetworkData/vervet/Kinx2Apr2012.txt -l LDPrunedMerged.genome  -O Kinx2Apr2012VsIBDPI -u yh
		-D IBDCheckPIHat -x 2xKinship
	

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
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB
from pymodule.MatrixFile import MatrixFile
from pymodule import SNP
import numpy, random, pylab
import numpy as np
from pymodule.plot.AbstractPlot import AbstractPlot
from vervet.src import VervetDB

class PlotPedigreeKinshipVsGeneticIBD(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict.copy()
	option_default_dict.pop(('xColumnHeader', 1, ))
	#option_default_dict.pop(('xColumnPlotLabel', 0, ))
	option_default_dict.update({
						('plinkIBDCheckOutputFname', 1, ): ["", 'l', 1, 'file that contains IBD check result'], \
						('drivername', 1,):['postgresql', '', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'q', 1, 'database password', ],\
						('port', 0, ):[None, '', 1, 'database port number'],\
						('doPairwiseLabelCheck', 0, int):[0, '', 0, 'toggle to output two more tsv files: \
			_SumAbsDelta.tsv, _pairwiseCorOfKinshipIBDDelta.tsv'],\
					})
	
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		if self.outputFname and not self.outputFnamePrefix:
			self.outputFnamePrefix = os.path.splitext(self.outputFname)[0]
		
		
	def getMonkeyKinshipData(self, inputFname=None):
		"""
		2012.8.22
			use SNP.readAdjacencyListDataIntoMatrix(), and defaultValue=0
		2012.2.10
		"""
		
		sys.stderr.write("Reading kinship from %s ... "%(inputFname))
		kinshipData = SNP.readAdjacencyListDataIntoMatrix(inputFname=inputFname, rowIDHeader=None, colIDHeader=None, rowIDIndex=0, colIDIndex=1, \
								dataHeader=None, dataIndex=2, hasHeader=False, defaultValue=0)
		#set kinshipData diagonal to 1
		for i in xrange(len(kinshipData.row_id_ls)):
			kinshipData.data_matrix[i][i] = 1
		return kinshipData
		"""
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey1_id_index = col_name2index.get("monkeyId1")
		monkey2_id_index = col_name2index.get("monkeyId2")
		kinship_index = col_name2index.get("kinship")
		"""
		"""
		#reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		monkey1_id_index = 0
		monkey2_id_index = 1
		kinship_index = 2
		monkey_id_pair2kinship = {}
		for row in reader:
			monkey1_id = row[monkey1_id_index]
			monkey2_id = row[monkey2_id_index]
			kinship = float(row[kinship_index])
			pair_in_ls = [monkey1_id, monkey2_id]
			pair_in_ls.sort()
			monkey_id_pair2kinship[tuple(pair_in_ls)] = kinship
		del reader
		sys.stderr.write("%s pairs of monkeys.\n"%(len(monkey_id_pair2kinship)))
		return monkey_id_pair2kinship
		"""
	
	def getMonkeyIDPair2Correlation(self, smartpcaCorrelationFname=None):
		"""
		2012.3.1
			smartpcaCorrelationFname is output from  PCAOnVCFWorkflow.py (with modified smartpca). tab-delimited.
				553_2_VRC_ref_GA_vs_524	555_15_1987079_GA_vs_524        Case    Case    0.025
				553_2_VRC_ref_GA_vs_524	556_16_1985088_GA_vs_524        Case    Case    -0.020
				553_2_VRC_ref_GA_vs_524	557_17_1986014_GA_vs_524        Case    Case    -0.106
				553_2_VRC_ref_GA_vs_524	558_18_1988009_GA_vs_524        Case    Case    -0.059
		
		"""
		sys.stderr.write("Reading correlation from %s ... "%(smartpcaCorrelationFname))
		monkey_id_pair2genotype_correlation = {}
		import csv
		reader = csv.reader(open(smartpcaCorrelationFname), delimiter=figureOutDelimiter(smartpcaCorrelationFname))
		monkey_id_extract = lambda x: x.split('_')[2]
		for row in reader:
			monkey1 = row[0]
			monkey2 = row[1]
			cor = float(row[4])
			pair_in_ls = [monkey_id_extract(monkey1), monkey_id_extract(monkey2)]
			pair_in_ls.sort()
			pair_key = tuple(pair_in_ls)
			monkey_id_pair2genotype_correlation[pair_key] = cor
		sys.stderr.write("%s pairs .\n"%(len(monkey_id_pair2genotype_correlation)))
		return monkey_id_pair2genotype_correlation
	
		
	def getMonkeyIBDCheckData(self, inputFname=None):
		"""
		2012.8.21
			inputFname is output of plink ibd check.
 FID1     IID1 FID2     IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO
   1  1996093   1  1995025 OT     0  1.0000  0.0000  0.0000  0.0000  -1  0.654218  0.3630  1.9764
   1  1996093   1  2001039 OT     0  0.9832  0.0000  0.0168  0.0168  -1  0.653608  0.0318  1.8792
   1  1996093   1  1984011 OT     0  1.0000  0.0000  0.0000  0.0000  -1  0.645011  0.0168  1.8624
   1  1996093   1  1987004 OT     0  0.9260  0.0628  0.0113  0.0427  -1  0.660490  0.9999  2.2805
   		
		"""
		sys.stderr.write("Reading PI_hat from %s ... "%(inputFname))
		ibdData = SNP.readAdjacencyListDataIntoMatrix(inputFname=inputFname, rowIDHeader="IID1", colIDHeader="IID2", rowIDIndex=None, colIDIndex=None, \
								dataHeader="PI_HAT", dataIndex=None, hasHeader=True)
		return ibdData
		"""
		monkey_id_pair2genotype_correlation = {}
		
		reader = MatrixFile(inputFname=inputFname)
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		monkey1Index = col_name2index.get('IID1')
		monkey2Index = col_name2index.get('IID2')
		PI_HATIndex = col_name2index.get("PI_HAT")
		
		for row in reader:
			monkey1 = row[monkey1Index]
			monkey2 = row[monkey2Index]
			PI_HAT = float(row[PI_HATIndex])
			pair_in_ls = [monkey1, monkey2]
			pair_in_ls.sort()
			pair_key = tuple(pair_in_ls)
			monkey_id_pair2genotype_correlation[pair_key] = PI_HAT
		sys.stderr.write("%s pairs .\n"%(len(monkey_id_pair2genotype_correlation)))
		return monkey_id_pair2genotype_correlation
		"""
	
	def getMonkeyDBEntry(self, db_vervet=None, ucla_id=None):
		"""
		2012.8.21
		"""
		return db_vervet.getIndividualDBEntry(ucla_id=ucla_id)
		
	
	def plotPairwiseKinshipFromPedigreeVsGenotype(self, db_vervet=None, kinshipFname=None, plinkIBDCheckOutputFname=None, \
									outputFnamePrefix=None, doPairwiseLabelCheck=False):
		"""
		2012.8.17
			do full join, output pairs with data in only one source
		2012.3.1
			smartpcaCorrelationFname is output from  PCAOnVCFWorkflow.py (with modified smartpca). tab-delimited.
			kinshipFname is from Sue (estimated by SOLAR based on pedigree)
			
		"""
		ibdData = self.getMonkeyIBDCheckData(plinkIBDCheckOutputFname)
		kinshipData = self.getMonkeyKinshipData(kinshipFname)
		
		tableOutputWriter = csv.writer(open("%s_table.tsv"%(outputFnamePrefix), 'w'), delimiter='\t')
		header = ['monkey_pair', 'pedigree_kinship', 'IBD_PI_HAT', 'kinshipIBDDelta','age_difference']
		tableOutputWriter.writerow(header)
		
		monkey_id_pair_ls = []
		x_ls = []
		y_ls = []
		
		monkey_id2index = {}	#2012.8.21
		monkey_id_pair2kinship_ibd_delta = {}
		shared_monkey_id_set = set(ibdData.row_id_ls)&set(kinshipData.row_id_ls)
		shared_monkey_id_ls = list(shared_monkey_id_set)
		no_of_monkeys = len(shared_monkey_id_ls)
		for i in xrange(no_of_monkeys):
			monkey1_id = shared_monkey_id_ls[i]
			
			for j in xrange(i+1, no_of_monkeys):
				monkey2_id = shared_monkey_id_ls[j]
				
				ibd = ibdData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				kinship = kinshipData.getCellDataGivenRowColID(monkey1_id, monkey2_id)
				if (not hasattr(ibd, 'mask')) and not numpy.isnan(ibd) and (not hasattr(kinship, 'mask'))\
						and (not numpy.isnan(kinship)):
					x_ls.append(kinship)
					y_ls.append(ibd)
					monkey_id_pair = [monkey1_id, monkey2_id]
					monkey_id_pair.sort()
					monkey_id_pair  = tuple(monkey_id_pair)
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
					data_row = ['_'.join(monkey_id_pair), kinship, ibd, kinship-ibd, ageDelta]
					
					tableOutputWriter.writerow(data_row)
					monkey_id_pair2kinship_ibd_delta[monkey_id_pair] = kinship-ibd
		del tableOutputWriter
		
		sys.stderr.write("%s pairs overlap between pedigree kinship and genotype correlation, among %s monkeys.\n"%\
						(len(x_ls), len(monkey_id2index)))
		from pymodule import yh_matplotlib
		fig_fname = '%s_scatter.png'%(outputFnamePrefix)
		yh_matplotlib.drawScatter(x_ls, y_ls, fig_fname=fig_fname, title=None, xlabel=self.xColumnPlotLabel, \
								ylabel=self.whichColumnPlotLabel, dpi=self.figureDPI)
		
		if doPairwiseLabelCheck:
			sys.stderr.write("Pairwise calculation using kinship-ibdcheck ...")
			no_of_monkeys = len(monkey_id2index)
			kinship_ibd_deltaMatrix = numpy.zeros([no_of_monkeys, no_of_monkeys], dtype=numpy.float32)
			kinship_ibd_deltaMatrix[:] = numpy.nan
			for monkey_id_pair, kinship_ibd_delta in  monkey_id_pair2kinship_ibd_delta.iteritems():
				monkey1_id = monkey_id_pair[0]
				monkey2_id = monkey_id_pair[1]
				monkey1_index = monkey_id2index.get(monkey1_id)
				monkey2_index = monkey_id2index.get(monkey2_id)
				kinship_ibd_deltaMatrix[monkey1_index][monkey2_index] = kinship_ibd_delta
				kinship_ibd_deltaMatrix[monkey2_index][monkey1_index] = kinship_ibd_delta
			
			maskedMatrix = np.ma.array(kinship_ibd_deltaMatrix, mask=np.isnan(kinship_ibd_deltaMatrix))
			outputFname = '%s_SumAbsDelta.tsv'%(outputFnamePrefix)
			writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
			header = ['monkey_pair', 'index', 'sumAbsDelta', 'noOfNonMissing', 'avgAbsDelta', 'medianAbsDelta']
			writer.writerow(header)
			monkey_id2medianAbsDelta = {}
			for monkey_id , index in monkey_id2index.iteritems():
				
				sumAbsDelta = numpy.ma.sum(numpy.ma.abs(maskedMatrix[index,:]))
				noOfNonMissing = len(maskedMatrix[index,:]) - numpy.ma.sum(maskedMatrix[index,:].mask)
				medianAbsDelta = numpy.ma.median(numpy.ma.abs(maskedMatrix[index,:]))
				data_row = [monkey_id, index, sumAbsDelta, noOfNonMissing, sumAbsDelta/float(noOfNonMissing), medianAbsDelta]
				monkey_id2medianAbsDelta[monkey_id] = medianAbsDelta
				writer.writerow(data_row)
			del writer
			
			outputFname = '%s_pairwiseCorOfKinshipIBDDelta.tsv'%(outputFnamePrefix)
			writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
			header = ['monkey_pair', 'noOfValidPairs', 'corr', 'monkey1_medianAbsDelta', 'monkey2_medianAbsDelta']
			writer.writerow(header)
			monkey_id_ls = monkey_id2index.keys()
			monkey_id_ls.sort()
			no_of_monkeys = len(monkey_id_ls)
			for i in xrange(no_of_monkeys):
				monkey1_id =  monkey_id_ls[i]
				monkey1_index = monkey_id2index.get(monkey1_id)
				for j in xrange(i+1, no_of_monkeys):
					monkey2_id =  monkey_id_ls[j]
					monkey2_index = monkey_id2index.get(monkey2_id)
					monkey1_vector = maskedMatrix[monkey1_index,:]
					monkey2_vector = maskedMatrix[monkey2_index,:]
					orMask = numpy.ma.mask_or(monkey1_vector.mask, monkey2_vector.mask)
					noOfValidPairs = len(monkey1_vector) - numpy.sum(orMask)
					corr = numpy.ma.corrcoef(monkey1_vector, monkey2_vector, \
									rowvar=True)[0,1]
					if not hasattr(corr, 'mask'):	#valid float correlation, not a masked array cell (which means not enough valid data) 
					#	if (abs(corr)>0.2):
							#output the plot
						data_row = ['%s_%s'%(monkey1_id, monkey2_id), noOfValidPairs, corr, monkey_id2medianAbsDelta.get(monkey1_id),\
							monkey_id2medianAbsDelta.get(monkey2_id)]
						writer.writerow(data_row)
			del writer
			sys.stderr.write(".\n")
		
		"""
		#2012.3.1
		kinshipFname = "/Network/Data/vervet/Kinx2Jan2012.txt"
		smartpcaCorrelationFname= "/Network/Data/vervet/vervetPipeline/PCAOnVCFWorkflow_VRC105_Top1000Contigs.2012.3.1T1505/pca/smartpca.cor"
		outputFnamePrefix = os.path.expanduser("~/script/vervet//data/VRC105_Top1000Contigs_pedigree_kinship_vs_genotype_correlation")
		DBVervet.plotPairwiseKinshipFromPedigreeVsGenotype(db_vervet, kinshipFname=kinshipFname, \
														smartpcaCorrelationFname=smartpcaCorrelationFname,\
														outputFnamePrefix=outputFnamePrefix)
		sys.exit(3)
		
		"""
	
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
		
		self.plotPairwiseKinshipFromPedigreeVsGenotype(db_vervet=db_vervet, kinshipFname=self.inputFname, \
									plinkIBDCheckOutputFname=self.plinkIBDCheckOutputFname, \
									outputFnamePrefix=self.outputFnamePrefix,\
									doPairwiseLabelCheck=self.doPairwiseLabelCheck)
	

if __name__ == '__main__':
	main_class = PlotPedigreeKinshipVsGeneticIBD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
