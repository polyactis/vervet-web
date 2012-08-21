#!/usr/bin/env python
"""
Examples:
	%s 
	
	# 2012.8.21
	%s -i ~/NetworkData/vervet/Kinx2Apr2012.txt -l LDPrunedMerged.genome  -o Kinx2Apr2012VsIBDPI.tsv -u yh 
	

Description:
	2012.8.21
		Program that plots kinship vs. PI_hat from plinkIBD result.
		inputFname is the kinship file from Sue.
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
import numpy, random, pylab
from pymodule.plot.AbstractPlot import AbstractPlot
from vervet.src import VervetDB

class PlotPedigreeKinshipVsGeneticIBD(AbstractPlot):
	__doc__ = __doc__
#						
	option_default_dict = AbstractPlot.option_default_dict
	option_default_dict.pop(('xColumnHeader', 1, ))
	option_default_dict.pop(('xColumnPlotLabel', 0, ))
	option_default_dict.update({
						('plinkIBDCheckFname', 1, ): ["", 'l', 1, 'file that contains IBD check result'], \
						('drivername', 1,):['postgresql', '', 1, 'which type of database? mysql or postgresql', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
							('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'q', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							('logFilename', 0, ): [None, '', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow \
		and has no output file'],\
							("dataDir", 0, ): ["", '', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractPlot.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		if self.outputFname and not self.outputFnamePrefix:
			self.outputFnamePrefix = os.path.splitext(self.outputFname)[0]
		
		self.uclaID2monkeyDBEntry = {}
		
	def getMonkeyIDPair2Kinship(self, inputFname):
		"""
		2012.2.10
		"""
		
		sys.stderr.write("Reading kinship from %s ... "%(inputFname))
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		"""
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		monkey1_id_index = col_name2index.get("monkeyId1")
		monkey2_id_index = col_name2index.get("monkeyId2")
		kinship_index = col_name2index.get("kinship")
		"""
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
	
		
	def getMonkeyIDPair2PIHAT(self, inputFname=None):
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
	
	def getMonkeyDBEntry(self, db_vervet=None, ucla_id=None):
		"""
		2012.8.21
		"""
		if ucla_id in self.uclaID2monkeyDBEntry:
			return self.uclaID2monkeyDBEntry.get(ucla_id)
		else:
			monkey = VervetDB.Individual.query.filter_by(code=ucla_id).first()
			self.uclaID2monkeyDBEntry[ucla_id] = monkey
			return monkey
			
	
	def plotPairwiseKinshipFromPedigreeVsGenotype(self, db_vervet=None, kinshipFname=None, plinkIBDCheckFname=None, outputFnamePrefix=None):
		"""
		2012.8.17
			do full join, output pairs with data in only one source
		2012.3.1
			smartpcaCorrelationFname is output from  PCAOnVCFWorkflow.py (with modified smartpca). tab-delimited.
			kinshipFname is from Sue (estimated by SOLAR based on pedigree)
			
		"""
		monkey_id_pair2genotype_correlation = self.getMonkeyIDPair2PIHAT(plinkIBDCheckFname)
		monkey_id_pair2pedigree_kinship = self.getMonkeyIDPair2Kinship(kinshipFname)
		
		tableOutputWriter = csv.writer(open("%s_table.tsv"%(outputFnamePrefix), 'w'), delimiter='\t')
		header = ['monkey_pair', 'pedigree_kinship', 'IBD_PI_HAT', 'age_difference']
		tableOutputWriter.writerow(header)
		
		monkey_id_pair_ls = []
		x_ls = []
		y_ls = []
		for monkey_id_pair, kinship in monkey_id_pair2pedigree_kinship.iteritems():
			if monkey_id_pair in monkey_id_pair2genotype_correlation:
				genotype_cor = monkey_id_pair2genotype_correlation.get(monkey_id_pair)
				x_ls.append(kinship)
				y_ls.append(genotype_cor)
				monkey1 = self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=monkey_id_pair[0])
				monkey2 =  self.getMonkeyDBEntry(db_vervet=db_vervet, ucla_id=monkey_id_pair[1])
				
				if monkey1 and monkey2 and monkey1.getCurrentAge() and monkey2.getCurrentAge():
					ageDelta = abs(monkey1.getCurrentAge() - monkey2.getCurrentAge())
				else:
					ageDelta = ''
				data_row = ['_'.join(monkey_id_pair), kinship, genotype_cor, ageDelta]
				tableOutputWriter.writerow(data_row)
		del tableOutputWriter
		
		sys.stderr.write("%s pairs overlap between pedigree kinship and genotype correlation.\n"%(len(x_ls)))
		from pymodule import yh_matplotlib
		fig_fname = '%s_scatter.png'%(outputFnamePrefix)
		yh_matplotlib.drawScatter(x_ls, y_ls, fig_fname=fig_fname, title='%s pairs'%(len(x_ls)), xlabel='pedigree kinship', \
								ylabel='IBD PI_HAT', dpi=300)
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
		
		self.plotPairwiseKinshipFromPedigreeVsGenotype(db_vervet=db_vervet, kinshipFname=self.inputFname, plinkIBDCheckFname=self.plinkIBDCheckFname, \
										outputFnamePrefix=self.outputFnamePrefix)
	

if __name__ == '__main__':
	main_class = PlotPedigreeKinshipVsGeneticIBD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
