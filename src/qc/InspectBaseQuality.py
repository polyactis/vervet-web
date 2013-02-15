#!/usr/bin/env python
"""
Examples:
	%s -i /Network/Data/vervet/db/individual_sequence/3_Barbados_GA/gerald_62FGFAAXX_3_1.fastq.gz -u yh -r  -e 0.0001
	
	%s 

Description:
	2011-8-15 Five types of figures will be generated.
		qualityHist.png
		quality_per_position.png
		no_of_bases_per_position.png
		diNuc_count.png
		diNuc_quality.png
	2011-8-17 the db saving part is not implemented, you can supply any random password.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import random
import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
from pymodule import PassingData, ProcessOptions, utils, yh_matplotlib
from vervet.src import VervetDB

class InspectBaseQuality(object):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('ind_sequence_id', 0, int): [None, 'n', 1, 'The individual_sequence id corresponding to the inputFname, if not given, no db commit action.', ],\
						('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('quality_score_format', 1, ): ['Sanger', 'q', 1, 'could be Sanger (phred+33), or Illumina1.3 (~phred+64). Illumina1.8+ (after 2011-02) is Sanger.', ],\
						('outputFnamePrefix', 0, ): [None, 'o', 1, 'output filename prefix. if not given, derive prefix from inputFname'],\
						('read_sampling_rate', 0, float): [0.05, 'e', 1, 'the probability that a read will be included.', ],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path) and qsub jobs'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		if not self.outputFnamePrefix:
			self.outputFnamePrefix = utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(self.inputFname)[0]
	
	def getQualityData(self, inputFname, read_sampling_rate=0.05, quality_score_format='Sanger'):
		"""
		2011-8-15
		"""
		sys.stderr.write("Getting base quality data from %s ...\n"%(inputFname))
		
		quality_ls_per_position = []
		quality_ls = []
		no_of_bases_per_position = []
		diNuc2count = {}
		diNuc2quality_ls = {}
		
		
		fname_prefix, fname_suffix = os.path.splitext(inputFname)
		if fname_suffix=='.gz':	#the input file is gzipped. get the new prefix
			import gzip
			inf = gzip.open(inputFname, 'rb')
		else:
			inf = open(inputFname, 'r')
		counter = 0
		real_counter = 0
		for line in inf:
			if line[0]=='@':	#a new read
				counter += 1
				coin_toss = random.random()
				base_string = inf.next().strip()
				inf.next()
				quality_string = inf.next().strip()
				if coin_toss<=read_sampling_rate:
					real_counter += 1
					read_length = len(base_string)
					if len(quality_ls_per_position)<read_length:	# extend quality_ls_per_position to house more data
						extraNoOfBases = read_length-len(quality_ls_per_position)
						for j in xrange(extraNoOfBases):
							quality_ls_per_position.append([])
							no_of_bases_per_position.append(0)
						
					for i in range(read_length):
						base = base_string[i]
						base_quality = quality_string[i]
						if quality_score_format=='Illumina1.3':
							phredScore = utils.getPhredScoreOutOfSolexaScore(base_quality)
						else:
							phredScore = ord(base_quality)-33
						quality_ls_per_position[i].append(phredScore)
						quality_ls.append(phredScore)
						if base!='N':
							no_of_bases_per_position[i] += 1
							if i<read_length-1:
								nextBase = base_string[i+1]
								if nextBase!='N':
									diNuc = base + nextBase
									if diNuc not in diNuc2quality_ls:
										diNuc2quality_ls[diNuc] = []
										diNuc2count[diNuc] = 0
									diNuc2quality_ls[diNuc].append(phredScore)
									diNuc2count[diNuc] += 1
			if counter%5000==0 and self.report:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, real_counter, counter))
			
			#if baseCount>10000:	#temporary, for testing
			#	break
		del inf
		sys.stderr.write("%s/%s reads selected. Done.\n"%(real_counter, counter))
		return PassingData(quality_ls_per_position=quality_ls_per_position, quality_ls=quality_ls, \
					no_of_bases_per_position=no_of_bases_per_position, diNuc2quality_ls=diNuc2quality_ls, diNuc2count=diNuc2count)
	
	def drawQualityData(self, qualityDataStructure, outputFnamePrefix, sequence_id=''):
		"""
		2011-8-15
		"""
		sys.stderr.write("Making plots on quality data ...")
		
		yh_matplotlib.drawHist(qualityDataStructure.quality_ls, title='histogram of phredScore from %s'%(sequence_id), xlabel_1D=None, \
							outputFname='%s_qualityHist.png'%(outputFnamePrefix), \
							min_no_of_data_points=50, needLog=False, dpi=200)
		
		yh_matplotlib.drawBoxPlot(qualityDataStructure.quality_ls_per_position, title='quality box plot from %s'%(sequence_id), \
								xlabel_1D='base position in read', xticks=None, outputFname='%s_quality_per_position.png'%(outputFnamePrefix), \
								dpi=200)
		
		no_of_bases_per_position = qualityDataStructure.no_of_bases_per_position
		readLength = len(no_of_bases_per_position)
		yh_matplotlib.drawBarChart(range(1, readLength+1), no_of_bases_per_position, title='no of base calls from %s'%(sequence_id),\
						xlabel_1D='base position in read', xticks=None, outputFname='%s_no_of_bases_per_position.png'%(outputFnamePrefix), \
						bottom=0, needLog=False, dpi=200)
		
		diNuc2count = qualityDataStructure.diNuc2count
		diNuc2quality_ls = qualityDataStructure.diNuc2quality_ls
		
		diNuc_key_ls = diNuc2count.keys()
		diNuc_key_ls.sort()
		diNuc_count_ls = []
		diNuc_quality_ls_ls = []
		for diNuc in diNuc_key_ls:
			diNuc_count_ls.append(diNuc2count.get(diNuc))
			diNuc_quality_ls_ls.append(diNuc2quality_ls.get(diNuc))
		
		yh_matplotlib.drawBarChart(range(1, len(diNuc_count_ls)+1), diNuc_count_ls, title='di-nucleotide counts from %s'%(sequence_id),\
						xlabel_1D=None, xticks=diNuc_key_ls, outputFname='%s_diNuc_count.png'%(outputFnamePrefix), \
						bottom=0, needLog=False, dpi=200)
		
		yh_matplotlib.drawBoxPlot(diNuc_quality_ls_ls, title='di-Nucleotide quality box plot from %s'%(sequence_id), \
								xlabel_1D=None, xticks=diNuc_key_ls, outputFname='%s_diNuc_quality.png'%(outputFnamePrefix), \
								dpi=200)
		
		sys.stderr.write("Done.\n")
	
	def saveDataIntoDB(self, db_vervet, ind_sequence_id=None):
		"""
		2011-8-15
		"""
		pass
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		session = db_vervet.session
		session.begin()	#no transaction for input node as there is no data insertion
		
		qualityDataStructure = self.getQualityData(self.inputFname, read_sampling_rate=self.read_sampling_rate,\
									quality_score_format=self.quality_score_format)
		
		sequence_id = os.path.split(self.outputFnamePrefix)[1]	#to be part of title in each figure
		self.drawQualityData(qualityDataStructure, self.outputFnamePrefix, sequence_id=sequence_id)
		
		self.saveDataIntoDB(db_vervet, ind_sequence_id=self.ind_sequence_id)
			
		if self.commit:
			session.commit()

if __name__ == '__main__':
	main_class = InspectBaseQuality
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
