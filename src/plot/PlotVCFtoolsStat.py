#!/usr/bin/env python
"""
Examples:
	%s -o outputFname input1 input2 ...
	
	%s 
	

Description:
	2011-11-28
		this program draws a manhattan plot (gwas plot) and a histogram for some vcftools outputted windowed statistics.
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

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB
from PlotTrioInconsistencyOverFrequency import PlotTrioInconsistencyOverFrequency
import numpy

class PlotVCFtoolsStat(PlotTrioInconsistencyOverFrequency):
	__doc__ = __doc__
	option_default_dict = PlotTrioInconsistencyOverFrequency.option_default_dict
	option_default_dict.update({('whichColumn', 1, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
							('need_svg', 0, ): [0, '', 0, 'whether need svg output', ]})
	option_for_DB_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ]}
	option_default_dict.update(option_for_DB_dict)
	
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		PlotTrioInconsistencyOverFrequency.__init__(self, inputFnameLs, **keywords)
		#Super(PlotVCFtoolsStat, self).__init__(inputFnameLs, **keywords)
		self.chr2xy_ls = {}
		
	def vcftoolsOutputStatFileWalker(self, inputFname, processFunc=None, run_type=1, minChrLength=1000000):
		"""
		2011-11-2
			remove the maxDepth filter. apply afterwards through filterDataByDepth().
		2011-9-30
		
		"""
		chr2xy_ls = self.chr2xy_ls
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		chr_id_index = col_name2index.get("CHROM")
		bin_start_index = col_name2index.get("BIN_START")
		chrLength_index = col_name2index.get("chrLength")
		for row in reader:
			chr_id = row[chr_id_index]
			bin_start = int(float(row[bin_start_index]))
			chrLength = int(row[chrLength_index])
			if chrLength<minChrLength:
				continue
			yValue = float(row[self.whichColumn])
			if chr_id not in chr2xy_ls:
				chr2xy_ls[chr_id] = [[],[]]
			chr_cumu_start = self.chr_id2cumu_start.get(chr_id)
			
			chr2xy_ls[chr_id][0].append(chr_cumu_start + bin_start + 1)
			chr2xy_ls[chr_id][1].append(yValue)
		del reader
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#without commenting out db_vervet connection code. schema "genome" wont' be default path.
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		#chrOrder=2 means chromosomes are not ordered alphabetically but by their sizes (descendingly)
		oneGenomeData = db_genome.getOneGenomeData(tax_id=60711, chr_gap=0, chrOrder=2, sequence_type_id=9)
		chr2size = db_genome.getTopNumberOfChomosomes(contigMaxRankBySize=80000, contigMinRankBySize=1, tax_id=60711, \
											sequence_type_id=9)
		
		self.chr_id2cumu_start = oneGenomeData.chr_id2cumu_start
		size_chr_id_ls = [(value, key) for key, value in chr2size.iteritems()]
		size_chr_id_ls.sort()
		size_chr_id_ls.reverse()
		
		sys.stderr.write("Reading in data ...")
		for inputFname in self.inputFnameLs:
			if not os.path.isfile(inputFname):
				continue
			try:
				self.vcftoolsOutputStatFileWalker(inputFname, processFunc=None)
			except:	#in case something wrong (i.e. file is empty)
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				print sys.exc_info()
		sys.stderr.write("Done.\n")
		
		import pylab
		pylab.clf()
		fig = pylab.figure(figsize=(20,2))
		#ax = pylab.axes()
		ax = fig.gca()
		import numpy
		max_y = None
		min_y = None
		value_ls = []
		for size, chr in size_chr_id_ls:
			xy_ls = self.chr2xy_ls.get(chr)
			if xy_ls:
				if max_y is None:
					max_y = max(xy_ls[1])
				else:
					max_y = max(max_y, max(xy_ls[1]))
				if min_y is None:
					min_y = min(xy_ls[1])
				else:
					min_y = min(min_y, min(xy_ls[1]))
				ax.plot(xy_ls[0], xy_ls[1], '.', markeredgewidth=0, markersize=4, alpha=0.8)
				value_ls += xy_ls[1]
		#separate each chromosome
		#for chr in chr_ls[:-1]:
		#	print chr
		#	ax.axvline(chr_id2cumu_size[chr], linestyle='--', color='k', linewidth=0.8)
		
		
		#draw the bonferroni line
		#bonferroni_value = -math.log10(0.01/len(genome_wide_result.data_obj_ls))
		#ax.axhline(bonferroni_value, linestyle='--', color='k', linewidth=0.8)
		
		#ax.set_ylabel("-log(P-value)")
		#ax.set_xlabel('Chromosomal Position')
		#ax.set_xlim([0, chr_id2cumu_size[chr_ls[-1]]])
		ylim_type=1
		if ylim_type==1:
			ylim = ax.get_ylim()
			ax.set_ylim([0, ylim[1]])
		elif ylim_type==2:
			ax.set_ylim([min_y, max_y])
		
		outputFnamePrefix = os.path.splitext(self.outputFname)[0]
		pylab.savefig('%s.png'%outputFnamePrefix, dpi=self.figureDPI)
		if self.need_svg:
			pylab.savefig('%s.svg'%outputFnamePrefix, dpi=self.figureDPI)
		outputFname = '%s_hist.png'%(outputFnamePrefix)
		yh_matplotlib.drawHist(value_ls, title='', \
				xlabel_1D="#SNPs in 100kb window", xticks=None, outputFname=outputFname, min_no_of_data_points=20, needLog=True, \
				dpi=self.figureDPI, min_no_of_bins=200)

if __name__ == '__main__':
	main_class = PlotVCFtoolsStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()