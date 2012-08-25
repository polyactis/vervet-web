#!/usr/bin/env python
"""
Examples:
	%s -o outputFname input1 input2 ...
	
	%s  -s 0.1 -O hweplot -W "P" -l POS -g -C CHR -D "hweLogPvalue" -x position 
		-z dl324b-1.cmb.usc.edu -u yh /tmp/5988_VCF_Contig966.hwe
	

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
from pymodule import yh_matplotlib, GenomeDB, utils
from PlotTrioInconsistencyOverFrequency import PlotTrioInconsistencyOverFrequency
import numpy, random, pylab

class PlotVCFtoolsStat(PlotTrioInconsistencyOverFrequency):
	__doc__ = __doc__
	option_default_dict = PlotTrioInconsistencyOverFrequency.option_default_dict
	option_default_dict.pop(('outputFname', 1, ))
	option_default_dict.update({
			('whichColumn', 0, int): [3, 'w', 1, 'data from this column (index starting from 0) is plotted as y-axis value'],\
			('whichColumnLabel', 0, ): ["", 'W', 1, 'column label (in the header) for the data to be plotted as y-axis value, substitute whichColumn'],\
			('logWhichColumn', 0, int): [0, 'g', 0, 'whether to take -log of the whichColumn'],\
			('need_svg', 0, ): [0, 'n', 0, 'whether need svg output', ],\
			('whichColumnPlotLabel', 1, ): ['#SNPs in 100kb window', 'D', 1, 'plot label for data of the whichColumn', ],\
			('posColumnPlotLabel', 1, ): ['position', 'x', 1, 'x-axis label (posColumn) in manhattan plot', ],\
			('chrLengthColumnLabel', 1, ): ['chrLength', 'c', 1, 'label of the chromosome length column', ],\
			('chrColumnLabel', 1, ): ['CHR', 'C', 1, 'label of the chromosome column', ],\
			('minChrLength', 1, int): [1000000, 'm', 1, 'minimum chromosome length for one chromosome to be included', ],\
			('posColumnLabel', 1, ): ['BIN_START', 'l', 1, 'label of the position column, BIN_START for binned vcftools output. POS for others.', ],\
			('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional).'],\
			('logCount', 0, int): [0, '', 0, 'whether to take log on the y-axis of the histogram, the raw count'], \
			})
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
		
	def vcftoolsOutputStatFileWalker(self, inputFname, processFunc=None, run_type=1, \
									chrColumnLabel='CHR', minChrLength=1000000, chrLengthColumnLabel='chrLength',\
									posColumnLabel="BIN_START"):
		"""
		2012.8.13 bugfix. pass inf to figureOutDelimiter
		2012.8.1
		2011-11-2
			remove the maxDepth filter. apply afterwards through filterDataByDepth().
		2011-9-30
		
		"""
		sys.stderr.write("walking through %s ..."%(inputFname))
		counter =0
		chr2xy_ls = self.chr2xy_ls
		inf = utils.openGzipFile(inputFname)
		delimiter=figureOutDelimiter(inf)	#2012.8.13 bugfix. pass inf to figureOutDelimiter
		sys.stderr.write(" delimiter is '%s'  "%(delimiter))
		reader = csv.reader(inf, delimiter=delimiter)
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		chr_id_index = col_name2index.get(chrColumnLabel, None)
		if chr_id_index is None:
			chr_id_index = col_name2index.get("CHROM", None)
		if chr_id_index is None:
			chr_id_index = col_name2index.get("CHR", None)
		if chr_id_index is None:
			sys.stderr.write("Error chr_id_index is None.\n")
			sys.exit(3)
		bin_start_index = col_name2index.get(posColumnLabel, None)
		chrLength_index = col_name2index.get(chrLengthColumnLabel, None)
		if self.whichColumnLabel:
			whichColumn = col_name2index.get(self.whichColumnLabel, None)
		else:
			whichColumn = self.whichColumn
		
		for row in reader:
			if self.samplingRate<1 and self.samplingRate>=0:
				r = random.random()
				if r>self.samplingRate:
					continue
			if chrLength_index:
				chrLength = int(row[chrLength_index])
				if chrLength<minChrLength:
					continue
			chr_id = row[chr_id_index]
			bin_start = int(float(row[bin_start_index]))
			
			yValue = float(row[whichColumn])
			if self.logWhichColumn:
				if yValue>0:
					yValue = -math.log10(yValue)
				else:
					yValue = 50
			if chr_id not in chr2xy_ls:
				chr2xy_ls[chr_id] = [[],[]]
			chr_cumu_start = self.chr_id2cumu_start.get(chr_id)
			
			chr2xy_ls[chr_id][0].append(chr_cumu_start + bin_start + 1)
			chr2xy_ls[chr_id][1].append(yValue)
			counter += 1
		del reader
		inf.close()
		sys.stderr.write("%s data.\n"%(counter))
		
	
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
				self.vcftoolsOutputStatFileWalker(inputFname, processFunc=None, chrColumnLabel=self.chrColumnLabel,\
									minChrLength=self.minChrLength, chrLengthColumnLabel=self.chrLengthColumnLabel,\
									posColumnLabel=self.posColumnLabel)
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
		
		ax.set_xlabel(self.posColumnPlotLabel)
		ax.set_ylabel(self.whichColumnPlotLabel)
		#ax.set_xlim([0, chr_id2cumu_size[chr_ls[-1]]])
		if self.ylim_type==1:
			ylim = ax.get_ylim()
			ax.set_ylim([0, ylim[1]])
		elif self.ylim_type==2:
			if max_y is not None and min_y is not None:
				delta = abs(max_y-min_y)/12.0
				ax.set_ylim([min_y-delta, max_y+delta])
		
		#outputFnamePrefix = os.path.splitext(self.outputFname)[0]
		outputFnamePrefix = self.outputFnamePrefix
		pylab.savefig('%s.png'%outputFnamePrefix, dpi=self.figureDPI)
		if self.need_svg:
			pylab.savefig('%s.svg'%outputFnamePrefix, dpi=self.figureDPI)
		outputFname = '%s_hist.png'%(outputFnamePrefix)
		yh_matplotlib.drawHist(value_ls, title='', \
				xlabel_1D=self.whichColumnPlotLabel, xticks=None, outputFname=outputFname, min_no_of_data_points=self.minNoOfTotal, \
				needLog=self.logCount, \
				dpi=self.figureDPI, min_no_of_bins=40)

if __name__ == '__main__':
	main_class = PlotVCFtoolsStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()