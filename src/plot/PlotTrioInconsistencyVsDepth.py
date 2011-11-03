#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011-11-2
		this program ignores loci whose depth is >4*Median for each member of the trio.
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
from pymodule import yh_matplotlib
from PlotTrioInconsistencyOverFrequency import PlotTrioInconsistencyOverFrequency
import numpy

class PlotTrioInconsistencyVsDepth(PlotTrioInconsistencyOverFrequency):
	__doc__ = __doc__
	option_default_dict = PlotTrioInconsistencyOverFrequency.option_default_dict
	#option_default_dict.update({('maxDepth', 1, int): [50, 'm', 1, 'genotype depth of every member of the trio should be below this number']})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		PlotTrioInconsistencyOverFrequency.__init__(self, inputFnameLs, **keywords)
		#Super(PlotTrioInconsistencyVsDepth, self).__init__(inputFnameLs, **keywords)
		self.fa_depth_ls = []
		self.mo_depth_ls = []
		self.child_depth_ls = []
		self.inconsistent_ls = []
		
	def trioInconsistentRateFileWalker(self, inputFname, processFunc=None, minNoOfTotal=100, run_type=1):
		"""
		2011-11-2
			remove the maxDepth filter. apply afterwards through filterDataByDepth().
		2011-9-30
		
		"""
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header, skipEmptyColumn=True)
		isInconsistent_index = col_name2index.get("isInconsistent")
		index_of_fa_depth = col_name2index.get("depthOfFather")
		index_of_mo_depth = col_name2index.get('depthOfMother')
		index_of_child_depth = col_name2index.get('depthOfChild')
		for row in reader:
			fa_depth = int(float(row[index_of_fa_depth]))
			mo_depth = int(float(row[index_of_mo_depth]))
			child_depth = int(float(row[index_of_child_depth]))
			isInconsistent = float(float(row[isInconsistent_index]))
			#if fa_depth<=self.maxDepth and mo_depth <=self.maxDepth and child_depth<=self.maxDepth:
			self.fa_depth_ls.append(fa_depth)
			self.mo_depth_ls.append(mo_depth)
			self.child_depth_ls.append(child_depth)
			self.inconsistent_ls.append(isInconsistent)
		del reader
	
	def filterDataByDepth(self, maxFaDepth=50, maxMoDepth=50, maxChildDepth=50):
		"""
		2011-11-2
			
		"""
		sys.stderr.write("Filter loci by father-depth<=%s, mother depth<=%s, child depth<=%s ..."%\
						(maxFaDepth, maxMoDepth, maxChildDepth))
		new_fa_depth_ls = []
		new_mo_depth_ls = []
		new_child_depth_ls = []
		new_inconsistent_ls = []
		no_of_total_loci = len(self.fa_depth_ls)
		for i in xrange(no_of_total_loci):
			fa_depth = self.fa_depth_ls[i]
			mo_depth = self.mo_depth_ls[i]
			child_depth = self.child_depth_ls[i]
			isInconsistent = self.inconsistent_ls[i]
			if fa_depth<=maxFaDepth and mo_depth <=maxMoDepth and child_depth<=maxChildDepth:
				new_fa_depth_ls.append(fa_depth)
				new_mo_depth_ls.append(mo_depth)
				new_child_depth_ls.append(child_depth)
				new_inconsistent_ls.append(isInconsistent)
		
		self.fa_depth_ls = new_fa_depth_ls
		self.mo_depth_ls = new_mo_depth_ls
		self.child_depth_ls = new_child_depth_ls
		self.inconsistent_ls = new_inconsistent_ls
		sys.stderr.write("%s/%s loci retained.\n"%(len(new_inconsistent_ls), no_of_total_loci))
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		import pylab
		pylab.clf()
		
		if len(self.inputFnameLs)>1:	#multiple input file (multiple contigs). use lines.
			self.formatString = '-'
		else:	#only one input file (one contig). use dot.
			self.formatString = '.'
		
		sys.stderr.write("Reading in data ...")
		for inputFname in self.inputFnameLs:
			self.trioInconsistentRateFileWalker(inputFname, processFunc=None)
		sys.stderr.write("Done.\n")
		
		medianFaDepth = numpy.median(self.fa_depth_ls)
		medianMoDepth = numpy.median(self.mo_depth_ls)
		medianChildDepth = numpy.median(self.child_depth_ls)
		self.filterDataByDepth(maxFaDepth=4*medianFaDepth, maxMoDepth=4*medianMoDepth, maxChildDepth=4*medianChildDepth)
		
		if self.title is None:
			title = " %s refs"%(len(self.inputFnameLs))
		else:
			title = self.title
		
		outputFnamePrefix = self.outputFname
		
		loci_count_ls = [1]*len(self.fa_depth_ls)
		colorBarLabelForInconsistency = 'inconsistent rate'
		colorBarLabelForLociCount = 'log10(no. of loci)'
		
		#if uniformly distributed over 2D, each hexagon has ~50 points.
		# the maximum gridsize is 30.
		gridsize = min(30, int(math.sqrt(len(loci_count_ls)/50.0)))
		
		fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.fa_depth_ls, self.mo_depth_ls, self.inconsistent_ls, \
								fig_fname=fa_mo_depth_Fname, gridsize=gridsize, title=title, xlabel='father depth', ylabel='mother depth',\
								colorBarLabel=colorBarLabelForInconsistency, reduce_C_function=numpy.mean, dpi=300)
		fa_mo_depth_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.fa_depth_ls, self.mo_depth_ls, loci_count_ls, \
								fig_fname=fa_mo_depth_Fname, gridsize=gridsize, title=title, xlabel='father depth', ylabel='mother depth',\
								colorBarLabel=colorBarLabelForLociCount, reduce_C_function=yh_matplotlib.logSum, dpi=300)
		
		fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.fa_depth_ls, self.child_depth_ls, self.inconsistent_ls, \
								fig_fname=fa_child_depth_Fname, gridsize=gridsize, title=title, xlabel='father depth', ylabel='child depth',\
								colorBarLabel=colorBarLabelForInconsistency, reduce_C_function=numpy.mean, dpi=300)
		fa_child_depth_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.fa_depth_ls, self.child_depth_ls, loci_count_ls, \
								fig_fname=fa_child_depth_Fname, gridsize=gridsize, title=title, xlabel='father depth', ylabel='child depth',\
								colorBarLabel=colorBarLabelForLociCount, reduce_C_function=yh_matplotlib.logSum, dpi=300)
		
		mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.mo_depth_ls, self.child_depth_ls, self.inconsistent_ls, \
								fig_fname=mo_child_depth_Fname, gridsize=gridsize, title=title, xlabel='mother depth', ylabel='child depth',\
								colorBarLabel=colorBarLabelForInconsistency, reduce_C_function=numpy.mean, dpi=300)
		mo_child_depth_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
		yh_matplotlib.drawHexbin(self.mo_depth_ls, self.child_depth_ls, loci_count_ls, \
								fig_fname=mo_child_depth_Fname, gridsize=gridsize, title=title, xlabel='mother depth', ylabel='child depth',\
								colorBarLabel=colorBarLabelForLociCount, reduce_C_function=yh_matplotlib.logSum, dpi=300)
		
if __name__ == '__main__':
	main_class = PlotTrioInconsistencyVsDepth
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
