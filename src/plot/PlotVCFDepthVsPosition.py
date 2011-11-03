#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 
	

Description:
	2011-10-29
		this program plots depth of genotype call v.s. its position.
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
import csv, re
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter
from pymodule import yh_matplotlib
from PlotVCFAAF_vs_Position import PlotVCFAAF_vs_Position
from pymodule.VCFFile import VCFFile


class PlotVCFDepthVsPosition(PlotVCFAAF_vs_Position):
	__doc__ = __doc__
	option_default_dict = PlotVCFAAF_vs_Position.option_default_dict

	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		PlotVCFAAF_vs_Position.__init__(self, inputFnameLs, **keywords)
		#Super(PlotVCFDepthVsPosition, self).__init__(inputFnameLs, **keywords)
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
		"""
	
	def getLocusAndData(self, inputFname, VCFOutputType=2):
		"""
		
		"""
		contig_id_pattern = re.compile(r'Contig(\d+).*')
		contig2locus2frequency = {}
		fname = inputFname
		if fname[-6:]!='vcf.gz' and fname[-3:]!='vcf':
			return None
		sys.stderr.write("%s ..."%fname)
		contig_id_pattern_sr = contig_id_pattern.search(inputFname)
		if contig_id_pattern_sr:
			contig_id = contig_id_pattern_sr.group(1)
		else:
			contig_id = os.path.splitext(os.path.split(inputFname)[1])[0]
		
		vcfFile = VCFFile(inputFname=self.inputFname)
		counter = 0
		real_counter = 0
		
		locus_ls = []
		xData_ls = []
		yData_ls = []
		
		for vcfRecord in vcfFile.parseIter():
			locus = vcfRecord.locus
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			
			yData = vcfRecord.info_tag2value.get("DP", None)
			
			if yData:
				yData = float(yData)
				locus_ls.append(locus)
				xData_ls.append(pos)
				yData_ls.append(yData)
		
		sys.stderr.write("%s loci. Done.\n"%(len(yData_ls)))
		return PassingData(contig_id=contig_id, locus_ls=locus_ls, yData_ls=yData_ls, xData_ls=xData_ls)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		self.plotVCFAlternativeAlleleFrequencyOverPosition(self.inputFname, self.outputFname, min_no_of_data_points=10, \
													xlabel='contig position', ylabel='depth', dpi=300)

if __name__ == '__main__':
	main_class = PlotVCFDepthVsPosition
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
