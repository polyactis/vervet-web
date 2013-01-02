#!/usr/bin/env python
"""
Examples:
	%s -i Contig149.vcf.gz -o 5NevisContig149_AAF_vs_position.png
	
	%s -i Contig83.vcf.gz -o 5NevisContig83_AAF_vs_position.png
	

Description:
	2011-10-17
		a program extracts the AF or AF1 (alternative allele frequency estimated by GATK/samtools) from vcf files
			and plot it along the contig position.
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
from pymodule import SNPData
import numpy, csv, os, sys, re, gzip
from pymodule.utils import runLocalCommand, getColName2IndexFromHeader
from pymodule.VCFFile import VCFFile

class PlotVCFAAF_vs_Position(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): [None, 'i', 1, 'input VCF file'],\
						('outputFname', 1, ): [None, 'o', 1, 'output file for the figure.'],\
						('minNoOfTotal', 1, int): [100, '', 1, 'minimum no of total variants (denominator of inconsistent rate)'],\
						('figureDPI', 1, int): [200, '', 1, 'dpi for the output figures (png)'],\
						('title', 0, ): [None, 't', 1, 'title for the figure.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def getLocusAndData(self, inputFname, VCFOutputType=2):
		"""
		2011-9-21
		
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
			locus_id = vcfRecord.locus_id
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			
			AF1 = vcfRecord.info_tag2value.get("AF", vcfRecord.info_tag2value.get("AF1", None))
			
			if AF1:
				AF1 = float(AF1)
				locus_ls.append(locus_id)
				xData_ls.append(pos)
				yData_ls.append(AF1)
		
		sys.stderr.write("%s loci. Done.\n"%(len(yData_ls)))
		return PassingData(contig_id=contig_id, locus_ls=locus_ls, yData_ls=yData_ls, xData_ls=xData_ls)
	
	@classmethod
	def isLocusPolymorphicBasedOnAAF(cls, frequency, frequencyDelta=0.001):
		"""
		2011-9-28
		"""
		distanceTo1 = abs(frequency-1.0)
		distanceTo0 = abs(frequency-0.0)
		if distanceTo0<frequencyDelta or distanceTo1<frequencyDelta:
			return False
		else:
			return True
	
	def plotVCFAlternativeAlleleFrequencyOverPosition(self, inputFname, outputFname, min_no_of_data_points=10, \
													xlabel='contig position', ylabel='AAF', dpi=300, frequencyDelta=0.001):
		"""
		2011-10-17
			purpose is to check alternative allele frequency along the contig position
		
			argument frequencyDelta is used to judge whether one frequency is close to 0 or 1, which essentially means
				these loci are not polymorphic.
		"""
		
		passingData = self.getLocusAndData(inputFname, VCFOutputType=2)
		
		if passingData is None:
			return
		
		locus_ls = passingData.locus_ls
		xData_ls = passingData.xData_ls
		yData_ls = passingData.yData_ls
		
		import os, sys
		from pymodule import yh_matplotlib
		
		
		import pylab
		pylab.clf()
		no_of_data_points = len(xData_ls)
		if no_of_data_points>=min_no_of_data_points:
			pylab.plot(xData_ls, yData_ls, ".",)
			pylab.title('%s sites in Contig%s'%(no_of_data_points, passingData.contig_id))
			if xlabel:
				pylab.xlabel(xlabel)
			if ylabel:
				pylab.ylabel(ylabel)
			pylab.savefig(outputFname, dpi=dpi)
		
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		self.plotVCFAlternativeAlleleFrequencyOverPosition(self.inputFname, self.outputFname, min_no_of_data_points=10, \
													xlabel='contig position', ylabel='AAF', dpi=self.figureDPI, frequencyDelta=0.001)
		

if __name__ == '__main__':
	main_class = PlotVCFAAF_vs_Position
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
