#!/usr/bin/env python
"""
Examples:
	%s -r -E ./*
	
	%s -c -E -r ./*
	

Description:
	2011-12.19
		program that deletes all empty VCF files.
		If -E is on, "empty" includes ones with VCF header but no site.
		Default is a dry-run. "-c" would really delete files.
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
import csv, numpy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, NextGenSeq
from pymodule import yh_matplotlib
from pymodule.utils import runLocalCommand


class DeleteEmptyVCF(object):
	__doc__ = __doc__
	option_default_dict = {
				('checkEmptyVCFByReading', 0, int):[0, 'E', 0, 'toggle to check if a vcf file is empty by reading its content'],\
				('commit', 0, int):[0, 'c', 0, 'toggle to really delete the empty vcf files.'],\
				('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
				('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
				}


	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		counter = 0
		no_of_vcf = 0
		real_counter = 0
		for inputFname in self.inputFnameLs:
			counter += 1
			if os.path.isfile(inputFname):
				try:
					if NextGenSeq.isFileNameVCF(inputFname, includeIndelVCF=False):
						no_of_vcf += 1
						if NextGenSeq.isVCFFileEmpty(inputFname, checkContent=self.checkEmptyVCFByReading):
							if self.commit:
								commandline = 'rm %s'%(inputFname)
								return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
							real_counter += 1
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
			if self.report and counter%500==0:
				sys.stderr.write("%s%s\t%s\t%s"%('\x08'*80, counter, no_of_vcf, real_counter))
		sys.stderr.write("%s%s\t%s\t%s\n"%('\x08'*80, counter, no_of_vcf, real_counter))
		sys.stderr.write("%s files in total.\n"%(counter))
		sys.stderr.write("Out of %s VCF files, %s are empty and were deleted.\n"%(no_of_vcf, real_counter))
		
		
if __name__ == '__main__':
	main_class = DeleteEmptyVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
