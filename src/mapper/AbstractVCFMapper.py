#!/usr/bin/env python
"""

Description:
	2012.1.17 an abstract class for vcf mapper
"""

import sys, os, math
__doc__ = __doc__

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile
from pymodule import SNP

class AbstractVCFMapper(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("chromosome", 1, ): [None, 'c', 1, 'chromosome name for these two VCF.'],\
						("chrLength", 1, int): [0, 'l', 1, 'length of the reference used for the input VCF file.'],\
						('minDepth', 1, float): [1, 'm', 1, 'minimum depth for a call to regarded as non-missing', ],\
						('outputFname', 1, ): [None, 'o', 1, 'output file'],\
						('outputFnamePrefix', 0, ): [None, 'O', 1, 'output filename prefix (optional). *_overlapSitePos.tsv lists position of overlap sites.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	