#!/usr/bin/env python
"""
Examples:
	#testing merge three identical genotype files
	%s -o /tmp/ccc.tsv /tmp/call_1.tsv /tmp/call_1.tsv /tmp/call_1.tsv
	
	%s 
	
Description:
	2011-7-12
		this program doesn't check whether multiple input genotype files contain same SNPs or not.
		
		The output format is VariantDiscovery.discoverHetsFromBAM() from vervet/src/misc.py.
			locus_id        locus_id        Barbados_GA_vs_top156Contigs    sabaeus_GA_vs_top156Contigs     VRC_ref_GA_vs_top156Contigs
			0_37    0_37    GG      CC      CC
			0_279   0_279   AA      AA      AC
			0_327   0_327   GG      GG      AG

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions

class MergeGenotypeMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, inputFnameLs, **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.inputFnameLs = inputFnameLs
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		header = None
		outf = open(self.outputFname, 'w')
		for inputFname in self.inputFnameLs:
			inf = open(inputFname)
			if header is None:
				header = inf.readline()
				outf.write(header)
			else:	#skip the header for other input files
				inf.readline()
			for line in inf:
				outf.write(line)
		

if __name__ == '__main__':
	main_class = MergeGenotypeMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	import copy
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
