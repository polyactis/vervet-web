#!/usr/bin/env python
"""
Examples:
	%s -a 454_vs_Contig0 -e Contig0,Contig1 -i /home/vervetData/subspecies/Barbados/vs_top150Contigs_by_aln.bam -o /tmp/
	
	%s -a Barbados_vs_Contig0 -p ILLUMINA -e Contig0,Contig1 -i /home/vervetData/subspecies/Barbados/vs_top150Contigs_by_aln.bam -o /tmp/
	
Description:
	2011-7-11
		select and split an input bam file into several small ones.
		each contains aligned reads based on one chosen reference name. It also attaches read group to each aligned read.
		
		It doesn't scale very well. Quickly running out of memory for large bam files.
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
import VervetDB
from pymodule.BamFile import BamFile
from pymodule import ProcessOptions

class SelectAndSplitAlignment(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('refNameLs', 1, ): ['', 'e', 1, 'a comma separated list of reference names', ],\
						("readGroup", 1, ): ["", 'a', 1, 'read-group to be added to each small bam file'],\
						("platform", 1, ): ['ILLUMINA', 'p', 1, 'or LS454, depends on the reads in the bam file. used in the readGroup info.'],\
						('outputDir', 1, ): [None, 'o', 1, 'directory to contain small bam files.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.refNameLs = self.refNameLs.split(',')
	
	class FilterAlignmentByReferenceIDs(object):
		"""
		2011-7-8
			This class filters bwa output bam file based on two criteria.
				read.mapq>=minMapQ
				score >= minPerBaseAS
			only for bwasw, which has AS (alignment score) tag.
			
		"""
		def __init__(self, **keywords):
			"""
			keywords shall include bamOutputF, scoreType, minMapQ, minPerBaseAS
			"""
			self.real_counter = 0
			
			# 2010-7-29
			for keyword, value in keywords.iteritems():
				setattr(self, keyword, value)
		
		def run(self, read, param_obj=None):
			"""
			2011-7-8
			"""
			param_obj.counter += 1
			if read.is_unmapped:
				return
			if read.qname not in param_obj.qname2count:
				param_obj.qname2count[read.qname] = 0
			param_obj.qname2count[read.qname] += 1
			score = None
			refName = param_obj.samfile.getrname(read.tid)
			if refName in self.refNameSet:
				param_obj.real_counter += 1
				if self.readGroup:	#2011-7-11 add a read group
					read.tags = read.tags + [("RG", self.readGroup)]
				bamOutputF = self.refName2bamOutputF.get(refName)	#find the corresponding bam output file Handler
				bamOutputF.write(read)
				return False	#like None, break the bam file reading
			else:
				return False	#if return True, means breaking once the target is no longer in refNameSet.
			
			
	@classmethod
	def filterAlignmentByReferenceIDs(cls, inputFname, outputDir=None, refNameSet = set(['Contig0', 'Contig1', 'Contig2']), \
									readGroup="", platform='LS454'):
		"""
		2011-7-8
		
		"""
		import os, sys
		import pysam, copy
		samfile = BamFile(inputFname, 'rb')
		header = copy.deepcopy(samfile.header)
		if readGroup:
			if "RG" not in header:
				header['RG'] = []
			header['RG'].append({'ID':readGroup, 'PL':platform, 'LB':platform, 'SM':readGroup})
		refName2bamOutputF = {}
		for refName in refNameSet:
			inputFileBaseNamePrefix = os.path.splitext(os.path.basename(inputFname))[0]
			outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix,refName))
			bamOutputF = pysam.Samfile(outputFname, 'wb', header=header)	# template=samfile)
			refName2bamOutputF[refName] = bamOutputF
		sys.stderr.write("Retain reads from %s only from these references %s ...\n"%(inputFname, refNameSet))
		processor = cls.FilterAlignmentByReferenceIDs(refName2bamOutputF=refName2bamOutputF, refNameSet=refNameSet, readGroup=readGroup)
		samfile.traverseBamByRead(processor=processor)
		for refName, bamOutputF in refName2bamOutputF.iteritems():
			bamOutputF.close()
	
	
	
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
		
		if not os.path.isdir(self.outputDir):
			os.makedirs(self.outputDir)
		self.filterAlignmentByReferenceIDs(self.inputFname, outputDir=self.outputDir, refNameSet = set(self.refNameLs), \
									readGroup=self.readGroup, platform=self.platform)

if __name__ == '__main__':
	main_class = SelectAndSplitAlignment
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
