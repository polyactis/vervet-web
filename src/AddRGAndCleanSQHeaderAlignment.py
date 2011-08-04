#!/usr/bin/env python
"""
Examples:
	%s -a 454_vs_Contig0 -p LS454 -e Contig0,Contig1 -i /home/vervetData/subspecies/Barbados/vs_top150Contigs_by_aln.bam
		-o /tmp/Contig_RG.bam
	
	%s -a Barbados_vs_Contig0 -p ILLUMINA -e Contig0,Contig1 -i /home/vervetData/subspecies/Barbados/vs_top150Contigs_by_aln.bam
		-o /tmp/contigs_by_aln_RG.bam
	
Description:
	2011-7-26
	This program does:
		1. retain reads that are only within refNameLs, from an input bam file.
		2. add read group information to the bam header
		3. remove "SQ"s from the bam header whose name is not in refNameLs
		4. add read group to each read.
		
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

class AddRGAndCleanSQHeaderAlignment(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('refNameLs', 1, ): ['', 'e', 1, 'a comma separated list of reference names. reads not aligned to these references will be removed', ],\
						("readGroup", 1, ): ["", 'a', 1, 'read-group to be added to each small bam file'],\
						("platform", 1, ): ['ILLUMINA', 'p', 1, 'or LS454, depends on the reads in the bam file. used in the readGroup info.'],\
						('outputFname', 1, ): [None, 'o', 1, 'file to contain the final output.'],\
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
				bamOutputF = self.bamOutputF
				bamOutputF.write(read)
				return False	#like None, break the bam file reading
			else:
				return False	#if return True, means breaking once the target is no longer in refNameSet.
			
			
	@classmethod
	def filterAlignmentByReferenceIDs(cls, inputFname, outputFname=None, refNameSet = set(['Contig0', 'Contig1', 'Contig2']), \
									readGroup="", platform='LS454'):
		"""
		2011-7-8
		
		"""
		import os, sys
		import pysam, copy
		samfile = BamFile(inputFname, 'rb')
		header = copy.deepcopy(samfile.header)
		if readGroup:	#add read group if it's missing
			if "RG" not in header:
				header['RG'] = []
			header['RG'].append({'ID':readGroup, 'PL':platform, 'LB':platform, 'SM':readGroup})
		
		# remove SQ entries that are not in refNameSet, from the header
		newSQList = []
		for SQ_entry in header["SQ"]:
			if SQ_entry['SN'] in refNameSet:
				newSQList.append(SQ_entry)
		header["SQ"] = newSQList
		
		refName2bamOutputF = {}
		bamOutputF = pysam.Samfile(outputFname, 'wb', header=header)	# template=samfile)
		sys.stderr.write("Retain reads from %s only from these references %s ...\n"%(inputFname, refNameSet))
		processor = cls.FilterAlignmentByReferenceIDs(refNameSet=refNameSet, readGroup=readGroup, bamOutputF=bamOutputF)
		samfile.traverseBamByRead(processor=processor)
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
		
		self.filterAlignmentByReferenceIDs(self.inputFname, outputFname=self.outputFname, refNameSet = set(self.refNameLs), \
									readGroup=self.readGroup, platform=self.platform)

if __name__ == '__main__':
	main_class = AddRGAndCleanSQHeaderAlignment
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
