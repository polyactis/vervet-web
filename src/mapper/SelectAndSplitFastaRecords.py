#!/usr/bin/env python
"""
Examples:
	%s -i /home/vervetData/subspecies/Barbados/vs_top150Contigs_by_aln.bam -o /tmp/ Contig0 Contig1

Description:
	2011-7-26
		select and split an input fasta file into several small ones. The input arguments are fasta record titles.
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0])

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions

class SelectAndSplitFastaRecords(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input fasta file.', ],\
						('outputDir', 1, ): [None, 'o', 1, 'directory to contain single-record fasta files.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self, refNameLs, **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		self.refNameLs = refNameLs
	
	
	@classmethod
	def splitMultiRecordFastaFileIntoSingleRecordFastaFiles(cls, inputFname, outputDir, refNameSet=None,
														chunkSize=70):
		"""
		2011-10-18
			this function gets a speed-up.
				writing the sequences in fixed-chunk through a for-loop is the cause of slow-speed.
				stop using Bio.SeqIO and write the input line out immediately.
		2011-7-7
			use Bio.SeqIO and add argument fastaRecordFilterHandler
		2010-12-3
			One fasta file contains >1 sequences. This function writes each sequence to an individual file in outputDir.
			The sequence title would be the filename in outputDir.
			
			SyMAP requires this kind of splitting for it to work on fast files.
		"""
		import os,sys
		sys.stderr.write("Outputting each sequence in %s into single file ...\n"%(inputFname))
		inf = open(inputFname, 'rU')
		outf = None
		counter = 0
		real_counter = 0
		for line in inf:
			if line[0]=='>':
				if real_counter >= len(refNameSet):	#exit if all required refs have been selected.
					#this break has to be put here, after the last fasta block has been out.
					break
				counter += 1
				title = line.strip()[1:].split()[0]
				if title not in refNameSet:
					outf = None
				else:
					if outf is not None:
						outf.close()
						del outf
					output_fname = os.path.join(outputDir, '%s.fasta'%title)
					outf = open(output_fname, 'w')
					outf.write(">%s\n"%title)
					real_counter += 1
			else:
				if outf is not None:
					outf.write(line)
			
			if counter%1000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		del inf
		sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
	
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
		self.splitMultiRecordFastaFileIntoSingleRecordFastaFiles(self.inputFname, self.outputDir, \
												refNameSet=set(self.refNameLs), chunkSize=70)

if __name__ == '__main__':
	main_class = SelectAndSplitFastaRecords
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
