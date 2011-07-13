#!/usr/bin/env python
"""
Examples:
	# call one genome (nothing because no polymorphic sites)
	%s -i /tmp/vs_top150Contigs_by_aln_Contig0.bam -n 1 -o /tmp/vs_top150Contigs_by_aln_Contig0.call
	
	# call genotype from 8 genomes
	%s -n 8 -i script/vervet/data/8_genome_vs_Contig0.RG.bam -o script/vervet/data/8_genome_vs_Contig0.RG.call
	
Description:
	2011-7-12
		A multi-sample genotype caller based entirely on coverage of reads.
		It invokes VariantDiscovery.discoverHetsFromBAM() from vervet/src/misc.py
		sam/bam file has to be indexed beforehand.
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
from pymodule import ProcessOptions


class GenotypeCallByCoverage(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('numberOfReadGroups', 1, int): [None, 'n', 1, 'number of read groups/genomes in the inputFname', ],\
						('minMinorAlleleCoverage', 1, int): [3, '', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, int): [7, '', 1, 'maximum read depth for the minor allele of a heterozygous call', ],\
						('maxNoOfReadsForGenotypingError', 1, int): [1, '', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, int): [20, '', 1, 'maximum read depth for one base to be considered'],\
						('maxMajorAlleleCoverage', 1, int): [10, '', 1, 'maximum read depth'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, int): [3, '', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-12
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		outputDir = os.path.split(self.outputFname)[0]
		if outputDir and not os.path.isdir(outputDir):
			os.makedirs(outputDir)
		
		from vervet.src.misc import VariantDiscovery
		maxNoOfReadsForAllSamples = self.numberOfReadGroups*self.maxNoOfReads*self.maxNoOfReadsMultiSampleMultiplier
		VariantDiscovery.discoverHetsFromBAM(self.inputFname, self.outputFname, \
						maxNoOfReads=self.maxNoOfReads, minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
						maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
						maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
						maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
						maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples)

if __name__ == '__main__':
	main_class = GenotypeCallByCoverage
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()