#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i ~namtran/panasas/Experiment/RNA-seq/Freimer/Developmental/ASE/Variant/MultipleSampleCalling/genome.algn.split.part17/samtools.var.filt.vcf.gz 
		-o /tmp/ -m 2
	
	#run in a for loop. "ls -d" only lists the folders and doesn't list contents of those folders.
	for i in `ls -d ~namtran/panasas/Experiment/RNA-seq/Freimer/Developmental/ASE/Variant/MultipleSampleCalling/genome*`; 
		do echo $i; 
		%s -i  $i/samtools.var.filt.vcf.gz  -o VariantsOf36RNA-SeqMonkeysFromNam_minDepth5/ -m 5;
	done
	
Description:
	2012.5.10
		
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import VCFFile
from AbstractVCFMapper import AbstractVCFMapper

class SplitNamVCFIntoMultipleSingleChrVCF(AbstractVCFMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVCFMapper.option_default_dict.copy()
	option_default_dict.pop(("chromosome", 0, ))
	option_default_dict.pop(("chrLength", 1, int))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('includeIndels', 0, ): [0, 'n', 0, 'By default, indels are excluded. INFO tag contains INDEL.', ],\
						('outputDir', 1, ): [None, 'o', 1, 'The folder to contain single-chromosome VCF files.', ],\
						('maxContigNumber', 1, int): [1000, 'a', 1, 'The maximum contig number', ],\
						})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractVCFMapper.__init__(self, **keywords)
		import re
		self.chr_pattern = re.compile(r'(\w+\d+).*')
		self.UCLAID_Pattern = re.compile(r'\w+/[\w\.]+/(?P<UCLAID>\d+)')
		self.contig_id_pattern = re.compile(r'vervet1_scaffolds_(?P<contigID>Contig\d+)')
		self.contig_number_pattern = re.compile(r'Contig(?P<contigNumber>\d+)')

	def splitNamVCFIntoMultipleSingleChrVCF(self, inputFname, outputDir, minDepth=1, includeIndels=False, maxContigNumber=1000):
		"""
		2012.5.10
			Two things in Nam's VCF file are to be modified. 
				1. extract VRC UCLAID from its sample ID
				2. replace vervet1_scaffolds_Contig137 with simply "Contig137"
		"""
		sys.stderr.write("Converting %s from VCF to EigenStrat ...\n"%(inputFname))
		from pymodule.VCFFile import VCFFile
		
		vcfFile = VCFFile(inputFname=inputFname, minDepth=minDepth)
		#replace Variant/PooledTissues/2002053/genome.algn.split.part17/5tissues.pooled.rmdup.bam with just monkey ID
		import re
		
		newSampleIDHeader = []
		for sampleID in vcfFile.sampleIDHeader:
			search_result = self.UCLAID_Pattern.search(sampleID)
			UCLAID = search_result.group('UCLAID')
			newSampleIDHeader.append(UCLAID)
		#new header for every output contig
		newHeader = vcfFile.header[:vcfFile.sampleStartingColumn] + newSampleIDHeader
		
		
		chr2outVCFFile = {}
		counter = 0
		real_counter = 0
		for vcfRecord in vcfFile.parseIter():
			counter += 1
			if not includeIndels and (len(vcfRecord.refBase)!=1 or len(vcfRecord.altBase)!=1):
				#it's an indel if refBase or altBase is not just one base
				continue
			
			contig_id_pattern_result = self.contig_id_pattern.search(vcfRecord.chr)
			chr = contig_id_pattern_result.group('contigID')
			if maxContigNumber:
				contigNumber = int(self.contig_number_pattern.search(chr).group('contigNumber'))
				if contigNumber>maxContigNumber:
					continue
			real_counter += 1
			vcfRecord.chr = chr
			pos = vcfRecord.pos
			if chr not in chr2outVCFFile:
				outputFname = os.path.join(outputDir, '%s.vcf'%(chr))
				outVCFFile = VCFFile(outputFname=outputFname)
				outVCFFile.metaInfoLs = vcfFile.metaInfoLs
				outVCFFile.header = newHeader
				outVCFFile.writeMetaAndHeader()
				chr2outVCFFile[chr] = outVCFFile
			outVCFFile = chr2outVCFFile.get(chr)
			
			# set genotype whose depth is below minDepth to ./. (=missing)
			for i in xrange(1, len(vcfRecord.data_row)):	#[0] is the ref base
				callData = vcfRecord.data_row[i]
				if callData is None or callData.get('DP',0)<minDepth:
					sampleColumnIndex = i+vcfFile.sampleStartingColumn-1
					vcfRecord.row[sampleColumnIndex] = './.'
			outVCFFile.writeVCFRecord(vcfRecord)
		
		vcfFile.close()
		#close all output files
		for chr, outVCFFile in chr2outVCFFile.iteritems():
			outVCFFile.close()
		
		sys.stderr.write("%s (out of %s) loci from %s chromosomes.\n"%(real_counter, counter, len(chr2outVCFFile)))
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.splitNamVCFIntoMultipleSingleChrVCF(self.inputFname, self.outputDir, minDepth=self.minDepth,\
												includeIndels=self.includeIndels, maxContigNumber=self.maxContigNumber)

if __name__ == '__main__':
	main_class = SplitNamVCFIntoMultipleSingleChrVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()