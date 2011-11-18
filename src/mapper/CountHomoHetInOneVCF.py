#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i gatk/Contig799.vcf.gz -j samtools/Contig799.vcf.gz -l 1000000 -c Contig799 -o /tmp/output

Description:
	2011-11-7
		count the number of homo-ref/homo-alt/het calls from one vcf
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

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.VCFFile import VCFFile


class CountHomoHetInOneVCF(object):
	__doc__ = __doc__
	option_default_dict = {('inputFname', 1, ): ['', 'i', 1, 'VCF input file. either plain vcf or gzipped is ok. could be unsorted.', ],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("chromosome", 1, ): [None, 'c', 1, 'chromosome name for these two VCF.'],\
						("chrLength", 0, int): [0, 'l', 1, 'length of the reference used for the input VCF file.'],\
						('outputFname', 1, ): [None, 'o', 1, 'output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	
	def countHomoHetCallsForEachSampleFromVCF(self, inputFname, outputFname, chromosome=None, chrLength=None):
		"""
		2011-11-2
			given a VCF file, count the number of homo-ref, homo-alt, het calls
			
		"""
		sys.stderr.write("Count the number of homozygous-ref/alt & het from %s .\n"%(inputFname))
		from pymodule.VCFFile import VCFFile
		vcfFile = VCFFile(inputFname=inputFname)
		
		sampleID2data = {}	#key is sampleID, value is a list of 3 numbers. 'NoOfHomoRef', 'NoOfHomoAlt', 'NoOfHet'
		
		no_of_total = 0.
		minStart = None
		for vcfRecord in vcfFile.parseIter():
			locus = vcfRecord.locus
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			refBase = vcfRecord.data_row[0].get("GT")[0]
			
			for sample_id, sample_index in vcfFile.sample_id2index.iteritems():
				if sample_id=='ref':	#ignore the reference
					continue
				if sample_id not in sampleID2data:
					sampleID2data[sample_id] = [0, 0, 0]
				if not vcfRecord.data_row[sample_index]:	#None for this sample
					continue
				callForThisSample = vcfRecord.data_row[sample_index].get('GT')
				if not callForThisSample or callForThisSample=='NA':
					continue
				if callForThisSample[0]==refBase and callForThisSample[1]==refBase:
					#homozygous reference allele
					sampleID2data[sample_id][0]+=1
				elif callForThisSample[0]==callForThisSample[1] and callForThisSample[0]!=refBase:
					#homozygous alternative allele
					sampleID2data[sample_id][1]+=1
				elif callForThisSample[0]!=callForThisSample[1]:
					sampleID2data[sample_id][2]+=1
			
		import csv
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['sampleID', 'chromosome', 'length', "NoOfTotal", 'NoOfHomoRef', 'NoOfHomoAlt', "FractionOfHomoAlt", 'NoOfHet', "FractionOfHet"])
		sampleIDLs = sampleID2data.keys()
		sampleIDLs.sort()
		for sampleID in sampleIDLs:
			count_data = sampleID2data.get(sampleID)
			noOfHomoRef, noOfHomoAlt, noOfHet = count_data[:3]
			no_of_calls = float(sum(count_data))
			if no_of_calls>0:
				fractionOfHomoAlt = noOfHomoAlt/no_of_calls
				fractionOfHet = noOfHet/no_of_calls
			else:
				fractionOfHomoAlt = -1
				fractionOfHet = -1
			writer.writerow([sampleID, self.chromosome, self.chrLength, int(no_of_calls), noOfHomoRef, noOfHomoAlt, \
							fractionOfHomoAlt, noOfHet, fractionOfHet])
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		#outputFname = "%s.homoHetCountPerSample.tsv"%(outputFnamePrefix)
		self.countHomoHetCallsForEachSampleFromVCF(self.inputFname, self.outputFname)

if __name__ == '__main__':
	main_class = CountHomoHetInOneVCF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()