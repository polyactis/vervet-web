#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i /tmp/outputStat.txt.gz -a 1043 -o /tmp/outputStat.tsv.gz

Description:
	2012.4.3
		reformat output of samtools flagstat so that it could be inserted into db.
		both input and output files could be either gzipped or not.

	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.yhio.VCFFile import VCFFile
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper

class ReformatFlagstatOutput(AbstractMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('alignmentID', 1, ): [None, 'a', 1, 'ID of this alignment from which all the stats are extracted.'],\
						
							})
	def __init__(self,  **keywords):
		"""
		"""
		AbstractMapper.__init__(self, **keywords)
		import re;
		self.numberGrabPattern = re.compile(r'^(\d+) \+ (\d+)')
	
	def getNumberOutOfFlagStatLine(self, line):
		"""
		2012.4.3
			the output of samtools flagstat looks like:

2725130 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
2725130 + 0 mapped (100.00%:-nan%)
2725130 + 0 paired in sequencing
1360823 + 0 read1
1364307 + 0 read2
576948 + 0 properly paired (21.17%:-nan%)
609252 + 0 with itself and mate mapped
2115878 + 0 singletons (77.64%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
		"""
		searchResult = self.numberGrabPattern.search(line)
		if searchResult:
			n1 = int(searchResult.group(1))
			n2 = int(searchResult.group(2))
			return n1+n2
		else:
			sys.stderr.writer("Error: could not parse numbers out of this line (%s).\n"%(line))
			sys.exit(2)
			return None
	
	def run(self):
		"""
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		"""
		2012.4.3
		the output of samtools flagstat looks like:

2725130 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
2725130 + 0 mapped (100.00%:-nan%)
2725130 + 0 paired in sequencing
1360823 + 0 read1
1364307 + 0 read2
576948 + 0 properly paired (21.17%:-nan%)
609252 + 0 with itself and mate mapped
2115878 + 0 singletons (77.64%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

		"""
		
		inf = utils.openGzipFile(self.inputFname, openMode='r')
		writer = csv.writer(utils.openGzipFile(self.outputFname, openMode='w'), delimiter='\t')
		header = ['alignmentID', 'total_no_of_reads', 'perc_reads_mapped', 'perc_duplicates', 'perc_paired', 'perc_properly_paired', \
				'perc_both_mates_mapped', 'perc_singletons',\
				'perc_mapped_to_diff_chrs',
				'perc_mapq5_mapped_to_diff_chrs']
		writer.writerow(header)
		total_no_of_reads = float(self.getNumberOutOfFlagStatLine(inf.next()))	#float it now so that no "float" upon division
		no_of_duplicates = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_mapped = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_paired = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_read1 = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_read2 = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_properly_paired = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_both_mates_mapped = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_singletons = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_mates_mapped_to_diff_chrs = self.getNumberOutOfFlagStatLine(inf.next())
		no_of_mates_mapped_to_diff_chrs_mapQAbove5 = self.getNumberOutOfFlagStatLine(inf.next())
		#
		del inf
		
		data_row = [self.alignmentID, total_no_of_reads, no_of_mapped/total_no_of_reads*100, no_of_duplicates/total_no_of_reads*100,\
				no_of_paired/total_no_of_reads*100, no_of_properly_paired/total_no_of_reads*100,\
				no_of_both_mates_mapped/total_no_of_reads*100, no_of_singletons/total_no_of_reads*100,\
				no_of_mates_mapped_to_diff_chrs/total_no_of_reads*100,
				no_of_mates_mapped_to_diff_chrs_mapQAbove5/total_no_of_reads*100]
		writer.writerow(data_row)
		del writer
		

if __name__ == '__main__':
	main_class = ReformatFlagstatOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()