#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice /tmp/outputStat.tsv

Description:
	2012.4.3
		Put output of "samtools flagstat" into db. part of InspectAlignmentPipeline.py.
		Input files are added after all the arguments on the commandline.
	
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class PutFlagstatOutput2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs, **keywords)
	
	
	def run(self):
		"""
		2012.4.3
			each input has this as its header:
			
			['alignmentID', 'total_no_of_reads', 'perc_reads_mapped', 'perc_duplicates', 'perc_paired', 'perc_properly_paired', \
				'perc_both_mates_mapped', 'perc_singletons',\
				'perc_mapped_to_diff_chrs']
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		session.begin()
		
		no_of_total_lines = 0
		for inputFname in self.inputFnameLs:
			reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
			header = reader.next()
			colName2Index = utils.getColName2IndexFromHeader(header, skipEmptyColumn=True)
			alignment_id_index = colName2Index.get('alignmentID')
			total_no_of_reads_index = colName2Index.get('total_no_of_reads')
			perc_reads_mapped_index = colName2Index.get("perc_reads_mapped")
			perc_duplicates_index = colName2Index.get("perc_duplicates")
			perc_paired_index = colName2Index.get("perc_paired")
			perc_properly_paired_index = colName2Index.get("perc_properly_paired")
			perc_both_mates_mapped_index = colName2Index.get("perc_both_mates_mapped")
			perc_singletons_index = colName2Index.get("perc_singletons")
			perc_mapped_to_diff_chrs_index = colName2Index.get("perc_mapped_to_diff_chrs")
			perc_mapq5_mapped_to_diff_chrs_index = colName2Index.get("perc_mapq5_mapped_to_diff_chrs")
			for row in reader:
				alignmentID = int(row[alignment_id_index])
				alignment = VervetDB.IndividualAlignment.get(alignmentID)
				alignment.perc_reads_mapped = float(row[perc_reads_mapped_index])
				alignment.perc_duplicates = float(row[perc_duplicates_index])
				alignment.perc_paired = float(row[perc_paired_index])
				alignment.perc_properly_paired = float(row[perc_properly_paired_index])
				alignment.perc_both_mates_mapped = float(row[perc_both_mates_mapped_index])
				alignment.perc_singletons = float(row[perc_singletons_index])
				alignment.perc_mapped_to_diff_chrs = float(row[perc_mapped_to_diff_chrs_index])
				alignment.perc_mapq5_mapped_to_diff_chrs = float(row[perc_mapq5_mapped_to_diff_chrs_index])
				alignment.total_no_of_reads = int(float(row[total_no_of_reads_index]))
				session.add(alignment)
				no_of_total_lines += 1
			del reader
		sys.stderr.write("%s alignments in total.\n"%(no_of_total_lines))
		
		if self.logFilename:
			logF = open(self.logFilename, 'w')
			logF.write("%s alignments in total.\n"%(no_of_total_lines))
			del logF
		
		if self.commit:
			self.db_vervet.session.flush()
			self.db_vervet.session.commit()

if __name__ == '__main__':
	main_class = PutFlagstatOutput2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()