#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c /tmp/outputStat.tsv /tmp/depth.mode.tsv

Description:
	2012.4.3
		Put output of gatk's DepthOfCoverage walker into db. part of InspectAlignmentPipeline.py.
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


class PutDOCOutput2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
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
		2012.5.7 new input looks like this (tab-delimited):
			alignmentID     total_base_count        sampled_base_count      meanDepth       medianDepth     modeDepth
			100     1005506 301614  70.0441756682   9.0     8.0
			27     1005506 301614  70.0441756682   9.0     8.0

		2012.4.3
			each input looks like this:
			
sample_id       total   mean    granular_third_quartile granular_median granular_first_quartile %_bases_above_15
553_2_VRC_ref_GA_vs_524 2434923137      8.25    11      9       6       4.4
Total   2434923137      8.25    N/A     N/A     N/A
554_3_Barbados_GA_vs_524        2136011136      7.23    11      8       6       3.5
Total   2136011136      7.23    N/A     N/A     N/A
...

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
			col_name2index = utils.getColName2IndexFromHeader(header, skipEmptyColumn=True)
			
			sample_id_index = col_name2index.get("alignmentID")
			total_base_count_index = col_name2index.get('total_base_count')
			mean_depth_index = col_name2index.get("meanDepth")
			median_depth_index = col_name2index.get("medianDepth")
			mode_depth_index = col_name2index.get("modeDepth")
			for row in reader:
				sample_id = row[sample_id_index]
				if sample_id=='Total':	#ignore rows with this as sample id
					continue
				alignment_id = int(sample_id.split("_")[0])
				total_base_count = int(row[total_base_count_index])
				mean_depth = float(row[mean_depth_index])
				median_depth = float(row[median_depth_index])
				mode_depth = float(row[mode_depth_index])
				individual_alignment = VervetDB.IndividualAlignment.get(alignment_id)
				individual_alignment.pass_qc_read_base_count = total_base_count	#2012.9.17 no longer trustworthy because CalculateMedianMeanOfInputColumn skips data.
				individual_alignment.mean_depth = mean_depth
				individual_alignment.median_depth = median_depth
				individual_alignment.mode_depth = mode_depth
				session.add(individual_alignment)
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
	main_class = PutDOCOutput2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()