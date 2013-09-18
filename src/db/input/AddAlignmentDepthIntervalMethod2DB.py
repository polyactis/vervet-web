#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s  --db_user yh --schema public --methodShortName PopSabaeusCoverageOn3488
		--alignmentIDList 6115-6177 --data_dir /u/home/p/polyacti/NetworkData/vervet/db/
		--commit --logFilename  folderLog/AddAlignmentDepthIntervalMethod2DB.log

Description:
	2013.08.29 Add a new AlignmentDepthIntervalMethod into db
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, utils
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
parentClass = AbstractVervetMapper

class AddAlignmentDepthIntervalMethod2DB(parentClass):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('min_segment_length', 0, int): [None, '', 1, 'a parameter of segmentation algorithm used in segmenting the depth file', ],\
						('alignmentIDList', 1, ): ['', '', 1, 'coma/dash-separated list of alignment IDs.', ],\
						('methodShortName', 1, ):[None, 's', 1, 'column short_name of AlignmentDepthIntervalMethod table,\
			will be created if not present in db.'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		self.alignmentIDList = utils.getListOutOfStr(list_in_str=self.alignmentIDList, data_type=int)
		
	
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		
		session.begin()
		if not self.data_dir:
			self.data_dir = self.db_vervet.data_dir
		
		
		alignmentList = self.db_vervet.getAlignmentsFromAlignmentIDList(self.alignmentIDList)
		
		method = self.db_vervet.getAlignmentDepthIntervalMethod(short_name=self.methodShortName, description=None, ref_ind_seq_id=None, \
						individualAlignmentLs=alignmentList, parent_db_entry=None, parent_id=None, \
						no_of_alignments=None, no_of_intervals=None, \
						sum_median_depth=None, sum_mean_depth=None,\
						min_segment_length=self.min_segment_length, \
						data_dir=self.data_dir)
		self.checkIfAlignmentListMatchMethodDBEntry(individualAlignmentLs=alignmentList, methodDBEntry=method, session=session)
		
		self.outputLogMessage(logMessage="AlignmentDepthIntervalMethod %s (%s) has been added into db.\n"%\
							(method.id, self.methodShortName))
		if self.commit:
			try:
				session.flush()
				session.commit()
			except:
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
				import traceback
				traceback.print_exc()
				self.cleanUpAndExitOnFailure(exitCode=3)
		else:
			session.rollback()
			#delete all target files
			self.cleanUpAndExitOnFailure(exitCode=0)

if __name__ == '__main__':
	main_class = AddAlignmentDepthIntervalMethod2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()