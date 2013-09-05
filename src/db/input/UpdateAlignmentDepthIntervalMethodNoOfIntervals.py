#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice --methodShortName 56PygerythrusCoverageOn3488 --data_dir ~/NetworkData/vervet/db/
		--commit --logFilename  folderLog/updateMethodNoOfIntervals.log

Description:
	2013.08
		Update the number of intervals in one AlignmentDepthIntervalMethod by summing the no_of_intervals of its associated AlignmentDepthIntervalFile entries.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB

class UpdateAlignmentDepthIntervalMethodNoOfIntervals(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('methodID', 0, int): [None, 'i', 1, 'AlignmentDepthIntervalMethod.id, used to fetch db entry, non-zero exit if not present in db', ],\
						('methodShortName', 0, ):[None, 's', 1, 'column short_name of AlignmentDepthIntervalMethod table,\
			non-zero exit if not present in db.'],\
						})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
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
		
		if self.methodID:
			method_db_entry = VervetDB.AlignmentDepthIntervalMethod.get(self.methodID)
		elif self.methodShortName:
			method_db_entry = VervetDB.AlignmentDepthIntervalMethod.query.filter_by(short_name=self.methodShortName).first()
		else:
			sys.stderr.write("ERROR: Both methodID (%s) and methodShortName (%s) are null.\n"%\
							(self.methodID, self.methodShortName))
			sys.exit(2)
		
		if method_db_entry is None:
			sys.stderr.write("ERROR: method_db_entry with methodID (%s) or methodShortName (%s) doesn't exist in db.\n"%\
							(self.methodID, self.methodShortName))
			sys.exit(4)
		logMessage = "method_db_entry %s (%s) has %s intervals.\n"%(method_db_entry.id, self.methodShortName, method_db_entry.no_of_intervals)
		self.db_vervet.updateAlignmentDepthIntervalMethodNoOfIntervals(db_entry=method_db_entry)
		
		logMessage += "method_db_entry %s (%s) updated with %s intervals.\n"%\
							(method_db_entry.id, self.methodShortName, method_db_entry.no_of_intervals)
		self.outputLogMessage(logMessage=logMessage)
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
	main_class = UpdateAlignmentDepthIntervalMethodNoOfIntervals
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
