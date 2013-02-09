#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s 

Description:
	2012.3.13
		abstract mapper for vervet mappers/reducers.

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from pymodule.AbstractDBInteractingClass import AbstractDBInteractingClass
from vervet.src import VervetDB

class AbstractVervetMapper(AbstractDBInteractingClass):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	
	option_default_dict.update({
							('logFilename', 0, ): [None, '', 1, 'file to contain logs. use it only if this program is at the end of pegasus workflow \
		and has no output file'],\
							})
	option_default_dict.update(AbstractMapper.db_option_dict)
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractDBInteractingClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
	
	def connectDB(self):
		"""
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, db_user=self.db_user, db_passwd=self.db_passwd, \
									hostname=self.hostname, dbname=self.dbname, schema=self.schema, port=self.port)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
	
