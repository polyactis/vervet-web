#!/usr/bin/env python
"""
Examples:
	%s -i
	
	%s -i smartpca.evec -o smartpca_withMetaInfo.tsv
		

Description:
	2012.9.5
		a program that adds country, site, latitute, longitude, ucla_id, species-name to smartpca .evec output.
		Its output looks like this:
		
#eigvals:       36.853  19.837  6.413   3.946 
1034_743_VWP00312_GA_vs_524     0.1642  0.0318  -0.0340 -0.0118 Case
...

The extra column(s) (beyond what the header indicates) in the data matrix portion will be discarded.
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib; matplotlib.use("Agg")	#to disable pop-up requirement
import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, getColName2IndexFromHeader, figureOutDelimiter, SNPData
from pymodule.utils import getColName2IndexFromHeader, getListOutOfStr, figureOutDelimiter
from pymodule import yh_matplotlib, GenomeDB, utils
from pymodule import MatrixFile
from pymodule import SNP
import numpy, random, pylab
import numpy as np
from vervet.src import VervetDB
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
import networkx as nx
import matplotlib as mpl
from pymodule import TaxonomyDB


class AppendInfo2SmartPCAOutput(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict
	option_default_dict.update({
						('inversePCValue', 0, int):[0, '', 0, 'toggle to inverse PC values (match geography better usually)'],\
						
					})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	
	def connectDB(self):
		"""
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		AbstractVervetMapper.connectDB(self)
		
		db_taxonomy = TaxonomyDB.TaxonomyDB(drivername=self.drivername, db_user=self.db_user, db_passwd=self.db_passwd, \
									hostname=self.hostname, dbname=self.dbname, schema="taxonomy", port=self.port)
		db_taxonomy.setup(create_tables=False)
		self.db_taxonomy = db_taxonomy
	
	
	def appendInfo(self, inputFname=None, db_vervet=None, outputFname=None,\
				inversePCValue=True):
		"""
		#2012.9.25 skip samples whose individual_alignment entry could not be parsed.
		2012.9.5
		"""
		sys.stderr.write("Appending info to %s ..."%(inputFname))
		reader = MatrixFile(inputFname)
		header = reader.next()
		newHeader = ['individualID']
		for i in xrange(1, len(header)):
			newHeader.append('PC%s'%(i))
		newHeader.extend(['sex|string', 'country|string', 'site-id', 'site-name|string', 'latitude', 'longitude', 'ucla_id|string', \
						'tax_id|string',\
						'species|string', 'collectionYear', 'medianDepth'])
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(newHeader)
		counter = 0
		for row in reader:
			row = row[:len(header)]	#don't take extra columns
			sampleID = row[0]
			individualAlignment = db_vervet.parseAlignmentReadGroup(sampleID).individualAlignment
			if individualAlignment is None:
				#2012.9.25
				#sampleID is not beginned with alignment ID, probably "ref" but could be something , skip them
				sys.stderr.write("Warning: sampleID %s is not parsable to get alignment out of it. Skip.\n"%(sampleID))
				continue
			individual = individualAlignment.individual_sequence.individual
			data_row = ['%s_%s'%(individual.code, individualAlignment.id)]
			
			floatValue_row = row[1:]
			if inversePCValue:
				floatValue_row = map(float, floatValue_row)
				floatValue_row = numpy.array(floatValue_row)
				floatValue_row = -floatValue_row
			data_row.extend(list(floatValue_row))
			scientifcName = self.db_taxonomy.returnScientificNameGivenTaxID(individual.tax_id)
			if scientifcName is None:
				scientifcName = ""
			if individual.collection_date:
				collectionYear = individual.collection_date.year
			else:
				collectionYear = ''
			data_row.extend([individual.sex, individual.site.country.name, individual.site.id, individual.site.short_name, \
							individual.latitude, individual.longitude, individual.ucla_id, \
							individual.tax_id, scientifcName, collectionYear, individualAlignment.median_depth])
			writer.writerow(data_row)
			counter += 1
		del writer
		sys.stderr.write("%s rows outputted.\n"%(counter))
		
	
	def run(self):
		"""
		2012.9.5
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.appendInfo(inputFname=self.inputFname, db_vervet=self.db_vervet, outputFname=self.outputFname, \
				inversePCValue=self.inversePCValue)
		

if __name__ == '__main__':
	main_class = AppendInfo2SmartPCAOutput
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()