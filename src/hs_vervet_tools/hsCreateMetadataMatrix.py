#!/usr/bin/env python
"""
Examples:
	%s -m 4 -o /tmp/Method4_Contig917.tsv -c Contig917 -u yh
	
	# add -b to enable debug session.
	%s -m 4 -o /tmp/Method4_Contig917.tsv -c Contig917 -b

Description:
	2012.8.27
		example to get VCF file, parse it, and output it.
"""

import sys, os, math, types
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB

class ExampleToFetchVCFFromDB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
							('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
							('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
							("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. \
									If not given, use the default stored in db.'],\
							('genotypeMethodID', 1, int):[4, 'm', 1, 'genotype method ID'],\
							("chromosome", 1, ): [None, 'c', 1, 'chromosome for the genotype_file to be fetched'],\
							("format", 1, ): ['VCF', 'f', 1, 'format of genotype_file entry'],\
							})
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		
	
	def connectDB(self):
		"""
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname, schema=self.schema, port=self.port)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
	
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		
		session.begin()
		if not self.dataDir:
			self.dataDir = self.db_vervet.data_dir
		dataDir = self.dataDir
		
		genotypeFile = self.db_vervet.getGenotypeFile(genotype_method_id=self.genotypeMethodID, chromosome=self.chromosome, format=self.format)
		
		if not genotypeFile:
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
			sys.exit(2)
		filename = os.path.join(dataDir, genotypeFile.path)
		if os.path.isfile(filename):
			counter= 0
			from pymodule.VCFFile import VCFFile
			
			#allow 0 depth-> no missing data
			vcfFile = VCFFile(inputFname=filename,minDepth=0)
			sampleIDList = vcfFile.getSampleIDList()
			writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			sampleIDlist = ['sampleID']
			columnIndexList = []
			countryid_row=['country_id']
			uclaIDList=['ucla_id']
			speciesid_row=['tax_id']
			longitudeList=['longitude'];
			latitudeList=['latitude'];
			for i in xrange(len(sampleIDList)):
				sampleID = sampleIDList[i]
				individualAlignment = self.db_vervet.parseAlignmentReadGroup(sampleID).individualAlignment
				site = individualAlignment.ind_sequence.individual.site				
				sampleIDlist.append(sampleID)
				columnIndexList.append(i)
				uclaIDList.append(individualAlignment.ind_sequence.individual.ucla_id);
				countryid_row.append(individualAlignment.ind_sequence.individual.site.country_id)
				speciesid_row.append(individualAlignment.ind_sequence.individual.tax_id)
				longitudeList.append(individualAlignment.ind_sequence.individual.longitude);
				latitudeList.append(individualAlignment.ind_sequence.individual.latitude);
			writer.writerow(sampleIDlist)
			writer.writerow(uclaIDList)
			writer.writerow(speciesid_row)
			writer.writerow(countryid_row)
			writer.writerow(longitudeList)
			writer.writerow(latitudeList)
			del writer
		


if __name__ == '__main__':
	main_class = ExampleToFetchVCFFromDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()