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

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule import TaxonomyDB
from pymodule import AbstractMapper
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class ExampleToFetchVCFFromDB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.update({
							('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
							('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['vervetdb', 'd', 1, 'database name', ],\
							('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, '', 1, 'database port number'],\
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
		AbstractVervetMapper.connectDB(self)
		
		db_taxonomy = TaxonomyDB.TaxonomyDB(drivername=self.drivername, db_user=self.db_user, db_passwd=self.db_passwd, \
									hostname=self.hostname, dbname=self.dbname, schema="taxonomy", port=self.port)
		db_taxonomy.setup(create_tables=False)
		self.db_taxonomy = db_taxonomy
	
	
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
		data_dir = self.data_dir
		
		genotypeFile = self.db_vervet.getGenotypeFile(genotype_method_id=self.genotypeMethodID, chromosome=self.chromosome, format=self.format)
		#query = VervetDB.GenotypeFile.query.filter_by(genotype_method_id=self.genotypeMethodID).filter_by(format=self.format)
		#for genotypeFile in query:
		if not genotypeFile:
			
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
			sys.exit(2)
		filename = os.path.join(data_dir, genotypeFile.path)
		if os.path.isfile(filename):
			counter= 0
			from pymodule import VCFFile
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			sampleIDList = vcfFile.getSampleIDList()
			writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			header = ['Chromosome', 'position', 'ref']
			columnIndexList = []
			for i in xrange(len(sampleIDList)):
				sampleID = sampleIDList[i]
				individualAlignment = self.db_vervet.parseAlignmentReadGroup(sampleID).individualAlignment
				site = individualAlignment.individual_sequence.individual.site
				#2012.8.29 get scientific name from the taxonomy db
				scientifcName = self.db_taxonomy.returnScientificNameGivenTaxID(individualAlignment.individual_sequence.individual.tax_id)
				#if individualAlignment.individual_sequence.individual.tax_id==60711 and (site.country_id!=144 and site.country_id!=135 \
				#																and site.country_id!=136 and site.country_id!=148): 
				header.append('%s %s'%(sampleID, scientifcName))
				columnIndexList.append(i)
			writer.writerow(header)
			for vcfRecord in vcfFile:
				data_row = [vcfRecord.chr, vcfRecord.pos]
				refCall = vcfRecord.data_row[0]
				data_row.append(refCall['GT'])
				#get alternative allele frequency
				AF_list = vcfRecord.info_tag2value.get('AF')	#info_tag2value['AF']
				AF_list = AF_list.split(',')
				AF_list = map(float, AF_list)
				for columnIndex in columnIndexList:
					#for vcfCall in vcfRecord.data_row[1:]: #data_row is a list of dictionary {'GT': base-call, 'DP': depth}, or None if missing.
					#it includes one extra sample in index 0, which is the reference individual (from the ref column of VCF).
					vcfCall = vcfRecord.data_row[columnIndex+1]
					if vcfCall:
						data_row.append(vcfCall['GT'])
					else:
						data_row.append('NA')
						
				writer.writerow(data_row)
				counter += 1
			sys.stderr.write("%s loci outputted.\n"%(counter))
			del writer
	


if __name__ == '__main__':
	main_class = ExampleToFetchVCFFromDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()