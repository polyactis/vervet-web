#!/usr/bin/env python
"""
Examples:
	# enter debug mode
	misc.py -b
	
	# change the db connection setting 
	misc.py -z localhost -u yh
	
Description:
	It's a file with all sorts of classes that haven't gone on to be independent.
	
	All those classes are manually chosen to be run in Main.run().
"""
import os, sys, numpy

class VariantDiscovery(object):
	"""
	2011-1-6
	"""
	def __init__(self):
		pass
	
	@classmethod
	def discoverHets(cls, input_fname, output_fname, minNumberOfReads=4):
		import csv
		reader =csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['id', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		writer.writerow(header)
		reader.next()
		reader.next()
		for row in reader:
			chr = row[0][8:]
			pos = row[1]
			quality = row[5]
			info = row[7]
			info_ls = info.split(';')
			info_tag2value = {}
			outputHet= False
			for info in info_ls:
				try:
					tag, value = info.split('=')
				except:
					sys.stderr.write("Error %s.\n"%info)
					continue
				if tag=='DP4':
					tag = 'DP4_ratio'
					value = value.split(',')
					value = map(int, value)
					no_of_ref_allele = sum(value[0:2])
					no_of_non_ref_allele = sum(value[2:])
					if no_of_ref_allele>=minNumberOfReads and no_of_non_ref_allele>=minNumberOfReads:
						outputHet = True
						value = float(no_of_ref_allele)/no_of_non_ref_allele
						info_tag2value['minDP4'] = min(no_of_ref_allele, no_of_non_ref_allele)
					else:
						value = None
				info_tag2value[tag] = value
			if outputHet:
				output_row = ['chr%s:%s'%(chr, pos), chr, pos, quality, info_tag2value.get('DP'), \
							info_tag2value.get('minDP4'), info_tag2value.get('DP4_ratio'), info_tag2value.get('MQ')]
				writer.writerow(output_row)
		del reader, writer
	"""
		#2011-1-6
		input_fname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.vcf'
		output_fname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		VariantDiscovery.discoverHets(input_fname, output_fname, minNumberOfReads=4)
		
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_p1'
		input_fname = '%s.raw.vcf'%(common_prefix)
		minNumberOfReads=2
		output_fname = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHets(input_fname, output_fname, minNumberOfReads=minNumberOfReads)
		sys.exit(0)
	"""
	
	@classmethod
	def checkOverlapping(cls, myVariantFile, jessicaVariantFname, output_fname):
		"""
		2011-1-6
		"""
		import csv
		sys.stderr.write("Reading in Jessica's variants ... ")
		reader = csv.reader(open(jessicaVariantFname), delimiter='\t')
		reader.next()	#skip header
		snp_id2work_status = {}	
		chr_pos2data_ls = {}	#map chr_pos to [snp_id, work_status]
		for row in reader:
			if row[3]=='':	# the genotyping work/failure section. it comes first
				snp_id = row[0]
				genotyped = int(row[1])
				status = int(row[2])
				if genotyped==1:
					snp_id2work_status[snp_id] = status
				
			else:	# the snp_id, chromosome, position part. 2nd part. at this stage, snp_id2work_status should be ready.
				snp_id = row[1]
				chr_pos = row[2]
				if snp_id in snp_id2work_status:
					chr_pos2data_ls[chr_pos] = [snp_id, snp_id2work_status.get(snp_id)]
		del reader
		sys.stderr.write("%s genotyped jessica variants. Done.\n"%(len(chr_pos2data_ls)))
		
		sys.stderr.write("Reading in my variants ... ")
		reader = csv.reader(open(myVariantFile), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['snp_id', 'chr_pos', 'status', 'chr', 'pos', 'qual', 'DP', 'minDP4', 'DP4_ratio', 'MQ']
		writer.writerow(header)
		
		from pymodule.utils import getColName2IndexFromHeader
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter = 0
		for row in reader:
			chr_pos = row[col_name2index['id']]
			pos = row[col_name2index['pos']]
			if chr_pos in chr_pos2data_ls:
				snp_id, work_status = chr_pos2data_ls[chr_pos]
				data_row = [snp_id, chr_pos, work_status] + row[1:]
				writer.writerow(data_row)
				real_counter += 1
			
			counter += 1
		del writer
		sys.stderr.write("%s/%s of mine overlapped with jessica's. Done.\n"%(real_counter, counter))
	"""
		#2011-1-6
		myVariantFile = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		output_fname = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.overlap.tsv'
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, output_fname)
		sys.exit(0)
	"""

class Main(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 0, ): ['', 'i', 1, 'common input file.', ],\
							('output_fname', 0, ): ['', 'o', 1, 'common output file', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		This class is the entry to all others.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self,):
		if self.debug:	# 2010-4-18 enter debug mode "~/.../variation/misc.py -b"
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		#import Stock_250kDB
		#db_250k = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
		#				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		#db_250k.setup(create_tables=False)
		#self.db_250k = db_250k
		
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		#curs = conn.cursor()
		
		#2011-1-6
		common_prefix = '/Network/Data/vervet/ref/454_vs_hg19_20101230_ANXA_p1'
		input_fname = '%s.raw.vcf'%(common_prefix)
		minNumberOfReads=3
		output_fname = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		VariantDiscovery.discoverHets(input_fname, output_fname, minNumberOfReads=minNumberOfReads)
		sys.exit(0)
		
		#2011-1-6
		myVariantFile = '/Network/Data/vervet/ref/454_vs_hg19_20101230_3eQTL_D100.raw.hets'
		myVariantFile = '%s_min%s.raw.hets'%(common_prefix, minNumberOfReads)
		jessicaVariantFname = os.path.expanduser('~/script/vervet-web/data/eQTL summary.txt')
		output_fname = '%s_min%s.overlap.tsv'%(common_prefix, minNumberOfReads)
		VariantDiscovery.checkOverlapping(myVariantFile, jessicaVariantFname, output_fname)
		sys.exit(0)
		
		
		
		
		#2010-11-24
		CNV.updateCNVQCCallProbeInfo(db_250k, cnv_method_id_ls=[9])
		sys.exit(0)
		


#2007-03-05 common codes to initiate database connection
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	#sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	#sys.path.insert(0, os.path.expanduser('~/lib/python'))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	#sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment

"""
from codense.common import db_connect, form_schema_tables
hostname='dl324b-1'
dbname='yhdb'
schema = 'dbsnp'
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
#conn, curs = db_connect(hostname, dbname, schema)

hostname='localhost'
dbname='stock20070829'
import MySQLdb
conn0 = MySQLdb.connect(db=dbname,host=hostname)
curs0 = conn0.cursor()

hostname='papaya.usc.edu'
dbname='stock_250k'
db_user='yh'
db_passwd = ''
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname, user=db_user, passwd=db_passwd)
curs = conn.cursor()

drivername='mysql'
schema = None
import Stock_250kDB
db_250k = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_250k.setup(create_tables=False)

drivername='mysql'
hostname='papaya.usc.edu'
dbname='stock'
db_user='yh'
db_passwd = ''
schema = None
from StockDB import StockDB
db_149 = StockDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_149.setup(create_tables=False)

"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Main
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
