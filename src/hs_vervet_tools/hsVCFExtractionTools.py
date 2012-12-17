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
__doc__ = __doc__ % (sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
import numpy as np
from subprocess import call
from vervet.src.hs_vervet_tools import hs_support_tools as hs


class hsVCFExtractionTools(AbstractVervetMapper):
	__doc__ = __doc__
	"""
	option_default_dict = AbstractMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.update({
							('statname', 1,):['aaf', '', 1, 'which statistic to calculate? aaf... list of alternative allele frequencies', ],\
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
	"""						
							
	def __init__(self, projectName, dataDir='/net/gmi.oeaw.ac.at/nordborg/lab/vervetpopgen'):
		"""
		projectName... name of the directory in which all the output of this calculation is stored
		 and where previous results are searched for
		"""
		self.projectName = projectName
		self.projectDir = os.path.join(dataDir, 'analyses', projectName)
		self.dataDir = dataDir
		# create project dir if it doesn't exist:
		import errno
		try:
			os.makedirs(self.projectDir)
		except OSError as exc:  # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(self.projectDir):
				pass
			else: raise
        
		# AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		# to use this class without calling it, we should run selectSubPop here
	
	# def manual_init(self):
	# 	{'outputFname': '/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/testplink1.tsv', 'hostname': 'dl324b-1.cmb.usc.edu', 'dataDir': '/media/Data/Akademisches_data/vervetpopgen', 'genotypeMethodID': '27', 'db_passwd': '1o9p2a4', 'db_user': 'hannes', 'statname': 'aaf', 'chromosome': 'Contig917'}
	
	
	def createProjectDescription(self,description_str,exist_action='a'):
		"""
		exist_action...
		'w'... overwrite existing file
		'a'... append to existing file
		'p'... pass if the file exists
		
		description_str... string describing the aim of this analysis, with trailing \n
		"""
		fname = os.path.join(self.projectDir,'_ProjectDescription.txt')		
		hs.safeWrite(self,fname,description_str,exist_action)
		
	
	def DBobject(self, drivername='postgresql', hostname='crocea.mednet.ucla.edu', \
				db_user='hannes', db_passwd='1o9p2a4', schema='public', dbname='vervetdb', port=None):
		"""
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		db_vervet = VervetDB.VervetDB(drivername=drivername, username=db_user, password=db_passwd, \
									hostname=hostname, database=dbname, schema=schema, port=port)
		db_vervet.setup(create_tables=False)		
		return db_vervet
	
	def saveUCLA_id_ls(self,UCLA_id_ls,exist_action):
		fname = os.path.join(self.projectDir,'UCLA_ids.tsv')
		if exist_action=='check':
			if os.path.isfile(fname):
				return True
			else:
				return False
		elif exist_action=='p':
			if os.path.isfile(fname):
				pass
			else:
				with open(fname, 'w') as f:
					writ=csv.writer(f, delimiter='\t')
					writ.writerow(UCLA_id_ls)
					del writ		
		elif exist_action=='w':
			with open(fname, 'w') as f:
				writ=csv.writer(f, delimiter='\t')
				writ.writerow(UCLA_id_ls)
				del writ
		else:
			raise hs.hsError('variable existAction must be w (overwrite) or p (pass) or check')
	
	def loadUCLA_id_ls(self):
		fname = os.path.join(self.projectDir,'UCLA_ids.tsv')
		with open(fname, 'r') as f:
			UCLA_id_ls=f.readline().split()
		return 	UCLA_id_ls
	
		
	def getVCFInd(self, uclaidlist, genotypeMethodID, chromosome, format1='VCF'):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
			and return genotype matrix
		"""
		db_vervet = self.DBobject()
		session = db_vervet.session
		
		session.begin()
		# if not self.dataDir:
		# 	self.dataDir = self.db_vervet.data_dir
		dataDir = self.dataDir
		
		genotypeFile = db_vervet.getGenotypeFile(genotype_method_id=genotypeMethodID, chromosome=chromosome, format=format1)
		
		if not genotypeFile:
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n" % (genotypeMethodID, chromosome))
			sys.exit(2)
		filename = os.path.join(dataDir, genotypeFile.path)
		if os.path.isfile(filename):
			from pymodule.yhio.VCFFile import VCFFile
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			# this is a list with the read-group names
			readgroupIDList = vcfFile.getSampleIDList()
			
			new_ucla_id_ls = []
			columnIndexList = []
	
			for i in xrange(len(readgroupIDList)):
				readgroupID = readgroupIDList[i]
				# this is the first part of the read group
				individualAlignment = db_vervet.parseAlignmentReadGroup(readgroupID).individualAlignment
				uclaid = individualAlignment.individual_sequence.individual.ucla_id
				if uclaid in uclaidlist:			
					# header.append(readgroupID)
					columnIndexList.append(i)
					new_ucla_id_ls.append(str(uclaid))
			session.close()		
			return (columnIndexList, new_ucla_id_ls)		
					
	def saveVCFInd(self,VCF_index_ls,genotype_method_id,exist_action):
		"""
		exist_action...
		'w'... overwrite existing file
		'p'... pass if the file exists
		
		description_str... string describing the aim of this analysis, with trailing \n
		"""
		fname = os.path.join(self.projectDir,'VCF_indices_Meth'+str(genotype_method_id)+'.tsv')		
		if exist_action=='p':
			if os.path.isfile(fname):
				pass
			else:
				with open(fname, 'w') as f:
					writ=csv.writer(f, delimiter='\t')
					writ.writerow(VCF_index_ls)
					del writ		
		elif exist_action=='w':
			with open(fname, 'w') as f:
				writ=csv.writer(f, delimiter='\t')
				writ.writerow(VCF_index_ls)
				del writ
		else:
			raise hs.hsError('variable existAction must be w (overwrite) or p (pass)')	
	
	def createVCFfilename_ls(self):
		pass
	
	def createGenotypeData(self,vcf_indx_ls,ucla_id_ls,vcffilename,genotype_method_id,contig=None):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
			and return genotype matrix
		"""
		#import pdb
		filename = vcffilename
		if os.path.isfile(filename):
			counter= 0
			from pymodule.yhio.VCFFile import VCFFile
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			#this is a list with the read-group names
			#readgroupIDList = vcfFile.getSampleIDList()
			#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			#header = ['Chromosome', 'position', 'ref','alt']
			chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]
			columnIndexList = vcf_indx_ls
			datalist=[]
			for vcfRecord in vcfFile:
				data_row=[]
				chrom_ls.append(vcfRecord.chr)
				snp_pos_ls.append(vcfRecord.pos)
				refBase = vcfRecord.refBase
				nonRefBase = vcfRecord.altBase
				ref_ls.append(refBase)
				alt_ls.append(nonRefBase)
				for columnIndex in columnIndexList:
					#for vcfCall in vcfRecord.data_row[1:]: #data_row is a list of dictionary {'GT': base-call, 'DP': depth}, or None if missing.
					#it includes one extra sample in index 0, which is the reference individual (from the ref column of VCF).
					vcfCall = vcfRecord.data_row[columnIndex+1]
					if vcfCall:
						if vcfCall['GT'][0]==refBase and vcfCall['GT'][1]==refBase:
							gt=0
						elif vcfCall['GT'][0]==refBase or vcfCall['GT'][1]==refBase:
							gt=1
						else:
							gt=2
						data_row.append(gt)
					else:
						data_row.append(-9)
				counter += 1
				datalist.append(data_row)
			sys.stderr.write("%s loci in %i individuals outputted.\n"%(counter,len(columnIndexList)))
			#pdb.set_trace()
			data=np.array(datalist,dtype=np.float)
			datastruct=hsContigDataStruct(ucla_id_ls=np.array(ucla_id_ls), chrom_ls=np.array(chrom_ls),\
										ref_ls=np.array(ref_ls),snp_pos_ls=np.array(snp_pos_ls),alt_ls=np.array(alt_ls),\
										data=data,genotype_method=genotype_method_id,contig=contig)
			return datastruct		
	
	def saveContigGenotypeData(self,datastruct):
		if datastruct.contig==False:
			print 'datastruct has no contig specified, can\'t save with saveContigGenotypeData'
		else:
			fnameDat = os.path.join(self.projectDir,'data_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			fnamePos = os.path.join(self.projectDir,'SNPpos_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			fnameRefAlt = os.path.join(self.projectDir,'RefAlt_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			writerDat=csv.writer(open(fnameDat, 'w'), delimiter='\t')
			writerDat.writerows(datastruct.data)
			del writerDat
			writerPos=csv.writer(open(fnamePos, 'w'), delimiter='\t')
			writerPos.writerow(datastruct.chrom_ls)
			writerPos.writerow(datastruct.snp_pos_ls)
			del writerPos
			writerRefAlt=csv.writer(open(fnameRefAlt, 'w'), delimiter='\t')
			writerRefAlt.writerow(datastruct.ref_ls)
			writerRefAlt.writerow(datastruct.alt_ls)
			del writerRefAlt
			
					
	
	def run(self):
		"""
		2012.7.13
		"""
		pass
			
class hsContigDataStruct(object):
	# var1 = 'hi'
	def __init__(self, ind_id_ls=None, chrom_ls=None, ref_ls=None, snp_pos_ls=None, alt_ls=None, data=None,genotype_method=None,contig=None):
		"""
		row_id_ls ... individual ids
		col_info_ls ...  
		
		make a length check:(adapt)
		if not genotypeFile:
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
			sys.exit(2)
		"""
		self.ind_id_ls = ind_id_ls
		self.chrom_ls = chrom_ls
		self.ref_ls = ref_ls
		self.snp_pos_ls = snp_pos_ls
		self.alt_ls = alt_ls
		self.data = data
		self.genotype_method=genotype_method
		self.contig=contig
		#self.chrom_length= figure out length of the contigs!!!
	

if __name__ == '__main__':
	"""
	main_class = CalculateStatsForSubPop
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	print po.arguments
	print po.long_option2value
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
	"""
	print "main"
