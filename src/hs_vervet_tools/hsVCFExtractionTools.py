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
#from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils
#from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
import numpy as np
#from subprocess import call
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
	
	def loadVCFInd(self,genotype_method_id):
		fname = os.path.join(self.projectDir,'VCF_indices_Meth'+str(genotype_method_id)+'.tsv')
		with open(fname, 'r') as f:
			columnIndexList=f.readline().split()
		return columnIndexList	
	
	def createVCFfilename_ls(self,genotypeMethodID, format1='VCF'):
		"""
		create a list of VCF Filenames for each contig
		"""
		db_vervet = self.DBobject()
		session = db_vervet.session
		
		#session.begin()
		# if not self.dataDir:
		# 	self.dataDir = self.db_vervet.data_dir
		#dataDir = self.dataDir
		
		query = VervetDB.GenotypeFile.query.filter_by(genotype_method_id=genotypeMethodID)
		contig_ls=[]
		VCFfilename_ls=[]
		for entry in query:
			VCFfilename_ls.append(str(entry.path))
			contig_ls.append(str(entry.chromosome))
		session.close()
			
		return (VCFfilename_ls,contig_ls)	
	
	def saveVCFfilename_ls(self,VCFfilename_ls,contig_ls,genotypeMethodID,exist_action):
		fname = os.path.join(self.projectDir,'VCF_filename_ls_Meth'+str(genotypeMethodID)+'.tsv')
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
					writer=csv.writer(f, delimiter='\t')
					writer.writerow(VCFfilename_ls)
					writer.writerow(contig_ls)
					del writer		
		elif exist_action=='w':
			with open(fname, 'w') as f:
				writer=csv.writer(f, delimiter='\t')
				writer.writerow(VCFfilename_ls)
				writer.writerow(contig_ls)
				del writer	
		else:
			raise hs.hsError('variable existAction must be w (overwrite) or p (pass) or check')
			
	def loadVCFfilename_ls(self,genotypeMethodID):
		"""
		load a list of VCF Filenames for each contig
		"""
		fname = os.path.join(self.projectDir,'VCF_filename_ls_Meth'+str(genotypeMethodID)+'.tsv')
		with open(fname, 'r') as f:
			VCFfilename_ls=f.readline().split()
			contig_ls=f.readline().split()
		return (VCFfilename_ls,contig_ls)

				
	
	def createGenotypeData(self,vcf_indx_ls,ucla_id_ls,vcffilename,genotype_method_id,contig=None):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with vcf file index in vcf_indx_ls
			and return genotype matrix
		"""
		#import pdb
		filename = os.path.join(self.dataDir,vcffilename)
		if os.path.isfile(filename):
			counter= 0
			from pymodule.yhio.VCFFile import VCFFile
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			
				
			#this is a list with the read-group names
			#readgroupIDList = vcfFile.getSampleIDList()
			#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			#header = ['Chromosome', 'position', 'ref','alt']
			chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]; quality_ls=[];  info_ls=[]
			columnIndexList = map(int,vcf_indx_ls)
			datalist=[]
			for vcfRecord in vcfFile:
				#objects in VCFRecord: [vcfRecord.chr, vcfRecord.pos, vcfRecord.vcf_locus_id, vcfRecord.refBase, vcfRecord.altBase, \
#					vcfRecord.quality, vcfRecord.filter,\
#					vcfRecord.info, vcfRecord.format] + vcfRecord.row[self.sampleStartingColumn:]				
				data_row=[]
				chrom_ls.append(vcfRecord.chr)
				snp_pos_ls.append(vcfRecord.pos)
				refBase = vcfRecord.refBase
				nonRefBase = vcfRecord.altBase
				ref_ls.append(refBase)
				alt_ls.append(nonRefBase)
				quality_ls.append(vcfRecord.quality)
				info_ls.append(vcfRecord.info)
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
										quality_ls=np.array(map(float,quality_ls)),info_ls=np.array(info_ls),\
										data=data,genotype_method=genotype_method_id,contig=contig)
			return datastruct		
	
	def calculateAAF(self,data):
		"""
		Calculate alternative allele frequency
		"""
		return data.sum(axis=1)/(2.*data.shape[1])
	
	def nonFixedGenotypeData(self,datastruct):
		aaf=self.calculateAAF(datastruct.data)
		nonfixedpos=np.nonzero((aaf>0)*(aaf<1))
		nonfixedpos=nonfixedpos[0]
		datastructNF=hsContigDataStruct(ucla_id_ls=datastruct.ucla_id_ls, chrom_ls=datastruct.chrom_ls[nonfixedpos],
									ref_ls=datastruct.ref_ls[nonfixedpos], snp_pos_ls=datastruct.snp_pos_ls[nonfixedpos],
									alt_ls=datastruct.alt_ls[nonfixedpos], data=datastruct.data[nonfixedpos,:],\
									quality_ls=datastruct.quality_ls[nonfixedpos],info_ls=datastruct.info_ls[nonfixedpos],\
									genotype_method=datastruct.genotype_method,contig=datastruct.contig)
		return datastructNF
	
	
	def saveContigGenotypeData(self,datastruct):
		if datastruct.contig==False:
			print 'datastruct has no contig specified, can\'t save with saveContigGenotypeData'
		else:
			fnameDat = os.path.join(self.projectDir,'data_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			
			writerDat=csv.writer(open(fnameDat, 'w'), delimiter='\t')
			writerDat.writerows(datastruct.data)
			del writerDat
			
			fnamePos = os.path.join(self.projectDir,'SNPpos_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			writerPos=csv.writer(open(fnamePos, 'w'), delimiter='\t')
			writerPos.writerow(datastruct.chrom_ls)
			writerPos.writerow(datastruct.snp_pos_ls)
			del writerPos
			
			fnameRefAlt = os.path.join(self.projectDir,'RefAlt_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			writerRefAlt=csv.writer(open(fnameRefAlt, 'w'), delimiter='\t')
			writerRefAlt.writerow(datastruct.ref_ls)
			writerRefAlt.writerow(datastruct.alt_ls)
			del writerRefAlt
			
			fnameQualMeta = os.path.join(self.projectDir,'QualMeta_Meth'+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
			writerQualMeta=csv.writer(open(fnameQualMeta, 'w'), delimiter='\t')
			writerQualMeta.writerow(datastruct.quality_ls)
			writerQualMeta.writerow(datastruct.info_ls)
			del writerQualMeta
	
	
	"""			
	def saveContigGenotypeData2(self,datastruct):
		if datastruct.contig==False:
			print 'datastruct has no contig specified, can\'t save with saveContigGenotypeData'
		else:
			for attribute, value in vars(datastruct).iteritems():
				if (value =! None) and attribute =! "contig"
				fname = os.path.join(self.projectDir,str(attribute)+str(datastruct.genotype_method)+\
								'_'+datastruct.contig+'.tsv')
				with open(fnameDat, 'w') as f:
					writer=csv.writer(f, delimiter='\t')
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
		"""			
	
	def parseVariantInfo(self,info_ls):
		"""
		this function takes the np.array of the VCF info column and returns a dictionary of 
		VCF info lists, where each entry corresponds to a variant
		"""
		dict_ls=[]
		for variant_entry in info_ls:
			variant_entry2=variant_entry.split(';')
			dict_creator=[]
			for entry in variant_entry2:
				entry2=entry.split("=")
				try:
					entry2[1]=int(entry2[1])
				except:
					try:
						entry2[1]=float(entry2[1])
					except:
						pass
				dict_creator.append(tuple(entry2))
			dict_ls.append(dict(dict_creator))
#			dict_ls=map(lambda s: s.split("="),dict_ls)
		
		#convert list of dicts into dict of lists:
		#numpy array has problems with mixed data types...
		key_ls=dict_ls[0].keys()
		empt_array_ls=[np.empty(len(dict_ls),dtype=type(dict_ls[0][key])) for key in key_ls]
		for el in empt_array_ls:
			el[:]=np.NAN
#		empt_array[:]=np.NAN
		info_dic=dict([(key,empt_array) for key,empt_array in zip(key_ls,empt_array_ls)])
		for (indx,dic) in enumerate(dict_ls):
			for key in key_ls:
#				print key,indx,dic[key],info_dic[key][indx]
				try:
					info_dic[key][indx]=dic[key]
				except:
					pass		
		return info_dic
	
	def parseVariantInfoOld(self,info_ls):
		"""
		this function takes the np.array of the VCF info column and returns a dictionary of 
		VCF info lists, where each entry corresponds to a variant
		"""
		dict_ls=[]
		for variant_entry in info_ls:
			variant_entry2=variant_entry.split(';')
			dict_creator=[]
			for entry in variant_entry2:
				entry2=entry.split("=")
				try:
					entry2[1]=int(entry2[1])
				except:
					try:
						entry2[1]=float(entry2[1])
					except:
						pass
				dict_creator.append(tuple(entry2))
			dict_ls.append(dict(dict_creator))
#			dict_ls=map(lambda s: s.split("="),dict_ls)
		return dict_ls
	
	def run(self):
		"""
		2012.7.13
		"""
		pass
			
class hsContigDataStruct(object):
	# var1 = 'hi'
	def __init__(self, ucla_id_ls=None, chrom_ls=None, ref_ls=None, snp_pos_ls=None, alt_ls=None, \
				quality_ls=None,info_ls=None,\
				data=None,genotype_method=None,contig=None):
		#row_id_ls ... individual idsq
		"""
		col_info_ls ...  
		
		make a length check:(adapt)
		if not genotypeFile:
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
			sys.exit(2)
		"""
		self.ucla_id_ls = ucla_id_ls
		self.chrom_ls = chrom_ls
		self.ref_ls = ref_ls
		self.snp_pos_ls = snp_pos_ls
		self.alt_ls = alt_ls
		self.quality_ls=quality_ls
		self.info_ls=info_ls
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
