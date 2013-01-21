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
import numpy as np
from subprocess import call


class CalculateStatsForSubPop(AbstractVervetMapper):
	__doc__ = __doc__
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
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)	#self.connectDB() called within its __init__()
		#to use this class without calling it, we should run selectSubPop here
	
	#def manual_init(self):
	#	{'outputFname': '/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/testplink1.tsv', 'hostname': 'dl324b-1.cmb.usc.edu', 'dataDir': '/media/Data/Akademisches_data/vervetpopgen', 'genotypeMethodID': '27', 'db_passwd': '1o9p2a4', 'db_user': 'hannes', 'statname': 'aaf', 'chromosome': 'Contig917'}
	
	def connectDB(self):
		"""
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname, schema=self.schema, port=self.port)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
	def loadMetadataMat(self,Fname):
		metadata=[]
		meta=open(Fname).readlines()
		for line in meta:
			metadata.append(line.split())
		self.metadata=metadata
		#self.metadata=np.fromfile(file, dtype=float, count=-1, sep='')
	
	def createMetadataMat(self):
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
			self.metadata=[uclaIDList,countryid_row,speciesid_row,longitudeList,latitudeList]
			session.close()
		
		
	
	def selectUCLAids(self,metadata):#cond_dict
		#for now, just take all individuals
		return metadata[1][1:]	
	
		
	def getVCFInd(self,uclaidlist):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
			and return genotype matrix
		"""
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
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			#this is a list with the read-group names
			readgroupIDList = vcfFile.getSampleIDList()
			#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			#header = ['Chromosome', 'position', 'ref','alt']
			ind_id_ls=[]; chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]
			columnIndexList = []
			datalist=[]
			for i in xrange(len(readgroupIDList)):
				readgroupID = readgroupIDList[i]
				#this is the first part of the read group
				individualAlignment = self.db_vervet.parseAlignmentReadGroup(readgroupID).individualAlignment
				uclaid=individualAlignment.ind_sequence.individual.ucla_id
				if uclaid in uclaidlist:			
					#header.append(readgroupID)
					columnIndexList.append(i)
					ind_id_ls.append(uclaid)
			session.close()		
			return (columnIndexList,ind_id_ls)		
					
	def selectSubPopNoDB(self,columnindexlist,ind_id_ls,vcffilename):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
			and return genotype matrix
		"""
		#import pdb
		filename = vcffilename
		if os.path.isfile(filename):
			counter= 0
			from pymodule.VCFFile import VCFFile
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			#this is a list with the read-group names
			readgroupIDList = vcfFile.getSampleIDList()
			#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			#header = ['Chromosome', 'position', 'ref','alt']
			chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]
			columnIndexList = columnindexlist
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
			datastruct=hsContigDataStruct(ind_id_ls=np.array(ind_id_ls), chrom_ls=np.array(chrom_ls),ref_ls=np.array(ref_ls),snp_pos_ls=np.array(snp_pos_ls),alt_ls=np.array(alt_ls), data=data)
			return datastruct
	
	def produceContigFileNameLS(self,contig_ls):
		session = self.db_vervet.session
		session.begin()
		contFName_ls=[]
		if not self.dataDir:
			self.dataDir = self.db_vervet.data_dir
		dataDir = self.dataDir
		
		for chrom in contig_ls:
			genotypeFile = self.db_vervet.getGenotypeFile(genotype_method_id=self.genotypeMethodID, chromosome=chrom, format=self.format)
		
			if not genotypeFile:
				sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
				sys.exit(2)
			filename = os.path.join(dataDir, genotypeFile.path) 
			contFName_ls.append(filename)
		session.close()
		return contFName_ls
					
	def selectSubPop(self,uclaidlist):
		"""
		2012.9.19
			get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
			and return genotype matrix
		"""
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
			
			vcfFile = VCFFile(inputFname=filename, minDepth=0)
			#this is a list with the read-group names
			readgroupIDList = vcfFile.getSampleIDList()
			#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			#header = ['Chromosome', 'position', 'ref','alt']
			ind_id_ls=[]; chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]
			columnIndexList = []
			datalist=[]
			for i in xrange(len(readgroupIDList)):
				readgroupID = readgroupIDList[i]
				#this is the first part of the read group
				individualAlignment = self.db_vervet.parseAlignmentReadGroup(readgroupID).individualAlignment
				uclaid=individualAlignment.ind_sequence.individual.ucla_id
				if uclaid in uclaidlist:			
					#header.append(readgroupID)
					columnIndexList.append(i)
					ind_id_ls.append(uclaid)
			#writer.writerow(header)
			#datalist.append(header)
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
						data_row.append(-9)#missing data
				counter += 1
				datalist.append(data_row)
			sys.stderr.write("%s loci in %i individuals outputted.\n"%(counter,len(columnIndexList)))
			#pdb.set_trace()
			data=np.array(datalist,dtype=np.float)
			datastruct=hsContigDataStruct(ind_id_ls=np.array(ind_id_ls), chrom_ls=np.array(chrom_ls),ref_ls=np.array(ref_ls),snp_pos_ls=np.array(snp_pos_ls),alt_ls=np.array(alt_ls), data=data)
			session.close()
			return datastruct


		
	def thinSNPs(self,datastruct,snpbin):#,r2limit #use snpbin=1000 first
		#first, we only take a snip each snpbin
		selectind=[i for i in range(0,len(datastruct.snp_pos_ls),snpbin)]
		datastructTh=hsContigDataStruct(ind_id_ls=datastruct.ind_id_ls, chrom_ls=datastruct.chrom_ls[selectind],
									ref_ls=datastruct.ref_ls[selectind], snp_pos_ls=datastruct.snp_pos_ls[selectind],
									alt_ls=datastruct.alt_ls[selectind], data=datastruct.data[selectind,:])
		return datastructTh
			
	
		
	def producePLINK(self,datastruct):
		pedFname='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/test2.ped'
		#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		writerPED = csv.writer(open(pedFname, 'w'), delimiter='\t')
		data=datastruct.data
		#---plink .ped file---
		#family id
		col1=np.tile(self.genotypeMethodID,(len(datastruct.ind_id_ls),1))
		col2=datastruct.ind_id_ls
		col3_5=np.tile(0,(len(datastruct.ind_id_ls),3))
		col6=np.tile(-9,(len(datastruct.ind_id_ls),1))
		#colums 7+
		refmat=np.tile(datastruct.ref_ls.reshape((-1,1)),(1,data.shape[1]))
		altmat=np.tile(datastruct.alt_ls.reshape((-1,1)),(1,data.shape[1]))
		gamete1mat=np.where(data-2,refmat,altmat)
		gamete2mat=np.where(data,altmat,refmat)
		pairlist=range(0,2*data.shape[0],2)
		impairlist=range(1,2*data.shape[0],2)
		zygotemat=np.zeros((2*data.shape[0],data.shape[1]),dtype=str)
		zygotemat[pairlist,:]=gamete1mat
		zygotemat[impairlist,:]=gamete2mat
		col7=zygotemat.T
		#stack ped file
		pedmat=np.column_stack((col1,col2,col3_5,col6,col7))
		writerPED.writerows(pedmat)
		#sys.stderr.write("%s loci outputted.\n"%(counter))
		del writerPED
		#---plink .map file---
		mapFname='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/test2.map'
		writerMAP=csv.writer(open(mapFname, 'w'), delimiter='\t')
		col1=datastruct.chrom_ls
		col2=np.core.defchararray.add(datastruct.chrom_ls,np.tile('p',(len(datastruct.chrom_ls),)))
		col2=np.core.defchararray.add(col2,datastruct.snp_pos_ls)
		col3=np.tile(0,(len(datastruct.chrom_ls),))
		col4=datastruct.snp_pos_ls
		#print col1.shape,col2.shape,col3.shape,col4.shape
		mapmat=np.column_stack((col1,col2,col3,col4))
		writerMAP.writerows(mapmat)
		#sys.stderr.write("%s loci outputted.\n"%(counter))
		del writerMAP
	
	def produceTFAM(self,datastruct):
		tfamFname='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Meth40AllCont1.tfam'
		#writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		writerTFAM = csv.writer(open(tfamFname, 'w'), delimiter='\t')
		data=datastruct.data
		col1=np.tile(self.genotypeMethodID,(len(datastruct.ind_id_ls),1))
		col2=datastruct.ind_id_ls
		col3_5=np.tile(0,(len(datastruct.ind_id_ls),3))
		col6=np.tile(-9,(len(datastruct.ind_id_ls),1))
		TFAMmat=np.column_stack((col1,col2,col3_5,col6))
		writerTFAM.writerows(TFAMmat)
		#sys.stderr.write("%s loci outputted.\n"%(counter))
		del writerTFAM
	
	def appendTransPLINK(self,datastruct):
		#---plink  file---
		#import pdb
		data=datastruct.data
		tpedFname='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Meth26AllCont1.tped'
		writerTPED=csv.writer(open(tpedFname, 'a'), delimiter='\t')
		col1=datastruct.chrom_ls
		col2=np.core.defchararray.add(datastruct.chrom_ls,np.tile('p',(len(datastruct.chrom_ls),)))
		col2=np.core.defchararray.add(col2,datastruct.snp_pos_ls)
		col3=np.tile(0,(len(datastruct.chrom_ls),))
		col4=datastruct.snp_pos_ls
		#colums 5+
		refmat=np.tile(datastruct.ref_ls.reshape((-1,1)),(1,data.shape[1]))
		altmat=np.tile(datastruct.alt_ls.reshape((-1,1)),(1,data.shape[1]))
		zeromat=np.tile('0',(data.shape[0],data.shape[1]))
		#pdb.set_trace()
		gamete1mat=np.where(data-2,refmat,altmat)
		gamete1mat=np.where((0==(data+9)),zeromat,gamete1mat)
		gamete2mat=np.where(data,altmat,refmat)
		gamete2mat=np.where((0==(data+9)),zeromat,gamete2mat)
		pairlist=range(0,2*data.shape[1],2)
		impairlist=range(1,2*data.shape[1],2)
		zygotemat=np.zeros((data.shape[0],2*data.shape[1]),dtype=str)
		zygotemat[:,pairlist]=gamete1mat
		zygotemat[:,impairlist]=gamete2mat		
		tpedmat=np.column_stack((col1,col2,col3,col4,zygotemat))
		writerTPED.writerows(tpedmat)
		#sys.stderr.write("%s loci outputted.\n"%(counter))
		del writerTPED
	
	def saveVCF_ids(self,vcfid_ls):
		outFnameID='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Meth26AllCont_VCFind_ls.tsv'
		writerID =csv.writer(open(outFnameID, 'w'), delimiter='\t')
		writerID.writerow(vcfid_ls)
		del writerID
		
	def append01file(self,datastruct):
		data=datastruct.data
		#write 0 1 matrix (e.g. for t-sne)
		Fname01='/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Meth26AllCont_each1000snp_01_mat.tsv'
		writer01mat=csv.writer(open(Fname01, 'a'), delimiter='\t')
		np.where
		#set missing data to reference allele
		refmat=np.tile(datastruct.ref_ls.reshape((-1,1)),(1,data.shape[1]))
		altmat=np.tile(datastruct.alt_ls.reshape((-1,1)),(1,data.shape[1]))
		zeromat=np.tile('0',(data.shape[0],data.shape[1]))
		onemat=np.tile('1',(data.shape[0],data.shape[1]))
		#pdb.set_trace()
		gamete1mat=np.where(data-2,zeromat,onemat)
		gamete1mat=np.where((0==(data+9)),zeromat,gamete1mat)
		gamete2mat=np.where(data,onemat,zeromat)
		gamete2mat=np.where((0==(data+9)),zeromat,gamete2mat)
		pairlist=range(0,2*data.shape[1],2)
		impairlist=range(1,2*data.shape[1],2)
		zygotemat=np.zeros((data.shape[0],2*data.shape[1]),dtype=str)
		zygotemat[:,pairlist]=gamete1mat
		zygotemat[:,impairlist]=gamete2mat		
		writer01mat.writerows(zygotemat)
		#sys.stderr.write("%s loci outputted.\n"%(counter))
		del writer01mat		
	
	def produceBinaryPLINK(self,filename):	
		call(["p-link", "--file",filename,"--make-bed"])
		
	def calculateAAF(self,data):
		"""
		Calculate alternative allele frequency
		"""
		return data.sum(axis=1)/(2.*data.shape[1])
		
	def nonFixedDataList(self,datastruct):
		aaf=self.calculateAAF(datastruct.data)
		nonfixedpos=np.nonzero((aaf>0)*(aaf<1))
		nonfixedpos=nonfixedpos[0]
		datastructNF=hsContigDataStruct(ind_id_ls=datastruct.ind_id_ls, chrom_ls=datastruct.chrom_ls[nonfixedpos],
									ref_ls=datastruct.ref_ls[nonfixedpos], snp_pos_ls=datastruct.snp_pos_ls[nonfixedpos],
									alt_ls=datastruct.alt_ls[nonfixedpos], data=datastruct.data[nonfixedpos,:])
		return datastructNF
		
	def calculateMAF(self,aaflist):
		pass
	
	def run(self):
		"""
		2012.7.13
		"""
		#import pdb
		print self.chromosome
		if self.debug:
			import pdb
			pdb.set_trace()
		self.loadMetadataMat("/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/AfriCaribMetaData.tsv")
		idls=self.selectUCLAids(self.metadata)	
		Fname="/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Method27ContigNames.tsv"
		contiglist=open(Fname).readlines()
		contiglist = map(lambda s: s.strip(), contiglist)
		genotypeFname_ls=self.produceContigFileNameLS(contiglist)
		#pdb.set_trace()
		(indx,ind_id_ls)=self.getVCFInd(idls)
		self.saveVCF_ids(ind_id_ls)
		tok=1
		ind=0 
		for contig in contiglist:
			self.chromosome=contig		
			datastruct=self.selectSubPopNoDB(indx,ind_id_ls,genotypeFname_ls[ind])
			datastruct=self.thinSNPs(datastruct,1000)
			if tok==1:
				self.produceTFAM(datastruct)
			tok=0
			self.appendTransPLINK(datastruct)
			self.append01file(datastruct)
			ind=ind+1
			print ind,"), plink Meth40AllCont1.tped appended with ", len(datastruct.snp_pos_ls), "snps from", contig 
		#self.nonFixedDataList()
		#if self.statname=='aaf':
		#	t123=self.nonFixedDataList()
			
class hsContigDataStruct(object):
	#var1 = 'hi'
	def __init__(self, ind_id_ls=None, chrom_ls=None,ref_ls=None,snp_pos_ls=None,alt_ls=None, data=None):
		"""
		row_id_ls ... individual ids
		col_info_ls ...  
		
		make a length check:(adapt)
		if not genotypeFile:
			sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n"%(self.genotypeMethodID, self.chromosome))
			sys.exit(2)
		"""
		self.ind_id_ls=ind_id_ls
		self.chrom_ls=chrom_ls
		self.ref_ls=ref_ls
		self.snp_pos_ls=snp_pos_ls
		self.alt_ls=alt_ls
		self.data=data
	

if __name__ == '__main__':
	main_class = CalculateStatsForSubPop
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	print po.arguments
	print po.long_option2value
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()