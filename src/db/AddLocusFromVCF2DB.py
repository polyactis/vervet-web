#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice -a 524 -s 1
		-i FilterVCF_trioCallerTop7559Contigs.n10.x1.2FoldMedianDepth.2012.5.1T1748/trioCaller_vcftoolsFilter/Contig103.filter_by_vcftools.recode.vcf.gz

Description:
	2012.5.2
		Add locus from one VCF file into database. 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB


class AddLocusFromVCF2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('ref_ind_seq_id', 1, int):[524, 'a', 1, 'IndividualSequence.id for the reference sequence in generating this VCF'],\
							('locus_type_id', 1, int):[1, 's', 1, 'LocusType.id'],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs, **keywords)
		self.sequence2DBentry = {}	#2012.5.2 a cache to store the relevant db entries for sequence
	
	def getSequenceDBEntry(self, db_vervet, sequence=None, comment=None):
		"""
		2012.5.2
		"""
		if sequence in self.sequence2DBentry:
			return self.sequence2DBentry.get(sequence)
		else:
			dbEntry = db_vervet.getSequence(sequence=sequence, comment=comment)
			self.sequence2DBentry[sequence] = dbEntry
			return dbEntry
	
	def addLocusFromVCF2DB(self, db_vervet, inputFname=None, ref_ind_seq_id=None, locus_type_id=None, minDepth=0):
		"""
		2012-5.2
			given a VCF file, find all the loci and submit them into db
		"""
		sys.stderr.write("Adding loci from %s into db ... "%(inputFname))
		from pymodule.VCFFile import VCFFile
		vcfFile = VCFFile(inputFname=inputFname, minDepth=minDepth)
		
		counter = 0
		previous_reported_counter = ''
		for vcfRecord in vcfFile.parseIter():
			chr = vcfRecord.chr
			pos = vcfRecord.pos
			pos = int(pos)
			refBase = vcfRecord.data_row[0].get("GT")[0]
			refBaseDBEntry = self.getSequenceDBEntry(db_vervet, sequence=refBase, comment=None)
			altBase = vcfRecord.altBase
			altBaseDBEntry = self.getSequenceDBEntry(db_vervet, sequence=altBase, comment=None)
			locus = db_vervet.getLocus(chr=chr, start=pos, stop=pos, ref_seq=refBaseDBEntry, alt_seq=altBaseDBEntry, \
							ref_ind_seq_id=ref_ind_seq_id, \
							locus_type_id=locus_type_id)
			counter += 1
			if counter%500==0:
				sys.stderr.write("%s%s"%('\x08'*len(previous_reported_counter), counter))
				previous_reported_counter = repr(counter)
		sys.stderr.write("%s%s"%(len(previous_reported_counter), counter))
		sys.stderr.write(" Done.\n")
	
	def run(self):
		"""
		2012.4.25
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		self.session = self.db_vervet.session
		
		if not self.debug:	#in debug mode, no transaction, auto-commit
			self.session.begin()
		self.addLocusFromVCF2DB(self.db_vervet, inputFname=self.inputFname, ref_ind_seq_id=self.ref_ind_seq_id, locus_type_id=self.locus_type_id)
		
		if not self.debug:
			if self.commit:
				self.session.flush()
				self.session.commit()
			else:
				self.session.rollback()
		

if __name__ == '__main__':
	main_class = AddLocusFromVCF2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()