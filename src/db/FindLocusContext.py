#!/usr/bin/env python
"""

Examples:
	%s -m 1 -u yh -n /tmp/genomeRBDictTaxID60711.pickle  -c
	
	%s 
Description:
	2012.5.2
	program to find the context (nearby genes) of a locus. It fills results into db table locus_context.
	It also annotates structural variation or SNPs (syn, non-syn, intron, etc.) db table locus_annotation.
	

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData, GenomeDB
from pymodule.db import formReadmeObj
from pymodule.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, get_overlap_ratio
from pymodule.RBTree import RBDict
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
from Bio.Seq import Seq	#translate cds sequence, decide whether synonymous/non SNP
from Bio.Alphabet import IUPAC
from pymodule.SNP import nt2complement	#to complement single nucleotide

class FindLocusContext(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('genome_db_schema', 1, ): ['genome', 'g', 1, 'genome schema name within db (postgresql)', ],\
							('locus_type_id', 1, int): [None, 'm', 1, 'construct contexts for loci from this locus_type_id'],\
							('genomeRBDictPickleFname', 1, ): ['', 'n', 1, 'The file to contain pickled genomeRBDict.'],\
							('max_distance', 0, int): [20000, 'x', 1, "maximum distance allowed between a CNV and a gene"],\
							('tax_id', 0, int): [60711, 'a', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']
							})
	
	def __init__(self,  **keywords):
		"""
		2012.5.05
		"""
		
		inputFnameLs = []
		AbstractVervetMapper.__init__(self, inputFnameLs, **keywords)
	
	def getLocusAnnotationShortName2dbEntry(self, db_vervet):
		"""
		2012.5.14
		"""
		sys.stderr.write("Getting locus_annotation_short_name2db_entry ...")
		locus_annotation_short_name2db_entry = {}
		rows = VervetDB.LocusAnnotationType.query.all()
		for row in rows:
			locus_annotation_short_name2db_entry[row.short_name] = row
		sys.stderr.write("Done.\n")
		return locus_annotation_short_name2db_entry
	
	def _constructSNPAnnotation(self, db_vervet, locus=None, oneGeneData=None, locus_context=None, locus_annotation_short_name2db_entry=None,\
							geneCommentaryRBDict=None, geneSegmentKey = None, compareIns=None, param_obj=None):
		"""
		2012.5.18
			add which_codon into the db_vervet.getLocusAnnotation()
		2012.5.14
			adapted from variation/src/ConstructSNPAnnotation.py
			
			locus is of type VervetDB.Locus
			oneGeneData = PassingData(strand = row.strand, gene_id = row.id, gene_start = row.start, \
										gene_stop = row.stop, geneCommentaryRBDictLs=[],\
										ncbi_gene_id=row.ncbi_gene_id)
		2009-2-5
			bug fixed. when adding a box as UTR, make sure after the forloop above, the current box is still the UTR.
		2009-1-5
		"""
		sys.stderr.write("Constructing LocusAnnotation for SNP ...\n")
		#standard_translator = Translate.unambiguous_dna_by_id[1]	#2012.5.23 Translate is to be deprecated.
		counter = 0
		real_counter = 0
		counter += 1
		chr=locus.chromosome
		pos = locus.start
		allele1 = locus.ref_seq.sequence.encode('ascii')	#2012.5.17 unicode is not accepted by cds_seq.tomutable()
		allele2 = locus.alt_seq.sequence.encode('ascii')
		snp_annotation_type_short_name_ls = []	#each element is (snp_annotation_type_short_name, gene_id, gene_commentary_id, 
				# which_exon_or_intron, pos_within_codon)
		disp_pos = locus_context.disp_pos
		gene_id = locus_context.gene_id
		#gene_box_node_ls = []
		#geneCommentaryRBDict.findNodes(segmentKey, node_ls=gene_box_node_ls, compareIns=compareIns)
		gene_commentary_id = geneCommentaryRBDict.gene_commentary_id
		box_ls = geneCommentaryRBDict.box_ls
		protein_box_ls = geneCommentaryRBDict.protein_box_ls
		#gene_commentary = GeneCommentary.get(gene_commentary_id)
		#box_ls = gene_commentary.construct_annotated_box()
		cds_sequence = geneCommentaryRBDict.cds_sequence
		geneCommentaryRBDict.gene_commentary_type_name
		CDS_5_end_pos = geneCommentaryRBDict.CDS_5_end_pos
		no_of_introns = geneCommentaryRBDict.no_of_introns
		
		detailed_box_type = geneSegmentKey.label
		is_translated = geneSegmentKey.is_translated
		which_intron = geneSegmentKey.intron_number	#which intron the SNP resides, starting from 1
		which_coding_exon = geneSegmentKey.cds_number	#which exon the SNP resides in terms of the CDS sequence, starting from 1
		cumulativeWithinCDSUTRAndIntronLen = geneSegmentKey.cumulativeWithinCDSUTRAndIntronLen
		gene_segment_id = geneSegmentKey.gene_segment_id
		
		if protein_box_ls:	#this is a protein coding gene
			if detailed_box_type.find('UTR')>=0 and detailed_box_type=='exon' and is_translated==0:	#it's UTR. bug fixed. make sure after the forloop above, 
					# the current box is still the UTR.
				snp_annotation_type_short_name_ls.append((detailed_box_type, gene_id, gene_commentary_id, None , None, None))
			else:
				if oneGeneData.strand=='-1':	#reverse the order of exon/intron
					no_of_coding_exons = len(protein_box_ls)
					#no_of_introns = len(gene_commentary.mrna_box_ls)-no_of_coding_exons	#not right
					which_coding_exon = no_of_coding_exons-which_coding_exon+1
					which_intron = no_of_introns-which_intron + 1
				if detailed_box_type=='intron':
					snp_annotation_type_short_name_ls.append(('intron', gene_id, gene_commentary_id, which_intron, None, None))
					if pos-geneSegmentKey.start<=1:	#within the donor/acceptor two-nucleotide
						if oneGeneData.strand=='-1':
							snp_annotation_type_short_name = 'splice-acceptor'	#on the 3' of this intron
						else:
							snp_annotation_type_short_name = 'splice-donor'	#on the 5' of this intron
						snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, \
																gene_commentary_id, which_intron, None, None))
					elif geneSegmentKey.stop-pos<=1:
						if oneGeneData.strand=='-1':
							snp_annotation_type_short_name = 'splice-donor'	#on the 5' of this intron
						else:
							snp_annotation_type_short_name = 'splice-acceptor'	#on the 3' of this intron
						snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, \
																gene_commentary_id, which_intron, None, None))
				elif detailed_box_type=='CDS':	#must be translated
					SNP_index_in_CDS = pos - CDS_5_end_pos - cumulativeWithinCDSUTRAndIntronLen
					if oneGeneData.strand=='-1':	#reverse
						SNP_index_in_CDS = len(cds_sequence)-SNP_index_in_CDS-1	#-1 because SNP_index_in_CDS starts from 0
						gene_allele1 = nt2complement[allele1]
						gene_allele2 = nt2complement[allele2]
					else:
						gene_allele1 = allele1
						gene_allele2 = allele2
					SNP_index_in_CDS = int(SNP_index_in_CDS)	
					# SNP_index_in_CDS is type long. without int(), cds_seq[SNP_index_in_CDS] returns a Bio.Seq with one nucleotide, 
					# rather than a single-char string
					SNP_index_in_peptide = SNP_index_in_CDS/3
					
					SNP_index_in_peptide = int(SNP_index_in_peptide)	#ditto as int(SNP_index_in_CDS), not necessary
					
					pos_within_codon = SNP_index_in_CDS%3+1	#pos_within_codon starts from 1
					cds_seq = Seq(cds_sequence, IUPAC.unambiguous_dna)
					if SNP_index_in_CDS>=len(cds_seq):
						sys.stderr.write("Error: SNP (%s, %s), SNP_index_in_CDS=%s, is beyond any of the boxes from gene %s (chr=%s, %s-%s), \
								gene_commentary_id %s (%s-%s), box_ls=%s, cds-length=%s. counted as intergenic.\n"%\
								(chr, pos, SNP_index_in_CDS, gene_id, oneGeneData.chromosome, oneGeneData.gene_start, oneGeneData.gene_stop, \
								gene_commentary_id, geneCommentaryRBDict.start, geneCommentaryRBDict.stop, repr(box_ls), \
								len(cds_seq)))
						sys.exit(3)
						snp_annotation_type_short_name_ls.append(['intergenic'])
					if cds_seq[SNP_index_in_CDS]!=gene_allele1 and cds_seq[SNP_index_in_CDS]!=gene_allele2:
						sys.stderr.write("Error: Neither allele (%s, %s) from SNP (%s,%s) matches the nucleotide, %s, from the cds seq of gene %s \
							(gene_commentary_id=%s).\n"%\
							(gene_allele1, gene_allele2, chr, pos, cds_seq[SNP_index_in_CDS], gene_id, gene_commentary_id))
						sys.exit(3)
					cds_mut_ar = cds_seq.tomutable()
					cds_mut_ar[SNP_index_in_CDS] = gene_allele1
					peptide = cds_mut_ar.toseq().translate()	#2012.5.23 no more translator. table=1
					
					alt_cds_mut_ar = cds_seq.tomutable()
					alt_cds_mut_ar[SNP_index_in_CDS] = gene_allele2
					alt_peptide = alt_cds_mut_ar.toseq().translate()
					aa = peptide[SNP_index_in_peptide]
					alt_aa = alt_peptide[SNP_index_in_peptide]
					if aa != alt_aa:
						snp_annotation_type_short_name = 'non-synonymous'
						comment = '%s->%s'%(aa, alt_aa)
					else:
						snp_annotation_type_short_name = 'synonymous'
						comment = None
					snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary_id, \
															which_coding_exon, pos_within_codon, comment, SNP_index_in_peptide))
					
					if aa != alt_aa:
						if aa=='*' or alt_aa=='*':
							snp_annotation_type_short_name = 'premature-stop-codon'	#could also be the last stop codon changing to something else 
							# and thereby extending the cds
							snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary_id, \
															which_coding_exon, pos_within_codon, comment, SNP_index_in_peptide))
						if SNP_index_in_peptide==0:
							snp_annotation_type_short_name = 'init-Met'
							snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary_id, \
															which_coding_exon, pos_within_codon, comment, SNP_index_in_peptide))
					"""
					except:
						traceback.print_exc()
						sys.stderr.write('%s.\n'%repr(sys.exc_info()))
						sys.stderr.write("Except encountered for SNP (%s, %s), gene %s (chr=%s, %s-%s), gene_commentary_id %s (%s-%s), box_ls=%s.\n"%\
							(chr, pos, gene_id, locus.chromosome, oneGeneData.start, oneGeneData.stop, \
							gene_commentary_id, geneCommentaryRBDict.start, geneCommentaryRBDict.stop, repr(box_ls)))
					"""
		else:
			if oneGeneData.type_of_gene=='pseudo':
				snp_annotation_type_short_name = oneGeneData.type_of_gene
			else:
				snp_annotation_type_short_name = geneCommentaryRBDict.gene_commentary_type_name
			snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary_id, \
													None, None, None))
		#else:	#integenic
		#	snp_annotation_type_short_name_ls.append(['intergenic'])
			
		#now save everything into db
		locus_annotation_ls = []
		for snp_annotation_type_tup in snp_annotation_type_short_name_ls:
			snp_annotation_type_short_name = snp_annotation_type_tup[0]
			if snp_annotation_type_short_name not in locus_annotation_short_name2db_entry:
				ty = db_vervet.getLocusAnnotationType(locus_annotation_type_short_name=snp_annotation_type_short_name)
				locus_annotation_short_name2db_entry[snp_annotation_type_short_name] = ty
			if len(snp_annotation_type_tup)>=3:
				gene_id = snp_annotation_type_tup[1]
				gene_commentary_id = snp_annotation_type_tup[2]
			else:
				gene_id = None
				gene_commentary_id = None
			if len(snp_annotation_type_tup)>=4:
				which_exon_or_intron = snp_annotation_type_tup[3]
			else:
				which_exon_or_intron = None
			if len(snp_annotation_type_tup)>=5:
				pos_within_codon = snp_annotation_type_tup[4]
			else:
				pos_within_codon = None
			if len(snp_annotation_type_tup)>=6:
				comment = snp_annotation_type_tup[5]
			else:
				comment = None
			if len(snp_annotation_type_tup)>=7:
				which_codon = snp_annotation_type_tup[6] +1	#[6] is the SNP_index_in_peptide
			else:
				which_codon = None
			locus_annotation_type = locus_annotation_short_name2db_entry.get(snp_annotation_type_short_name)
			
			locus_annotation = db_vervet.getLocusAnnotation(locus_id=locus.id, locus_context_id=locus_context.id, \
													locus_context=locus_context, gene_id=gene_id, \
													gene_commentary_id=gene_commentary_id, \
													gene_segment_id=gene_segment_id, \
													locus_annotation_type=locus_annotation_type, locus_annotation_type_id=locus_annotation_type.id,\
													which_exon_or_intron=which_exon_or_intron, pos_within_codon=pos_within_codon, \
													which_codon=which_codon, label=geneSegmentKey.label, \
													utr_number=geneSegmentKey.utr_number, cds_number=geneSegmentKey.cds_number, \
													intron_number=geneSegmentKey.intron_number,\
													exon_number=geneSegmentKey.exon_number, overlap_length=None, \
													overlap_fraction_in_locus=None, overlap_fraction_in_gene=None,\
													comment=comment)
			if locus_annotation:
				param_obj.no_of_locus_annotations_already_in_db += 1
			else:
				param_obj.no_of_into_db += 1
			param_obj.no_of_total_annotations += 1
			locus_annotation_ls.append(locus_annotation)
			real_counter += 1
		"""
		if self.report and counter%2000==0:
			sys.stderr.write("%s%s\t%s"%('\x08'*40, counter, real_counter))
			session.flush()
		if self.report:
			sys.stderr.write("%s%s\t%s\n"%('\x08'*40, counter, real_counter))
		"""
		sys.stderr.write("Done.\n")
		return locus_annotation_ls
	
	
	def findLocusContext(self, db_vervet=None, genomeRBDict=None, locus_type_id=None, compareIns=None, max_distance=50000, debug=0,
					param_obj=None, locus_annotation_short_name2db_entry=None):
		"""
		2012.5.14
		2011-3-25
			cast row.chromosome (from db) into str type.
		2010-10-3
			bug fixed: (chr, start, stop) is not unique. There are genes with the same coordinates.
		2010-8-18
		"""
		sys.stderr.write("Finding Locus context ... \n")
		session = db_vervet.session
		TableClass = VervetDB.Locus
		query = TableClass.query.filter_by(locus_type_id=locus_type_id)
		for row in query:
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=str(row.chromosome), \
							span_ls=[row.start, row.stop], \
							min_reciprocal_overlap=0.0000001, )	#min_reciprocal_overlap doesn't matter here.
				# it's decided by compareIns.
			node_ls = []
			genomeRBDict.findNodes(segmentKey, node_ls=node_ls, compareIns=compareIns)
			for node in node_ls:
				geneSegKey = node.key
				for oneGeneData in node.value:
					# geneSegKey.span_ls expands 20kb upstream or downstream of the gene.
					overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
															[oneGeneData.gene_start, oneGeneData.gene_stop])[:5]
					if overlap_length>0:	#use fraction of length as coordinates.
						gene_length = oneGeneData.gene_stop - oneGeneData.gene_start + 1
						try:
							if oneGeneData.strand == '+1':
								term5_disp_pos = abs(overlap_start_pos - oneGeneData.gene_start)/float(gene_length)
								term3_disp_pos = abs(overlap_stop_pos - oneGeneData.gene_start + 1)/float(gene_length)
							else:
								term5_disp_pos = abs(oneGeneData.gene_stop - overlap_stop_pos)/float(gene_length)
								term3_disp_pos = abs(oneGeneData.gene_stop - overlap_start_pos + 1)/float(gene_length)
						except:
							import pdb
							pdb.set_trace()
					else:	#no overlap at the gene itself, but within the neighborhood (<=max_distance)
						term3_disp_pos = None
						if oneGeneData.strand == '+1':
							if row.stop<=oneGeneData.gene_start:	#upstream
								term5_disp_pos = row.stop - oneGeneData.gene_start
							elif row.start>=oneGeneData.gene_stop:	# downstream
								term5_disp_pos = row.start - oneGeneData.gene_stop
						else:
							if row.stop<=oneGeneData.gene_start:	#downstream
								term5_disp_pos = oneGeneData.gene_start - row.stop
							elif row.start>=oneGeneData.gene_stop:	# upstream
								term5_disp_pos = oneGeneData.gene_stop - row.start
					locus_context = db_vervet.getLocusContext(locus_id=row.id, gene_id=oneGeneData.gene_id, \
											gene_strand=oneGeneData.strand, disp_pos=term5_disp_pos, \
											overlap_length=overlap_length,\
											overlap_fraction_in_locus=overlap1, overlap_fraction_in_gene=overlap2)
					
					if locus_context.id:
						param_obj.no_of_locus_contexts_already_in_db += 1
					else:
						param_obj.no_of_into_db += 1
					param_obj.no_of_total_contexts += 1
					
					for geneCommentaryRBDict in oneGeneData.geneCommentaryRBDictLs:
						gene_box_node_ls = []
						geneCommentaryRBDict.findNodes(segmentKey, node_ls=gene_box_node_ls, compareIns=compareIns)
						for geneSegmentNode in gene_box_node_ls:
							geneSegmentKey = geneSegmentNode.key
							if row.stop!=row.start: 	#2012.5.14 structural variation, not single-base locus 
								overlap1, overlap2, overlap_length, overlap_start_pos, overlap_stop_pos = get_overlap_ratio(segmentKey.span_ls, \
														geneSegmentKey.span_ls)[:5]
								locus_annotation = db_vervet.getLocusAnnotation(locus_id=row.id, locus_context_id=locus_context.id, locus_context=locus_context, \
														gene_id=oneGeneData.gene_id,
														gene_commentary_id=geneCommentaryRBDict.gene_commentary_id, \
														gene_segment_id=geneSegmentKey.gene_segment_id, \
														locus_annotation_type =None, locus_annotation_type_id=None,\
														which_exon_or_intron = None, pos_within_codon=None, \
														label=geneSegmentKey.label, \
														utr_number=geneSegmentKey.utr_number, cds_number=geneSegmentKey.cds_number, \
														intron_number=geneSegmentKey.intron_number,\
														exon_number=geneSegmentKey.exon_number, overlap_length=overlap_length, \
														overlap_fraction_in_locus=overlap1, overlap_fraction_in_gene=overlap2,\
														comment=None)
								if locus_annotation:
									param_obj.no_of_locus_annotations_already_in_db += 1
								else:
									param_obj.no_of_into_db += 1
								param_obj.no_of_total_annotations += 1
							else:	#single-nucleotide 
								self._constructSNPAnnotation(db_vervet, locus=row, locus_context=locus_context, \
													oneGeneData=oneGeneData, \
													geneCommentaryRBDict=geneCommentaryRBDict, geneSegmentKey=geneSegmentKey,\
													locus_annotation_short_name2db_entry=locus_annotation_short_name2db_entry, param_obj=param_obj)
							
					if param_obj.no_of_into_db>2000:
						session.flush()
						param_obj.no_of_into_db = 0
						sys.stderr.write("\t %s/%s LocusContext(s) & %s/%s LocusAnnotation(s) already in db.\n"%(\
									param_obj.no_of_locus_contexts_already_in_db, param_obj.no_of_total_contexts, \
									param_obj.no_of_locus_annotations_already_in_db, param_obj.no_of_total_annotations))
		session.flush()
		session.expunge_all()
		sys.stderr.write("\t %s/%s LocusContext(s) & %s/%s LocusAnnotation(s) already in db.\n"%(\
								param_obj.no_of_locus_contexts_already_in_db, param_obj.no_of_total_contexts, \
								param_obj.no_of_locus_annotations_already_in_db, param_obj.no_of_total_annotations))
	
		
	def connectDB(self):
		"""
		2012.5.14
		2012.4.29
			split out of __init__() so that derived classes could overwrite this function
		"""
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, \
									hostname=self.hostname, database=self.dbname, schema=self.schema, port=self.port)
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, \
						schema=self.genome_db_schema, port=self.port)
		db_genome.setup(create_tables=False)	#do the setup() after both db have been instantiated.
		self.db_genome = db_genome
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
	
	def run(self):
		"""
		2012.5.14
			
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
			
		
		#self.db_genome.session.begin()
		#self.db_genome.metadata.bind.execute('set search_path to %s'%(self.genome_db_schema))
		
		genomeRBDict = self.db_genome.dealWithGenomeRBDict(self.genomeRBDictPickleFname, tax_id=self.tax_id, \
													max_distance=self.max_distance, debug=self.debug)
		
		self.db_vervet.session.begin()
		param_obj = PassingData(no_of_total_annotations=0, session=self.db_vervet.session, \
					locus_type_id=self.locus_type_id, no_of_total_contexts=0, no_of_into_db=0, report=self.report,\
					no_of_locus_contexts_already_in_db=0, no_of_locus_annotations_already_in_db=0)
		compareIns = CNVCompare(min_reciprocal_overlap=0.0000001)	#any overlap is an overlap
		
		locus_annotation_short_name2db_entry = self.getLocusAnnotationShortName2dbEntry(self.db_vervet)
		
		self.findLocusContext(self.db_vervet, genomeRBDict, locus_type_id=self.locus_type_id, compareIns=compareIns, \
					max_distance=self.max_distance, debug=self.debug, param_obj=param_obj, \
					locus_annotation_short_name2db_entry=locus_annotation_short_name2db_entry)
		if self.commit:
			self.db_vervet.session.flush()
			self.db_vervet.session.commit()
		else:
			self.db_vervet.session.rollback()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = FindLocusContext
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
