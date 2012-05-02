#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -u yh -c -z uclaOffice -k genome -t 60711 -s 9 -i ../data/IRF7/Contig103_protein_msa_cufflinksIntronHints_g129_IRF7.gff3

Description:
	2012.4.25
		Put gene models predicted by augustus into database. 
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
from pymodule import GenomeDB


class AllGenes(object):
	def __init__(self, tax_id=None, sequence_type_id=None):
		"""
		"""
		self.tax_id = tax_id
		self.sequence_type_id = sequence_type_id
		
		self.geneID2GeneModel = {}
		self.transcriptID2geneID = {}
	
	def getNumberOfGenes(self):
		return len(self.geneID2GeneModel)
	
	def getNumberOfTranscripts(self):
		return len(self.transcriptID2geneID)
	
	def addOneGene(self, chromosome=None, start=None, stop=None, strand=None, geneID=None):
		"""
		2012.4.25
		"""
		if geneID not in self.geneID2GeneModel:
			geneModel = GeneModel(chromosome=chromosome, start=start, stop=stop, strand=strand, geneID=geneID)
			self.geneID2GeneModel[geneID] = geneModel
		else:
			sys.stderr.write("Error: geneID %s already in geneID2GeneModel.\n"%(geneID))
			sys.exit(3)
	
	def getGeneByGeneID(self, geneID):
		"""
		2012.4.25
		"""
		return self.geneID2GeneModel.get(geneID)
	
	def addTranscript(self, geneID=None, transcriptID=None, chromosome=None, start=None, stop=None, strand=None):
		"""
		2012.4.25
		"""
		if geneID not in self.geneID2GeneModel:
			self.geneID2GeneModel[geneID] = GeneModel(chromosome=chromosome, start=start, stop=stop, strand=strand, geneID=geneID)
		geneModel = self.geneID2GeneModel[geneID]
		geneModel.addTranscript(transcriptID=transcriptID, start=start, stop=stop, strand=strand)
		if transcriptID in self.transcriptID2geneID:
			sys.stderr.write("Warning: transcriptID %s already in transcriptID2geneID.\n"%(transcriptID))
		else:
			self.transcriptID2geneID[transcriptID] = geneID
	
	def getGeneByTranscriptID(self, transcriptID):
		geneToReturn = None
		geneID = self.transcriptID2geneID.get(transcriptID)
		if geneID and geneID in self.geneID2GeneModel:
			geneToReturn = self.geneID2GeneModel.get(geneID)
		return geneToReturn
	
	def getGeneModelLs(self):
		"""
		2012.4.26
		"""
		return self.geneID2GeneModel.values()

class GeneModel(object):
	def __init__(self, chromosome=None, start=None, stop=None, strand=None, geneID=None):
		"""
		2012.4.25
			to hold a geneModel structure from GFF3 file
		"""
		self.chromosome = chromosome
		self.start = start
		self.stop = stop
		self.strand = strand
		self.geneID = geneID
		self.transcriptID2Transcript = {}
	
	def addTranscript(self, transcriptID=None, start=None, stop=None, strand=None):
		"""
		2012.4.25
		"""
		if self.start is None or self.start>start:
			self.start = start
		if self.stop is None or self.stop<stop:
			self.stop = stop
		if self.strand is None:
			self.strand = strand
		if transcriptID not in self.transcriptID2Transcript:
			self.transcriptID2Transcript[transcriptID] = Transcript(transcriptID=transcriptID, \
													start=start, stop=stop, strand=strand)
	
	def addExon(self, transcriptID=None, start=None, stop=None, strand=None):
		"""
		2012.4.25
		"""
		if transcriptID not in self.transcriptID2Transcript:
			self.transcriptID2Transcript[transcriptID] = Transcript(transcriptID=transcriptID, start=start,\
													stop=stop, strand=strand)
		transcript = self.transcriptID2Transcript.get(transcriptID)
		transcript.addExon(start=start, stop=stop, strand=strand)
	
	def addCDS(self, transcriptID=None, start=None, stop=None, strand=None):
		"""
		2012.4.25
		"""
		if transcriptID not in self.transcriptID2Transcript:
			self.transcriptID2Transcript[transcriptID] = Transcript(transcriptID=transcriptID, start=start,\
													stop=stop, strand=strand)
		transcript = self.transcriptID2Transcript.get(transcriptID)
		transcript.addCDS(start=start, stop=stop, strand=strand)
	
	def getTranscriptLs(self):
		"""
		2012.4.25
		"""
		return self.transcriptID2Transcript.values()

class Transcript(object):
	def __init__(self, transcriptID=None, start=None, stop=None, strand=None):
		self.transcriptID = transcriptID
		self.start = start
		self.stop = stop
		self.strand = strand
		self.exon_ls = []
		self.protein = None
	
	def addCDS(self, start=None, stop=None, strand=None):
		"""
		"""
		if self.protein is None:
			self.protein = Protein(start=start, stop=stop, strand=strand)
		self.protein.addCDS(start=start, stop=stop, strand=strand)
		
		
	def addExon(self, start=None, stop=None, strand=None):
		"""
		"""
		if self.start is None or self.start>start:
			self.start = start
		if self.stop is None or self.stop<stop:
			self.stop = stop
		if self.strand is None:
			self.strand = strand
		self.exon_ls.append([start, stop])
	
	def getExonLs(self):
		return self.exon_ls
	
	def getProtein(self):
		return self.protein
	
	def getCDSLs(self):
		if self.protein:
			return self.protein.CDS_ls
		else:
			return None
	
class Protein(object):
	def __init__(self, start=None, stop=None, strand=None):
		self.start = start
		self.stop = stop
		self.strand = strand
		self.CDS_ls = []
	
	def addCDS(self, start=None, stop=None, strand=None):
		"""
		2012.4.25
		"""
		if self.start is None or self.start>start:
			self.start = start
		if self.stop is None or self.stop<stop:
			self.stop = stop
		if self.strand is None:
			self.strand = strand
		self.CDS_ls.append([start, stop])
	
	def getCDSLs(self):
		return self.CDS_ls

class PutAugustusGenePrediction2DB(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 1, ))
	option_default_dict.pop(('outputFname', 0, ))
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
							('tax_id', 1, int):[60711, 't', 1, 'taxonomy ID for the reference species'],\
							('sequence_type_id', 1, int):[9, 's', 1, 'sequence type ID of the reference sequence (scaffold, BAC, chromosome?)'],\
							})
	def __init__(self, inputFnameLs, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs, **keywords)
	
	def parseGFF3Attributes(self, attributes):
		"""
		2012.4.25
		
		"""
		attributeLs = attributes.split(';')
		attributeName2Value = {}
		for attribute in attributeLs:
			if attribute:
				name, value = attribute.split('=')
				attributeName2Value[name] = value
		return attributeName2Value
	
	def constructGeneCommentaryFromProtein(self, db_genome, gene=None, transcriptGeneCommentary=None, protein=None, \
											gene_commentary_type_id=8):
		"""
		2012.4.26
			gene is an instance of GenomeDB.Gene, not class GeneModel in this file.
			gene_commentary_type_id 3 is mRNA. 8 is peptide.
		"""
		gene_commentary = db_genome.getGeneCommentary(start=protein.start, stop=protein.stop, \
						gene_commentary_type_id=gene_commentary_type_id, gene_commentary_id=None,\
						accession=None, version=None, gi=None, label=None, text=None, comment=None,
						gene=gene, gene_commentary=transcriptGeneCommentary)
		return gene_commentary
	
	def constructGeneCommentaryFromTranscript(self, db_genome, gene=None, transcript=None, gene_commentary_type_id=3):
		"""
		2012.4.26
			gene is an instance of GenomeDB.Gene, not class GeneModel in this file.
			gene_commentary_type_id 3 is mRNA. 8 is peptide.
		"""
		transcriptGeneCommentary = db_genome.getGeneCommentary(start=transcript.start, stop=transcript.stop, \
						gene_commentary_type_id=gene_commentary_type_id, gene_commentary_id=None,\
						accession=None, version=None, gi=None, label=None, text=None, comment=None,
						gene=gene)
		transcriptGeneCommentary.gene_segments = db_genome.constructGeneSegments(box_ls=transcript.getExonLs(), \
																			gene_commentary=transcriptGeneCommentary)
		
		protein_gene_commentary = self.constructGeneCommentaryFromProtein(db_genome, gene=gene, \
								transcriptGeneCommentary=transcriptGeneCommentary, protein=transcript.protein, gene_commentary_type_id=8)
		protein_gene_commentary.gene_segments = db_genome.constructGeneSegments(box_ls=transcript.getCDSLs(), \
																		gene_commentary=protein_gene_commentary)
		transcriptGeneCommentary.gene_commentaries.append(protein_gene_commentary)
		db_genome.session.add(transcriptGeneCommentary)
		db_genome.session.flush()
		return transcriptGeneCommentary
		
	
	def parseAugustusGFF3File(self, db_genome=None, inputFname=None, tax_id=None, sequence_type_id=None):
		"""
		2012.4.25

		Augustus commandline:
		# /home/crocea/bin/augustus/bin/augustus --species=human --proteinprofile=IRF7_multi.prfl --UTR=on --alternatives-from-evidence=true 
			--hintsfile=hintsFromCufflinks.gff --gff3=on --extrinsicCfgFile=~/bin/augustus/config/extrinsic/extrinsic.E.cfg Contig103.fa

			Each bock of augustus output looks like this:
			
# start gene g129
Contig103	AUGUSTUS        gene    4837992 4841222 1       +       .       ID=g129
Contig103       AUGUSTUS        transcript      4837992 4841222 .       +       .       ID=g129.t1;Parent=g129
Contig103       AUGUSTUS        transcription_start_site        4837992 4837992 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        exon    4837992 4838133 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        exon    4838386 4838723 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        start_codon     4838704 4838706 .       +       0       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4838704 4838723 .       +       0       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        CDS     4838809 4838971 .       +       1       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4838809 4838971 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4839055 4839268 .       +       0       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4839055 4839268 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4839513 4839571 .       +       2       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4839513 4839571 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4839648 4839873 .       +       0       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4839648 4839873 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4840022 4840108 .       +       2       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4840022 4840108 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4840194 4840274 .       +       2       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4840194 4840274 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4840347 4840736 .       +       2       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4840347 4840736 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4840825 4840943 .       +       2       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4840825 4840943 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        CDS     4841029 4841184 .       +       0       ID=g129.t1.cds;Parent=g129.t1
Contig103       AUGUSTUS        exon    4841029 4841222 .       +       .       Parent=g129.t1
Contig103       AUGUSTUS        stop_codon      4841182 4841184 .       +       0       Parent=g129.t1
Contig103       AUGUSTUS        transcription_end_site  4841222 4841222 .       +       .       Parent=g129.t1
....
# end gene g129
###

		"""
		sys.stderr.write("Parsing %s ..."%(inputFname))
		reader = csv.reader(open(inputFname), delimiter=figureOutDelimiter(inputFname))
		counter = 0
		allGenes = AllGenes(tax_id=tax_id, sequence_type_id=sequence_type_id)
		
		for row in reader:
			if row[0][0]=='#':	#skip all comments
				continue
			try:
				chromosome, source, span_type, start, stop, score, strand, phase, attributes = row
			except:
				sys.stderr.write("Except while reading this row: %s.\n"%(repr(row)))
				sys.exit(4)
			start = int(start)
			stop = int(stop)
			try:
				attributeName2Value = self.parseGFF3Attributes(attributes)
			except:
				sys.stderr.write("Except while parsing attributes of this row: %s.\n"%(repr(row)))
				sys.exit(4)
			if span_type=='gene':
				geneID = attributeName2Value.get("ID")
				allGenes.addOneGene(chromosome=chromosome, start=start, stop=stop, strand=strand, geneID=geneID)
			elif span_type == "transcript":
				transcriptID = attributeName2Value.get("ID")
				parentGeneID = attributeName2Value.get("Parent")
				allGenes.addTranscript(geneID=parentGeneID, transcriptID=transcriptID, chromosome=chromosome, \
									start=start, stop=stop, strand=strand)
			elif span_type == 'exon':
				transcriptID = attributeName2Value.get("Parent")
				geneModel = allGenes.getGeneByTranscriptID(transcriptID)
				geneModel.addExon(transcriptID=transcriptID, start=start, stop=stop, strand=strand)
				
			elif span_type =='CDS':
				transcriptID = attributeName2Value.get("Parent")
				geneModel = allGenes.getGeneByTranscriptID(transcriptID)
				geneModel.addCDS(transcriptID=transcriptID, start=start, stop=stop, strand=strand)
		sys.stderr.write("%s genes, %s transcripts.\n"%(allGenes.getNumberOfGenes(), allGenes.getNumberOfTranscripts()))
		
		for geneModel in allGenes.getGeneModelLs():
			annot_assembly = db_genome.getAnnotAssemblyFromTaxIDSequenceTypeChromosome(tax_id=tax_id, sequence_type_id=sequence_type_id, \
														chromosome=geneModel.chromosome)
			
			gene = db_genome.getGene(annot_assembly_id=annot_assembly.id, locustag=geneModel.geneID, tax_id=tax_id, \
									chromosome=geneModel.chromosome, strand="%s1"%(geneModel.strand),\
									start=geneModel.start, stop=geneModel.stop, type_of_gene="protein_coding", entrezgene_type_id=1006, \
									chromosome_sequence_type_id=sequence_type_id,
									annot_assembly=annot_assembly)
			for transcript in geneModel.getTranscriptLs():
				transcriptGeneCommentary = self.constructGeneCommentaryFromTranscript(db_genome, gene=gene, transcript=transcript, \
																					gene_commentary_type_id=3)
	def connectDB(self):
		"""
		2012.4.29
			overwrite the ancestral function
		"""
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_genome.setup(create_tables=False)	#2010-6-22
		self.session = db_genome.session
		self.db_genome = db_genome
		
	
	def run(self):
		"""
		2012.4.25
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		if not self.debug:	#in debug mode, no transaction, auto-commit
			self.session.begin()
		
		self.parseAugustusGFF3File(self.db_genome, self.inputFname, tax_id=self.tax_id, sequence_type_id=self.sequence_type_id)
		if not self.debug:
			if self.commit:
				self.session.flush()
				self.session.commit()
			else:
				self.session.rollback()
		

if __name__ == '__main__':
	main_class = PutAugustusGenePrediction2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()