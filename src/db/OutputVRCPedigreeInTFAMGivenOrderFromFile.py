#!/usr/bin/env python
"""
Examples:
	%s 
	
	%s -i AlignmentToCall_AllVRC_vs_524_AllContig195Sites.2012.8.5T1549/samtools/Contig195.vcf.gz
		-o /tmp/723VRC.tfam -u yh

Description:
	2012.8.9
		outputs VRC pedigrees (only the monkeys in inputFname/vcf format) in plink tfam format.
		parse the order of monkeys from VCF file (inputFname) and preserve it in TFAM file.
	six columns:
	 Family ID
	 Individual ID
	 Paternal ID
	 Maternal ID
	 Sex (1=male; 2=female; other=unknown)
	 Phenotype
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, NextGenSeq, Genome
from pymodule.VCFFile import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB

class OutputVRCPedigreeInTFAMGivenOrderFromFile(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						})
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def getIndividual(self, db_vervet=None, individual_id=None, individual_id2individual=None):
		"""
		2012.8.14
		"""
		if individual_id2individual:
			individual = individual_id2individual.get(individual_id)
			if individual:
				return individual
			else:
				individual = VervetDB.Individual.get(individual_id)
				individual_id2individual[individual_id] = individual
				return individual
		else:
			individual = VervetDB.Individual.get(individual_id)
			return individual
		
		
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		db_vervet = self.db_vervet
		session.begin()
		
		vcfFile = VCFFile(inputFname=self.inputFname)
		alignmentLs = []
		for read_group in vcfFile.getSampleIDList():	#sample_id_ls is not good because its index 0 is "ref"
			individualAlignment = db_vervet.parseAlignmentReadGroup(read_group).individualAlignment
			alignmentLs.append(individualAlignment)
		"""
		pedigreeGraphData = db_vervet.constructPedgreeGraphOutOfAlignments(alignmentLs)
		DG = pedigreeGraphData.DG
		individual_id2alignmentLs = pedigreeGraphData.individual_id2alignmentLs
		"""
		DG = db_vervet.constructPedgree()
		individual_id2individual = {}
		
		writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
		counter = 0
		family_id= 1	#all in one family
		for alignment in alignmentLs:
			node_id = alignment.ind_sequence.individual_id
			individual = self.getIndividual(db_vervet=db_vervet, individual_id=node_id, \
										individual_id2individual=individual_id2individual)
			if node_id in DG:
				parents = DG.predecessors(node_id)
				if len(parents)==2:
					parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
												individual_id2individual=individual_id2individual)
					parent2 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[1], \
												individual_id2individual=individual_id2individual)
					parent1Sex = parent1.codeSexInNumber()
					if parent1Sex==2:
						#swap the father and mother row
						tmp = parent1
						parent1 = parent2
						parent2 = tmp
					
					father_id = parent1.ucla_id
					mother_id = parent2.ucla_id
				elif len(parents)==1:
					parent1 = self.getIndividual(db_vervet=db_vervet, individual_id=parents[0], \
										individual_id2individual=individual_id2individual)
					parent1Sex = parent1.codeSexInNumber()
					if parent1Sex==2:
						father_id = 0
						mother_id = parent1.ucla_id
					else:
						father_id = parent1.ucla_id
						mother_id = 0
				elif len(parents)==0:
					father_id = 0
					mother_id = 0
				else:
					sys.stderr.write("Error: number of parents (%s) for %s is %s.\n"%(repr(parents), node_id, len(parents)))
					sys.exit(3)
			else:
				father_id = 0
				mother_id = 0
			data_row = [family_id, individual.ucla_id, father_id, mother_id, \
					individual.codeSexInNumber(), 1]
			writer.writerow(data_row)
			counter += 1
		sys.stderr.write("%s alignments outputted.\n"%(counter))
		del writer
		


if __name__ == '__main__':
	main_class = OutputVRCPedigreeInTFAMGivenOrderFromFile
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()