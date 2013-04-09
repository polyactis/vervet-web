#!/usr/bin/env python
"""
Examples:
	#for plink
	%s -i AlignmentToCall_AllVRC_vs_524_AllContig195Sites.2012.8.5T1549/samtools/Contig195.vcf.gz
		-o /tmp/723VRC.tfam -u yh

	#2013.1.3 for TrioCaller
	%s -i /Network/Data/vervet/db/genotype_file/method_25/25999_VCF_24115_VCF_Contig520.filterByMaxSNPMissingRate.recode.vcf.gz
		--outputFileFormat 2 -u yh -o /tmp/trioCaller_VRC.merlin --sampleID2FamilyCountFname=/tmp/sampleID2familyCount_y2.tsv
	
	#2013.1.3 output for polymutt
	%s -i /Network/Data/vervet/db/genotype_file/method_25/25999_VCF_24115_VCF_Contig520.filterByMaxSNPMissingRate.recode.vcf.gz
		--outputFileFormat 3 -u yh -o /tmp/polymutt_VRC.merlin --polymuttDatFname /tmp/datfile
		--sampleID2FamilyCountFname=/tmp/sampleID2familyCount_y3.tsv --db_passwd secret
	
Description:
	2013.04.08 this program takes a VCF file as input,
		matches the sample IDs with IndividualAlignment entries inside db,
		outputs minDepth & maxDepth & GQ of that sample, to be used by workflows such as FilterVCFPipeline.py
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import heapq, copy
from pymodule import ProcessOptions, getListOutOfStr, PassingData, utils, figureOutDelimiter, NextGenSeq, Genome
from pymodule import VCFFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB

class OutputVCFAlignmentDepthRange(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						('depthFoldChange', 0, int): [2, '', 1, 'fold-change of median depth, determining minDepth & maxDepth in output'],\
						("minGQ", 0, int): [30, '', 0, 'minimum genotype quality to output'],\
						})
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def outputAlignmentDepthAndOthersForFilter(self, db_vervet=None, outputFname=None, ref_ind_seq_id=524, \
											foldChange=2, minGQ=30):
		"""
		2012.6.12
			added argument db_vervet, moved from FilterVCFPipeline.py
		2011-9-2
		"""
		sys.stderr.write("Outputting sequence coverage to %s ..."%outputFname)
		import csv
		counter = 0
		writer = csv.writer(open(outputFname, 'w'), delimiter='\t')
		writer.writerow(['aln.id', 'minDepth', 'maxDepth', 'minGQ'])
		TableClass = VervetDB.IndividualAlignment
		query = TableClass.query.filter(TableClass.median_depth!=None)
		if ref_ind_seq_id:
			query = query.filter(TableClass.ref_ind_seq_id==ref_ind_seq_id)
		query = query.order_by(TableClass.id)
		for row in query:
			minDepth = row.median_depth/float(foldChange)
			if abs(minDepth-0)<=0.001:	#if it's too close to 0, regard it as 0.
				minDepth = 0
			read_group
			writer.writerow([read_group, minDepth, \
							row.median_depth*float(foldChange), minGQ])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
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
		
	
if __name__ == '__main__':
	main_class = OutputVCFAlignmentDepthRange
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()