#!/usr/bin/env python
"""
Examples:
	%s -z localhost -u yh -o /tmp/alignmentDepth.tsv

	%s 
	
	%s 
Description:
	2013.06.13 this program outputs alignment depth, part of workflows such as FilterVCFPipeline.py
		If inputFname (VCF file) is given, it takes a VCF file as input 
			and matches the sample IDs with IndividualAlignment entries inside db.
		Otherwise, it outputs all alignments (--ref_ind_seq_id).
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import copy
from pymodule import ProcessOptions
from pymodule import MatrixFile
from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
parentClass = AbstractVervetMapper

class OutputVCFAlignmentDepthRange(parentClass):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						('depthFoldChange', 0, int): [2, '', 1, 'fold-change of median depth, determining minDepth & maxDepth for outputFileFormat 2'],\
						("minGQ", 0, int): [30, '', 1, 'minimum genotype quality. for outputFileFormat 2'],\
						('outputFileFormat', 0, int): [1, '', 1, 'tab-delimited. 1: read_group, median_depth, individual.id; \
	2: read_group, minDepth, maxDepth, minGQ.'],\
						("ref_ind_seq_id", 0, int): [None, '', 1, 'ID of the reference sequence used in alignment. used to filter alignments.'],\
						
						("reduce_reads", 0, int): [None, '', 1, 'To filter which input alignments to fetch from db'],\
						('excludeContaminant', 0, int):[0, '', 0, 'toggle this to exclude alignments or sequences that are from contaminated individuals, \n\
		(IndividualSequence.is_contaminated=1)'],\
						("sequence_filtered", 0, int): [None, 'Q', 1, 'to filter alignments/individual_sequences. None: whatever; 0: unfiltered sequences, 1: filtered sequences: 2: ...'],\
						('completedAlignment', 0, int):[None, '', 1, 'a flag requiring whether user chooses alignment that has been completed or not.\n\
	--completedAlignment 0 is same as --skipDoneAlignment. --completedAlignment 1 gets you only the alignments that has been completed. Default (None) has no effect.'],\
						('alignment_outdated_index', 0, int): [0, '', 1, 'filter based on value of IndividualAlignment.outdated_index.', ],\
						("alignment_method_id", 0, int): [None, 'G', 1, 'To filter alignments. None: whatever; integer: AlignmentMethod.id'],\
						("local_realigned", 0, int): [None, '', 1, 'To filter which input alignments to fetch from db (i.e. AlignmentReadBaseQualityRecalibrationWorkflow.py)\
	OR to instruct whether local_realigned should be applied (i.e. ShortRead2AlignmentWorkflow.py)'],\
						
						})
	#add these arguments to filter alignments. to make them unique (one individual has one alignment in the output)
	#sequence_filtered 1 --local_realigned 0 --reduce_reads 0 --completedAlignment 1 --excludeContaminant --alignment_method_id
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		parentClass.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
	def outputAlignmentDepthAndOthersForFilter(self, db_vervet=None, inputFname=None, \
						ref_ind_seq_id=None, depthFoldChange=2, minGQ=30, \
						outputFname=None, outputFileFormat=1):
		"""
		2012.6.12
			added argument db_vervet, moved from FilterVCFPipeline.py
		2011-9-2
		"""
		sys.stderr.write("Outputting alignment (from %s) coverage to %s ..."%(inputFname, outputFname))
		if inputFname:
			alignmentLs = db_vervet.getAlignmentsFromVCFFile(inputFname=inputFname)
		else:
			alignmentLs = db_vervet.getAlignments(ref_ind_seq_id=self.ref_ind_seq_id, \
										alignment_method_id=self.alignment_method_id, data_dir=self.data_dir,\
										local_realigned=self.local_realigned, outdated_index=self.alignment_outdated_index,\
										completedAlignment=self.completedAlignment, \
										reduce_reads=self.reduce_reads)
			"""
			TableClass = VervetDB.IndividualAlignment
			query = TableClass.query.filter(TableClass.median_depth!=None)
			if ref_ind_seq_id:
				query = query.filter(TableClass.ref_ind_seq_id==ref_ind_seq_id)
			alignmentLs = query.order_by(TableClass.id)
			"""
			
		alignmentLs = db_vervet.filterAlignments(data_dir=self.data_dir, alignmentLs=alignmentLs, sequence_filtered=self.sequence_filtered, \
						mask_genotype_method_id=None, parent_individual_alignment_id=None,\
						excludeContaminant=self.excludeContaminant,local_realigned=self.local_realigned,\
						reduce_reads=self.reduce_reads,\
						completedAlignment=self.completedAlignment,\
						alignment_method_id=self.alignment_method_id, \
						outdated_index=self.alignment_outdated_index)
		writer = MatrixFile(inputFname=outputFname, openMode='w', delimiter='\t')
		if outputFileFormat==1:
			header = ['alignmentID', 'medianDepth', "individualID"]
		else:
			header = ['alignmentID', 'minDepth', 'maxDepth', 'minGQ']
		writer.writeHeader(header)
		
		counter = 0
		for row in alignmentLs:
			read_group = row.read_group
			if outputFileFormat==1:
				data_row = [read_group, row.median_depth, row.individual_sequence.individual.id]
			else:
				minDepth = row.median_depth/float(depthFoldChange)
				if abs(minDepth-0)<=0.001:	#if it's too close to 0, assign 0.
					minDepth = 0
				data_row = [read_group, minDepth, row.median_depth*float(depthFoldChange), minGQ]
			writer.writerow(data_row)
			counter += 1
		writer.close()
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
	def run(self):
		"""
		2012.7.13
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		session = self.db_vervet.session
		session.begin()
		
		self.outputAlignmentDepthAndOthersForFilter(db_vervet=self.db_vervet, \
							inputFname=self.inputFname, ref_ind_seq_id=self.ref_ind_seq_id, \
							depthFoldChange=self.depthFoldChange, minGQ=self.minGQ, \
							outputFname=self.outputFname, outputFileFormat=self.outputFileFormat)
	
if __name__ == '__main__':
	main_class = OutputVCFAlignmentDepthRange
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()