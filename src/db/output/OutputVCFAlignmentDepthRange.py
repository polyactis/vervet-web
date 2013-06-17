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

class OutputVCFAlignmentDepthRange(AbstractVervetMapper):
	__doc__ = __doc__
	option_default_dict = AbstractVervetMapper.option_default_dict.copy()
	option_default_dict.pop(('outputFnamePrefix', 0, ))
	option_default_dict.update({
						('outputFname', 1, ):option_default_dict.get(('outputFname', 0, )),\
						('depthFoldChange', 0, int): [2, '', 1, 'fold-change of median depth, determining minDepth & maxDepth for outputFileFormat 2'],\
						("minGQ", 0, int): [30, '', 1, 'minimum genotype quality. for outputFileFormat 2'],\
						('outputFileFormat', 0, int): [1, '', 1, 'tab-delimited. 1: read_group, median_depth; 2: read_group, minDepth, maxDepth, minGQ.'],\
						("ref_ind_seq_id", 0, int): [None, '', 1, 'ID of the reference sequence used in alignment. used to filter alignments.'],\
						})
	option_default_dict.pop(('outputFname', 0, ))	#pop after its value has been used above
	def __init__(self, inputFnameLs=None, **keywords):
		"""
		"""
		AbstractVervetMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
	
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
			TableClass = VervetDB.IndividualAlignment
			query = TableClass.query.filter(TableClass.median_depth!=None)
			if ref_ind_seq_id:
				query = query.filter(TableClass.ref_ind_seq_id==ref_ind_seq_id)
			alignmentLs = query.order_by(TableClass.id)
			
		writer = MatrixFile(inputFname=outputFname, openMode='w', delimiter='\t')
		if outputFileFormat==1:
			header = ['alignmentID', 'medianDepth']
		else:
			header = ['alignmentID', 'minDepth', 'maxDepth', 'minGQ']
		writer.writeHeader(header)
		
		counter = 0
		for row in alignmentLs:
			read_group = row.read_group
			if outputFileFormat==1:
				data_row = [read_group, row.median_depth]
			else:
				minDepth = row.median_depth/float(depthFoldChange)
				if abs(minDepth-0)<=0.001:	#if it's too close to 0, assign 0.
					minDepth = 0
				data_row = [read_group, minDepth, \
							row.median_depth*float(depthFoldChange), minGQ]
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