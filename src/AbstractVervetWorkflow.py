#!/usr/bin/env python
"""
2012.6.12
	a NGS-workflow that derives from AbstractVCFWorkflow and specific for vervet repository
"""
import sys, os, math

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, NextGenSeq
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractVCFWorkflow import AbstractVCFWorkflow
import VervetDB

class AbstractVervetWorkflow(AbstractVCFWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractVCFWorkflow.option_default_dict.copy()
	option_default_dict.update(AbstractVCFWorkflow.db_option_dict)
	
	option_default_dict.update({
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.6.12
		"""
		AbstractVCFWorkflow.__init__(self, **keywords)
	
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
			writer.writerow([row.getReadGroup(), minDepth, \
							row.median_depth*float(foldChange), minGQ])
			counter += 1
		del writer
		sys.stderr.write("%s entries fetched.\n"%(counter))
	
