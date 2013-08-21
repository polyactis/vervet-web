#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.0 -c 1
		-o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	%s  -a 525 -f 9 -I 8GenomeVsTop156Contigs_GATK_all_bases/call/ -N 0.8 -M0
		-c 1 -o 8GenomeVsTop156Contigs_GATK_all_bases_maxNA0.8_minMAF0_het2NA.xml -j condorpool -l condorpool  -u yh -z uclaOffice
	
	#2012.5.11 on hoffman condor, no job clustering (-C1), always need db connection on hcondor (-H)
	# set minDepth=1 (-m1)
	# add -U 0 -Z 3000 if u want to change the interval configuration
	%s ~/NetworkData/vervet/db/genotype_file/method_109/
		--oldRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3280_vervet_ref_6.0.3.fasta
		--newRefFastaFname ~/NetworkData/vervet/db/individual_sequence/3488_indID1_codeVRC_ref_sequencer3_seqType1_filtered0_version3.fasta
		--maxNoOfMismatches 2 -H -C 4 --no_of_aln_threads 1
		-j hcondor -l hcondor -D ~/NetworkData/vervet/db/ -t ~/NetworkData/vervet/db/
		-u yh -z localhost -o dags/LiftPolymorphismCoordinates/FindNewRefCoordinates_Method109_vs_3488_BWA.xml
		--alignmentMethodType 2 --intervalSize 4000
		#-U 0 -Z 3000

Description:
	2012-10-03 unfinished child of pymodule/polymorphism/FindNewRefCoordinatesGivenVCFFolderWorkflow()
	arguments for fetching the new reference chromosome info from db (GenomeDB)
	
    #. mapEachVCF:
    		split VCF into several small ones, N-SNP each
      #. mapEachInterval
         #. extract flanking sequences from the input VCF ():
         	ExtractFlankingSequenceForVCFLoci.py -i vcf -a reference fasta -o output.fasta -f flankingLength (=24)
         #. blast them
         #. run FindSNPPositionOnNewRefFromFlankingBlastOutput.py
             #. where hit length match query length, and no of mismatches <=2 => good => infer new coordinates
         #. output a mapping file between old SNP and new SNP coordinates.
         		#. reduce this thing by combining everything
         #. make a new VCF file based on the input split VCF file
            (replace contig ID , position with the new one's, remove the header part regarding chromosomes or replace it)
	(#. merge same contig small VCF files into one big file (ligateVCF job))
    #. reduce()
         #. merge all small VCF files into one big one
         #. fetch chromosome info of the new reference
         #. sort it by chromosome position
         #. split the big one into individual chromosomes
         	#. or reverse, first split , then sort
         #. add the contig-only file into db, with new reference and new genotype method ID (parent=genotype-method-id) (AddVCF2DB.py)
         
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, PassingData
from pymodule.polymorphism.FindNewRefCoordinatesGivenVCFFolderWorkflow import FindNewRefCoordinatesGivenVCFFolderWorkflow
from vervet.src.pegasus.AbstractVervetWorkflow import AbstractVervetWorkflow

parentClass = FindNewRefCoordinatesGivenVCFFolderWorkflow
class FindNewRefCoordinatesGivenGenotypeMethodWorkflow(parentClass, AbstractVervetWorkflow):
	__doc__ = __doc__
	option_default_dict = parentClass.option_default_dict.copy()
	option_default_dict.update(AbstractVervetWorkflow.db_option_dict.copy())
	option_default_dict.update({
					})
	
	#2012.9.25 no overlap and make the interval a lot smaller, (for VCF file)
	option_default_dict[('intervalOverlapSize', 1, int)][0] = 0
	option_default_dict[('intervalSize', 1, int)][0] = 5000
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		parentClass.__init__(self, **keywords)
	
	connectDB=AbstractVervetWorkflow.connectDB
	
	getReferenceSequence=AbstractVervetWorkflow.getReferenceSequence
	
if __name__ == '__main__':
	main_class = FindNewRefCoordinatesGivenGenotypeMethodWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()