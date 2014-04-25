#!/usr/bin/env python
"""
Examples:
	%s  -o /tmp/inconsistentParentsOfSharedKids.tsv -y 2
		~/../sservice/Vervet/SeqData2012/Linda.fam 
		PlinkIBDCheck/PlinkIBDCheck_Method14_W100Z10R0.4.2012.8.13T1413/vcf2plinkVCF2PlinkMerged/pedigree.2.tfam
		/u/home/eeskin/sservice/Vervet/Expression/SolarPed.txt

	%s -o /tmp/inconsistentParentsOfSharedIndividuals.tsv
		~/../sservice/Vervet/SeqData2012/Linda.fam 
		PlinkIBDCheck/PlinkIBDCheck_Method14_W100Z10R0.4.2012.8.13T1413/vcf2plinkVCF2PlinkMerged/pedigree.2.tfam
		/u/home/eeskin/sservice/Vervet/Expression/SolarPed.txt

Description:
	2013.3.5 program to simulate genotype of a whole pedigree based on an input pool of haplotypes (founders).
"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, numpy
import networkx as nx
from pymodule import ProcessOptions, getListOutOfStr, utils, figureOutDelimiter
from pymodule import SNP
from pymodule import MatrixFile
from pymodule.utils import PassingData
from pymodule.yhio.PolymorphismTableFile import PolymorphismTableFile, OneIndividualPolymorphismData
from pymodule.pegasus.mapper.AbstractMapper import AbstractMapper
from ComparePedigreeFromMultipleInput import ComparePedigreeFromMultipleInput

class SimulatePedigreeHaplotype(ComparePedigreeFromMultipleInput):
	__doc__ = __doc__
	option_default_dict = ComparePedigreeFromMultipleInput.option_default_dict.copy()
	#option_default_dict.pop(('inputFname', 0, ))
	option_default_dict.update({
						('recombinationRate', 1, float):[1.5e-8, '', 1, 'recombination rate per meiosis'],\
						('enableIntraLocusRecombination', 0, int):[0, '', 1, 'toggle to enable intra-locus recombination, \n\
	only applicable to non-single-nucleotide loci'],\
						('inputPedigreeFname', 1, ):['', '', 1, "pedigree file, LINKAGE format, Column 1: pedigree identifier.\n\
Column 2: individual's ID\n\
Column 3: the ID of the individual's father\n\
Column 4: the ID of the individual's mother\n\
Column 5: sex (1=male, 2=female)\n\
Column 6-N: affection status (optional)"],\
						('mutationRate', 0, float):[0, '', 1, 'mutation rate per nucleotide per generation'],\
						
						('ploidy', 0, int):[2, '', 1, 'how many sets of chromosomes one individual carries (for output). ploidy=2 means diploid.'],\
						('speciesName', 0, ):['doggie', '', 1, 'a phantom species name (for output).'],\
						
						})
	
	def __init__(self, inputFnameLs=None, **keywords):
		AbstractMapper.__init__(self, inputFnameLs=inputFnameLs, **keywords)
		
	def sampleRecombinantHaplotype(self, parentalHaplotypeList=None, locusPositionList=None, recombinationRate=1.5e-8, enableIntraLocusRecombination=False,\
								chromosomeLength=None):
		"""
		2013.3.6
			argument enableIntraLocusRecombination is not used at this moment.
			
			#. figure out how many recombinations (=n) for this chromosome, from binomial distribution
			#. uniformly choose n positions from the chromosome
			#. there are N-1 possible locations for recombination. N=chromosome length
		"""
		
		
		noOfRecombinationEvents = numpy.random.binomial(chromosomeLength, recombinationRate)
		recombinationLocationList = None
		ploidy = len(parentalHaplotypeList)
		if noOfRecombinationEvents>0:
			recombinationLocationList = numpy.random.randint(1, chromosomeLength+1, size=noOfRecombinationEvents)
				#chromosomeLength+1 is excluded
			recombinationLocationList.sort()
			
			#reverse it so that the first recombination event is now the last, which first pops out upon pop(). 
			recombinationLocationList.reverse()
			nextRecombinationLocation = recombinationLocationList.pop()
			locusIndexSpanAndHaplotypeIndexList = []
			startLocusIndex =0
			currentHaplotypeIndex = 0	#start from the first haplotype
			noOfPolymorphicLoci = len(locusPositionList)
			for i in xrange(noOfPolymorphicLoci):
				locusPosition = locusPositionList[i]
				if locusPosition>=nextRecombinationLocation:
					locusIndexSpanAndHaplotypeIndexList.append((startLocusIndex, i+1, currentHaplotypeIndex))
						#i is included in current haplotype block, next one starts from i+1
					startLocusIndex = i+1
					#a recombination, switches to a different haplotype for the next block
					currentHaplotypeIndex = (currentHaplotypeIndex+1)%(ploidy)	#alternate to the next parental haplotype index
					if len(recombinationLocationList)>0:
						nextRecombinationLocation = recombinationLocationList.pop()
					else:	#last recombination event has just been visited
						#add the last block and exit
						if startLocusIndex<noOfPolymorphicLoci:	#make sure not beyond
							locusIndexSpanAndHaplotypeIndexList.append((startLocusIndex, noOfPolymorphicLoci, currentHaplotypeIndex))
						break
				elif i==noOfPolymorphicLoci-1:	#reach the last locus but there are still recombination to its right
					#add the last block
					locusIndexSpanAndHaplotypeIndexList.append((startLocusIndex, noOfPolymorphicLoci, currentHaplotypeIndex))
			
			recombinantHaplotype = parentalHaplotypeList[0]
			for locusIndexSpanAndHaplotypeIndex in locusIndexSpanAndHaplotypeIndexList:
				startLocusIndex, stopLocusIndex, currentHaplotypeIndex = locusIndexSpanAndHaplotypeIndex
				currentHaplotype = parentalHaplotypeList[currentHaplotypeIndex]
				recombinantHaplotype[startLocusIndex:stopLocusIndex] = currentHaplotype[startLocusIndex:stopLocusIndex]
		else:	#no recombination, just pick one haplotype
			i = numpy.random.randint(0, len(parentalHaplotypeList))
			recombinantHaplotype = parentalHaplotypeList[i]
		return PassingData(recombinantHaplotype=recombinantHaplotype, recombinationLocationList=recombinationLocationList)
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		#read in the pedigree graph
		graphData = self.constructPedigreeGraphFromOneFile(inputFname=self.inputPedigreeFname)
		
		# read in the haplotype pool for founders
		self.inputPolymorphismTableFile = PolymorphismTableFile(self.inputFname, openMode='r', constructSNPData=False)
		if not self.inputPolymorphismTableFile.isPhased:
			sys.stderr.write("Error input file %s has unphased polymorphism data, can't sample haplotypes.\n"%(self.inputFname))
			sys.exit(4)
		
		# set up output
		self.outputPolymorphismFile = PolymorphismTableFile(self.outputFname, openMode='w', isPhased=1, ploidy=self.ploidy)
		self.outputPolymorphismFile.writePedigreeDiGraph2IndividualTable(diGraph=graphData.DG, \
											populationName="", speciesName=self.speciesName, \
											ploidy=self.ploidy)
		self.outputPolymorphismFile.writeChrStartStopTupleList2LocusTable(chr_start_stop_list=\
										self.inputPolymorphismTableFile.locusChrStartStopList,\
										speciesName=self.speciesName, ploidy=self.ploidy)
		
		# order the pedigree members based on their distance to the founders
		founderDistance2NodeList = graphData.DG.orderMembersByDistanceToFounders()
		#2013.10.16 YH: bug needs to be fixed here. orderMembersByDistanceToFounders() does not return a data structure like founderDistance2NodeList.
		
		individualName2polymorphismData = {}
		chromosomeLength = self.inputPolymorphismTableFile.snpData.col_id_ls[-1][2]	#stop of the last locus is chromosomeLength
		#. sample haplotypes for founders / their descendents
		for founderDistance, individualNameList in founderDistance2NodeList.iteritems():
			for individualName in individualNameList:
				polymorphismData = OneIndividualPolymorphismData(isPhased=True, ploidy=self.ploidy)
				if founderDistance==0:	#sample haplotypes for founders
					if self.inputPolymorphismTableFile.ploidy==2 and self.inputPolymorphismTableFile.isPhased:
						#input is diploid, phased haplotype data
						polymorphismData = self.inputPolymorphismTableFile.sampleOneIndividualPolymorphismWithReplacement()
					else:
						# self.inputPolymorphismTableFile.ploidy==1 or self.inputPolymorphismTableFile.ploidy is None:
						#input is haplotype or unknown ploidy
						for i in xrange(self.ploidy):
							haplotype = self.inputPolymorphismTableFile.sampleOneRandomHaplotypeWithReplacement()
							polymorphismData.addHaplotype(haplotype)
				else:	#sample recombinant haplotype based on two parents' four haplotypes
					polymorphismData = OneIndividualPolymorphismData()
					parents = graphData.DG.predecessors(individualName)
					for parentName in parents:
						parentalHaplotypeList = individualName2polymorphismData.get(parentName).haplotypeList
						returnData = self.sampleRecombinantHaplotype(parentalHaplotypeList=parentalHaplotypeList, \
											locusPositionList=self.inputPolymorphismTableFile.locusStartPositionList,\
											recombinationRate=self.recombinationRate, \
											enableIntraLocusRecombination=self.enableIntraLocusRecombination, \
											chromosomeLength=chromosomeLength)
						#output recombination events
						if returnData.recombinationLocationList:
							self.outputPolymorphismFile.writeRecombinationEvents(parentName=parentName, childName=individualName, \
																recombinationLocationList=returnData.recombinationLocationList)
						polymorphismData.addHaplotype(returnData.haplotype)
				
				individualName2polymorphismData[individualName] = polymorphismData
		#. output
		self.outputPolymorphismFile.writeIndividualName2PolymorphismData(\
								individualName2polymorphismData=individualName2polymorphismData,\
								speciesName=self.speciesName, ploidy=self.ploidy)
		self.outputPolymorphismFile.close()


if __name__ == '__main__':
	main_class = SimulatePedigreeHaplotype
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(po.arguments, **po.long_option2value)
	instance.run()
