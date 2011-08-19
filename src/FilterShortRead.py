#!/usr/bin/env python
"""
Examples:
	#both input and output are gzipped
	%s -i /tmp/3_Barbados_GA/gerald_62FGFAAXX_3_1.fastq.gz -o /tmp/gerald_62FGFAAXX_3_1_filtered.fastq.gz
	
	#open an unzipped fastq and output to a zipped fastq
	%s -i /tmp/gerald_62FGFAAXX_3_1.fastq  -o /tmp/gerald_62FGFAAXX_3_1_filtered.fastq.gz

Description:
	2011-8-18
		deserted now cuz it's too slow. Use picard/dist/FilterRead.jar instead.
	2011-8-15
		this filter program does following:
			0. get smoothed phred score
			1. head trimming
			2. tail trimming
			3. determine if the read after trimming is still long enough (>=minFinalReadLength)
			4. convert low-quality base to N
			5. if percentage of Ns is <=maxNPercentage

"""

import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import numpy
import VervetDB
from pymodule import PassingData, ProcessOptions, utils, yh_matplotlib

class FilterShortRead(object):
	__doc__ = __doc__
	option_default_dict = {
						('inputFname', 1, ): ['', 'i', 1, 'The input Bam file.', ],\
						('quality_score_format', 1, ): ['Sanger', 'q', 1, 'could be Standard (phred+33, Sanger), or Illumina (~phred+64). Illumina1.8+ (after 2011-02) is Sanger.', ],\
						('outputFname', 1, ): [None, 'o', 1, 'output filename'],\
						('minValidPhredScore', 0, int): [20, 'm', 1, 'any base with phred score below this number will be turned into N', ],\
						('maxNPercentage', 0, float): [0.30, 'a', 1, 'any read with Ns more than this percentage will be discarded', ],\
						('minFinalReadLength', 0, int): [30, 'n', 1, 'final read must have length>=this number', ],\
						('halfWindowSize', 0, int): [1, 'l', 1, 'during head/tail trimming, smoothing phred score is applied. amounts to how many flanking bases used', ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def readFilter(self, readData, pmData=None):
		"""
		2011-8-17
			0. get smoothed phred score
			1. head trimming
			2. tail trimming
			3. determine if the read after trimming is still long enough (>=minFinalReadLength)
			4. convert low-quality base to N
			5. if percentage of Ns is <=maxNPercentage
			
		"""
		base_string = readData.base_string
		quality_string = readData.quality_string
		
		halfWindowSize = getattr(pmData, 'halfWindowSize', 1)
		
		read_length = len(base_string)
		smooth_phred_score_ls = []
		phred_score_ls = []
		for i in range(read_length):
			base = base_string[i]
			base_quality = quality_string[i]
			if pmData.quality_score_format=='Illumina':
				phredScore = utils.getPhredScoreOutOfSolexaScore(base_quality)
			else:
				phredScore = ord(base_quality)-33
			phred_score_ls.append(phredScore)
			if i>=halfWindowSize:	#calculate the smooth quality for base at i-halfWindowSize
				smooth_start_index = max(i-halfWindowSize*2, 0)
				smooth_stop_index = i+1
				smooth_quality = numpy.median(phred_score_ls[smooth_start_index:smooth_stop_index])
				smooth_phred_score_ls.append(smooth_quality)
		
		
		leftOver = len(phred_score_ls)-len(smooth_phred_score_ls)
		for i in range(leftOver):
			indexOfInterest = len(smooth_phred_score_ls)
			smooth_start_index = max(indexOfInterest-halfWindowSize, 0)
			smooth_stop_index = min(indexOfInterest+halfWindowSize+1, len(phred_score_ls))
			smooth_quality = numpy.median(phred_score_ls[smooth_start_index:smooth_stop_index])
			smooth_phred_score_ls.append(smooth_quality)
		
		badHeadStopIndex = -1	#start from the one before 0
		foundBadHeadStop = False
		badTailStart = len(base_string)	#start from the one after the final base
		foundBadTailStart = False
		for i in range(read_length):
			if (foundBadHeadStop and foundBadTailStart) or badHeadStopIndex>=badTailStart:	#stop right here
				break
			if foundBadHeadStop is False:	#check starts from the beginning
				phredScore = smooth_phred_score_ls[i]
				if phredScore<pmData.minValidPhredScore:
					badHeadStopIndex += 1
				elif phredScore>=pmData.minValidPhredScore:
					foundBadHeadStop = True
			if foundBadTailStart is False:	#check starts from the tail
				phredScore = smooth_phred_score_ls[-(i+1)]
				if phredScore<pmData.minValidPhredScore:
					badTailStart -= 1
				elif phredScore>=pmData.minValidPhredScore:
					foundBadTailStart = True
		
		filtered_read_length = badTailStart-badHeadStopIndex-1
		returnData = None
		if filtered_read_length>=pmData.minFinalReadLength:
			filtered_old_base_string = base_string[badHeadStopIndex+1:badTailStart]
			filtered_base_string = ''
			filtered_quality_string = quality_string[badHeadStopIndex+1:badTailStart]
			filtered_phred_score_ls = phred_score_ls[badHeadStopIndex+1:badTailStart]
			
			no_of_Ns = 0.
			for i in range(filtered_read_length):
				if filtered_phred_score_ls[i]<pmData.minValidPhredScore:
					filtered_base_string += 'N'
					no_of_Ns += 1
				else:
					filtered_base_string += filtered_old_base_string[i]
			N_percentage = no_of_Ns/filtered_read_length
			if N_percentage<=pmData.maxNPercentage:
				returnData=PassingData(base_string=filtered_base_string, quality_string=filtered_quality_string)
		return returnData
		
	
	def walkFastq(self, inputFname, outputFname = None, filterParamData=None):
		"""
		2011-8-15
		"""
		sys.stderr.write("Walking through %s ...\n"%(inputFname))
		
		fname_prefix, fname_suffix = os.path.splitext(inputFname)
		if fname_suffix=='.gz':	#the file is gzipped
			import gzip
			inf = gzip.open(inputFname, 'rb')
		else:
			inf = open(inputFname, 'r')
		
		fname_prefix, fname_suffix = os.path.splitext(outputFname)
		if fname_suffix=='.gz':	#the file is gzipped
			import gzip
			outf = gzip.open(outputFname, 'wb')
		else:
			outf = open(outputFname, 'w')
		
		counter = 0
		real_counter = 0
		rawBaseCount = 0
		filteredBaseCount = 0
		for line in inf:
			if line[0]=='@':	#a new read
				counter += 1
				fastaTitle = line #no need for stripping.
				base_string = inf.next().strip()
				qualityTitle = inf.next()	#no need for stripping.
				quality_string = inf.next().strip()
				rawBaseCount += len(base_string)
				
				readData = PassingData(base_string=base_string, quality_string=quality_string)
				filteredReadData = self.readFilter(readData, filterParamData)
				if filteredReadData:
					real_counter += 1
					filteredBaseCount  += len(filteredReadData.base_string)
					outf.write(fastaTitle)
					outf.write("%s\n"%filteredReadData.base_string)
					outf.write(qualityTitle)
					outf.write("%s\n"%filteredReadData.quality_string)
			
			if counter%5000==0 and self.report:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, real_counter, counter))
			
		del inf, outf
		sys.stderr.write("%s/%s reads selected; %s/%s bases left.\n"%(real_counter, counter, filteredBaseCount, rawBaseCount))
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		filterParamData = PassingData(minValidPhredScore=self.minValidPhredScore, maxNPercentage=self.maxNPercentage, \
							minFinalReadLength=self.minFinalReadLength, quality_score_format=self.quality_score_format,\
							halfWindowSize=self.halfWindowSize)
		
		qualityDataStructure = self.walkFastq(self.inputFname, outputFname=self.outputFname,\
									filterParamData=filterParamData)
		
		

if __name__ == '__main__':
	main_class = FilterShortRead
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
