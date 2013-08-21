#!/usr/bin/env python
"""
Examples:
	# 2012.8.15
	%s -f ~/NetworkData/vervet/db/individual_sequence/10_hs_genome.fasta
		-q ~/NetworkData/vervet/db/individual_sequence/12_mm_genome.fasta
		-o workflow/MummerBetweenHumanAndMacaque.xml -l hcondor -j hcondor  -C 1
	
	# 2011-9-29
	%s -f ... -q ...
		-o workflow/MummerBetweenTwoGenomes.xml -l condorpool -j condorpool -C 1
	
Description:
	2011-10-05
		grab two sequences (fasta format)
		split one or both into single-fasta entries
		run nucmer between the two (show-coords, show-aligns, draw plot)
		run PostNucmer (requiring gnumplot to generat png synteny figures)
	
	2012.8.15 ToDo (error in gnuplot script output by PostNucmer, remove "tiny", )
	
	polyacti@n6203:/u/scratch/p/polyacti/pegasus/MummerTwoGenomes/MummerBetweenHumanAndMacaque.2012.8.15T1531$ gnuplot plot/chr10_vs_10_hs_genome_chr10_plot.gp
	
	set terminal png tiny size 800,800
	                 ^
	"plot/chr10_vs_10_hs_genome_chr10_plot.gp", line 1: unrecognized terminal option
	
	polyacti@n6203:/u/scratch/p/polyacti/pegasus/MummerTwoGenomes/MummerBetweenHumanAndMacaque.2012.8.15T1531$ gnuplot plot/chr10_vs_10_hs_genome_chr10_plot.gp
	
	set ticscale 0 0
	    ^
	"plot/chr10_vs_10_hs_genome_chr10_plot.gp", line 7: warning: Deprecated syntax - please use 'set tics scale' keyword
	
		put the filter a bit more stringent so that only large synteny is displayed
			run mummerplot only after delta-filter returns something
			
		group different results into different folders to reduce the number of files in one folder
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])


#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *
from pymodule.pegasus.AbstractWorkflow import AbstractWorkflow

class MummerTwoGenomesPipeline(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	option_default_dict.update({
						('query_seq_fname', 1, ): ['', 'q', 1, 'the query sequence fasta file', ],\
						('ref_seq_fname', 1, ): ['', 'f', 1, 'the reference sequence fasta file', ],\
						("mummer_path", 1, ): ["%s/bin/MUMmer", '', 1, 'path to mummer binary programs'],\
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
						#('ref_ind_seq_id', 1, int): [120, '', 1, 'IndividualSequence.id. the reference sequence in fasta', ],\
						#('query_ind_seq_id', 1, int): [120, 'q', 1, 'IndividualSequence.id. the query sequence in fasta', ],\
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		self.pathToInsertHomePathList.append('mummer_path')	#inserted before AbstractWorkflow.__init__()
		AbstractWorkflow.__init__(self, **keywords)
	
	def getFastaRecordTitleLs(self, fastaFname=None):
		"""
		2011-10-5
			
		"""
		sys.stderr.write("Getting fasta record titles from %s ..."%(fastaFname))
		fastaTitleLs = []
		handle = open(fastaFname, "rU")
		for line in handle:
			if line[0]=='>':
				title = line.strip()[1:].split()[0]
				fastaTitleLs.append(title)
		sys.stderr.write("%s records. Done.\n"%(len(fastaTitleLs)))
		return fastaTitleLs
	
	def addSplitFastaFileJobs(self, workflow, refFastaF, SelectAndSplitFastaRecords, fastaTitleLs, mkdirWrap=None, 
								site_handler=None, namespace='workflow', version='1.0', fastaOutputDir = "fasta"):
		"""
		2011-7-25
			split the whole fasta file into files, each containing one fasta record (from fastaTitleLs)
			return the data
		"""
		sys.stderr.write("Adding job to split %s into %s records ..."%(refFastaF.name, len(fastaTitleLs)))
		# Add a mkdir job
		mkDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=fastaOutputDir, namespace=namespace, version=version)
		
		
		selectAndSplitFastaJob = Job(namespace=namespace, name=SelectAndSplitFastaRecords.name, version=version)
		selectAndSplitFastaJob.addArguments('-i', refFastaF, "-o", fastaOutputDir)
		selectAndSplitFastaJob.uses(refFastaF, transfer=True, register=True, link=Link.INPUT)
		workflow.addJob(selectAndSplitFastaJob)
		workflow.depends(parent=mkDirJob, child=selectAndSplitFastaJob)
		
		refName2jobDataLs = {}
		for refName in fastaTitleLs:
			if refName not in refName2jobDataLs:
				refName2jobDataLs[refName] = []
			selectAndSplitFastaJob.addArguments(refName)
			
			fastaFname = os.path.join(fastaOutputDir, '%s.fasta'%(refName))
			fastaFile = File(fastaFname)
			selectAndSplitFastaJob.uses(fastaFile, transfer=False, register=True, link=Link.OUTPUT)
			
			refName2jobDataLs[refName] = [selectAndSplitFastaJob, fastaFile]
		sys.stderr.write("Done.\n")
		return PassingData(refName2jobDataLs=refName2jobDataLs, workflow=workflow)
	
	def registerCustomExecutables(self, workflow=None):
		"""
		2011-11-28
		"""
		AbstractWorkflow.registerCustomExecutables(self, workflow)
		
		namespace = workflow.namespace
		version = workflow.version
		operatingSystem = workflow.operatingSystem
		architecture = workflow.architecture
		clusters_size = workflow.clusters_size
		site_handler = workflow.site_handler
		vervetSrcPath = self.vervetSrcPath
		
		executableClusterSizeMultiplierList = []	#2012.8.7 each cell is a tuple of (executable, clusterSizeMultipler (0 if u do not need clustering)
		
		
		nucmer = Executable(namespace=namespace, name="nucmer", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		nucmer.addPFN(PFN("file://" + os.path.join(self.mummer_path, "nucmer"), site_handler))
		executableClusterSizeMultiplierList.append((nucmer, 0))
		
		
		PostNucmer = Executable(namespace=namespace, name="PostNucmer", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PostNucmer.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "shell/PostNucmer.sh"), site_handler))
		executableClusterSizeMultiplierList.append((PostNucmer, 0))
		
		self.addExecutableAndAssignProperClusterSize(executableClusterSizeMultiplierList, defaultClustersSize=self.clusters_size)
		
		#2013.07.26
		self.addOneExecutableFromPathAndAssignProperClusterSize(path= os.path.join(vervetSrcPath, "mapper/SelectAndSplitFastaRecords.py"), \
												name='SelectAndSplitFastaRecords', clusterSizeMultipler=0)
		
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		
		workflow = self.initiateWorkflow()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables(workflow)
		site_handler =self.site_handler
		input_site_handler = self.input_site_handler
		
		ref_seq_f = self.registerOneInputFile(workflow, self.ref_seq_fname, folderName=self.pegasusFolderName)
		
		query_seq_f = self.registerOneInputFile(workflow, self.query_seq_fname, folderName=self.pegasusFolderName)
		
		# Add a mkdir job
		deltaOutputDir = "delta"
		deltaOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=self.mkdirWrap, outputDir=deltaOutputDir)
		
		coordsOutputDir = "coords"
		coordsOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=self.mkdirWrap, outputDir=coordsOutputDir)
		
		filterOutputDir = "filter"
		filterOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=self.mkdirWrap, outputDir=filterOutputDir)
		
		plotOutputDir = "plot"
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=self.mkdirWrap, outputDir=plotOutputDir)
		
		#plotScriptOutputDir = "plotScript"
		#plotScriptOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=self.mkdirWrap, outputDir=plotScriptOutputDir)
		
		refNameLs = self.getFastaRecordTitleLs(self.ref_seq_fname)
		returnData3 = self.addSplitFastaFileJobs(workflow, ref_seq_f, self.SelectAndSplitFastaRecords, refNameLs, mkdirWrap=self.mkdirWrap,\
						site_handler=site_handler, namespace=self.namespace, version=self.version, fastaOutputDir='refFasta')
		refName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		queryNameLs = self.getFastaRecordTitleLs(self.query_seq_fname)
		returnData3 = self.addSplitFastaFileJobs(workflow, query_seq_f, self.SelectAndSplitFastaRecords, queryNameLs, mkdirWrap=self.mkdirWrap,\
						site_handler=site_handler, namespace=self.namespace, version=self.version, fastaOutputDir='queryFasta')
		queryName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		noOfJobs = len(refName2splitFastaJobDataLs) + len(queryName2splitFastaJobDataLs)
		ref_seq_prefix = os.path.splitext(os.path.basename(ref_seq_f.name))[0]
		for queryName, jobDataLs in queryName2splitFastaJobDataLs.iteritems():
			for refName, refJobDataLs in refName2splitFastaJobDataLs.iteritems():
				refSelectAndSplitFastaJob, refFastaFile = refJobDataLs[:2]
				selectAndSplitFastaJob, fastaFile = jobDataLs[:2]
				nucmerJob = Job(namespace=self.namespace, name=self.nucmer.name, version=self.version)
				outputPrefix = "%s_vs_%s_%s"%(queryName, ref_seq_prefix, refName)
				deltaFnamePrefix = os.path.join(deltaOutputDir, outputPrefix)
				nucmerJob.addArguments("--maxgap=500", "--mincluster=100", "--prefix", deltaFnamePrefix, \
									refFastaFile, fastaFile)
				nucmerJob.uses(refFastaFile, transfer=False, register=True, link=Link.INPUT)
				nucmerJob.uses(fastaFile, transfer=False, register=True, link=Link.INPUT)
				deltaFname = "%s.delta"%(deltaFnamePrefix)
				deltaF = File(deltaFname)
				nucmerJob.uses(deltaFname, transfer=True, register=True, link=Link.OUTPUT)
				#3000M for one nucmer job with human as ref
				job_max_memory = 5000	#in MB
				yh_pegasus.setJobProperRequirement(nucmerJob, job_max_memory=job_max_memory)
				workflow.addJob(nucmerJob)
				
				workflow.depends(parent=refSelectAndSplitFastaJob, child=nucmerJob)
				workflow.depends(parent=selectAndSplitFastaJob, child=nucmerJob)
				workflow.depends(parent=deltaOutputDirJob, child=nucmerJob)
				
				coordsFname = os.path.join(coordsOutputDir, "%s.coords"%(outputPrefix))
				coordsF = File(coordsFname)
				filterFname = os.path.join(filterOutputDir, "%s.filter"%(outputPrefix))
				filterF = File(filterFname)
				plotPrefix = os.path.join(plotOutputDir, "%s_plot"%(outputPrefix))
				png_plotF = File("%s.png"%plotPrefix)
				gp_plotF = File("%s.gp"%plotPrefix)
				fplot_plotF = File("%s.fplot"%plotPrefix)
				rplot_plotF = File("%s.rplot"%plotPrefix)
				postNucJob = Job(namespace=self.namespace, name=self.PostNucmer.name, version=self.version)
				postNucJob.addArguments(deltaF, coordsF, filterF, refFastaFile, fastaFile, plotPrefix)
				postNucJob.uses(deltaF, transfer=True, register=True, link=Link.INPUT)
				postNucJob.uses(refFastaFile, transfer=False, register=True, link=Link.INPUT)
				postNucJob.uses(fastaFile, transfer=False, register=True, link=Link.INPUT)
				
				postNucJob.uses(coordsF, transfer=True, register=True, link=Link.OUTPUT)
				postNucJob.uses(filterF, transfer=True, register=True, link=Link.OUTPUT)
				postNucJob.uses(png_plotF, transfer=True, register=True, link=Link.OUTPUT)
				#leave files below behind
				#postNucJob.uses(gp_plotF, transfer=True, register=True, link=Link.OUTPUT)
				#postNucJob.uses(fplot_plotF, transfer=True, register=True, link=Link.OUTPUT)
				#postNucJob.uses(rplot_plotF, transfer=True, register=True, link=Link.OUTPUT)
				
				yh_pegasus.setJobProperRequirement(postNucJob, job_max_memory=2000)
				workflow.addJob(postNucJob)
				workflow.depends(parent=nucmerJob, child=postNucJob)
				workflow.depends(parent=coordsOutputDirJob, child=postNucJob)
				workflow.depends(parent=filterOutputDirJob, child=postNucJob)
				workflow.depends(parent=plotOutputDirJob, child=postNucJob)
				#workflow.depends(parent=plotScriptOutputDirJob, child=postNucJob)
				noOfJobs += 2
		sys.stderr.write(" %s jobs. \n"%(noOfJobs))
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = MummerTwoGenomesPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
