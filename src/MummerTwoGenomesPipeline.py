#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s  -a 524 -i ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o TrioInconsistency92VRC_20110922T2216.xml -l condorpool -j condorpool  -u yh -z uclaOffice
	
Description:
	2011-10-05
		grab two sequences (fasta format)
		split one or both into single-fasta entries
		run nucmer between the two (show-coords, show-aligns, draw plot)
		
		put the filter a bit more stringent so that only large synteny is displayed
			run mummerplot only after delta-filter returns something
			
		group different results into different folders to reduce the number of files in one folder
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

from sqlalchemy.types import LargeBinary

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import subprocess, cStringIO
import VervetDB
from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus
from Pegasus.DAX3 import *


class MummerTwoGenomesPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('query_seq_fname', 1, ): ['', 'q', 1, 'the query sequence fasta file', ],\
						('ref_seq_fname', 1, ): ['', 'f', 1, 'the reference sequence fasta file', ],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["condorpool", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\
						#('ref_ind_seq_id', 1, int): [120, '', 1, 'IndividualSequence.id. the reference sequence in fasta', ],\
						#('query_ind_seq_id', 1, int): [120, 'q', 1, 'IndividualSequence.id. the query sequence in fasta', ],\
	
	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
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
	
	def addSplitFastaFileJobs(self, workflow, refFastaF, selectAndSplitFasta, fastaTitleLs, mkdirWrap=None, 
								site_handler=None, namespace='workflow', version='1.0', fastaOutputDir = "fasta"):
		"""
		2011-7-25
			split the whole fasta file into files, each containing one fasta record (from fastaTitleLs)
			return the data
		"""
		sys.stderr.write("Adding job to split %s into %s records ..."%(refFastaF.name, len(fastaTitleLs)))
		# Add a mkdir job
		mkDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=fastaOutputDir, namespace=namespace, version=version)
		
		
		selectAndSplitFastaJob = Job(namespace=namespace, name=selectAndSplitFasta.name, version=version)
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
	
	def run(self):
		"""
		2011-9-28
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		
		# Create a abstract dag
		workflowName = os.path.splitext(os.path.basename(self.outputFname))[0]
		workflow = ADAG(workflowName)
		vervetSrcPath = self.vervetSrcPath
		site_handler = self.site_handler
		
		
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		namespace = "workflow"
		version="1.0"
		#clusters_size controls how many jobs will be aggregated as a single job.
		clusters_size = 20
		
		#mkdirWrap is better than mkdir that it doesn't report error when the directory is already there.
		mkdirWrap = Executable(namespace=namespace, name="mkdirWrap", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		mkdirWrap.addPFN(PFN("file://" + os.path.join(self.vervetSrcPath, "mkdirWrap.sh"), site_handler))
		mkdirWrap.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mkdirWrap)
		
		#mv to rename files and move them
		mv = Executable(namespace=namespace, name="mv", version=version, os=operatingSystem, arch=architecture, installed=True)
		mv.addPFN(PFN("file://" + "/bin/mv", site_handler))
		mv.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(mv)
		
		nucmer = Executable(namespace=namespace, name="nucmer", version=version, \
						os=operatingSystem, arch=architecture, installed=True)
		nucmer.addPFN(PFN("file://" + os.path.join(self.home_path, "bin/nucmer"), site_handler))
		#nucmer.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(nucmer)
		
		
		PostNucmer = Executable(namespace=namespace, name="PostNucmer", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PostNucmer.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "PostNucmer.sh"), site_handler))
		#PostNucmer.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PostNucmer)
		
		selectAndSplitFasta = Executable(namespace=namespace, name="SelectAndSplitFastaRecords", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		selectAndSplitFasta.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "SelectAndSplitFastaRecords.py"), site_handler))
		workflow.addExecutable(selectAndSplitFasta)
		
		ref_seq_f = File(os.path.basename(self.ref_seq_fname))
		ref_seq_f.addPFN(PFN("file://" + self.ref_seq_fname, self.input_site_handler))
		workflow.addFile(ref_seq_f)
		
		query_seq_f = File(os.path.basename(self.query_seq_fname))
		query_seq_f.addPFN(PFN("file://" + self.query_seq_fname, self.input_site_handler))
		workflow.addFile(query_seq_f)
		
		# Add a mkdir job
		deltaOutputDir = "delta"
		deltaOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=deltaOutputDir, namespace=namespace, version=version)
		
		coordsOutputDir = "coords"
		coordsOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=coordsOutputDir, namespace=namespace, version=version)
		
		filterOutputDir = "filter"
		filterOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=filterOutputDir, namespace=namespace, version=version)
		
		plotOutputDir = "plot"
		plotOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=plotOutputDir, namespace=namespace, version=version)
		
		#plotScriptOutputDir = "plotScript"
		#plotScriptOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=plotScriptOutputDir, namespace=namespace, version=version)
		
		refNameLs = self.getFastaRecordTitleLs(self.ref_seq_fname)
		returnData3 = self.addSplitFastaFileJobs(workflow, ref_seq_f, selectAndSplitFasta, refNameLs, mkdirWrap=mkdirWrap,\
						site_handler=site_handler, namespace=namespace, version=version, fastaOutputDir='refFasta')
		refName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		queryNameLs = self.getFastaRecordTitleLs(self.query_seq_fname)
		returnData3 = self.addSplitFastaFileJobs(workflow, query_seq_f, selectAndSplitFasta, queryNameLs, mkdirWrap=mkdirWrap,\
						site_handler=site_handler, namespace=namespace, version=version, fastaOutputDir='queryFasta')
		queryName2splitFastaJobDataLs = returnData3.refName2jobDataLs
		
		ref_seq_prefix = os.path.splitext(ref_seq_f.name)[0]
		for queryName, jobDataLs in queryName2splitFastaJobDataLs.iteritems():
			for refName, refJobDataLs in refName2splitFastaJobDataLs.iteritems():
				refSelectAndSplitFastaJob, refFastaFile = refJobDataLs[:2]
				selectAndSplitFastaJob, fastaFile = jobDataLs[:2]
				nucmerJob = Job(namespace=namespace, name=nucmer.name, version=version)
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
				nucmerJob.addProfile(Profile(Namespace.GLOBUS, key="maxmemory", value="%s"%job_max_memory))
				nucmerJob.addProfile(Profile(Namespace.CONDOR, key="requirements", value="(memory>=%s)"%job_max_memory))
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
				postNucJob = Job(namespace=namespace, name=PostNucmer.name, version=version)
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
				
				workflow.addJob(postNucJob)
				workflow.depends(parent=nucmerJob, child=postNucJob)
				workflow.depends(parent=coordsOutputDirJob, child=postNucJob)
				workflow.depends(parent=filterOutputDirJob, child=postNucJob)
				workflow.depends(parent=plotOutputDirJob, child=postNucJob)
				#workflow.depends(parent=plotScriptOutputDirJob, child=postNucJob)
				
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = MummerTwoGenomesPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
