#!/usr/bin/env python
"""
Examples:
	# 2011-9-29
	%s  -a 524 -i ./AlignmentToCallPipeline_AllVRC_Barbados_552_554_555_626_630_649_vs_524_top_156Contigs_condor_20110922T2216/call/ 
		-o TrioInconsistency92VRC_20110922T2216.xml -l condorpool -j condorpool  -u yh -z uclaOffice
	
Description:
	2011-9-29
		a program which generates a pegasus workflow dag (xml file) to 
		
		1. distribution of inconsistent rate
		2. inconsistent rate v.s. AAF
		3. inconsistent rate v.s. position of contig
		4. distribution of AAF (todo) 
				
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


class CalculateTrioInconsistencyPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. Choose trios from alignments with this sequence as reference', ],\
						('inputDir', 0, ): ['', 'i', 1, 'input folder that contains vcf or vcf.gz files', ],\
						("vervetSrcPath", 1, ): ["%s/script/vervet/src", '', 1, 'vervet source code folder'],\
						("home_path", 1, ): [os.path.expanduser("~"), 'e', 1, 'path to the home directory on the working nodes'],\
						("site_handler", 1, ): ["condorpool", 'l', 1, 'which site to run the jobs: condorpool, hoffman2'],\
						("input_site_handler", 1, ): ["local", 'j', 1, 'which site has all the input files: local, condorpool, hoffman2. \
							If site_handler is condorpool, this must be condorpool and files will be symlinked. \
							If site_handler is hoffman2, input_site_handler=local induces file transfer and input_site_handler=hoffman2 induces symlink.'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		self.inputDir = os.path.abspath(self.inputDir)
		self.vervetSrcPath = self.vervetSrcPath%self.home_path
	
	def getAllTrios(self, db_vervet, aln_ref_ind_seq_id):
		"""
		2011-9-28
			find all trios, a list of "father,mother,child", all of which are identified by IndividualSequence.id.
			
			for each individual_alignment, -> ind_seq -> individual
				get its parents from ind2ind
				check if both parents have also been sequenced and aligned.
		"""
		sys.stderr.write("Finding trios, sequenced and aligned against IndividualSequence %s ..."%(aln_ref_ind_seq_id))
		#find out all available alignments.
		alignmentQuery = VervetDB.IndividualAlignment.query.filter_by(ref_ind_seq_id=aln_ref_ind_seq_id)
		individual_with_alignment_id2alignment = {}
		
		for individual_alignment in alignmentQuery:
			ind = individual_alignment.ind_sequence.individual
			individual_with_alignment_id2alignment[ind.id] = individual_alignment
		
		trioLs = []
		for ind_id, individual_alignment in individual_with_alignment_id2alignment.iteritems():
			#find the father first
			ind2ind = VervetDB.Ind2Ind.query.filter_by(individual2_id=ind_id).filter_by(relationship_type_id=1).first()
			if ind2ind:
				fatherIndividualID = ind2ind.individual1_id
			else:
				fatherIndividualID = None
			#find mother
			ind2ind = VervetDB.Ind2Ind.query.filter_by(individual2_id=ind_id).filter_by(relationship_type_id=2).first()
			if ind2ind:
				motherIndividualID = ind2ind.individual1_id
			else:
				motherIndividualID = None
			motherAlignment = individual_with_alignment_id2alignment.get(motherIndividualID)
			fatherAlignment = individual_with_alignment_id2alignment.get(fatherIndividualID)
			if fatherAlignment and motherAlignment:
				trio = '%s,%s,%s'%(fatherAlignment.ind_seq_id, motherAlignment.ind_seq_id, individual_alignment.ind_seq_id)
				trioLs.append(trio)
		sys.stderr.write(" %s trios fetched.\n"%(len(trioLs)))
		return trioLs
	
	def addTrioInconsistencyCalculationJob(self, workflow, executable=None, inputF=None, outputFnamePrefix=None, \
						namespace=None, version=None,\
						parentJob=None, additionalArgumentLs=[], windowSize=200000):
		"""
		2011-9-29
		"""
		job = Job(namespace=namespace, name=executable.name, version=version)
		
		job.addArguments("-i", inputF)
		job.uses(inputF, transfer=True, register=True, link=Link.INPUT)
		job.addArguments("-o", outputFnamePrefix)
		
		summaryOutputF = File("%s.summary.tsv"%outputFnamePrefix)
		job.uses(summaryOutputF, transfer=True, register=True, link=Link.OUTPUT)
		job.summaryOutputF = summaryOutputF
		
		windowOutputF = File("%s.window.%s.tsv"%(outputFnamePrefix, windowSize))
		job.uses(windowOutputF, transfer=True, register=True, link=Link.OUTPUT)
		job.windowOutputF = windowOutputF
		
		frequencyOutputF = File("%s.frequency.tsv"%(outputFnamePrefix))
		job.uses(frequencyOutputF, transfer=True, register=True, link=Link.OUTPUT)
		job.frequencyOutputF = frequencyOutputF
		
		depthOutputF = File("%s.vs.depth.tsv"%(outputFnamePrefix))
		job.uses(depthOutputF, transfer=True, register=True, link=Link.OUTPUT)
		job.depthOutputF = depthOutputF
				
		for argument in additionalArgumentLs:
			job.addArguments(argument)
		workflow.addJob(job)
		if parentJob:
			workflow.depends(parent=parentJob, child=job)
		return job
	
	def addPlotJob(self, workflow, PlotExecutable=None, outputFname=None, outputFnameLs=[], namespace=None, version=None,\
						parentJob=None, title=None):
		"""
		2011-10-20
			outputFname is the argument passed to plotJob.
			outputFnameLs contains a list of all output from this job that need to be registered and transferred.
		2011-9-29
		"""
		plotJob = Job(namespace=namespace, name=PlotExecutable.name, version=version)
		plotJob.addArguments("-o", outputFname, "-t", title)
		for outputFname in outputFnameLs:
			outputF = File(outputFname)
			plotJob.uses(outputF, transfer=True, register=True, link=Link.OUTPUT)
		workflow.addJob(plotJob)
		if parentJob:
			workflow.depends(parent=parentJob, child=plotJob)
		return plotJob
	
	def addParentToPlotJob(self, workflow, parentJob=None, parentOutputF=None, plotJob=None):
		"""
		2011-9-29
		"""
		plotJob.addArguments(parentOutputF)
		plotJob.uses(parentOutputF, transfer=True, register=True, link=Link.INPUT)
		workflow.depends(parent=parentJob, child=plotJob)
	
	def registerAllInputFiles(self, workflow, inputDir, input_site_handler=None):
		"""
		2011-9-29
			vcf files only
		"""
		sys.stderr.write("Registering input files from %s ..."%(inputDir))
		inputFLs = []
		fnameLs = os.listdir(inputDir)
		for fname in fnameLs:
			if fname[-6:]!='vcf.gz' and fname[-3:]!='vcf':	#ignore non-vcf files
				continue
			inputFname = os.path.join(inputDir, fname)
			
			inputF = File(os.path.basename(inputFname))
			inputF.addPFN(PFN("file://" + inputFname, input_site_handler))
			workflow.addFile(inputF)
			inputFLs.append(inputF)
		sys.stderr.write("%s files.\n"%(len(inputFLs)))
		return inputFLs
	
	@classmethod
	def getContigIDFromFname(cls, filename):
		"""
		2011-10-20
		"""
		import re
		contig_id_pattern = re.compile(r'Contig(\d+).*')
		contig_id_pattern_sr = contig_id_pattern.search(filename)
		if contig_id_pattern_sr:
			contig_id = contig_id_pattern_sr.group(1)
		else:
			contig_id = os.path.splitext(os.path.split(filename)[1])[0]
		return contig_id
	
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
		
		trioLs = self.getAllTrios(self.db_vervet, aln_ref_ind_seq_id=self.aln_ref_ind_seq_id)
		
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
		
		CalculateTrioInconsistency = Executable(namespace=namespace, name="CalculateTrioInconsistency", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		CalculateTrioInconsistency.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "CalculateTrioInconsistency.py"), site_handler))
		CalculateTrioInconsistency.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(CalculateTrioInconsistency)
		
		
		PlotTrioInconsistencySummaryHist = Executable(namespace=namespace, name="PlotTrioInconsistencySummaryHist", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PlotTrioInconsistencySummaryHist.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencySummaryHist.py"), site_handler))
		#PlotTrioInconsistencySummaryHist.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PlotTrioInconsistencySummaryHist)
		
		PlotTrioInconsistencyOverPosition = Executable(namespace=namespace, name="PlotTrioInconsistencyOverPosition", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PlotTrioInconsistencyOverPosition.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyOverPosition.py"), site_handler))
		#PlotTrioInconsistencyOverPosition.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PlotTrioInconsistencyOverPosition)
		
		PlotTrioInconsistencyOverFrequency = Executable(namespace=namespace, name="PlotTrioInconsistencyOverFrequency", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PlotTrioInconsistencyOverFrequency.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyOverFrequency.py"), site_handler))
		#PlotTrioInconsistencyOverFrequency.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PlotTrioInconsistencyOverFrequency)
		
		PlotTrioInconsistencyVsDepth = Executable(namespace=namespace, name="PlotTrioInconsistencyVsDepth", version=version, \
										os=operatingSystem, arch=architecture, installed=True)
		PlotTrioInconsistencyVsDepth.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "plot/PlotTrioInconsistencyVsDepth.py"), site_handler))
		#PlotTrioInconsistencyVsDepth.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%clusters_size))
		workflow.addExecutable(PlotTrioInconsistencyVsDepth)
		
		inputFLs = self.registerAllInputFiles(workflow, self.inputDir, input_site_handler=self.input_site_handler)
		
		for trio in trioLs:
			# Add a mkdir job for the call directory.
			#letting numerou genotype call jobs detect&create this directory runs into race condition.
			trioDir = trio.replace(",", "_")
			trioDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=trioDir, namespace=namespace, version=version)
			homoOnlyDir = os.path.join(trioDir, "homoOnly")
			homoOnlyDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=homoOnlyDir, namespace=namespace, version=version)
			workflow.depends(parent=trioDirJob, child=homoOnlyDirJob)
			
			allSitesDir = os.path.join(trioDir, "homoHet")
			allSitesDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=allSitesDir, namespace=namespace, version=version)
			workflow.depends(parent=trioDirJob, child=allSitesDirJob)
			
			depthOutputDir = os.path.join(trioDir, "depth")
			depthOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=depthOutputDir, namespace=namespace, version=version)
			
			positionOutputDir = os.path.join(trioDir, "position")
			positionOutputDirJob = yh_pegasus.addMkDirJob(workflow, mkdir=mkdirWrap, outputDir=positionOutputDir, namespace=namespace, version=version)
			
			# no space in a single argument
			title = "trio-%s-%s-contigs"%(trio, len(inputFLs))
			outputFname = '%s_inconsistency_summary_hist_homo_only.png'%(trioDir)
			summaryHomoOnlyPlotJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencySummaryHist, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFname = '%s_inconsistency_summary_hist_homo_het.png'%(trioDir)
			summaryAllSitesPlotJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencySummaryHist, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFname = '%s_inconsistency_over_position_homo_only.png'%(trioDir)
			inconsistencyOverPositionHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverPosition, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFname = '%s_inconsistency_over_position_homo_het.png'%(trioDir)
			inconsistencyOverPositionAllSitesJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverPosition, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFname = '%s_inconsistency_over_frequency_homo_only.png'%(trioDir)
			inconsistencyOverFrequencyHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverFrequency, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFname = '%s_inconsistency_over_frequency_homo_het.png'%(trioDir)
			inconsistencyOverFrequencyAllSitesJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverFrequency, \
						outputFname=outputFname, outputFnameLs=[outputFname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFnamePrefix = '%s_homo_only_inconsistency_over'%(trioDir)
			fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
			fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
			fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
			fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
			mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
			mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)

			inconsistencyOverDepthHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyVsDepth, \
						outputFname=outputFnamePrefix, outputFnameLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
								fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			outputFnamePrefix = '%s_homo_het_inconsistency_over'%(trioDir)
			fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
			fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
			fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
			fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
			mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
			mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
			inconsistencyOverDepthAllSitesJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyVsDepth, \
						outputFname=outputFnamePrefix, outputFnameLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
								fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
						namespace=namespace, version=version, parentJob=trioDirJob,\
						title=title)
			
			for inputF in inputFLs:
				contig_id = self.getContigIDFromFname(inputF.name)
				windowSize = 200000
				calculateTrioInconsistencyCommonArgumentLs = ["-t", trio, "-w", repr(windowSize)]
				
				outputFnamePrefix  = os.path.join(homoOnlyDir, '%s.inconsistency.homoOnly'%(os.path.basename(inputF.name)))
				trioInconsistencyCaculationHomoOnlyJob = self.addTrioInconsistencyCalculationJob(workflow, \
							executable=CalculateTrioInconsistency, \
							inputF=inputF, outputFnamePrefix=outputFnamePrefix, \
							namespace=namespace, version=version,\
							parentJob=homoOnlyDirJob, additionalArgumentLs=calculateTrioInconsistencyCommonArgumentLs + ['-m'], \
							windowSize=windowSize)	#homoOnly
				
				outputFnamePrefix  = os.path.join(allSitesDir, '%s.inconsistency.homoHet'%(os.path.basename(inputF.name)))
				trioInconsistencyCaculationAllSitesJob = self.addTrioInconsistencyCalculationJob(workflow, \
							executable=CalculateTrioInconsistency, \
							inputF=inputF, outputFnamePrefix=outputFnamePrefix, \
							namespace=namespace, version=version,\
							parentJob=allSitesDirJob, additionalArgumentLs=calculateTrioInconsistencyCommonArgumentLs,\
							windowSize=windowSize)	#all sites
				
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.summaryOutputF, \
								plotJob=summaryHomoOnlyPlotJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.windowOutputF, \
								plotJob=inconsistencyOverPositionHomoOnlyJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.frequencyOutputF, \
								plotJob=inconsistencyOverFrequencyHomoOnlyJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.depthOutputF, \
								plotJob=inconsistencyOverDepthHomoOnlyJob)
				
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.summaryOutputF, \
								plotJob=summaryAllSitesPlotJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.windowOutputF, \
								plotJob=inconsistencyOverPositionAllSitesJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.frequencyOutputF, \
								plotJob=inconsistencyOverFrequencyAllSitesJob)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.depthOutputF, \
								plotJob=inconsistencyOverDepthAllSitesJob)
				
				newTitle = "trio-%s-contig-%s"%(trio, contig_id)
				outputFname = os.path.join(positionOutputDir, '%s_contig_%s_inconsistency_over_position_homo_only.png'%(trioDir, contig_id))
				contigInconsistencyOverPositionHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverPosition, \
							outputFname=outputFname, outputFnameLs=[outputFname], \
							namespace=namespace, version=version, parentJob=positionOutputDirJob,\
							title=newTitle)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.windowOutputF, \
								plotJob=contigInconsistencyOverPositionHomoOnlyJob)
				
				outputFname = os.path.join(positionOutputDir, '%s_contig_%s_inconsistency_over_position_homo_het.png'%(trioDir, contig_id))
				contigInconsistencyOverPositionAllSitesJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyOverPosition, \
							outputFname=outputFname, outputFnameLs=[outputFname], \
							namespace=namespace, version=version, parentJob=positionOutputDirJob,\
							title=newTitle)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.windowOutputF, \
								plotJob=contigInconsistencyOverPositionAllSitesJob)
				
				outputFnamePrefix = os.path.join(depthOutputDir, '%s_contig_%s_homo_only_inconsistency_over'%(trioDir, contig_id))
				fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
				fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
				mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
				fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
				fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
				mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
				contigInconsistencyOverDepthHomoOnlyJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyVsDepth, \
							outputFname=outputFnamePrefix, outputFnameLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,\
								fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
							namespace=namespace, version=version, parentJob=depthOutputDirJob,\
							title=newTitle)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationHomoOnlyJob, \
								parentOutputF=trioInconsistencyCaculationHomoOnlyJob.depthOutputF, \
								plotJob=contigInconsistencyOverDepthHomoOnlyJob)
				
				outputFnamePrefix = os.path.join(depthOutputDir, '%s_contig_%s_homo_het_inconsistency_over'%(trioDir, contig_id))
				fa_mo_depth_Fname = '%s_fa_mo_depth.png'%(outputFnamePrefix)
				fa_child_depth_Fname = '%s_fa_child_depth.png'%(outputFnamePrefix)
				mo_child_depth_Fname = '%s_mo_child_depth.png'%(outputFnamePrefix)
				fa_mo_depth_loci_count_Fname = '%s_fa_mo_depth_loci_count.png'%(outputFnamePrefix)
				fa_child_depth_loci_count_Fname = '%s_fa_child_depth_loci_count.png'%(outputFnamePrefix)
				mo_child_depth_loci_count_Fname = '%s_mo_child_depth_loci_count.png'%(outputFnamePrefix)
				contigInconsistencyOverDepthAllSitesJob = self.addPlotJob(workflow, PlotExecutable=PlotTrioInconsistencyVsDepth, \
							outputFname=outputFnamePrefix, outputFnameLs=[fa_mo_depth_Fname, fa_child_depth_Fname, mo_child_depth_Fname,
									fa_mo_depth_loci_count_Fname, fa_child_depth_loci_count_Fname, mo_child_depth_loci_count_Fname], \
							namespace=namespace, version=version, parentJob=depthOutputDirJob,\
							title=newTitle)
				self.addParentToPlotJob(workflow, parentJob=trioInconsistencyCaculationAllSitesJob, \
								parentOutputF=trioInconsistencyCaculationAllSitesJob.depthOutputF, \
								plotJob=contigInconsistencyOverDepthAllSitesJob)
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


if __name__ == '__main__':
	main_class = CalculateTrioInconsistencyPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
