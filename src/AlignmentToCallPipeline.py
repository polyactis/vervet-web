#!/usr/bin/env python
"""
Examples:
	# scaffold (120) is used as reference in alignment. get genome sequence id from 1 to 8.
	%s -o workflow.xml -a 120 -i 1-8
	
	# 1Mb-BAC (9) is used as reference.
	%s -o workflow.xml -a 9 -i 1-4
	
Description:
	2011-7-12
		a program which generates a pegasus workflow dag (xml file) to call genotypes on all chosen alignments.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

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
from pymodule import ProcessOptions, getListOutOfStr
from Pegasus.DAX3 import *


class AlignmentToCallPipeline(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [120, 'a', 1, 'IndividualSequence.id. To pick alignments with this sequence as reference', ],\
						('ind_seq_id_ls', 1, ): ['1-3,4,5', 'i', 1, 'a comma/dash-separated list of IndividualSequence.id. alignments come from these', ],\
						('minMinorAlleleCoverage', 1, int): [3, '', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, int): [7, '', 1, 'maximum read depth for the minor allele of a heterozygous call', ],\
						('maxNoOfReadsForGenotypingError', 1, int): [1, '', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, int): [20, '', 1, 'maximum read depth for one base to be considered'],\
						('maxMajorAlleleCoverage', 1, int): [10, '', 1, 'maximum read depth'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, int): [3, '', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						("samtools_path", 1, ): [os.path.expanduser("~/bin/samtools"), '', 1, 'samtools binary'],\
						("picard_path", 1, ): [os.path.expanduser("~/script/vervet/bin/picard-tools-1.37"), '', 1, 'picard folder containing the jar binaries'],\
						("topNumberOfContigs", 1, int): [156, '', 1, 'number of contigs'],\
						('outputFname', 1, ): [None, 'o', 1, 'xml workflow output file'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-7-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if self.ind_seq_id_ls:
			self.ind_seq_id_ls = getListOutOfStr(self.ind_seq_id_ls, data_type=int)
	
	def getTopNumberOfContigs(self, topNumberOfContigs, tax_id=60711, sequence_type_id=9):
		"""
		2011-7-12
			get all the top contigs
		"""
		sys.stderr.write("Getting %s top big contigs ..."%(self.topNumberOfContigs))
		refNameSet = set([])
		
		from pymodule import GenomeDB
		db_genome = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema="genome")
		db_genome.setup(create_tables=False)
		from sqlalchemy import desc
		query = GenomeDB.AnnotAssembly.query.filter_by(tax_id=tax_id).filter_by(sequence_type_id=sequence_type_id).order_by(desc(GenomeDB.AnnotAssembly.stop))
		for row in query:
			refNameSet.add(row.chromosome)
			if len(refNameSet)>=topNumberOfContigs:
				break
		sys.stderr.write("Done.\n")
		return refNameSet
	
	def getAlignments(self, aln_ref_ind_seq_id, ind_seq_id_ls):
		"""
		2011-7-12
		
		"""
		sys.stderr.write("Getting all alignments for %s sequences with %s as reference ..."%(len(ind_seq_id_ls), aln_ref_ind_seq_id))
		alignmentLs = []
		TableClass = VervetDB.IndividualAlignment
		query = TableClass.query.filter_by(ref_ind_seq_id=aln_ref_ind_seq_id).filter(TableClass.ind_seq_id.in_(ind_seq_id_ls))
		for row in query:
			if row.path:	#it's not None
				alignmentLs.append(row)
		sys.stderr.write("%s alignments Done.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def run(self):
		"""
		2011-7-11
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		
		refNameSet = self.getTopNumberOfContigs(self.topNumberOfContigs)
		#refNameSet = set(['vervet_path2'])	#temporary when testing use the 1Mb-BAC
		refNameLs = list(refNameSet)
		refNameLsStr = ','.join(refNameLs)
		
		alignmentLs = self.getAlignments(self.aln_ref_ind_seq_id, self.ind_seq_id_ls)
		
		
		# Create a abstract dag
		workflow = ADAG("AlignmentToCallPipeline")
		
		vervetSrcPath = os.path.expanduser("~/script/vervet/src/")
		
		site_handler = "condorpool"
		
		#add the MergeSamFiles.jar file into workflow
		mergeSamFilesJar = File('MergeSamFiles.jar')
		mergeSamFilesJar.addPFN(PFN("file://" + os.path.join(self.picard_path, 'MergeSamFiles.jar'), site_handler))
		workflow.addFile(mergeSamFilesJar)
		
		"""
		
		"""
		# Add executables to the DAX-level replica catalog
		# In this case the binary is keg, which is shipped with Pegasus, so we use
		# the remote PEGASUS_HOME to build the path.
		architecture = "x86_64"
		operatingSystem = "linux"
		selectAndSplit = Executable(namespace="workflow", name="SelectAndSplitAlignment", version="1.0", os=operatingSystem, arch=architecture, installed=True)
		selectAndSplit.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "SelectAndSplitAlignment.py"), site_handler))
		workflow.addExecutable(selectAndSplit)
			
		samtools = Executable(namespace="workflow", name="samtools", version="1.0", os=operatingSystem, arch=architecture, installed=True)
		samtools.addPFN(PFN("file://" + self.samtools_path, site_handler))
		workflow.addExecutable(samtools)
			
		java = Executable(namespace="workflow", name="java", version="1.0", os=operatingSystem, arch=architecture, installed=True)
		java.addPFN(PFN("file://" + "/usr/bin/java", site_handler))
		workflow.addExecutable(java)
		
		genotypeCallByCoverage = Executable(namespace="workflow", name="GenotypeCallByCoverage", version="1.0", os=operatingSystem, arch=architecture, installed=True)
		genotypeCallByCoverage.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "GenotypeCallByCoverage.py"), site_handler))
		workflow.addExecutable(genotypeCallByCoverage)
		
		mergeGenotypeMatrix = Executable(namespace="workflow", name="MergeGenotypeMatrix", version="1.0", os=operatingSystem, arch=architecture, installed=True)
		mergeGenotypeMatrix.addPFN(PFN("file://" + os.path.join(vervetSrcPath, "MergeGenotypeMatrix.py"), site_handler))
		workflow.addExecutable(mergeGenotypeMatrix)
		
		refName2jobDataLs = {}
		for alignment in alignmentLs:
			# Add input file to the DAX-level replica catalog
			inputFname = os.path.join(self.db_vervet.data_dir, alignment.path)
			input = File(inputFname)
			input.addPFN(PFN("file://" + inputFname, site_handler))
			workflow.addFile(input)
			
			inputFileBaseNamePrefix = os.path.splitext(os.path.basename(alignment.path))[0]
			outputDir = inputFileBaseNamePrefix
			
			# Add a select & split job
			selectAndSplitJob = Job(namespace="workflow", name=selectAndSplit.name, version="1.0")
			
			# add RG to this bam
			sequencer = alignment.ind_sequence.sequencer
			read_group = '%s_%s_vs_top%sContigs'%(alignment.ind_sequence.individual.code, sequencer, self.topNumberOfContigs)
			if sequencer=='454':
				platform_id = 'LS454'
			elif sequencer=='GA':
				platform_id = 'ILLUMINA'
			else:
				platform_id = 'ILLUMINA'
			selectAndSplitJob.addArguments("-i", input, "-e", refNameLsStr, "-a", read_group, "-p", platform_id, "-o", outputDir)
			# selectAndSplitJob will create the directory itself if outputDir doesn't exist
			#selectAndSplitJob.uses(input, transfer=True, register=True, link=Link.INPUT)	#don't add this file. otherwise it'll be staged. transfer=False will result deletion of original file.
			workflow.addJob(selectAndSplitJob)
			
			for refName in refNameLs:
				if refName not in refName2jobDataLs:
					refName2jobDataLs[refName] = []
				
				outputFname = os.path.join(outputDir, '%s_%s.bam'%(inputFileBaseNamePrefix, refName))
				output = File(outputFname)
				
				
				#output.addPFN(PFN("file://" + outputFname, site_handler))
				#selectAndSplitJob.uses(output, transfer=True, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
				
				# add the index job
				samtools_index_job = Job(namespace="workflow", name=samtools.name, version="1.0")
				samtools_index_job.addArguments("index", output)
				#samtools_index_job.uses(output, transfer=True, register=False, link=Link.INPUT)	#this file is input & output
				#samtools_index_job.uses(output, transfer=True, register=True, link=Link.OUTPUT)	#registering here would result in their early deletion.
				bai_output = File('%s.bai'%outputFname)
				#samtools_index_job.uses(bai_output, transfer=True, register=True, link=Link.OUTPUT)
				workflow.addJob(samtools_index_job)
				
				workflow.addDependency(parent=selectAndSplitJob, child=samtools_index_job)
				refName2jobDataLs[refName].append((output, samtools_index_job, bai_output))
		
		# add merge jobs for every reference
		refName2mergedBamCallJob = {}
		for refName, jobDataLs in refName2jobDataLs.iteritems():
			# add the index job
			picard_job = Job(namespace="workflow", name=java.name, version="1.0")
			outputFname = '%s.bam'%(refName)
			picard_output = File(outputFname)
			picard_job.addArguments('-jar', os.path.join(self.picard_path, 'MergeSamFiles.jar'), \
				'USE_THREADING=true', 'ASSUME_SORTED=true', 'OUTPUT=%s'%outputFname)
			#picard_job.uses(picard_output, transfer=False, register=False, link=Link.OUTPUT)	#registering here would result in their early deletion.
			
			workflow.addJob(picard_job)
			
			input_list = []
			
			for jobData in jobDataLs:
				bamFile, samtools_index_job, bamFileBai = jobData[:3]
				input_list.append('INPUT=%s'%bamFile.name)
				picard_job.addArguments('INPUT=%s'%bamFile.name)
				picard_job.uses(bamFileBai, transfer=False, register=True, link=Link.INPUT)	#register them here to be deleted 
				picard_job.uses(bamFile, transfer=False, register=True, link=Link.INPUT)
				#this picard merge job depends on a bunch of prior samtools index jobs
				workflow.addDependency(parent=samtools_index_job, child=picard_job)
			
			# add the index job on the merged bam file
			samtools_index_job = Job(namespace="workflow", name=samtools.name, version="1.0")
			samtools_index_job.addArguments("index", picard_output)
			samtools_index_job.uses(picard_output, transfer=False, register=True, link=Link.OUTPUT)	#write this as OUTPUT, otherwise it'll be deleted and next program won't know 
			bai_output = File('%s.bai'%outputFname)
			samtools_index_job.uses(bai_output, transfer=False, register=True, link=Link.OUTPUT)
			workflow.addJob(samtools_index_job)
			workflow.addDependency(parent=picard_job, child=samtools_index_job)
			
			#add the genotype call job after index is done
			genotypeCallByCoverage_job = Job(namespace="workflow", name=genotypeCallByCoverage.name, version="1.0")
			genotypeCallOutputFname = 'call/%s.call'%(refName)	#genotypeCallByCoverage_job would create directory "call".
			genotypeCallOutput = File(genotypeCallOutputFname)
			genotypeCallByCoverage_job.addArguments("-i", picard_output, "-n", str(len(jobDataLs)), "-o", genotypeCallOutput)
			genotypeCallByCoverage_job.uses(bai_output, transfer=False, register=True, link=Link.INPUT)	#make sure the bai file is still there upon start of this job 
			genotypeCallByCoverage_job.uses(picard_output, transfer=False, register=True, link=Link.INPUT)
			genotypeCallByCoverage_job.uses(genotypeCallOutput, transfer=False, register=True, link=Link.OUTPUT)
			workflow.addJob(genotypeCallByCoverage_job)
			workflow.addDependency(parent=samtools_index_job, child=genotypeCallByCoverage_job)
			
			refName2mergedBamCallJob[refName] = [genotypeCallOutput, genotypeCallByCoverage_job]
		
		
		#merge all genotype call files
		mergeGenotypeMatrix_job = Job(namespace="workflow", name=mergeGenotypeMatrix.name, version="1.0")
		finalCallOutputFname = '%s_genomes_vs_top%sReferences_call.tsv'%(len(alignmentLs), len(refNameLs))
		finalCallOutput = File(finalCallOutputFname)
		mergeGenotypeMatrix_job.addArguments("-o", finalCallOutput)
		mergeGenotypeMatrix_job.uses(finalCallOutput, transfer=True, link=Link.OUTPUT, register=True)
		workflow.addJob(mergeGenotypeMatrix_job)
		inputFnameLs = []
		for refName, callJobData in refName2mergedBamCallJob.iteritems():
			callJobOutput, callJob = callJobData[:2]
			#callJobOutput = callJob.used_files[1]
			
			mergeGenotypeMatrix_job.uses(callJobOutput, transfer=False, register=True, link=Link.INPUT)
			inputFnameLs.append(callJobOutput.name)
			workflow.addDependency(parent=callJob, child=mergeGenotypeMatrix_job)
		
		mergeGenotypeMatrix_job.addArguments("-i", ','.join(inputFnameLs))
		
		# Write the DAX to stdout
		outf = open(self.outputFname, 'w')
		workflow.writeXML(outf)
		


	
if __name__ == '__main__':
	main_class = AlignmentToCallPipeline
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
