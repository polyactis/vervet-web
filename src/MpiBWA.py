#!/usr/bin/env python
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiBWA.py -t 8 -i ...

	#test parallel run on desktop
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython ~/script/vervet-web/src/MpiBWA.py
		-t 6 -f /Network/Data/NCBI/hs_genome.fasta -y2
		-i /Network/Data/vervet/ref/454/ -o /Network/Data/vervet/ref/454_vs_hg19_20101230
		-m /Network/Data/vervet/ref/454_vs_hg19_20101230.bam
	
	#pass "-c 30" to bwa bwasw
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython ~/script/vervet-web/src/MpiBWA.py
		-t 6 -a "-c 30" -f /Network/Data/NCBI/hs_genome.fasta -y2
		-i /Network/Data/vervet/ref/454/ -o /Network/Data/vervet/ref/454_vs_hg19_20101230
		-m /Network/Data/vervet/ref/454_vs_hg19_20101230.bam
	
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython ~/script/vervet-web/src/MpiBWA.py
		-f ~/script/vervet-web/data/ref/BAC/1Mb-BAC.fa -t 6 -i /Network/Data/vervet/subspecies/aethiops/fastq/
		-o /Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC_by_aln
		-m /Network/Data/vervet/subspecies/aethiops/aethiops_vs_1MbBAC.bam -p -y1
	
	# run under openmpi
	mpiexec -n 3 -machinefile ~/hostfile_4 ~/script/vervet-web/src/MpiBWA.py ...
	
Description:
	2011-1-18
		program to run BWA in multi-thread and multi-node. Each computing node would output on its own.
		Due to multi-threading, request one cpu per node for the parallel job. Otherwise,
			multiple mpi instances would run on one machine and each uses multiple cpus.
		The sorted individual bam files in the output directory do not contain unmapped reads.
		The finalBamFile is a merge of all bam files in output_dir.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, random
import cPickle, subprocess
from pymodule import PassingData, importNumericArray
from sets import Set
#from Scientific import MPI
from mpi4py import MPI
from pymodule.MPIwrapper import MPI4pywrapper

num = importNumericArray()

class MpiBWA(MPI4pywrapper):
	__doc__ = __doc__
	option_default_dict = {('input_dir', 1, ): [None, 'i', 1, 'input directory containing *.fastq or *.fastq.gz'],\
						('fasta_fname', 1, ): [None, 'f', 1, 'fasta file containing chromosome sequences of a ref genome'],\
						("output_dir", 1, ): [None, 'o', 1, 'directory to hold individual bam output before they are merged into finalBamFile'],\
						("bwa_path", 1, ): [os.path.expanduser("~/bin/bwa"), '', 1, 'bwa binary'],\
						("samtools_path", 1, ): [os.path.expanduser("~/bin/samtools"), '', 1, 'samtools binary'],\
						('no_of_threads', 1, int): [6, 't', 1, 'number of threads run on each node for bwa'],\
						('alnType', 1, int): [1, 'y', 1, 'Alignment type. 1: aln (short read), 2: bwasw (long read)'],\
						('additionalArguments', 0, ): ["", 'a', 1, 'space-separated list of additional arguments passed to aln or bwasw'],\
						('pairedEndInput', 0, int): [0, 'p', 0, 'input is paired-end (alnType=1) or not'],\
						('finalBamFile', 1, ): ["", 'm', 1, 'final bam output file'],\
						('debug', 0, int): [0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int): [0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-09-29
			add option min_LD_to_output and min_MAF
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		if self.additionalArguments:	
			self.additionalArguments = self.additionalArguments.split()
		else:
			self.additionalArguments = []
		# 2010-5-30
		#self.communicator = MPI.world.duplicate()
		self.communicator = MPI.COMM_WORLD
		MPI4pywrapper.__init__(self, self.communicator, debug=self.debug, report=self.report)
	
	@classmethod
	def getPEInputFiles(cls, input_dir, isPE=True):
		"""
		2011-8-5
			add argument isPE, which flags whether input_dir contains PE or single-end reads
			become a classmethod
		2011-2-7
			for paired-end files, sequence_628BWAAXX_1_1.fastq.gz and sequence_628BWAAXX_1_2.fastq.gz
				are regarded as one pair of two files.
		"""
		sys.stderr.write("Pair input files from %s ..."%input_dir)
		pairedEndPrefix2FileLs = {}
		files = os.listdir(input_dir)
		no_of_fastq_files = 0
		for fname in files:
			fname_prefix, fname_suffix = cls.seqInputFilenamePrefix(fname)
			if fname_suffix!='.fastq':		#skip non-fastq files
				continue
			no_of_fastq_files += 1
			if isPE==True:
				pairedEndPrefix = fname_prefix[:-2]
				pairedEndOrder = fname_prefix[-2:]
				
				if pairedEndPrefix not in pairedEndPrefix2FileLs:
					pairedEndPrefix2FileLs[pairedEndPrefix] = ['', '']
				
				if pairedEndOrder=='_1':	#the first file
					pairedEndPrefix2FileLs[pairedEndPrefix][0] = fname
				else:
					pairedEndPrefix2FileLs[pairedEndPrefix][1] = fname
			else:
				pairedEndPrefix2FileLs[fname_prefix] = [fname]	#single End
		no_of_files = len(files)
		no_of_pairedEndPrefix = len(pairedEndPrefix2FileLs)
		if no_of_pairedEndPrefix>0:
			avg_no_of_files_per_prefix = no_of_fastq_files/float(no_of_pairedEndPrefix)
		else:
			avg_no_of_files_per_prefix = 0.0
		sys.stderr.write("%.2f files per one pairedEnd prefix. %s fastq files. %s total files. Done.\n"%\
						(avg_no_of_files_per_prefix, no_of_fastq_files, no_of_files))
		return pairedEndPrefix2FileLs
	
	def generate_params(self, input_dir, \
					param_obj=None, debug=None, pairedEndPrefix2FileLs=None):
		"""
		2011-2-7
			add argument pairedEndPrefix2FileLs
		2010-12-21
		"""
		if pairedEndPrefix2FileLs is not None:
			no_of_pairs = len(pairedEndPrefix2FileLs)
			sys.stderr.write("%s paired-end pairs.\n"%(no_of_pairs))
			for pairedEndPrefix, fileLs in pairedEndPrefix2FileLs.iteritems():
				yield fileLs
		else:
			files = os.listdir(input_dir)
			sys.stderr.write("%s files.\n"%(len(files)))
			for fname in files:
				fname_prefix, fname_suffix = self.seqInputFilenamePrefix(fname)
				if fname_suffix!='.fastq':		#skip non-fastq files
					continue
				yield (fname,)
	
	@classmethod
	def seqInputFilenamePrefix(cls, fname):
		"""
		2011-2-7
			become an independent function because same functionality is used in multiple places.
			
			fname is either sequence_628BWAAXX_4_1.fastq.gz or sequence_628BWAAXX_4_1.fastq (without gz).
			Prefix is always sequence_628BWAAXX_4_1.fastq.
		"""
		fname_prefix, fname_suffix = os.path.splitext(fname)
		if fname_suffix=='.gz':	#the input file is gzipped. get the new prefix
			fname_prefix, fname_suffix = os.path.splitext(fname_prefix)
		return fname_prefix, fname_suffix
	
	def runAlignmentTillSamOutput(self, inputFnameLs=[], param_obj=None):
		"""
		2011-2-7
			function handles either bwasw or aln, and single-end or paired-end.
		"""
		if param_obj.alnType==2:
			alnCommand = 'bwasw'
		elif param_obj.alnType==1:
			alnCommand = 'aln'
		
		if len(inputFnameLs)==2:	# two paired-End files
			output_fname_ls = []
			input_full_path_ls = []
			for fname in inputFnameLs:
				fname_prefix = self.seqInputFilenamePrefix(fname)[0]
				input_fname = os.path.join(param_obj.input_dir, fname)
				input_full_path_ls.append(input_fname)
				output_fname = os.path.join(param_obj.input_dir, '%s.sai'%fname_prefix)
				output_fname_ls.append(output_fname)
				commandline = [param_obj.bwa_path, alnCommand, '-t', '%s'%param_obj.no_of_threads] + param_obj.additionalArguments + \
						['-f', output_fname, param_obj.fasta_fname, input_fname]
				sys.stderr.write("%s\n"%commandline)
				p0 = subprocess.Popen(commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=sys.stderr)
				stdout_content, stderr_content = p0.communicate()
			
			# remove the _1 or _2 in the end of the two paired-end filenames.
			fname_prefix = fname_prefix[:-2]
			### run sampe to combine two paired-end results into one sam file
			# -P of sampe speeds things up but requires 4-5G memory for a human-size genome
			sampe_commandline = [param_obj.bwa_path, "sampe", '-P', param_obj.fasta_fname] + output_fname_ls + input_full_path_ls
			"bwa sampe hsref.fa ga1.sai ga2.sai ga1.fq ga2.fq | gzip > ga.sam.gz"
			p1 = subprocess.Popen(sampe_commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=subprocess.PIPE)
		else:
			fname = inputFnameLs[0]
			fname_prefix = self.seqInputFilenamePrefix(fname)[0]
			input_fname = os.path.join(param_obj.input_dir, fname)
			commandline = [param_obj.bwa_path, alnCommand, '-t', '%s'%param_obj.no_of_threads] + param_obj.additionalArguments + \
					[param_obj.fasta_fname, input_fname]
			sys.stderr.write("%s\n"%commandline)
			p1 = subprocess.Popen(commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=subprocess.PIPE)
			
		
		return p1, fname_prefix
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2011-2-7
			split the part before sam-bam-convert into self.runAlignmentTillSamOutput()
			Now it can handle two paired-end input files.
		2010-12-21
			
		"""
		node_rank = communicator.rank
		#data = cPickle.loads(data)
		result_ls = []

		for one_data in data:
			sys.stderr.write("Node no.%s working on %s ...\n"%(node_rank, one_data))
			
			p1, fname_prefix = self.runAlignmentTillSamOutput(inputFnameLs=one_data, param_obj=param_obj)[:2]
			
			"""
			# 2010-2-4
				convert sam to bam
			"""
			## convert sam into bam and remove unmapped reads (or queries)
			#sam_convert_cmdline = "%s calmd -uS - %s"%(param_obj.samtools_path, param_obj.fasta_fname)
			sam_convert_cmdline = [param_obj.samtools_path, 'view',  '-F', '4', '-bSh', '-']
			#print sam_convert_cmdline
			p2 = subprocess.Popen(sam_convert_cmdline, shell=False, stdin=p1.stdout, stderr=sys.stderr, stdout=subprocess.PIPE)
			
			"""
			# 2010-2-4
				sort it so that it could be used for merge
			"""
			bam_output_fname_prefix = '%s.sorted'%(os.path.join(param_obj.output_dir, fname_prefix))
			sam_sort_cmdline = [param_obj.samtools_path, 'sort', '-m', '1000000000', '-', bam_output_fname_prefix]	#maxMemory is down to 1G
			#because bwa, samtools view, samtools sort are all running at the same time due to the unix pipe.
			p3 = subprocess.Popen(sam_sort_cmdline, shell=False, stdin=p2.stdout, stderr=sys.stderr, stdout=subprocess.PIPE)
			try:
				stdout_content, stderr_content = p3.communicate()
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write("Error on node %s in running alignment and etc.\n"%(node_rank))
				sys.stderr.write("%s\n"%sam_convert_cmdline)
				sys.stderr.write("%s\n"%sam_sort_cmdline)
				sys.stderr.write('Except error info: %s\n'%repr(sys.exc_info()))
				
			if stderr_content:
				sys.stderr.write('stderr of %s: %s \n'%(sam_sort_cmdline, stderr_content))
			if stdout_content:
				sys.stderr.write('stdout of %s: %s \n'%(sam_sort_cmdline, stdout_content))
			result_ls.append('%s.bam'%(bam_output_fname_prefix))
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2011-2-7
			collect bam output filename from each output node
		2010-12-21
		"""
		#result_ls = cPickle.loads(data)
		result_ls = data	#2011-2-10
		counter = 0
		for result in result_ls:
			output_param_obj.output_fname_ls.append(result)
			counter += 1
		sys.stderr.write("%s results were outputted.\n"%counter)
	
	
	def run(self):
		"""
		2010-12-21
		"""
		node_rank = self.communicator.rank
		output_node_rank = self.communicator.size-1
		free_computing_nodes = range(1, output_node_rank)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		
		if node_rank == 0:
			try:	#2011-1-18 directory creation now done by the master node (rather than each computing node)
				if not os.path.isdir(self.output_dir):
					os.makedirs(self.output_dir)
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			if self.pairedEndInput and self.alnType==1:	#2011-2-7	2011-6-28 alnType has to be 1. pair-end doesn't apply to 2(bwasw)
				pairedEndPrefix2FileLs = self.getPEInputFiles(self.input_dir)
			else:
				pairedEndPrefix2FileLs = None
			param_ls = self.generate_params(self.input_dir, debug=self.debug, pairedEndPrefix2FileLs=pairedEndPrefix2FileLs)
			
		elif node_rank in free_computing_node_set:
			pass
		else:
			pass
		
		#self.synchronize()
		if node_rank == 0:
			#parameter_list = [param_ls, array_id2col_index, data_matrix]
			#self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, block_size=self.block_size)
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0,
								pairedEndInput=self.pairedEndInput)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=1)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(input_dir=self.input_dir, fasta_fname=self.fasta_fname, output_dir=self.output_dir,\
												bwa_path=self.bwa_path, samtools_path=self.samtools_path, no_of_threads=self.no_of_threads,\
												additionalArguments = self.additionalArguments, alnType=self.alnType)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_param_obj = PassingData(output_fname_ls=[])
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
		self.synchronize()	#to avoid some node early exits
		
		if node_rank==output_node_rank:
			"""
			2011-2-7
				to merge all bam files into one and index it
			"""
			cmdline = [self.samtools_path, 'merge', self.finalBamFile] + output_param_obj.output_fname_ls
			p0 = subprocess.Popen(cmdline, shell=False, stdin=None, stderr=sys.stderr, stdout=sys.stderr)
			try:
				stdout_content, stderr_content = p0.communicate()
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write("Error on node %s in merging all bam files.\n"%(node_rank))
				sys.stderr.write("cmdline: %s\n"%cmdline)
				sys.stderr.write('Except error info: %s\n'%repr(sys.exc_info()))
				
			if stderr_content:
				sys.stderr.write('stderr of %s: %s \n'%(cmdline, stderr_content))
			if stdout_content:
				sys.stderr.write('stdout of %s: %s \n'%(cmdline, stdout_content))
				
			cmdline = [self.samtools_path, 'index', self.finalBamFile]
			p0 = subprocess.Popen(cmdline, shell=False, stdin=None, stderr=sys.stderr, stdout=sys.stderr)
			try:
				stdout_content, stderr_content = p0.communicate()
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write("Error on node %s in merging all bam files.\n"%(node_rank))
				sys.stderr.write("cmdline: %s\n"%cmdline)
				sys.stderr.write('Except error info: %s\n'%repr(sys.exc_info()))
			
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiBWA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
