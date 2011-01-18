#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiBWA.py -t 8 -i ...

	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython ~/script/vervet-web/src/MpiBWA.py
		-t 6 -f /Network/Data/NCBI/hs_genome.fasta
		-i /Network/Data/vervet/ref/454/ -o /Network/Data/vervet/ref/454_vs_hg19_20101230
	
	
Description:
	2011-1-18
		program to run BWA in multi-thread and multi-node. Each computing node would output on its own.
		Due to multi-threading, request one cpu per node for the parallel job. Otherwise,
		multiple mpi instances would run on one machine and each uses multiple cpus.
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
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper

num = importNumericArray()

class MpiBWA(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('input_dir', 1, ): [None, 'i', 1, 'a directory containing *.fastq files'],\
						('fasta_fname', 1, ): [None, 'f', 1, 'fasta file containing chromosome sequences of a genome'],\
						("output_dir", 1, ): [None, 'o', 1, 'Filename to store data matrix'],\
						("bwa_path", 1, ): [os.path.expanduser("~/bin/bwa"), '', 1, 'bwa binary'],\
						("sam_path", 1, ): [os.path.expanduser("~/bin/samtools"), '', 1, 'samtools binary'],\
						('no_of_threads', 1, int): [6, 't', 1, 'number of threads run on each node for bwa'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-09-29
			add option min_LD_to_output and min_MAF
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
		
		# 2010-5-30
		self.communicator = MPI.world.duplicate()
		MPIwrapper.__init__(self, self.communicator, debug=self.debug, report=self.report)
	
	def generate_params(cls, input_dir, \
					param_obj=None, debug=None,):
		"""
		2010-12-21
		"""
		files = os.listdir(input_dir)
		sys.stderr.write("%s files.\n"%(len(files)))
		
		for fname in files:
			yield (fname,)
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2010-12-21
			
		"""
		node_rank = communicator.rank
		data = cPickle.loads(data)
		result_ls = []

		for one_data in data:
			fname, = one_data[:1]
			sys.stderr.write("Node no.%s working on %s ...\n"%(node_rank, fname))
			fname_prefix = os.path.splitext(fname)[0]
			input_fname = os.path.join(param_obj.input_dir, fname)
			commandline = [param_obj.bwa_path, 'bwasw', '-t', param_obj.no_of_threads, param_obj.fasta_fname, input_fname]
			"""
			commandline = '%s bwasw %s %s | %s calmd -uS - %s | %s sort - %s.sorted'%\
				(param_obj.bwa_path, param_obj.fasta_fname, input_fname, \
				param_obj.sam_path, param_obj.fasta_fname, \
				param_obj.sam_path, os.path.join(param_obj.output_dir, fname_prefix))
			"""
			#calmd step is extremely slow compared with others.
			sys.stderr.write("%s\n"%commandline)
			#output_fname = os.path.join(param_obj.output_dir, '%s.sam'%fname)
			#outf = open(output_fname, 'w')
			p1 = subprocess.Popen(commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=subprocess.PIPE)
			#stdout_content, stderr_content = command_handler.communicate()
			
			#sam_convert_cmdline = "%s calmd -uS - %s"%(param_obj.sam_path, param_obj.fasta_fname)
			sam_convert_cmdline = [param_obj.sam_path, 'view', '-bS', '-']
			#print sam_convert_cmdline
			p2 = subprocess.Popen(sam_convert_cmdline, shell=False, stdin=p1.stdout, stderr=sys.stderr, stdout=subprocess.PIPE)
			sam_sort_cmdline = [param_obj.sam_path, 'sort', '-', '%s.sorted'%(os.path.join(param_obj.output_dir, fname_prefix))]
			p3 = subprocess.Popen(sam_sort_cmdline, shell=False, stdin=p2.stdout, stderr=sys.stderr, stdout=subprocess.PIPE)
			try:
				stdout_content, stderr_content = p3.communicate()
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write('Except in saving results_method_json (aborted): %s\n'%repr(sys.exc_info()))
				sys.stderr.write("Error on node %s: %s.\n"%node_rank)
				
			if stderr_content:
				sys.stderr.write('stderr of %s: %s \n'%(sam_sort_cmdline, stderr_content))
			if stdout_content:
				sys.stderr.write('stdout of %s: %s \n'%(sam_sort_cmdline, stdout_content))
			result_ls.append(1)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2010-12-21
		"""
		result_ls = cPickle.loads(data)
		counter = 0
		for result in result_ls:
			counter += 1
		sys.stderr.write("%s results were outputted.\n"%counter)
	
	def run(self):
		"""
		2010-12-21
		"""
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			try:	#2011-1-18 directory creation now done by the master node (rather than each computing node)
				if not os.path.isdir(self.output_dir):
					os.makedirs(self.output_dir)
			except:
				import traceback
				traceback.print_exc()
				sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			param_ls = self.generate_params(self.input_dir, debug=self.debug)
		elif node_rank in free_computing_node_set:
			pass
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			#parameter_list = [param_ls, array_id2col_index, data_matrix]
			#self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, block_size=self.block_size)
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=1)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(input_dir=self.input_dir, fasta_fname=self.fasta_fname, output_dir=self.output_dir,\
												bwa_path=self.bwa_path, sam_path=self.sam_path, no_of_threads=self.no_of_threads)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_param_obj = PassingData()
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiBWA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()