#!/usr/bin/env python
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec %s -u yh ...

	# 2010 test parallel mpich run on desktop (exclude uclaOffice as its ssh port is on 1999)
	mpirun -np 5 -machinefile ~/hostfile /usr/bin/mpipython %s
	 -u yh -o ~/sequence_baseCount.tsv  -i 3-5 -z dl324b-1.cmb.usc.edu -c -p secret
	
	# 2011-8-30 test parallel openmpi run on desktop (exclude uclaOffice as its ssh port is on 1999)
	mpiexec -n 5 -machinefile ~/hostfile %s
	 -u yh -o ~/sequence_baseCount.tsv  -i 3-5 -z dl324b-1.cmb.usc.edu -c -p secret
	
	# run under openmpi on the hoffman2 (tunnel first)
	ssh -N -f -L 5432:dbserver.com:5432 username@login4
	mpiexec -n 3 -machinefile ~/hostfile %s 
		 -u yh -o ~/sequence_baseCount.tsv  -i 3-5 -c -p secret -z localhost
		 -t ~/extraStorage/NetworkData/vervet/db/
	
Description:
	2011-8-5
		a program counts the number of reads for individual sequences recorded in db,
			individual file is in either fasta or fastq format.
			
		The output file has only two columns: individual_sequence.id baseCount.
			The output data will also be updated into db (individual_sequence.coverage & individual_sequence.base_count).
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, random
import cPickle, subprocess
from pymodule import PassingData, importNumericArray, getListOutOfStr
#from Scientific import MPI
from mpi4py import MPI
from pymodule.MPIwrapper import MPI4pywrapper
from MpiBWA import MpiBWA
from vervet.src import VervetDB

num = importNumericArray()

class MpiBaseCount(MPI4pywrapper):
	__doc__ = __doc__
	option_default_dict = {
						('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						("dataDir", 0, ): ["", 't', 1, 'the base directory where all db-affiliated files are stored. If not given, use the default stored in db.'],\
						("outputFname", 1, ): [None, 'o', 1, 'output file to with 2 columns: individual_sequence.id baseCount'],\
						('individualSequenceIDList', 1, ): ["", 'i', 1, 'comma,space-separated list of id in table individual_sequence'],\
						('genomeSize', 1, int): [3000000000, '', 1, 'the estimated genome size for calculating coverage=baseCount/genomeSize'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction (individual_sequence.path) and qsub jobs'],\
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
		if self.individualSequenceIDList:
			self.individualSequenceIDList = getListOutOfStr(self.individualSequenceIDList, data_type=int)	
		
		# 2010-5-30
		#self.communicator = MPI.world.duplicate()
		self.communicator = MPI.COMM_WORLD
		MPI4pywrapper.__init__(self, self.communicator, debug=self.debug, report=self.report)
	
	
	def generate_params(self, individualSequenceID2FilePairLs, \
					param_obj=None, dataDir=None, debug=None,):
		"""
		2011-8-30
			add argument dataDir
			pass dataDir to computing node as well
		"""
		for individualSequenceID, FilePairLs in individualSequenceID2FilePairLs.iteritems():
			for filePair in FilePairLs:
				#filename, format, sequence_type
				yield (individualSequenceID, filePair, dataDir)
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2011-8-30
			no more multiplying factor for Paired-end data. read every file
			
		"""
		node_rank = communicator.rank
		#data = cPickle.loads(data)
		result_ls = []
		ignore_set = set(['>', '+', '@'])
		
		individualSequenceID2baseCount = {}
		
		for one_data in data:
			sys.stderr.write("Node no.%s working on %s ...\n"%(node_rank, one_data))
			
			individualSequenceID, filePair, dataDir = one_data[:3]
			baseCount = 0
			for fileRecord in filePair:
				relativePath, format, sequence_type = fileRecord[:3]
				filename = os.path.join(dataDir, relativePath)
				fname_prefix, fname_suffix = os.path.splitext(filename)
				if fname_suffix=='.gz':	#the input file is gzipped. get the new prefix
					import gzip
					inf = gzip.open(filename, 'rb')
				else:
					inf = open(filename, 'r')
				
				for line in inf:
					if line[0] in ignore_set:
						if line[0]=='+':	#skip one more quality line
							inf.next()
						continue
					baseCount += len(line.strip())
					
					#if baseCount>10000:	#temporary, for testing
					#	break
				del inf
			individualSequenceID2baseCount[individualSequenceID] = baseCount
			
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(individualSequenceID2baseCount)))
		return individualSequenceID2baseCount
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		"""
		local_individualSequenceID2baseCount = data
		counter = 0
		for individualSequenceID, baseCount in local_individualSequenceID2baseCount.iteritems():
			if individualSequenceID not in output_param_obj.individualSequenceID2baseCount:
				output_param_obj.individualSequenceID2baseCount[individualSequenceID] = 0
			output_param_obj.individualSequenceID2baseCount[individualSequenceID] += baseCount
			counter += 1
		sys.stderr.write("%s results were outputted.\n"%counter)
	
	
	def run(self):
		"""
		2010-12-21
		"""
		node_rank = self.communicator.rank
		output_node_rank = self.communicator.size-1
		free_computing_nodes = range(1, output_node_rank)	#exclude the 1st and last node
		free_computing_node_set = set(free_computing_nodes)
		
		if node_rank == 0:
			db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db_vervet.setup(create_tables=False)
			session = db_vervet.session
			#session.begin()	#no transaction for input node as there is no data insertion
			
			if not self.dataDir:
				self.dataDir = db_vervet.data_dir
			
			individualSequenceID2FilePairLs = db_vervet.getIndividualSequenceID2FilePairLs(self.individualSequenceIDList, self.dataDir)
			param_ls = self.generate_params(individualSequenceID2FilePairLs, dataDir=self.dataDir, debug=self.debug)
			
		elif node_rank in free_computing_node_set:
			pass
		else:
			pass
		
		#self.synchronize()
		if node_rank == 0:
			#parameter_list = [param_ls, array_id2col_index, data_matrix]
			#self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, block_size=self.block_size)
			param_obj = PassingData(param_ls=param_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = param_ls, message_size=1)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData()
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_param_obj = PassingData(individualSequenceID2baseCount={})
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
		self.synchronize()	#to avoid some node early exits
		
		if node_rank==output_node_rank:
			# establish db connection and insert data into db
			db_vervet = VervetDB.VervetDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db_vervet.setup(create_tables=False)
			session = db_vervet.session
			session.begin()
			
			#output all the count collected
			import csv
			writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
			writer.writerow(['individualSequenceID', 'baseCount'])
			for individualSequenceID, baseCount in output_param_obj.individualSequenceID2baseCount.iteritems():
				writer.writerow([individualSequenceID, baseCount])
				individual_sequence = VervetDB.IndividualSequence.get(individualSequenceID)
				if individual_sequence.base_count is None or individual_sequence.base_count!=baseCount:
					individual_sequence.base_count = baseCount
					individual_sequence.coverage = baseCount/float(self.genomeSize)
					session.add(individual_sequence)
			del writer
			
			if self.commit:
				session.commit()
			else:
				session.rollback()
			
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiBaseCount
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
