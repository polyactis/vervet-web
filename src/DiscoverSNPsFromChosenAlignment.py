#!/usr/bin/env python
"""

Examples:
	%s -s 1 -u yh -o /tmp/topConservedHG19_vervetSNP.tsv
	
Description:
	2011-4-6
		This program pulls loci from the database, selects corresponding alignments from bam files and
			do SNP discovery purely based on coverage information.

	1. the program picks all loci based on (score_method_id, min_score, max_score)
	2. find all qualified alignment instances from db
	3. for each locus (parallelize)
		1. select region out of each alignment (samtools/pysam)
		2. sort it
		3. maybe add a custom RG to identify this genome
	4. picard merge all sub-alignments in this region
	5. snp-variant call
		
	6. output the variant matrix from all loci
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0])

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))

import subprocess, cStringIO
import VervetDB
from mpi4py import MPI
from pymodule.MPIwrapper import MPI4pywrapper


class DiscoverSNPsFromChosenAlignment(MPI4pywrapper):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['vervetdb', 'd', 1, 'stock_250k database name', ],\
						('schema', 0, ): ['public', 'k', 1, 'database schema name', ],\
						('genome_dbname', 1, ): ['genome', 'g', 1, 'genome database name', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('aln_ref_ind_seq_id', 1, int): [9, 'a', 1, 'IndividualSequence.id, Alignments using this sequence as reference will be included.', ],\
						('score_method_id', 1, int): [None, 's', 1, 'ScoreMethod.id from which the loci will be chosen', ],\
						('min_score', 1, float): [0.9, 'n', 1, 'minimum score for locus to be chosen out of Locus', ],\
						('max_score', 1, float): [1.0, 'm', 1, 'maximum score for locus to be chosen out of Locus', ],\
						('minMinorAlleleCoverage', 1, int): [3, '', 1, 'minimum read depth for an allele to be called (heterozygous or homozygous)', ],\
						('maxMinorAlleleCoverage', 1, int): [7, '', 1, 'maximum read depth for the minor allele of a heterozygous call', ],\
						('maxNoOfReadsForGenotypingError', 1, int): [1, '', 1, 'if read depth for one allele is below or equal to this number, regarded as genotyping error ', ],\
						('maxNoOfReads', 1, int): [20, '', 1, 'maximum read depth for one base to be considered'],\
						('maxMajorAlleleCoverage', 1, int): [10, '', 1, 'maximum read depth'],\
						('maxNoOfReadsMultiSampleMultiplier', 1, int): [3, '', 1, 'across n samples, ignore bases where read depth > n*maxNoOfReads*multiplier.'],\
						("tmp_dir", 1, ): ["/tmp", 'e', 1, 'temporary directory to hold temporary files'],\
						("samtools_path", 1, ): [os.path.expanduser("~/bin/samtools"), '', 1, 'samtools binary'],\
						("picard_path", 1, ): [os.path.expanduser("~/script/vervet-web/bin/picard-tools-1.37"), '', 1, 'picard folder containing the jar binaries'],\
						("addRGToBAM_path", 1, ): [os.path.expanduser("~/script/vervet-web/src/addRGToBAM.sh"), '', 1, 'script to add read-group to a bam file'],\
						('outputFname', 1, ): [None, 'o', 1, 'output the SNP data.'],\
						('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}

	def __init__(self,  **keywords):
		"""
		2011-4-6
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self)
	
	def selectLoci(self, score_method_id, min_score, max_score):
		"""
		2011-4-6
		"""
		sys.stderr.write("Selecting loci score_method_id=%s, min_score=%s, max_score=%s ... \n"%\
						(score_method_id, min_score, max_score))
		TableClass = VervetDB.LocusScore
		query = TableClass.query.filter_by(score_method_id=score_method_id).\
				filter(TableClass.score>=min_score).filter(TableClass.score<=max_score)
		locus_ls = []
		for row in query:
			locus_ls.append(row.locus)
		sys.stderr.write("%s loci. Done.\n"%(len(locus_ls)))
		return locus_ls
	
	def selectAlignments(self, aln_ref_ind_seq_id=None):
		"""
		2011-4-6
		"""
		sys.stderr.write("Selecting alignments aln_ref_ind_seq_id=%s ... \n"%\
						(aln_ref_ind_seq_id))
		alignment_ls = []
		query = VervetDB.IndividualAlignment.query.filter_by(ref_ind_seq_id=aln_ref_ind_seq_id)
		for row in query:
			alignment_ls.append(row)
		sys.stderr.write("%s alignments. Done.\n"%(len(alignment_ls)))
		return alignment_ls
	
	def mergeExistingBAMWithTmpPicardBAM(self, picard_input_bam_ls, previous_picard_output_fname=None, \
								picard_output_fname=None, picard_path=None):
		"""
		2011-4-7
			picard can't merge too many files altogether. so incrementally
		"""
		
		if previous_picard_output_fname and os.path.isfile(previous_picard_output_fname):
			picard_input_bam_ls.append(previous_picard_output_fname)
		
		input_list = ['INPUT=%s'%picard_input_bam for picard_input_bam in picard_input_bam_ls]
		
		picard_merge_commandline = ['java', '-jar', os.path.join(picard_path, 'MergeSamFiles.jar')] + \
				input_list + \
				['USE_THREADING=true', 'ASSUME_SORTED=true', 'OUTPUT=%s'%picard_output_fname]
		p4 = subprocess.Popen(picard_merge_commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=subprocess.PIPE)
		stdout_content, stderr_content = p4.communicate()
		
		#cleanup
		for picard_input_bam in picard_input_bam_ls:
			os.remove(picard_input_bam)
		return picard_output_fname
	
	def prepareMergedBAMGivenLoci(self, locus_ls, alignment_ls=None, tmp_dir='/tmp/', samtools_path=None, picard_path=None,\
								addRGToBAM_path=None, maxNumberOfBamsBeforeMerge=200):
		"""
		2011-4-8
			remove the intermediate RG-tagged sam file if it contains no reads and gets skipped.
		2011-4-6
			1. select region out of each alignment (samtools/pysam)
			2. sort it
			3. maybe add a custom RG to identify this genome
			4. picard merge all sub-alignments in this region
		"""
		picard_input_bam_ls = []
		locus_id_ls = []
		previous_picard_output_fname = None
		for locus in locus_ls:
			locus_id_ls.append(locus.id)
			interval = 'chr%s:%s-%s'%(locus.chromosome, locus.start, locus.stop)
			for alignment in alignment_ls:
				sam_select_cmdline = [samtools_path, 'view', '-h', alignment.path, interval]
				#print sam_select_cmdline
				p1 = subprocess.Popen(sam_select_cmdline, shell=False, stdin=None, stderr=sys.stderr, stdout=subprocess.PIPE)
				stdout_content, stderr_content = p1.communicate()
				
				
				# add RG to this bam
				sequencer = alignment.ind_sequence.sequencer
				read_group = '%s_%s'%(alignment.ind_sequence.individual.code, sequencer)
				if sequencer=='454':
					platform_id = 'LS454'
				elif sequencer=='GA':
					platform_id = 'ILLUMINA'
				else:
					platform_id = 'ILLUMINA'
				
				sam_output_fname = os.path.join(tmp_dir, '%s_%s.RG.sam'%(alignment.id, locus.id))
				outf = open(sam_output_fname, 'w')
				outf.write("@RG\tID:%s\tSM:%s\tLB:%s\tPL:%s\n"%(read_group, read_group, platform_id, platform_id))
				no_of_reads = 0
				for line in cStringIO.StringIO(stdout_content):
					if line[0]=='@':
						outf.write(line)
					else:
						no_of_reads += 1
						outf.write('%s\tRG:Z:%s\n'%(line.strip(), read_group))
				outf.close()
				if no_of_reads==0:	#skip this region for this alignment. nothing is there.
					os.remove(sam_output_fname)
					continue
				
				"""
				rgFname=os.path.join(tmp_dir, 'rg_%s_%s.txt'%(alignment.id, locus.id))
				rgF = open(rgFname, 'w')
				rgF.write("@RG\tID:%s\tSM:%s\tLB:%s\tPL:%s\n"%(read_group, read_group, platform_id, platform_id))
				rgF.close()
				
				
				cat_commandline = ['cat', rgFname, '-']
				cat_p = subprocess.Popen(cat_commandline, shell=False, stdin=p1.stdin, stderr=sys.stderr, \
									stdout=subprocess.PIPE)
				
				printf_str = '\tRG:Z:%s\n",$0'%(read_group)
				#awk_commandline = ['awk', "{ if (substr($1,1,1)=="@") print; else printf "%s\tRG:Z:"+ read_group + "\n",$0; }"]
				
				awk_p = subprocess.Popen(awk_commandline, shell=False, stdin=cat_p.stdin, stderr=sys.stderr, \
									stdout=subprocess.PIPE)
				"""
				sam2bam_convert_commandline = [samtools_path, 'view', '-bhSu', sam_output_fname]
				# '-u' is for uncompressed bam.
				sam2bam_p = subprocess.Popen(sam2bam_convert_commandline, shell=False, stdin=None, stderr=sys.stderr, \
									stdout=subprocess.PIPE)
				
				#addRGToBAM_commandline = ['-c', addRGToBAM_path, read_group, platform_id, '-']
				#p2 = subprocess.Popen(addRGToBAM_path, shell=True, stdin=p1.stdin, stderr=sys.stderr, \
				#					stdout=sys.stdout)
				"""
				# 2010-2-4
					sort it so that it could be used for merge
				"""
				bam_output_fname_prefix = os.path.join(tmp_dir, '%s_%s'%(alignment.id, locus.id))
				sam_sort_cmdline = [samtools_path, 'sort', '-m', '1000000000', '-', bam_output_fname_prefix]	#maxMemory is down to 1G
				#because bwa, samtools view, samtools sort are all running at the same time due to the unix pipe.
				p3 = subprocess.Popen(sam_sort_cmdline, shell=False, stdin=sam2bam_p.stdout, stderr=sys.stderr, stdout=subprocess.PIPE)
				stdout_content, stderr_content = p3.communicate()
				picard_input_bam_ls.append('%s.bam'%bam_output_fname_prefix)
				
				#cleanup
				os.remove(sam_output_fname)
				
				if len(picard_input_bam_ls)>maxNumberOfBamsBeforeMerge:	#too many files would break picard. merge now.
					next_picard_output_fname = os.path.join(tmp_dir, '%s.bam'%(os.path.basename(os.tmpnam())))
					previous_picard_output_fname = self.mergeExistingBAMWithTmpPicardBAM(picard_input_bam_ls, \
														previous_picard_output_fname=previous_picard_output_fname, \
														picard_output_fname=next_picard_output_fname,\
														picard_path=picard_path)
					#reset
					picard_input_bam_ls = []
		
		alignment_id_ls = [alignment.id for alignment in alignment_ls]
		alignment_id_ls = map(str, alignment_id_ls)
		locus_id_ls = map(str, locus_id_ls)
		locus_part = '%s_loci_%s'%(len(locus_id_ls), ''.join(locus_id_ls)[:10])
		alignment_part = '%s_alns_%s'%(len(alignment_id_ls), ''.join(alignment_id_ls)[:10])
		picard_output_fname = os.path.join(tmp_dir, '%s_%s.bam'%(alignment_part, locus_part))
	
		if len(picard_input_bam_ls)>0:	#it could be empty because of the last merging inside the loop
			# do the final picard merging
			self.mergeExistingBAMWithTmpPicardBAM(picard_input_bam_ls, previous_picard_output_fname=previous_picard_output_fname, \
											picard_output_fname=picard_output_fname,\
											picard_path=picard_path)
		elif previous_picard_output_fname is not None:
			os.rename(previous_picard_output_fname, picard_output_fname)
		else:
			picard_output_fname = None
		
		if picard_output_fname is not None:
			#without indexing, pysam won't work
			samtools_index_commandline = [samtools_path, 'index', picard_output_fname]
			p5 = subprocess.Popen(samtools_index_commandline, shell=False, stdin=None, stderr=sys.stderr, stdout=sys.stdout)
			stdout_content, stderr_content = p5.communicate()
		
		return picard_output_fname
	
	def run(self,):
		"""
		2011-4-6
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
			debug = True
		else:
			debug =False
		
		
		
		db_vervet = VervetDB.AutismDB(drivername=self.drivername, username=self.db_user,
					password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db_vervet.setup(create_tables=False)
		self.db_vervet = db_vervet
		
		locus_ls = self.selectLoci(self.score_method_id, min_score=self.min_score, max_score=self.max_score)
		alignment_ls = self.selectAlignments(self.aln_ref_ind_seq_id)
		if len(alignment_ls) ==0:
			sys.stderr.write("No alignments using individual_sequence %s as reference.\n"%(self.aln_ref_ind_seq_id))
			sys.exit(0)
		picard_output_fname = self.prepareMergedBAMGivenLoci(locus_ls, alignment_ls, tmp_dir=self.tmp_dir, \
									samtools_path=self.samtools_path, picard_path=self.picard_path,\
									addRGToBAM_path=self.addRGToBAM_path)
		if picard_output_fname is not None:
			from misc import VariantDiscovery
			maxNoOfReadsForAllSamples = len(alignment_ls)*self.maxNoOfReads*self.maxNoOfReadsMultiSampleMultiplier
			VariantDiscovery.discoverHetsFromBAM(picard_output_fname, self.outputFname, \
							maxNoOfReads=self.maxNoOfReads, minMinorAlleleCoverage=self.minMinorAlleleCoverage, \
							maxMinorAlleleCoverage=self.maxMinorAlleleCoverage,\
							maxNoOfReadsForGenotypingError=self.maxNoOfReadsForGenotypingError, \
							maxMajorAlleleCoverage=self.maxMajorAlleleCoverage, \
							maxNoOfReadsForAllSamples=maxNoOfReadsForAllSamples)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DiscoverSNPsFromChosenAlignment
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()