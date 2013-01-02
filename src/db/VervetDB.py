#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	%s -v postgresql -u yh -z localhost -d vervetdb -k public
	
	#setup database in mysql
	%s -u yh -z papaya.usc.edu
	
Description:
	2011-1-18
	This is a wrapper for the autism db database, build on top of elixir.
"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

from sqlalchemy.types import LargeBinary

#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:	   #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
#else:   #32bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, BigInteger, Integer, UnicodeText, Text, Boolean, Float, Binary, Enum
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import and_, or_, not_

from datetime import datetime

from pymodule.db import ElixirDB, TableClass
from pymodule import ProcessOptions, utils, NextGenSeq, PassingData
import hashlib


__session__ = scoped_session(sessionmaker(autoflush=False, autocommit=True))
#__metadata__ = ThreadLocalMetaData()	#2008-10 not good for pylon

__metadata__ = MetaData()
#used in getattr(individual_site_id_set, '__len__', returnZeroFunc)()
from pymodule.utils import returnZeroFunc
	
class AnalysisMethod(Entity):
	"""
	2011-4-5
		record the analysis method used in like ScoreMethod or others.
	"""
	short_name = Field(String(120))
	method_description = Field(String(8000))
	smaller_score_more_significant = Field(Integer)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='analysis_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	

class README(Entity, TableClass):
	title = Field(String(2000))
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Family(Entity, TableClass):
	short_name = Field(String(256))
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='family', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Country(Entity):
	"""
	2011-3-1
	"""
	name = Field(String(100), unique=True)
	abbr = Field(String(10))
	capital = Field(Text)
	latitude = Field(Float)
	longitude = Field(Float)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='country', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Site(Entity, TableClass):
	"""
	2012.12.6 add unique constraint on (latitude, longitude)
	2011-4-29
		add column altitude
	2011-3-1
	"""
	short_name = Field(String(256))
	description = Field(Text)
	latitude = Field(Float)
	longitude = Field(Float)
	altitude = Field(Float)	#2011-4-29
	city = Field(String(100))
	stateprovince = Field(String(100))
	region = Field(String(100))
	zippostal = Field(String(20))
	country = ManyToOne("Country", colname='country_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='site', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('short_name', 'latitude', 'longitude', 'city', 'stateprovince', 'country_id'))
	using_table_options(UniqueConstraint('short_name', 'latitude', 'longitude'))


class Group(Entity):
	"""
	2011-3-1
	"""
	id = Field(Integer,primary_key=True)
	name = Field(Unicode(512),required=True)
	user_ls = ManyToMany('User',tablename='user2group', local_colname='group_id')
	phenotype_method_ls = ManyToMany("PhenotypeMethod", tablename='group2phenotype_method', local_colname='group_id')
	individual_ls = ManyToMany("Individual",tablename='individual2group', local_colname='group_id')
	using_table_options(mysql_engine='InnoDB', useexisting=True)
	using_options(tablename='acl_group',metadata=__metadata__, session=__session__)
	#group is preserved keyword in postgresql (mysql likely)
	def __repr__(self):
		return (u'<Group: name=%s>' % self.name).encode('utf-8')

class User(Entity):
	"""
	2011-3-1
	"""
	title = Field(String(4))
	realname = Field(Unicode(512))
	email = Field(String(100))
	username = Field(String(10))
	_password = Field(String(40), colname='password', synonym='password')
	organisation = Field(Unicode(100))
	#isAdmin = Field(postgresql.Enum(('Y','N'), name=is_admin_enum_type), default='N', required=True,)
	isAdmin = Field(Enum("Y","N", name="is_admin_enum_type"), default='N', required=True,)
	group_ls = ManyToMany('Group', tablename='user2group', local_colname='user_id')
	phenotype_method_ls = ManyToMany("PhenotypeMethod",tablename='user2phenotype_method', local_colname='user_id')
	individual_ls = ManyToMany("Individual",tablename='individual2user', local_colname='user_id')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='acl_user', metadata=__metadata__, session=__session__)
	#user is preserved keyword in postgresql (mysql likely)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('realname', 'username'))
	
	def validate_password(self, password):
		"""
        Check the password against existing credentials.
        
        :param password: the password that was provided by the user to
            try and authenticate. This is the clear text version that we will
            need to match against the hashed one in the database.
        :type password: unicode object.
        :return: Whether the password is valid.
        :rtype: bool
        
        """
		hashed_pass = hashlib.sha1()
		hashed_pass.update(self._password[:8] + password)
		return self._password[8:] == hashed_pass.hexdigest()

	def _set_password(self, password):
		"""encrypts password on the fly using the encryption
        algo defined in the configuration
        """
		self._password = self.__encrypt_password(password)

	def _get_password(self):
		"""returns password
        """
		return self._password

	password = property(_get_password,_set_password)

	@classmethod
	def __encrypt_password(cls,  password):
		"""
		Hash the given password with the specified algorithm. Valid values
		for algorithm are 'md5' and 'sha1'. All other algorithm values will
		be essentially a no-op.
		"""
		hashed_password = password
		
		if isinstance(password, unicode):
			password_8bit = password.encode('UTF-8')
		else:
			password_8bit = password
		
		salt = hashlib.sha1()
		salt.update(os.urandom(60))
		salt_text = salt.hexdigest()
		hash = hashlib.sha1()
		hash.update(salt_text[:8] + password_8bit)
		hashed_password = salt_text[:8] + hash.hexdigest()
		print '*'*20,  salt_text[:8], " ", hashed_password[:8]
		
		if not isinstance(hashed_password, unicode):
			hashed_password = hashed_password.decode('UTF-8')

		return hashed_password

class GeographicIntegrity(Entity):
	"""
	2011-4-28
		a table storing different qualities of Geographic/GPS information associated with each ecotype
	"""
	short_name = Field(String(40), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='geographic_integrity', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Individual(Entity, TableClass):
	"""
	2012.12.6 get rid of latitude, longitude, altitude, it's now rolled into site.
	2012.9.27 add is_contaminated
	2012.7.5 added column sequence_batch_ls	
	2012.6.19 add column target_coverage
	2012.2.7
		make code unique
	2011-10-7
		add microchip ID
	2011-9-8
		add column comment
	2011-5-4
		add a unique constraint on ucla_id
	2011-4-8
		change (firstname, lastname) to one field, name.
	2011-3-1
		add tax_id, collector, site, group_ls, user_ls, latitude, longitude
	"""
	family = ManyToOne('Family', colname='family_id', ondelete='CASCADE', onupdate='CASCADE')
	code = Field(String(256), unique=True)
	name = Field(String(256))
	ucla_id = Field(String(256))	#2011-4-26
	sex = Field(String(256))
	birthdate = Field(DateTime)
	birthplace = Field(String(256))
	access = Field(Enum("public", "restricted", name="access_enum_type"), default='restricted')
	tax_id =Field(Integer)	#2011-3-1
	geographic_integrity = ManyToOne("GeographicIntegrity", colname='geographic_integrity_id', ondelete='CASCADE', onupdate='CASCADE')
	age = Field(Integer)	#2011-4-28
	age_cas = Field(Integer)	#2011-4-28 CAS stands for chris a. schmitt
	approx_age_group_at_collection = Field(String(256))	#2011-4-28
	collection_date = Field(DateTime)	#2011-4-28
	collector = ManyToOne("User", colname='collector_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-3-1
	site = ManyToOne("Site", colname='site_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-3-1
	group_ls = ManyToMany('Group',tablename='individual2group', local_colname='individual_id', remote_colname='group_id')	#2011-3-1
	user_ls = ManyToMany('User',tablename='individual2user', local_colname='individual_id', remote_colname='user_id')	#2011-3-1
	vrc_founder = Field(Boolean)
	comment = Field(String(4096))
	microchip_id = Field(String(4096))
	target_coverage = Field(Integer)	#2012.6.19
	is_contaminated = Field(Integer, default=0)	#2012.9.27 field to mark whether it's contaminated or not.
	sequence_batch_ls = ManyToMany('SequenceBatch', tablename='individual2batch', local_colname='individual_id')	#2012.7.5
	#ManyToOne('SequenceBatch', colname='sequence_batch_id', ondelete='CASCADE', onupdate='CASCADE')	#2012.7.5
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('family_id', 'code', 'tax_id'))
	using_table_options(UniqueConstraint('ucla_id', ))
	
	@classmethod
	def getIndividualsForACL(cls, user = None):
		"""
		2011-3-1
			get all individuals that could be accessed by this user
		"""
		TableClass = Individual
		query = TableClass.query
		clause = and_(TableClass.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause,TableClass.collector == user, TableClass.user_ls.any(User.id == user.id),
						TableClass.group_ls.any(Group.id.in_([group.id for group in user.group_ls])))
		query = query.filter(clause)
		return query
	
	def checkACL(self,user):
		"""
		2011-3-1
			check if the user could access this individual
		"""
		if self.access == 'public':
			return True
		if user is None:
			return False
		if user.isAdmin == 'Y':
			return True
		if self.collector == user: 
			return True
		if user in self.user_ls:
			return True
		if [group in self.group_ls for group in user.group_ls]: 
			return True
		return False
	
	def codeSexInNumber(self):
		"""
		2012.9.4 if sex is empty or None, return 0
		2011.12.4
			represent male as 1. represent female as 2.
		"""
		if self.sex:
			if self.sex[0]=='M':
				return 1
			else:
				return 2
		else:
			return 0
	
	def getCurrentAge(self):
		"""
		2012.2.22
			get the most current age of the monkey
		"""
		if self.birthdate:
			from datetime import datetime
			return datetime.now().year - self.birthdate.year
		elif self.collection_date and (self.age  or self.age_cas):
			extraYears = datetime.now().year - self.collection_date.year
			if self.age:
				ageInYears = self.age
			else:
				ageInYears = self.age_cas
			return ageInYears + extraYears
		else:
			return None
	
class Ind2Ind(Entity, TableClass):
	"""
	2011-5-5
		add a unique constraint
	"""
	individual1 = ManyToOne('Individual', colname='individual1_id', ondelete='CASCADE', onupdate='CASCADE')
	individual2 = ManyToOne('Individual', colname='individual2_id', ondelete='CASCADE', onupdate='CASCADE')
	relationship_type = ManyToOne('RelationshipType', colname='relationship_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ind2ind', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual1_id', 'individual2_id', 'relationship_type_id'))

class RelationshipType(Entity, TableClass):
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='relationship_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class AlignmentMethod(Entity, TableClass):
	"""
	2011-9-1
		add column command
	2011-3-3
	"""
	short_name = Field(String(256), unique=True)
	command = Field(String(256))	#sub-command of bwa
	description = Field(Text)
	individual_alignment_ls = OneToMany("IndividualAlignment")
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='alignment_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class IndividualAlignment(Entity, TableClass):
	"""
	2012.9.21 rename IndividualAlignment.ind_sequence to individual_sequence
		add file_size
	#2012.9.19 to distinguish alignments from different libraries/lanes/batches
		add individual_sequence_file_raw
	2012.7.26 add column parent_individual_alignment, mask_genotype_method
	2012.7.14 add md5sum
	2012.6.13
		add column outdated_index to accommodate old (incomplete or bad alignments) 
	2012.4.2
		add a couple of statistics, perc_reads_mapped, ... perc_mapped_to_diff_chrs.
	2011-11-28
		add pass_qc_read_base_count which counts the number of bases in all pass-QC reads (base quality>=20, read mapping quality >=30).
	2011-9-19
		add mean_depth, read_group_added
	2011-8-3
		add column median_depth, mode_depth
	2011-3-3
	"""
	individual_sequence = ManyToOne('IndividualSequence', colname='ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_sequence = ManyToOne('IndividualSequence', colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	alignment_method = ManyToOne('AlignmentMethod', colname='alignment_method_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_method_ls = ManyToMany("GenotypeMethod",tablename='genotype_method2individual_alignment', local_colname='individual_alignment_id')
	path = Field(Text)
	format = Field(String(512))
	median_depth = Field(Float)	#2011-8-2
	mode_depth = Field(Float)	#2011-8-2
	mean_depth = Field(Float)	#2011-9-12
	pass_qc_read_base_count = Field(BigInteger)	#2011-11-28	QC = (base quality>=20, read mapping quality >=30). 
	read_group_added = Field(Integer, default=0)	# 2011-9-15 0=No, 1=Yes
	perc_reads_mapped = Field(Float)	#2012.4.2
	perc_duplicates = Field(Float)	#2012.4.2
	perc_paired = Field(Float)	#2012.4.2
	perc_properly_paired = Field(Float)	#2012.4.2
	perc_both_mates_mapped = Field(Float)	#2012.4.2
	perc_singletons = Field(Float)	#2012.4.2
	perc_mapped_to_diff_chrs = Field(Float)	#2012.4.2
	perc_mapq5_mapped_to_diff_chrs = Field(Float)	#2012.4.2
	total_no_of_reads = Field(BigInteger)	#2012.4.2
	outdated_index = Field(Integer, default=0)	#2012.6.13 any non-zero means outdated. to allow multiple outdated alignments
	md5sum = Field(Text, unique=True)
	file_size = Field(BigInteger)	#2012.9.21
	#2012.7.26 the parent individual_alignment
	parent_individual_alignment = ManyToOne('IndividualAlignment', colname='parent_individual_alignment_id', ondelete='CASCADE', onupdate='CASCADE')
	#2012.7.26 mask loci of the alignment out for read-recalibration 
	mask_genotype_method = ManyToOne('GenotypeMethod', colname='mask_genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	#2012.9.19 to distinguish alignments from different libraries/lanes/batches
	individual_sequence_file_raw = ManyToOne('IndividualSequenceFileRaw', colname='individual_sequence_file_raw_id', \
											ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_alignment', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('ind_seq_id', 'ref_ind_seq_id', 'alignment_method_id', 'outdated_index', \
								'parent_individual_alignment_id', 'mask_genotype_method_id',\
								'individual_sequence_file_raw_id'))
	
	def getReadGroup(self):
		"""
		2012.9.19
			add three more optional additions to the read group
		2011-11-02
			read group is this alignment's unique identifier in a sam/bam file.
		"""
		sequencer = self.individual_sequence.sequencer
		read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
		if self.parent_individual_alignment_id:
			read_group += '_p%s'%(self.parent_individual_alignment_id)
		if self.mask_genotype_method_id:
			read_group += '_m%s'%(self.mask_genotype_method_id)
		if self.individual_sequence_file_raw_id:
			read_group += '_r%s'%(self.individual_sequence_file_raw_id)
		return read_group

	def getCompositeID(self):
		"""
		2012.1.25
			ID to be used to identify members of trios. almost same as getReadGroup()
		"""
		compositeID = '%s_%s_%s_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.id, \
								self.individual_sequence.individual.code)
		return compositeID
	
	def constructRelativePath(self, subFolder='individual_alignment', **keywords):
		"""
		2012.9.19
			add three more optional additions to the path
		2012.2.24
			fix a bug
		2012.2.10
			moved from VervetDB.constructRelativePathForIndividualAlignment
		2011-8-29
			called by getAlignment() and other programs
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		pathPrefix = '%s/%s_%s_vs_%s_by_%s'%(subFolder, self.id, self.ind_seq_id,\
											self.ref_ind_seq_id, self.alignment_method_id)
		if self.parent_individual_alignment_id:
			pathPrefix += '_p%s'%(self.parent_individual_alignment_id)
		if self.mask_genotype_method_id:
			pathPrefix += '_m%s'%(self.mask_genotype_method_id)
		if self.individual_sequence_file_raw_id:
			pathPrefix += '_r%s'%(self.individual_sequence_file_raw_id)
		dst_relative_path = '%s.%s'%(pathPrefix, self.format)
		
		return dst_relative_path
	
class IndividualSequence(Entity, TableClass):
	"""
	2012.2.27
		add column read_count
	2012.1.26
		add column individual_sequence_file_raw_ls
	2011-8-30
		add column chromosome
	2011-8-18
		add field original_path, quality_score_format, parent_individual_sequence, filtered
	2011-8-5
		change type of base_count to BigInteger
	2011-8-3
		add column base_count
	2011-5-8
		add column coverage
	2011-3-3
	"""
	individual = ManyToOne('Individual', colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	sequencer = Field(String(512))	# 454, GA, Sanger
	sequence_type = Field(String(512))	#genome, contig, SR (single-end read) or PE ...
	chromosome = Field(String(512))	#1,2,4,5,X,Y,etc
	tissue  = ManyToOne('Tissue', colname='tissue_id', ondelete='CASCADE', onupdate='CASCADE')	#2011-5-9
	coverage = Field(Float)	#2011-5-8
	read_count = Field(BigInteger)	#2012.2.27
	base_count = Field(BigInteger)	#2011-8-2
	path = Field(Text, unique=True)	#storage folder path
	format = Field(String(512))	#fasta, fastq
	original_path = Field(Text)	#the path to the original file
	quality_score_format = Field(String(512))	#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
	# Illumina1.8+ (after 2011-02) is Standard.
	parent_individual_sequence = ManyToOne('IndividualSequence', colname='parent_individual_sequence_id', ondelete='SET NULL', onupdate='CASCADE')
	filtered = Field(Integer, default=0)	#0 means not. 1 means yes.
	individual_sequence_file_ls = OneToMany("IndividualSequenceFile")
	individual_sequence_file_raw_ls = OneToMany("IndividualSequenceFileRaw")
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_sequence', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual_id', 'sequencer', 'sequence_type', 'tissue_id',\
		'filtered', "chromosome", 'parent_individual_sequence_id'))
	
	def constructRelativePath(self, subFolder='individual_sequence', **keywords):
		"""
		2012.7.13 link to constructRelativePathForIndividualSequence()
		"""
		return self.constructRelativePathForIndividualSequence(subFolder=subFolder)
	
	def constructRelativePathForIndividualSequence(self, subFolder='individual_sequence'):
		"""
		2012.2.10
			add "split" in the end of the path
		2011-8-3
			called by getIndividualSequence() and other outside programs
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		dst_relative_path = '%s/%s_%s_%s_%s_%s_%s'%(subFolder, self.id, self.individual_id, self.individual.code,\
										self.sequencer, getattr(self.tissue, 'id', 0), self.filtered)
		return dst_relative_path

class SequenceBatch(Entity, TableClass):
	"""
	2012.7.5 table to store batch names so that IndividualSequence could refer to them.
	"""
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	coverage = Field(Integer)
	individual_ls = ManyToMany('Individual', tablename='individual2batch', local_colname='sequence_batch_id')	#2012.7.5
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence_batch', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class IndividualSequenceFile(Entity, TableClass):
	"""
	2012.7.14 add md5sum
	2012.2.27
		add column read_count
	2012.2.10
		add column parent_individual_sequence_file_id
	2012.1.26
	2011-10-31
		a table recording a number of files (split fastq mostly) affiliated with one IndividualSequence entry
	"""
	individual_sequence = ManyToOne('IndividualSequence', colname='individual_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_file_raw = ManyToOne('IndividualSequenceFileRaw', colname='individual_sequence_file_raw_id', \
											ondelete='CASCADE', onupdate='CASCADE')
	library = Field(Text)	#id for the preparation library
	split_order = Field(Integer)	# the number that designates the order of this split fastq file within the large file
	mate_id = Field(Integer)	# id of the mate pair. 1 = 1st end. 2 = 2nd end. null = single-end.
	read_count = Field(BigInteger)	#2012.2.27
	base_count = Field(BigInteger)	#2011-8-2
	path = Field(Text, unique=True)	#path to the actual file
	format = Field(String(512))	#fasta, fastq
	quality_score_format = Field(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	filtered = Field(Integer, default=0)	#0 means not. 1 means yes.
	md5sum = Field(Text, unique=True)
	parent_individual_sequence_file = ManyToOne('IndividualSequenceFile', colname='parent_individual_sequence_file_id', \
											ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_sequence_file', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual_sequence_id', 'library', 'split_order', 'mate_id', 'filtered',\
										'parent_individual_sequence_file_id'))
	
	def constructRelativePath(self, subFolder='individual_sequence', sourceFilename="", **keywords):
		"""
		2012.7.13 
		"""
		folderRelativePath = self.individual_sequence.constructRelativePath(subFolder=subFolder)
		relativePath = os.path.join(folderRelativePath, '%s_%s'%(self.id, sourceFilename))
		return relativePath
	


class IndividualSequenceFileRaw(Entity, TableClass):
	"""
	2012.7.12 add column file_size
	2012.4.30
		add column mate_id
	2012.2.27
		add column read_count
	2012.1.26
		this table is used to store the bam files from WUSTL. The bam files hold either single-end or paired-end reads
			in one file.
		Technically, this table is parent of IndividualSequenceFile.
	"""
	individual_sequence = ManyToOne('IndividualSequence', colname='individual_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_file_ls = OneToMany("IndividualSequenceFile")
	library = Field(Text)	#id for the preparation library
	mate_id = Field(Integer)	#2012.4.30 to handle fastq raw sequence files.
	read_count = Field(BigInteger)	#2012.2.27
	base_count = Field(BigInteger)
	path = Field(Text)	#path to the actual file
	md5sum = Field(Text, unique=True)	#used to identify each raw file
	quality_score_format = Field(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	file_size = Field(BigInteger)	#2012.7.12
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_sequence_file_raw', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('library', 'md5sum', 'mate_id'))


class Tissue(Entity, TableClass):
	"""
	2011-5-9
	"""
	short_name = Field(String(512), unique=True)
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='tissue', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class IndividualSequence2Sequence(Entity, TableClass):
	"""
	2011-3-3
	"""
	individual1_sequence = ManyToOne('IndividualSequence', colname='individual1_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual1_chr = Field(String(512))
	individual1_start = Field(Integer)
	individual1_stop = Field(Integer)
	individual2_sequence = ManyToOne('IndividualSequence', colname='individual2_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual2_chr = Field(String(512))
	individual2_start = Field(Integer)
	individual2_stop = Field(Integer)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_seq2seq_map', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual1_sequence_id', 'individual1_chr', 'individual1_start', 'individual1_stop',\
										'individual2_sequence_id', 'individual2_chr', 'individual2_start', 'individual2_stop'))

class Locus(Entity, TableClass):
	"""
	2012.5.2
		make ref_seq_id, ref_ind_seq_id, locus_type_id required (not null)
	2012.4.29
		change locus_method_ls to locus_type
		include locus_type_id into the unique key constraint
	2011-4-5
		add table locus_method
	2011-2-3
		add ref_allele, alt_allele
	"""
	chromosome = Field(String(512))
	start = Field(Integer)
	stop = Field(Integer)
	ref_seq = ManyToOne('Sequence', colname='ref_seq_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	alt_seq = ManyToOne('Sequence', colname='alt_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_sequence = ManyToOne('IndividualSequence', colname='ref_ind_seq_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	locus_type = ManyToOne('LocusType', colname='locus_type_id', required=True, ondelete='CASCADE', onupdate='CASCADE')
	#locus_method_ls = ManyToMany('LocusMethod',tablename='locus2locus_method', local_colname='locus_id', \
	#							remote_colname='locus_method_id')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'ref_ind_seq_id', 'locus_type_id'))

class LocusScore(Entity, TableClass):
	"""
	2011-4-5
		score of locus
	"""
	locus = ManyToOne('Locus', colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	score_method = ManyToOne('ScoreMethod', colname='score_method_id', ondelete='CASCADE', onupdate='CASCADE')
	score = Field(Float)
	rank = Field(Integer)
	object = Field(LargeBinary(134217728), deferred=True)	#a python dictionary to store other attributes
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_score', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('locus_id', 'score_method_id'))

class LocusAnnotation(Entity):
	"""
	2012.5.18 add field which_codon
	2012.5.14 add fields for structural variation
	2012.4.25
		copied from SNPAnnotation of Stock_250kDB
	2009-01-05
		table to store annotation of SNPs, like synonymous, non-synonymous, ...
		information finer than SnpsContext
	"""
	locus = ManyToOne('%s.Locus'%__name__, colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_context = ManyToOne('%s.LocusContext'%__name__, colname='locus_context_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	gene_commentary_id = Field(Integer)
	gene_segment_id = Field(Integer)
	locus_annotation_type = ManyToOne('%s.LocusAnnotationType'%__name__, colname='locus_annotation_type_id', ondelete='CASCADE', onupdate='CASCADE')
	which_exon_or_intron = Field(Integer)
	pos_within_codon = Field(Integer)	#which position in the tri-nucleotide codon, this locus is at. for synonymou/non-syn nucleotide changes.
	which_codon = Field(Integer)	#2012.5.18 which AA this locus affects if synonymous, or non-synonymous, 
	label = Field(Text)	#what type of gene segment it is, exon, cds, intron, 3UTR, 5UTR
	utr_number = Field(Integer)	#2012.5.14 
	cds_number = Field(Integer)
	intron_number = Field(Integer)
	exon_number = Field(Integer)
	overlap_length = Field(Integer)	#for structural variation
	overlap_fraction_in_gene = Field(Float)		#for structural variation
	overlap_fraction_in_locus = Field(Float)	#for structural variation
	comment = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_annotation', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('locus_id', 'gene_commentary_id', 'locus_annotation_type_id'))

class LocusAnnotationType(Entity):
	"""
	2012.4.25
		copied from SNPAnnotationType of Stock_250kDB
	2009-01-05
		table to store types of SNP annotation, like synonymous, non-synonymous, ...
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_annotation_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class LocusContext(Entity):
	"""
	2012.5.14
		add overlap_length, overlap_fraction_in_locus, overlap_fraction_in_gene
	2012.4.25
		copied from SnpsContext of Stock_250kDB
	"""
	locus = ManyToOne('%s.Locus'%__name__, colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	disp_pos = Field(Float)	#[0,1) is for within gene fraction, <=-1 is upstream. >=1 is downstream
	gene_id = Field(Integer)
	gene_strand = Field(String(2))
	left_or_right = Field(String(200), deferred=True)
	disp_pos_comment = Field(String(2000), deferred=True)
	overlap_length = Field(Integer)	#for structural variation
	overlap_fraction_in_gene = Field(Float)		#for structural variation
	overlap_fraction_in_locus = Field(Float)	#for structural variation
	created_by = Field(String(200), deferred=True)
	updated_by = Field(String(200), deferred=True)
	date_created = Field(DateTime, default=datetime.now, deferred=True)
	date_updated = Field(DateTime, deferred=True)
	using_options(tablename='locus_context', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('locus_id', 'gene_id'))

class ScoreMethod(Entity, TableClass):
	"""
	2011-4-5
	"""
	short_name = Field(String(256), unique=True)
	filename = Field(String(1000), unique=True)
	original_filename = Field(Text)
	description = Field(Text)
	min_maf = Field(Float)
	genotype_method = ManyToOne('%s.GenotypeMethod'%__name__, colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	analysis_method = ManyToOne('%s.AnalysisMethod'%__name__, colname='analysis_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('%s.PhenotypeMethod'%__name__, colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	transformation_method = ManyToOne('%s.TransformationMethod'%__name__, colname='transformation_method_id', ondelete='CASCADE', onupdate='CASCADE')
	score_method_type = ManyToOne('%s.ScoreMethodType'%__name__, colname='score_method_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('genotype_method_id', 'analysis_method_id', 'phenotype_method_id', \
								'score_method_type_id', 'transformation_method_id'))

class TransformationMethod(Entity):
	"""
	2011-4-5
	"""
	name = Field(String(30))
	description = Field(Text)
	formular = Field(String(100))
	function = Field(String(20))
	using_options(tablename='transformation_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

	
class ScoreMethodType(Entity):
	"""
	2011-4-5
	"""
	short_name = Field(String(30), unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_method_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class LocusType(Entity, TableClass):
	"""
	2012.4.29
		renamed to LocusType
	2011-4-5
		to mark different sets of loci
	"""
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

"""
class Allele(Entity, TableClass):
	#locus = ManyToOne('Locus', colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_type = ManyToOne('AlleleType', colname='allele_type_id', ondelete='CASCADE', onupdate='CASCADE')
	sequence = Field(String(2048))
	sequence_length = Field(Integer)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='allele', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('allele_type_id'))
"""

class AlleleType(Entity, TableClass):
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='allele_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Sequence(Entity, TableClass):
	"""
	2011-2-4
		to store the base(s) of an allele
	"""
	sequence = Field(Text)
	sequence_length = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('sequence'))

class Genotype(Entity, TableClass):
	"""
	2011-2-3
		
	"""
	individual = ManyToOne('Individual', colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	locus = ManyToOne('Locus', colname='locus_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_order = Field(Integer)
	allele_type = ManyToOne('AlleleType', colname='allele_type_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_seq = ManyToOne('Sequence', colname='allele_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_seq_length = Field(Integer)
	score = Field(Float)
	target_locus = ManyToOne('Locus', colname='target_locus_id', ondelete='CASCADE', onupdate='CASCADE')	#for translocated allele only
	genotype_method = ManyToOne('GenotypeMethod', colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genotype', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual_id', 'locus_id', 'allele_order'))

class GenotypeMethod(Entity, TableClass):
	"""
	2012.9.6 add column min_neighbor_distance
	2012.7.12 add meta-columns: min_depth ... min_MAF
	2012.7.11
		get rid of all file-related columns, moved to GenotypeFile.
	2011-2-4
		file format:
						locus1	locus2
			individual1	allele1/allele2
			individual2
	"""
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	path = Field(Text, unique=True)
	individual_alignment_ls = ManyToMany("%s.IndividualAlignment"%(__name__), tablename='genotype_method2individual_alignment', \
							local_colname='genotype_method_id')
	ref_sequence = ManyToOne('%s.IndividualSequence'%(__name__), colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	parent = ManyToOne("%s.GenotypeMethod"%(__name__), colname='parent_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_file_ls = OneToMany('GenotypeFile')	#2012.7.17
	no_of_individuals = Field(Integer)
	no_of_loci = Field(BigInteger)
	min_depth = Field(Float)
	max_depth = Field(Float)
	max_missing_rate = Field(Float)
	min_maf = Field(Float)
	min_neighbor_distance = Field(Integer)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genotype_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('short_name', 'ref_ind_seq_id'))
	
	def constructRelativePath(self, subFolder='genotype_file', **keywords):
		"""
		2012.7.12
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		dst_relative_path = '%s/method_%s'%(subFolder, self.id)
		return dst_relative_path
	
	def getFileSize(self, data_dir=None):
		"""
		2012.7.12
		"""
		from pymodule import utils
		if data_dir is None:
			data_dir = VervetDB.get_data_dir()
		return utils.getFileOrFolderSize(os.path.join(data_dir, self.path))
	
class GenotypeFile(Entity, TableClass):
	"""
	2012.8.30 add column no_of_chromosomes
	2012.7.11
		modify it so that this holds files, each of which corresponds to GenotypeMethod.
	2011-2-4
		file format:
			locus.id	allele_order	allele_type	seq.id	score	target_locus
	"""
	#individual = ManyToOne('Individual', colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	path = Field(Text, unique=True)
	original_path = Field(Text)
	md5sum = Field(Text)	# unique=True
	file_size = Field(BigInteger)	#2012.7.12
	chromosome = Field(Text)
	no_of_chromosomes = Field(BigInteger, default=1)	#2012.8.30 BigInteger in case some huge number of contigs
	no_of_individuals = Field(Integer)
	no_of_loci = Field(BigInteger)
	format = Field(Text)
	genotype_method = ManyToOne('GenotypeMethod', colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(String(4096))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genotype_file', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('genotype_method_id', 'chromosome', 'format', 'no_of_chromosomes'))
	
	def constructRelativePath(self, subFolder='genotype_file', sourceFilename="", **keywords):
		"""
		2012.8.30
			name differently when the no_of_chromosomes is >1
		2012.7.12
			path relative to self.data_dir
		"""
		folderRelativePath = self.genotype_method.constructRelativePath(subFolder=subFolder)
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		folderRelativePath = folderRelativePath.lstrip('/')
		if self.no_of_chromosomes<=1:
			dst_relative_path = os.path.join(folderRelativePath, '%s_%s_%s'%(self.id, self.format, sourceFilename))
		else:	#files with more than 1 chromosome
			# 2012.8.30
			# it'll be something like genotype_file/method_6_$id_$format_#chromosomes_sourceFilename
			# do not put it in genotype_file/method_6/ folder.
			folderRelativePath= folderRelativePath.rstrip('/')
			dst_relative_path = '%s_%s_%s_%schromosomes_%s'%(folderRelativePath, self.id, self.format, self.no_of_chromosomes, sourceFilename)
		return dst_relative_path
	
	def getFileSize(self, data_dir=None):
		"""
		2012.7.12
		"""
		from pymodule import utils
		if data_dir is None:
			data_dir = VervetDB.get_data_dir()
		return utils.getFileOrFolderSize(os.path.join(data_dir, self.path))

class Phenotype(Entity, TableClass):
	individual = ManyToOne('Individual', colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	value = Field(Float)
	comment = Field(Text)
	replicate = Field(Integer, default=1)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual_id', 'replicate', 'phenotype_method_id'))

class BiologyCategory(Entity):
	#2011-5-4
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	phenotype_method_ls = OneToMany("PhenotypeMethod")
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='biology_category', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class PhenotypeMethod(Entity, TableClass):
	"""
	2011-5-4
		add biology_category, phenotype_scoring
	2011-3-1
		add user_ls/group_ls information
	"""
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	phenotype_scoring = Field(Text)
	biology_category = ManyToOne("BiologyCategory", colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	collector = ManyToOne("User",colname='collector_id')
	access = Field(Enum("public", "restricted", name="access_enum_type"), default='restricted')
	group_ls = ManyToMany('Group',tablename='group2phenotype_method', local_colname='phenotype_method_id', remote_colname='group_id')
	user_ls = ManyToMany('User',tablename='user2phenotype_method', local_colname='phenotype_method_id', remote_colname='user_id')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
	@classmethod
	def getPhenotypesForACL(cls, biology_category_id = None,user = None):
		#
		TableClass = PhenotypeMethod
		query = TableClass.query
		"""
		if biology_category_id is not None:
			query = query.filter(TableClass.biology_category_id==biology_category_id)
		"""
		clause = and_(TableClass.access == 'public')
		if user is not None:
			if user.isAdmin == 'Y':
				return query
			clause = or_(clause, TableClass.collector == user, TableClass.user_ls.any(User.id == user.id),
						TableClass.group_ls.any(Group.id.in_([group.id for group in user.group_ls])))
		query = query.filter(clause)
		return query
	
	def checkACL(self,user):
		if self.access == 'public':
			return True
		if user is None:
			return False
		if user.isAdmin == 'Y':
			return True
		if self.collector == user: 
			return True
		if user in self.user_ls:
			return True
		if [group in self.group_ls for group in user.group_ls]: 
			return True
		return False

class VervetDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['postgresql', 'v', 1, 'which type of database? mysql or postgresql', ],\
						('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
						('database', 1, ):['vervetdb', 'd', 1, '',],\
						('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
						('username', 1, ):[None, 'u', 1, 'database username',],\
						('password', 1, ):[None, 'p', 1, 'database password', ],\
						('port', 0, ):[None, 'o', 1, 'database port number'],\
						('pool_recycle', 0, int):[3600, 'l', 1, 'the length of time to keep connections open before recycling them.'],\
						('sql_echo', 0, bool):[False, 's', 0, 'display SQL Statements'],\
						('echo_pool', 0, bool):[False, 'e', 0, 'if True, the connection pool will log all checkouts/checkins to the logging stream, which defaults to sys.stdout.'],\
						('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-10-06
			add option 'pool_recycle' to recycle connection. MySQL typically close connections after 8 hours.
			__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle)
		2008-07-09
		"""
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
												class_to_have_attr=self)
		
		if self.echo_pool:	#2010-9-19 passing echo_pool to create_engine() causes error. all pool log disappeared.
			#2010-9-19 Set up a specific logger with our desired output level
			import logging
			#import logging.handlers
			logging.basicConfig()
			#LOG_FILENAME = '/tmp/sqlalchemy_pool_log.out'
			my_logger = logging.getLogger('sqlalchemy.pool')
	
			# Add the log message handler to the logger
			#handler = logging.handlers.RotatingFileHandler(LOG_FILENAME, maxBytes=1000000, backupCount=5)
			# create formatter
			#formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
			# add formatter to handler
			#handler.setFormatter(formatter)
			#my_logger.addHandler(handler)
			my_logger.setLevel(logging.DEBUG)
		"""
		ElixirDB.__init__(self, **keywords)
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
		self.uclaID2monkeyDBEntry = {}	#2012.9.4
		self.dbID2monkeyDBEntry = {}	#2012.9.6
		
		self.READMEClass = README	#2012.12.18 required to figure out data_dir
	
	def setup(self, create_tables=True):
		"""
		2008-09-07
			expose option create_tables, default=True. assign it to False if no new table is to be created.
		"""
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()
		#2008-08-26 setup_all() would setup other databases as well if they also appear in the program. Seperate this to be envoked after initialization
		# to ensure the metadata of other databases is setup properly.
	
	def getUniqueSequence(self, sequence):
		"""
		2011-2-11
		"""
		db_entry = Sequence.query.filter_by(sequence=sequence).first()
		if not db_entry:
			db_entry = Sequence(sequence=sequence, sequence_length=len(sequence))
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getAlignmentMethod(self, alignment_method_name=''):
		"""
		2011-5-6
		"""
		db_entry = AlignmentMethod.query.filter_by(short_name=alignment_method_name).first()
		if not db_entry:
			db_entry = AlignmentMethod(short_name=alignment_method_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getAlignment(self, individual_code=None, individual_id=None, individual=None, individual_sequence_id=None, \
					path_to_original_alignment=None, sequencer='GA', \
					sequence_type='SR', sequence_format='fastq', \
					ref_individual_sequence_id=10, \
					alignment_method_name='bwa-short-read', alignment_method_id=None, alignment_method=None,\
					alignment_format='bam', subFolder='individual_alignment', \
					createSymbolicLink=False, individual_sequence_filtered=0, read_group_added=None, dataDir=None, \
					outdated_index=0, mask_genotype_method_id=None, parent_individual_alignment_id=None,\
					individual_sequence_file_raw_id=None, md5sum=None):
		"""
		2012.9.19 add argument individual_sequence_file_raw_id, individual_id, individual,
			alignment_method_id, alignment_method, md5sum
		2012.7.26 add argument mask_genotype_method_id & parent_individual_alignment_id
		2012.6.13 add argument outdated_index
		2012.2.24
			add argument dataDir
		2011-9-15
			add argument read_group_added
		2011-8-30
			add argument individual_sequence_id
			use constructRelativePathForIndividualAlignment() to come up path
		2011-7-11
			add argument createSymbolicLink. default to False
				if True, create a symbolic from source file to target, instead of cp
		2011-5-6
			subFolder is the name of the folder in self.data_dir that is used to hold the alignment files.
		"""
		if dataDir is None:
			dataDir = self.data_dir
		if individual_code:
			individual = self.getIndividual(individual_code)
		elif individual_id:
			individual = Individual.get(individual_id)
		
		if individual_sequence_id:
			individual_sequence = IndividualSequence.get(individual_sequence_id)
			if individual is None:
				individual = individual_sequence.individual
		elif individual:
			individual_sequence = self.getIndividualSequence(individual_id=individual.id, sequencer=sequencer, \
								sequence_type=sequence_type,\
								sequence_format=sequence_format, filtered=individual_sequence_filtered)
		else:
			sys.stderr.write("Error: not able to get individual_sequence cuz individual_sequence_id=%s; individual_code=%s; individual_id=%s .\n"%\
							(individual_sequence_id, individual_code, individual_id))
			return None
		
		if alignment_method_name:
			alignment_method = self.getAlignmentMethod(alignment_method_name=alignment_method_name)
		elif alignment_method_id:
			alignment_method = AlignmentMethod.get(alignment_method_id)
		elif alignment_method is None:
			sys.stderr.write("Error: not able to get alignment_method cuz alignment_method_name=%s and alignment_method_id=%s.\n"%\
							(alignment_method_name, alignment_method_id))
			return None
		
		query = IndividualAlignment.query.filter_by(ind_seq_id=individual_sequence.id).\
				filter_by(ref_ind_seq_id=ref_individual_sequence_id).\
				filter_by(alignment_method_id=alignment_method.id).filter_by(outdated_index=outdated_index).\
				filter_by(mask_genotype_method_id=mask_genotype_method_id).\
				filter_by(parent_individual_alignment_id=parent_individual_alignment_id).\
				filter_by(individual_sequence_file_raw_id=individual_sequence_file_raw_id)
		
		if alignment_format:
			query = query.filter_by(format=alignment_format)
		db_entry = query.first()
		if not db_entry:
			db_entry = IndividualAlignment(ind_seq_id=individual_sequence.id, ref_ind_seq_id=ref_individual_sequence_id,\
								alignment_method_id=alignment_method.id, format=alignment_format, read_group_added=read_group_added,\
								outdated_index=outdated_index, mask_genotype_method_id=mask_genotype_method_id,\
								parent_individual_alignment_id=parent_individual_alignment_id,\
								individual_sequence_file_raw_id=individual_sequence_file_raw_id, md5sum=md5sum)
			self.session.add(db_entry)
			self.session.flush()
			#copy the file over
			
			if path_to_original_alignment and (os.path.isfile(path_to_original_alignment) or os.path.isdir(path_to_original_alignment)):
				from pymodule.utils import runLocalCommand
				#'/' must not be put in front of the relative path.
				# otherwise, os.path.join(dataDir, dst_relative_path) will only take the path of dst_relative_path.
				dst_relative_path = db_entry.constructRelativePath(subFolder=subFolder)
				#update its path in db to the relative path
				db_entry.path = dst_relative_path
				
				dst_pathname = os.path.join(dataDir, dst_relative_path)
				dst_dir = os.path.join(dataDir, subFolder)
				if not os.path.isdir(dst_dir):	#the upper directory has to be created at this moment.
					commandline = 'mkdir %s'%(dst_dir)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				if createSymbolicLink:
					commandline = 'ln -s %s %s'%(path_to_original_alignment, dst_pathname)
				else:
					commandline = 'cp -r %s %s'%(path_to_original_alignment, dst_pathname)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				
				#2011-7-11 copy the samtools index file as well
				indexFname = '%s.bai'%(path_to_original_alignment)
				if os.path.isfile(indexFname):
					dstIndexPathname = '%s.bai'%(dst_pathname)
					if createSymbolicLink:
						commandline = 'ln -s %s %s'%(indexFname, dstIndexPathname)
					else:
						commandline = 'cp -r %s %s'%(indexFname, dstIndexPathname)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				
				#since the .path has been updated, so add & flush
				self.session.add(db_entry)
				self.session.flush()
		return db_entry
	
	
	def getAlignments(self, ref_ind_seq_id=None, ind_seq_id_ls=None, ind_aln_id_ls=None, alignment_method_id=None, \
					dataDir=None, sequence_type=None, outdated_index=0, mask_genotype_method_id=None, \
					parent_individual_alignment_id=None, individual_sequence_file_raw_id_type=1,\
					country_id_ls=None, tax_id_ls=None, excludeAlignmentWithoutLocalFile=True):
		"""
		2012.11.29 added argument excludeAlignmentWithoutLocalFile, (exclude an alignment if it does not exist in local storage)
		2012.9.23 add argument country_id_ls and tax_id_ls
		2012.9.20 rename aln_method_id to alignment_method_id
		2012.9.19 & 2012.9.22 add argument individual_sequence_file_raw_id_type:
			1: only all-library-fused libraries, (individual_sequence_file_raw_id is null)
			2: only library-specific alignments, (individual_sequence_file_raw_id is non-null)
			3 or else: both all-library-fused and library-specific alignments, 
				(no filter based on individual_sequence_file_raw_id)
		2012.7.26 added argument mask_genotype_method_id & parent_individual_alignment_id
		2012.6.13 add argument outdated_index
		2012.4.13
			moved from AlignmentToCallPipeline.py
		2012.4.5
			select alignment using AND between all input arguments
		2011-11-27
			add argument sequence_type
		2011-9-16
			order each alignment by id. It is important because this is the order that gatk&samtools take input bams.
			#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
			# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		2011-9-13
			add argument dataDir, to filter out alignments that don't exist on file storage
		2011-8-31
			add argument aln_method_id
		2011-7-12
		
		"""
		if dataDir is None:
			dataDir = self.data_dir
		alignmentLs = []
		TableClass = IndividualAlignment
		query = TableClass.query
		if ind_aln_id_ls:
			sys.stderr.write("Adding filter via %s alignment IDs ... \n"%(len(ind_aln_id_ls)) )
			query = query.filter(TableClass.id.in_(ind_aln_id_ls))
		if ind_seq_id_ls:
			sys.stderr.write("Adding filter via %s sequence IDs  ... \n"%(len(ind_seq_id_ls)))
			query = query.filter(TableClass.ind_seq_id.in_(ind_seq_id_ls))
		if ref_ind_seq_id:
			sys.stderr.write("Adding filter via the reference sequence ID %s ... \n"%(ref_ind_seq_id))
			query = query.filter_by(ref_ind_seq_id=ref_ind_seq_id)
		if alignment_method_id:
			sys.stderr.write("Adding filter via the alignment_method_id=%s ... \n"%(alignment_method_id))
			query = query.filter_by(alignment_method_id=alignment_method_id)
		
		if country_id_ls:	#2012.9.23 not sure whether it'll work
			sys.stderr.write("Adding filter country_id_ls=%s ... \n"%(getattr(country_id_ls, '__len__', returnZeroFunc)()))
			query = query.filter(TableClass.individual_sequence.individual.site.has(Site.country_id.in_(country_id_ls)))
		if tax_id_ls:
			sys.stderr.write("Adding filter tax_id_ls=%s ... \n"%(getattr(tax_id_ls, '__len__', returnZeroFunc)()))
			query = query.filter(TableClass.individual_sequence.individual.has(Individual.tax_id.in_(tax_id_ls)))
		
		if not ind_aln_id_ls and not ind_seq_id_ls and not ref_ind_seq_id:
			sys.stderr.write("Both ind_seq_id_ls and ind_aln_id_ls are empty and ref_ind_seq_id is None. no alignment to be fetched.\n")
			sys.exit(3)
		
		if individual_sequence_file_raw_id_type==1:	#only all-library-fused alignments
			query = query.filter(TableClass.individual_sequence_file_raw_id==None)
		elif individual_sequence_file_raw_id_type==2:	#only library-specific alignments
			query = query.filter(TableClass.individual_sequence_file_raw_id!=None)
		else:
			#2012.9.22 do nothing = include both 
			pass
		
		#order by TableClass.id is important because this is the order that gatk&samtools take input bams.
		#Read group in each bam is beginned by alignment.id. GATK would arrange bams in the order of read groups.
		# while samtools doesn't do that and vcf-isect could combine two vcfs with columns in different order.
		
		query = query.filter_by(outdated_index=outdated_index).\
					filter_by(mask_genotype_method_id=mask_genotype_method_id).\
					filter_by(parent_individual_alignment_id=parent_individual_alignment_id).\
					order_by(TableClass.id)
		
		for row in query:
			if sequence_type is not None and row.individual_sequence.sequence_type!=sequence_type:
				continue
			if row.path:	#it's not None
				abs_path = os.path.join(dataDir, row.path)
				if excludeAlignmentWithoutLocalFile:
					if os.path.isfile(abs_path):
						alignmentLs.append(row)
				else:
					alignmentLs.append(row)
				
		sys.stderr.write("%s alignments Done.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def getProperAlignmentGivenIndividualID(self, ucla_id=None, individual_id=None, ref_ind_seq_id=524, alignment_method_id=2):
		"""
		2012.11.26
			Definition of proper alignment:
				1. not outdated
				2. from filtered reads.
				3. ref_ind_seq_id is the most recent reference (524 now)
				4. alignment_method_id is the consensus one(=2).
			
			query the view_alignment_with_country.
		"""
		query_string = "select * from view_alignment_with_country"
		where_condition_ls = ["filtered=1 and outdated_index=0 and ref_ind_seq_id=%s and alignment_method_id=%s "%\
							(ref_ind_seq_id, alignment_method_id)]
		if ucla_id:
			where_condition_ls.append("ucla_id='%s'"%(ucla_id))
		if individual_id:
			where_condition_ls.append("individual_id=%s"%(individual_id))
		query_string = "%s where %s "%(query_string, " and ".join(where_condition_ls))
		query = self.metadata.bind.execute(query_string)
		return query.fetchone()
	
	def getAlignmentsFromVCFSampleIDList(self, sampleIDList=None):
		"""
		2012.8.15
			each sample ID is constructed through IndividualAlignment.getReadGroup()
			use vcfFile.getSampleIDList() to generate sampleIDList. vcfFile is pymodule.VCFFile.
				vcfFile.sample_id_ls is not good because its index 0 is "ref".
		"""
		no_of_samples = 0
		if sampleIDList:
			no_of_samples = len(sampleIDList)
		sys.stderr.write("Getting alignments from %s samples ...\n"%(no_of_samples))
		alignmentLs = []
		if sampleIDList:
			for read_group in sampleIDList:
				individualAlignment = self.parseAlignmentReadGroup(read_group).individualAlignment
				alignmentLs.append(individualAlignment)
		sys.stderr.write("%s alignments.\n"%(len(alignmentLs)))
		return alignmentLs
	
	def getMonkeyID2ProperAlignment(self, ref_ind_seq_id=524, alignment_method_id=2, idType=1):
		"""
		2012.11.30
			Definition of proper alignment:
				1. not outdated
				2. from filtered reads.
				3. ref_ind_seq_id is the most recent reference (524 now)
				4. alignment_method_id is the consensus one(=2).
			
			query the view_alignment_with_country.
			
			idType
				1: ucla_id
				2: individual.id
				3: individual.code
		"""
		query_string = "select * from view_alignment_with_country"
		where_condition_ls = ["filtered=1 and outdated_index=0 and ref_ind_seq_id=%s and alignment_method_id=%s "%\
							(ref_ind_seq_id, alignment_method_id)]
		query_string = "%s where %s "%(query_string, " and ".join(where_condition_ls))
		query = self.metadata.bind.execute(query_string)
		monkeyID2ProperAlignment = {}
		for row in query:
			if idType ==1:
				monkeyID = row.ucla_id
			elif idType==2:
				monkeyID = row.individual_id
			else:
				monkeyID = row.code
			if monkeyID not in monkeyID2ProperAlignment:
				monkeyID2ProperAlignment[monkeyID] = row
			else:
				sys.stderr.write("Warning: monkey %s has >1 proper alignment. Only used the 1st one. (ref_ind_seq_id=%s, alignment_method_id=%s)\n"%\
								(ref_ind_seq_id, alignment_method_id))
		return monkeyID2ProperAlignment
	
	def filterAlignments(self, alignmentLs, max_coverage=None, individual_site_id=None, sequence_filtered=None,\
						individual_site_id_set=None, mask_genotype_method_id=None, parent_individual_alignment_id=None,\
						country_id_set=None, tax_id_set=None, excludeContaminant=False, excludeTissueIDSet=set([6]),\
						report=True):
		"""
		2012.10.2 add argument excludeTissueIDSet, default=6 (RNASASamples). Most alignments have either tissue_id=5 (ACD-blood) or null.
		2012.9.27 add argument excludeContaminant
		2012.9.22 add argument country_id_set, tax_id_set
		2012.7.26 added argument mask_genotype_method_id & parent_individual_alignment_id
		2012.5.8
			bugfix, individual_site_id_set could be none. so no_of_sites is un-defined.
		2012.4.13
			moved from AlignmentToCallPipeline.py
			become classmethod
			add argument individual_site_id_set
		2012.4.2
			add argument sequence_filtered
		2011-11-22
			447 in "individual_site_id=447" is VRC.
		"""
		if report:
			sys.stderr.write("Filter %s alignments to select individual_sequence.coverage <=%s & site-id=%s & sequence_filtered=%s & from %s sites & %s countries & %s taxonomies ..."%\
							(len(alignmentLs), max_coverage, individual_site_id, sequence_filtered, \
							getattr(individual_site_id_set, '__len__', returnZeroFunc)(),\
							getattr(country_id_set, '__len__', returnZeroFunc)(),\
							getattr(tax_id_set, '__len__', returnZeroFunc)(),\
							)
						)
		newAlignmentLs = []
		for alignment in alignmentLs:
			if not alignment:
				continue
			if max_coverage is not None and alignment.individual_sequence.coverage>max_coverage:
				continue
			if individual_site_id is not None and alignment.individual_sequence.individual.site_id!=individual_site_id:
				continue
			if sequence_filtered is not None and alignment.individual_sequence.filtered!=sequence_filtered:
				continue
			if individual_site_id_set and alignment.individual_sequence.individual.site_id not in individual_site_id_set:
				#2012.4.13
				continue
			if mask_genotype_method_id is not None and alignment.mask_genotype_method_id!=mask_genotype_method_id:
				continue
			if parent_individual_alignment_id is not None and alignment.parent_individual_alignment_id!=parent_individual_alignment_id:
				continue
			if country_id_set:
				if alignment.individual_sequence.individual.site is None:
					sys.stderr.write("Warning: alignment (id=%s, path=%s, %s) has no site.\n"%(alignment.id, alignment.path,\
																			alignment.individual_sequence.individual.code))
					continue
				elif (alignment.individual_sequence.individual.site is None or \
									alignment.individual_sequence.individual.site.country_id not in country_id_set):
					continue
			if tax_id_set:
				if alignment.individual_sequence.individual.tax_id is None:
					sys.stderr.write("Warning: alignment (id=%s, path=%s, %s) has no tax_id.\n"%(alignment.id, alignment.path,\
																	alignment.individual_sequence.individual.code))
					continue
				elif alignment.individual_sequence.individual.tax_id not in tax_id_set:
					continue
			if excludeContaminant and alignment.individual_sequence.individual.is_contaminated:	#2012.9.27
				continue
			if excludeTissueIDSet and alignment.individual_sequence.tissue_id in excludeTissueIDSet:	#2012.10.2
				continue
			newAlignmentLs.append(alignment)
		if report:
			sys.stderr.write(" kept %s alignments. Done.\n"%(len(newAlignmentLs)))
		return newAlignmentLs
	
	@classmethod
	def getCumulativeAlignmentMedianDepth(cls, alignmentLs=[], defaultSampleAlignmentDepth=10):
		"""
		2012.8.7
		"""
		sys.stderr.write("Getting cumulative median depth of %s alignments ..."%(len(alignmentLs)))
		cumulativeDepth = 0
		for alignment in alignmentLs:
			if alignment and alignment.median_depth is not None:
				medianDepth = alignment.median_depth
			else:
				medianDepth = defaultSampleAlignmentDepth
			cumulativeDepth += medianDepth
		sys.stderr.write("=%s.\n"%(cumulativeDepth))
		return cumulativeDepth
	
	def getAlleleType(self, allele_type_name):
		"""
		2011-2-11
		"""
		db_entry = AlleleType.query.filter_by(short_name=allele_type_name).first()
		if not db_entry:
			db_entry = AlleleType(short_name=allele_type_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividual(self, code=None, sex=None, age=None, age_cas=None, latitude=None, longitude=None, \
					altitude=None, ucla_id=None, site=None, site_name=None, city=None, stateprovince=None, country_name=None, \
					collection_date=None, collector=None, \
					approx_age_group_at_collection=None, tax_id=None, birthdate=None, vrc_founder=None, comment=None,\
					microchip_id=None):
		"""
		Examples:
			individual = db_vervet.getIndividual(ucla_id=monkey_id)
			individual = db_vervet.getIndividual(code=code)
			individual = db_vervet.getIndividual(code=code, ucla_id=ucla_id, site=None, site_name='YBK', country_name='Gambia')
			
		2012.12.6 in the case of creating a new individual, and site is None (to  be created on the fly),
			site_name and country_name must be not null.
			latitude, longitude, altitude is no longer part of table Individual.
		2012.5.29
			check Individual using code if ucla_id is missing
		2011-10-18
			add argument microchip_id
		2011-9-8
			add argument comment
		2011-5-6
			ucla_id should be unique.
			add birthdate
		2011-5-5
			take first letter of sex and upper case.
		2011-2-11
			code can't be None
		"""
		db_entry = None
		if ucla_id:	#ucla_id should be unique
			db_entry = Individual.query.filter_by(ucla_id=ucla_id).first()
		if not db_entry and code is not None:
			query = Individual.query.filter_by(code=code)
			if tax_id:
				query = query.filter_by(tax_id=tax_id)
			db_entry = query.first()
		if not db_entry:
			if sex is None or sex=='?' or sex=='':	#2011-4-29
				sex = None
			elif len(sex)>=1:	#2011-5-5 take the first letter
				sex = sex[0].upper()
			if not site and site_name and country_name:
				site = self.getSite(description=site_name, city=city, stateprovince=stateprovince, country_name=country_name,  \
								latitude=latitude,\
								longitude=longitude, altitude=altitude)
			db_entry = Individual(code=code, sex=sex, age=age, age_cas=age_cas, ucla_id=ucla_id, site=site, \
						collection_date=collection_date, collector=collector,\
						approx_age_group_at_collection=approx_age_group_at_collection, tax_id=tax_id, birthdate=birthdate,\
						vrc_founder=vrc_founder, comment=comment, microchip_id=microchip_id)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequence(self, individual_id=None, sequencer=None, sequence_type=None,\
						sequence_format=None, path_to_original_sequence=None, tissue_name=None, tissue_id=None, \
						coverage=None,\
						subFolder='individual_sequence', quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None, dataDir=None):
		"""
		2012.6.3
			columns that are None become part of the db query to see if entry is in db already
		2012.2.24
			add argument dataDir
		2012.2.10
			path_to_original_sequence is only given when you want to copy the file to db storage.
			add argument parent_individual_sequence_id
		2011-8-30
			add argument filtered
		2011-8-18
			add argument quality_score_format, default to "Standard"
		2011-8-3
			the path field is now considered a folder (rather than a file).
		2011-5-7
			subFolder is the name of the folder in self.data_dir that is used to hold the sequence files.
		"""
		if not dataDir:
			dataDir = self.data_dir
		
		query = IndividualSequence.query.filter_by(individual_id=individual_id)
		
		query = query.filter_by(sequencer=sequencer)
		query = query.filter_by(sequence_type=sequence_type)
		query = query.filter_by(format=sequence_format)
		
		if tissue_name:
			tissue = self.getTissue(short_name=tissue_name)
			query = query.filter_by(tissue_id=tissue.id)
		elif tissue_id:
			query = query.filter_by(tissue_id=tissue_id)
		else:
			query = query.filter_by(tissue_id=None)
		if parent_individual_sequence_id:
			query = query.filter_by(parent_individual_sequence_id=parent_individual_sequence_id)
		else:
			query = query.filter_by(parent_individual_sequence_id=None)
		query= query.filter_by(filtered=filtered)
		db_entry = query.first()
		if not db_entry:
			if tissue_name:
				tissue = self.getTissue(short_name=tissue_name)
			else:
				tissue = None
			individual = Individual.get(individual_id)
			db_entry = IndividualSequence(individual_id=individual_id, sequencer=sequencer, sequence_type=sequence_type,\
									format=sequence_format, tissue=tissue, coverage=coverage, \
									quality_score_format=quality_score_format, filtered=filtered,\
									parent_individual_sequence_id=parent_individual_sequence_id)
			#to make db_entry.id valid
			self.session.add(db_entry)
			self.session.flush()
			
			dst_relative_path = db_entry.constructRelativePathForIndividualSequence(subFolder=subFolder)
			
			#update its path in db to the relative path
			db_entry.path = dst_relative_path
			
			dst_abs_path = os.path.join(dataDir, dst_relative_path)
			from pymodule.utils import runLocalCommand
			if path_to_original_sequence and (os.path.isfile(path_to_original_sequence) or os.path.isdir(path_to_original_sequence)):
				dst_dir = os.path.join(dataDir, subFolder)
				if not os.path.isdir(dst_dir):	#the upper directory has to be created at this moment.
					commandline = 'mkdir -p %s'%(dst_dir)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				if not os.path.isdir(dst_abs_path):	#2011-8-3 create the directory to host all sequences.
					commandline = 'mkdir -p %s'%(dst_abs_path)
					return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
				commandline = 'cp -r %s %s'%(path_to_original_sequence, dst_abs_path)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualDBEntry(self, ucla_id=None):
		"""
		2012.9.4
			copied from vervet/src/plot/PlotPedigreeKinshipVsGeneticIBD.py
		2012.8.21
			has a cache inside
		"""
		if ucla_id in self.uclaID2monkeyDBEntry:
			return self.uclaID2monkeyDBEntry.get(ucla_id)
		else:
			monkey = Individual.query.filter_by(ucla_id=ucla_id).first()
			self.uclaID2monkeyDBEntry[ucla_id] = monkey
			return monkey

	def getIndividualDBEntryViaDBID(self, db_id=None):
		"""
		2012.9.6
			has a cache inside
		"""
		if db_id in self.dbID2monkeyDBEntry:
			return self.dbID2monkeyDBEntry.get(db_id)
		else:
			monkey = Individual.get(db_id)
			self.dbID2monkeyDBEntry[db_id] = monkey
			return monkey
	
	def copyParentIndividualSequence(self, parent_individual_sequence=None, parent_individual_sequence_id=None,\
									subFolder='individual_sequence', quality_score_format='Standard', filtered=1,\
									dataDir=None):
		"""
		2012.6.10
			call getIndividualSequence to construct individual_sequence, rather than construct it from here.
		2012.2.10
			constructRelativePathForIndividualSequence is now moved to Table IndividualSequence
			add argument filtered
		2011-9-1
			add quality_score_format
		2011-8-11
			make a copy of parent_individual_sequence_id individual_sequence entity
		"""
		if parent_individual_sequence is None:
			parent_individual_sequence = IndividualSequence.get(parent_individual_sequence_id)
		
		pis = parent_individual_sequence
		
		individual = Individual.get(pis.individual_id)
		
		individual_sequence = self.getIndividualSequence(individual_id=parent_individual_sequence.individual_id, \
						sequencer=parent_individual_sequence.sequencer, sequence_type=parent_individual_sequence.sequence_type,\
						sequence_format=parent_individual_sequence.format, path_to_original_sequence=None, \
						tissue_id=getattr(parent_individual_sequence.tissue, 'id', None), coverage=None,\
						quality_score_format=quality_score_format, filtered=filtered,\
						parent_individual_sequence_id=parent_individual_sequence.id, \
						dataDir=dataDir, subFolder=subFolder)
		"""
		db_entry = IndividualSequence(individual_id=pis.individual_id, sequencer=pis.sequencer, sequence_type=pis.sequence_type,\
							format=pis.format, tissue=pis.tissue, quality_score_format=quality_score_format)
		
		db_entry.filtered = filtered
		db_entry.parent_individual_sequence = pis
		#to make db_entry.id valid
		self.session.add(db_entry)
		self.session.flush()
		
		dst_relative_path = db_entry.constructRelativePathForIndividualSequence(subFolder=subFolder)
		
		#update its path in db to the relative path
		db_entry.path = dst_relative_path
		self.session.add(db_entry)
		self.session.flush()
		"""
		return individual_sequence
	
	def copyParentIndividualSequenceFile(self, parent_individual_sequence_file=None, parent_individual_sequence_file_id=None,\
									individual_sequence_id=None,\
									quality_score_format='Standard', filtered=1):
		"""
		2012.2.14
			call self.getIndividualSequenceFile() instead
		2012.2.10
		"""
		if parent_individual_sequence_file is None:
			if parent_individual_sequence_file_id is not None:
				parent_individual_sequence_file = IndividualSequenceFile.get(parent_individual_sequence_file_id)
			else:
				sys.stderr.write("Warning: parent individual_sequence_id is None. No IndividualSequenceFile instance to be created.\n")
				return None
		
		
		parent = parent_individual_sequence_file
		
		
		db_entry = self.getIndividualSequenceFile(individual_sequence_id, library=parent.library, mate_id=parent.mate_id, \
									split_order=parent.split_order, format=parent.format,\
									filtered=filtered, parent_individual_sequence_file_id=parent.id, \
									individual_sequence_file_raw_id=parent.individual_sequence_file_raw_id,\
									quality_score_format=quality_score_format)
		return db_entry
	
	def copyParentIndividualAlignment(self, parent_individual_alignment=None, \
									parent_individual_alignment_id=None, mask_genotype_method=None, \
									mask_genotype_method_id=None, dataDir=None):
		"""
		2012.7.28
		
		"""
		if parent_individual_alignment is None and parent_individual_alignment_id:
			parent_individual_alignment = IndividualAlignment.get(parent_individual_alignment_id)
		if mask_genotype_method is None and mask_genotype_method_id:
			mask_genotype_method = GenotypeMethod.get(mask_genotype_method_id)
		
		individual_sequence = parent_individual_alignment.individual_sequence
		ref_sequence = parent_individual_alignment.ref_sequence
		alignment_method = parent_individual_alignment.alignment_method
		individual_alignment = self.getAlignment(individual_code=individual_sequence.individual.code, \
								individual_sequence_id=individual_sequence.id,\
					path_to_original_alignment=None, sequencer=individual_sequence.sequencer, \
					sequence_type=individual_sequence.sequence_type, sequence_format=individual_sequence.format, \
					ref_individual_sequence_id=ref_sequence.id, \
					alignment_method_name=alignment_method.short_name, alignment_format=parent_individual_alignment.format,\
					individual_sequence_filtered=individual_sequence.filtered, read_group_added=1,
					dataDir=dataDir, outdated_index=0, mask_genotype_method_id=mask_genotype_method.id,\
					parent_individual_alignment_id=parent_individual_alignment.id,\
					individual_sequence_file_raw_id=parent_individual_alignment.individual_sequence_file_raw_id)
					#read-group addition is part of pipeline
		#if not individual_alignment.path:
		#	individual_alignment.path = individual_alignment.constructRelativePath()
		#	session.add(individual_alignment)
		#	session.flush()
		return individual_alignment
	
	def isPathInDBAffiliatedStorage(self, db_entry=None, relativePath=None, inputFileBasename=None, dataDir=None, \
								constructRelativePathFunction=None):
		"""
		2012.8.29
			check whether one relative path already exists in db-affliated storage.
				return -1 if db_entry.path is different from relativePath and could not be updated in db.
				return 1 if relativePath exists in db already
				return 0 if not.
			Example:
				# simplest
				db_vervet.isPathInDBAffiliatedStorage(relativePath=relativePath, dataDir=self.dataDir)
				# this will check if db_entry.path == relativePath, if not, update it in db with relativePath.
				db_vervet.isPathInDBAffiliatedStorage(db_entry=db_entry, relativePath=relativePath)
				# 
				db_vervet.isPathInDBAffiliatedStorage(db_entry=db_entry, inputFileBasename=inputFileBasename, dataDir=None, \
								constructRelativePathFunction=genotypeFile.constructRelativePath)
		"""
		if dataDir is None:
			dataDir = self.data_dir
		
		exitCode = 0
		if db_entry and (relativePath is None) and constructRelativePathFunction and inputFileBasename:
			relativePath = constructRelativePathFunction(db_entry=db_entry, sourceFilename=inputFileBasename)
		
		if db_entry and relativePath:
			if db_entry.path != relativePath:
				db_entry.path = relativePath
				try:
					self.session.add(db_entry)
					self.session.flush()
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					import traceback
					traceback.print_exc()
					exitCode = -1
					return exitCode
		
		dstFilename = os.path.join(dataDir, relativePath)
		if os.path.isfile(dstFilename):
			exitCode = 1
		else:
			exitCode = 0
		return exitCode
	
	def copyFileWithAnotherFilePrefix(self, inputFname=None, filenameWithPrefix=None, outputDir=None,\
									logMessage=None, srcFilenameLs=None, dstFilenameLs=None):
		"""
		2012.9.20
		"""
		srcFilename = inputFname
		prefix, suffix = os.path.splitext(os.path.basename(inputFname))
		newPrefix = os.path.splitext(filenameWithPrefix)[0]
		dstFilename = os.path.join(outputDir, '%s.%s'%(newPrefix, suffix))
		utils.copyFile(srcFilename=srcFilename, dstFilename=dstFilename)
		if logMessage:
			logMessage += "file %s has been copied to %s.\n"%(srcFilename, dstFilename)
		if srcFilenameLs:
			srcFilenameLs.append(srcFilename)
		if dstFilenameLs:
			dstFilenameLs.append(dstFilename)
		return logMessage
	
	def getIndividualSequenceFileRaw(self, individual_sequence_id=None, library=None, md5sum=None, path=None, mate_id=None,\
									file_size=None):
		"""
		2012.7.12 add argument file_size
		2012.4.30
			add argument mate_id
		2012.2.14
		"""
		#query first
		import os
		query = IndividualSequenceFileRaw.query.filter_by(individual_sequence_id=individual_sequence_id)
		if library:
			query = query.filter_by(library=library)
		if md5sum:
			query = query.filter_by(md5sum=md5sum)
		if path:
			path = os.path.realpath(path)	#get the realpath
			query = query.filter_by(path=path)
		if mate_id:
			query = query.filter_by(mate_id=mate_id)
		db_entry = query.first()
		if not db_entry:
			if path and os.path.isfile(path) and file_size is None:	#2012.7.12
				from pymodule import utils
				file_size = utils.getFileOrFolderSize(path)
			db_entry = IndividualSequenceFileRaw(individual_sequence_id=individual_sequence_id, library=library, md5sum=md5sum, \
										path=path, mate_id=mate_id, file_size=file_size)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequenceFile(self, individual_sequence_id, library=None, mate_id=None, split_order=None, format=None,\
							filtered=0, parent_individual_sequence_file_id=None, individual_sequence_file_raw_id=None,\
							quality_score_format='Standard'):
		"""
		2012.5.2
			add all filter_by() straight to IndividualSequenceFile.query
		2012.2.14
		"""
		#query first
		query = IndividualSequenceFile.query.filter_by(individual_sequence_id=individual_sequence_id).filter_by(library=library).\
					filter_by(mate_id=mate_id).filter_by(split_order=split_order).filter_by(format=format).filter_by(filtered=filtered).\
					filter_by(parent_individual_sequence_file_id=parent_individual_sequence_file_id).\
					filter_by(individual_sequence_file_raw_id=individual_sequence_file_raw_id)
		db_entry = query.first()
		if not db_entry:
			db_entry = IndividualSequenceFile(individual_sequence_id=individual_sequence_id,\
							library=library, mate_id=mate_id, split_order=split_order,\
							format=format, filtered=filtered, \
							parent_individual_sequence_file_id=parent_individual_sequence_file_id,\
							individual_sequence_file_raw_id=individual_sequence_file_raw_id,\
							quality_score_format=quality_score_format)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getRelationshipType(self, relationship_type_name=None):
		"""
		2011-5-5
			fetch (and add) a RelationshipType entry
		"""
		query = RelationshipType.query.filter_by(short_name=relationship_type_name)
		db_entry = query.first()
		if not db_entry:
			db_entry = RelationshipType(short_name=relationship_type_name)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def checkSpecificRelationOfIndividual2(self, individual2=None, relationship_type_name=None):
		"""
		2012.1.23
			check to see if specific relationship of individual2 exists in Ind2Ind already.
			This is to detect errors where one individual has >1 father/mother in the database.
		"""
		relationship_type = self.getRelationshipType(relationship_type_name)
		query = Ind2Ind.query.filter_by(individual2_id=individual2.id).filter_by(relationship_type_id=relationship_type.id)
		db_entry = query.first()
		if not db_entry:
			return None
		else:
			return db_entry
	
	def getInd2Ind(self, individual1=None, individual2=None, relationship_type_name=None):
		"""
		2011-5-5
			fetch (and add) a Ind2Ind entry
		"""
		relationship_type = self.getRelationshipType(relationship_type_name)
		query = Ind2Ind.query.filter_by(individual1_id=individual1.id).filter_by(individual2_id=individual2.id).\
			filter_by(relationship_type_id=relationship_type.id)
		db_entry = query.first()
		if not db_entry:
			db_entry = Ind2Ind(individual1=individual1, individual2=individual2, relationship_type=relationship_type)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry

	def getGenotypeMethod(self, short_name=None, description=None, ref_ind_seq_id=None, \
						individualAlignmentLs=None, parent_genotype_method=None, parent_id=None, \
						min_depth=None, max_depth=None, max_missing_rate=None, min_maf=None,\
						no_of_individuals=None, no_of_loci=None, dataDir=None, subFolder='genotype_file'):
		"""
		2012.7.12
			examples:
				genotypeMethod = self.db_vervet.getGenotypeMethod(short_name=self.genotypeMethodShortName, \
														individualAlignmentLs=individualAlignmentLs)
		"""
		if not short_name:
			sys.stderr.write("Error: short_name (%s) is empty.\n"%(short_name))
			return None
		
		if ref_ind_seq_id is None and individualAlignmentLs:
			firstAlignment = individualAlignmentLs[0]
			ref_ind_seq_id =firstAlignment.ref_ind_seq_id
			
		if parent_genotype_method and parent_id is None:
			parent_id = parent_genotype_method.id
		
		query = GenotypeMethod.query.filter_by(short_name=short_name)
		
		if ref_ind_seq_id:
			query = query.filter_by(ref_ind_seq_id=ref_ind_seq_id)
		if parent_id:
			query = query.filter_by(parent_id=parent_id)
		
		
		db_entry = query.first()
		if not db_entry:
			db_entry = GenotypeMethod(short_name=short_name, description=description, ref_ind_seq_id=ref_ind_seq_id,\
								parent_id=parent_id, min_depth=min_depth, max_depth=max_depth, \
								max_missing_rate=max_missing_rate, min_maf=min_maf,\
								no_of_individuals=no_of_individuals, no_of_loci=no_of_loci)
			db_entry.individual_alignment_ls = individualAlignmentLs
			self.session.add(db_entry)
			self.session.flush()
			from pymodule.utils import runLocalCommand
			if not dataDir:
				dataDir = self.data_dir
			
			db_entry.path = db_entry.constructRelativePath(subFolder=subFolder)
			genotypeFileUpperFolder = os.path.join(dataDir, db_entry.path)
			
			if not os.path.isdir(genotypeFileUpperFolder):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(genotypeFileUpperFolder)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
		
	def getGenotypeFile(self, genotype_method=None, genotype_method_id=None, chromosome=None, format=None, path=None, file_size=None, \
					no_of_individuals=None, no_of_loci=None, md5sum=None, original_path=None, dataDir=None, subFolder='genotype_file',\
					no_of_chromosomes=1):
		"""
		2012.8.30
			add argument no_of_chromosomes
		2012.7.13
		"""
		if genotype_method_id is None and genotype_method.id:
			genotype_method_id = genotype_method.id
			
		query = GenotypeFile.query.filter_by(genotype_method_id=genotype_method_id)
		if format:
			query = query.filter_by(format=format)
		if md5sum:	#2012.8.6 if format is not given, then use md5sum as the sole unique identifier 
			query = query.filter_by(md5sum=md5sum)
		if no_of_chromosomes:
			query = query.filter_by(no_of_chromosomes=no_of_chromosomes)
		if chromosome:
			query = query.filter_by(chromosome=chromosome)
		
		db_entry = query.first()
		if original_path:
			original_path = os.path.realpath(original_path)
		
		if not db_entry:
			
			db_entry = GenotypeFile(format=format, path=path, no_of_individuals=no_of_individuals, no_of_loci=no_of_loci, \
							chromosome=chromosome, md5sum=md5sum, file_size=file_size, \
							original_path=original_path, genotype_method_id=genotype_method_id, no_of_chromosomes=no_of_chromosomes)
			if not dataDir:
				dataDir = self.data_dir
			if db_entry.file_size is None and db_entry.path:
				db_entry.file_size = db_entry.getFileSize(data_dir=dataDir)
			self.session.add(db_entry)
			self.session.flush()
			
			from pymodule.utils import runLocalCommand
			
			genotypeFileUpperFolder = os.path.join(dataDir, db_entry.genotype_method.constructRelativePath(subFolder=subFolder))
			
			if not os.path.isdir(genotypeFileUpperFolder):	#the upper directory has to be created at this moment.
				commandline = 'mkdir -p %s'%(genotypeFileUpperFolder)
				return_data = runLocalCommand(commandline, report_stderr=True, report_stdout=True)
			
		return db_entry
	
	def getLocus(self, chr=None, start=None, stop=None, ref_seq=None, alt_seq=None, ref_ind_seq_id=None, \
				locus_type_id=None, locus_type=None, ref_seq_id=None, alt_seq_id=None):
		"""
		2012.5.2
			add ref_ind_seq_id & locus_type_id, locus_type
		2011-2-11
		
		"""
		if ref_seq_id is None and ref_seq is not None:
			ref_seq_id = ref_seq.id
		if alt_seq_id is None and alt_seq is not None:
			alt_seq_id = alt_seq.id
		if locus_type_id is None and locus_type is not None:
			locus_type_id = locus_type.id
		
		query = Locus.query.filter_by(chromosome=chr).filter_by(start=start).filter_by(stop=stop).filter_by(ref_seq_id=ref_seq_id).\
					filter_by(alt_seq_id=alt_seq_id).filter_by(locus_type_id=locus_type_id).\
					filter_by(ref_ind_seq_id=ref_ind_seq_id)
					
		db_entry = query.first()
		if not db_entry:
			db_entry = Locus(chromosome=chr, start=start, stop=stop, ref_ind_seq_id=ref_ind_seq_id, locus_type_id=locus_type_id,\
							ref_seq_id=ref_seq_id, alt_seq_id=alt_seq_id)
			if locus_type:
				db_entry.locus_type = locus_type
			if ref_seq:
				db_entry.ref_seq = ref_seq
			if alt_seq:
				db_entry.alt_seq = alt_seq
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getLocusContext(self, locus_id=None, gene_id=None, gene_strand=None, disp_pos=None, overlap_length=None,\
					overlap_fraction_in_locus=None, overlap_fraction_in_gene=None):
		"""
		2012.5.14
		"""
		locus_context = LocusContext.query.filter_by(locus_id=locus_id).filter_by(gene_id=gene_id).first()
		
		if locus_context:
			already_in_db = 1
		else:
			locus_context = LocusContext(locus_id=locus_id, gene_id = gene_id, \
								gene_strand=gene_strand, disp_pos=disp_pos, \
								overlap_length=overlap_length, \
								overlap_fraction_in_locus=overlap_fraction_in_locus, \
								overlap_fraction_in_gene=overlap_fraction_in_gene)
			self.session.add(locus_context)
			already_in_db = 0
		return locus_context
	
	def getLocusAnnotation(self, locus_id=None, locus_context_id=None, locus_context=None, gene_id=None, gene_commentary_id=None, \
						gene_segment_id=None, locus_annotation_type =None, locus_annotation_type_id=None,\
						which_exon_or_intron = None, pos_within_codon=None, which_codon=None, \
						label=None, utr_number=None, cds_number=None, intron_number=None,\
						exon_number=None, overlap_length=None, overlap_fraction_in_locus=None, overlap_fraction_in_gene=None,\
						comment=None):
		"""
		2012.5.14
		"""
		if locus_context_id is None and locus_context:
			locus_context_id = locus_context.id
		if locus_annotation_type_id is None and locus_annotation_type:
			locus_annotation_type_id = locus_annotation_type.id
		locus_annotation = LocusAnnotation.query.filter_by(locus_id=locus_id).filter_by(locus_context_id=locus_context.id).\
								filter_by(gene_id=gene_id).\
								filter_by(gene_commentary_id= gene_commentary_id).\
								filter_by(gene_segment_id= gene_segment_id).\
								filter_by(locus_annotation_type_id=locus_annotation_type_id).first()
		if not locus_annotation:
			locus_annotation = LocusAnnotation(locus_id=locus_id, locus_context_id=locus_context_id, gene_id=gene_id,\
							gene_commentary_id = gene_commentary_id, \
							gene_segment_id=gene_segment_id, \
							locus_annotation_type=locus_annotation_type, locus_annotation_type_id=locus_annotation_type_id,\
							which_exon_or_intron=which_exon_or_intron, pos_within_codon=pos_within_codon, \
							which_codon=which_codon, label=label, \
							utr_number = utr_number, cds_number = cds_number, \
							intron_number = intron_number, exon_number = exon_number,\
							overlap_length=overlap_length, \
							overlap_fraction_in_locus=overlap_fraction_in_locus, overlap_fraction_in_gene=overlap_fraction_in_gene,\
							comment=comment)
			locus_annotation.locus_context = locus_context
			self.session.add(locus_annotation)
		return locus_annotation
	
	def getLocusAnnotationType(self, locus_annotation_type_short_name=None):
		"""
		2012.5.16
		
		"""
		ty = LocusAnnotation.query.filter_by(short_name=locus_annotation_type_short_name).first()
		if not ty:
			ty = LocusAnnotationType(short_name=locus_annotation_type_short_name)
			self.session.add(ty)
			self.session.flush()
		return ty
	
	def getSequence(self, sequence=None, comment=None):
		"""
		2012.5.2
		
		"""
		db_entry = Sequence.query.filter_by(sequence=sequence).first()
		if not db_entry:
			db_entry = Sequence(sequence=sequence, comment=comment)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSequenceBatch(self, short_name=None, description=None):
		"""
		2012.7.5
		"""
		db_entry = SequenceBatch.query.filter_by(short_name=short_name).first()
		if not db_entry:
			db_entry = SequenceBatch(short_name=short_name, description=description)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSite(self, short_name=None, description=None, city=None, stateprovince=None, country_name=None, \
			latitude=None, longitude=None, altitude=None):
		"""
		2012.12.6 added argument short_name
		2011-4-28
		"""
		country = self.getCountry(country_name=country_name)
		if short_name :
			short_name = short_name
		elif description:
			short_name = description
		elif city:
			short_name = city
		else:
			short_name = None
		query = Site.query.filter_by(country_id=country.id)
		if short_name:
			query = query.filter_by(short_name=short_name)
		if city:
			query = query.filter_by(city=city)
		if stateprovince:
			query = query.filter_by(stateprovince=stateprovince)
		if latitude:
			query = query.filter_by(latitude=latitude)
		if longitude:
			query = query.filter_by(longitude=longitude)
		if altitude:
			query = query.filter_by(altitude=altitude)
		db_entry = query.first()
		if not db_entry:
			db_entry = Site(short_name=short_name, description=description, city=city, stateprovince=stateprovince, \
						country=country, latitude=latitude, longitude=longitude, altitude=altitude)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getTissue(self, short_name=None, description=None):
		"""
		2011-5-9
			fetch a tissue entry (create one if none existent)
		"""
		db_entry = Tissue.query.filter_by(short_name=short_name).first()
		if not db_entry:
			db_entry = Tissue(short_name=short_name, description=description)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getCountry(self, country_name=None, capital=None, abbr=None, latitude=None, longitude=None):
		"""
		2011-4-28
		
		"""
		db_entry = Country.query.filter_by(name=country_name).first()
		if not db_entry:
			db_entry = Country(name=country_name, capital=capital, abbr=abbr, latitude=latitude, longitude=longitude)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getUser(self, name=None, username=None):
		"""
		2011-4-28
		
		"""
		name = unicode(name)
		query = User.query.filter_by(realname=name)
		if username:
			query = query.filter_by(username=username)
		db_entry = query.first()
		if not db_entry:
			db_entry = User(realname=name, username=username, _password="123456")
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getPhenotype(self, phenotype_name=None, value=None, replicate=None, individual_id=None, comment=None, collector_name=None):
		"""
		2011-4-28
		"""
		phenotype_method = self.getPhenotypeMethod(phenotype_name=phenotype_name, collector_name=collector_name)
		if replicate is None:
			replicate = 1
		db_entry = Phenotype.query.filter_by(phenotype_method_id=phenotype_method.id).filter_by(individual_id=individual_id).\
			filter_by(replicate=replicate).first()
		if not db_entry:
			db_entry = Phenotype(phenotype_method=phenotype_method, value=value, replicate=replicate,\
								individual_id=individual_id, comment=comment)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
		
		
	def getPhenotypeMethod(self, phenotype_name=None, collector_name=None):
		"""
		2011-4-28
		"""
		db_entry = PhenotypeMethod.query.filter_by(short_name=phenotype_name).first()
		if not db_entry:
			if collector_name:
				collector = self.getUser(name=collector_name)
			else:
				collector = None
			db_entry = PhenotypeMethod(short_name=phenotype_name, collector=collector)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequenceID2FilePairLs(self, individualSequenceIDList=None, dataDir=None, needPair=True, checkOldPath=False):
		"""
		2012.2.10
			add argument checkOldPath.
				True: find files in IndividualSequence.path[:-6], (=path without the trailing '_split').
					This is the old format.
				False: find files in IndividualSequence.path 
		2011-8-30
			filename in individualSequenceID2FilePairLs is path relative to dataDir
		2011-8-28
			add argument needPair
			copied from MpiBaseCount.py
		2011-8-5
		"""
		sys.stderr.write("Getting individualSequenceID2FilePairLs ...")
		individualSequenceID2FilePairLs = {}
		if not dataDir:
			dataDir = self.data_dir
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = IndividualSequence.get(individualSequenceID)
			if individual_sequence and individual_sequence.path:
				if checkOldPath:
					path = individual_sequence.path[:-6]
				else:
					path = individual_sequence.path
				abs_path = os.path.join(dataDir, path)
				if individual_sequence.id not in individualSequenceID2FilePairLs:
					individualSequenceID2FilePairLs[individual_sequence.id] = []
				if os.path.isfile(abs_path):
					fileRecord = [path, individual_sequence.format, 'SR', individual_sequence.sequencer]
						#"SR" means it's single-end
					individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord])
				elif os.path.isdir(abs_path):	#it's a folder, sometimes it's nothing there
					if individual_sequence.sequence_type=='PE':
						isPE = True
					else:
						isPE = False
					pairedEndPrefix2FileLs = NextGenSeq.getPEInputFiles(abs_path, isPE=isPE)
					for pairedEndPrefix, fileLs in pairedEndPrefix2FileLs.iteritems():
						if isPE and len(fileLs)==2 and fileLs[0] and fileLs[1]:	#PE
							filename = os.path.join(path, fileLs[0])	#take one file only
							fileRecord = [filename, individual_sequence.format, 'PE', individual_sequence.sequencer]
							#"PE" means it's paired-end
							filename2 = os.path.join(path, fileLs[1])	#take one file only
							fileRecord2 = [filename2, individual_sequence.format, 'PE', individual_sequence.sequencer]
							#"PE" means it's paired-end
							individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord, fileRecord2])	#"PE" means it's paired-end
						else:
							for filename in fileLs:	#usually should be only one file
								if filename:
									filename = os.path.join(path, filename)
									fileRecord = [filename, individual_sequence.format, 'SR', individual_sequence.sequencer]
									#"SR" means it's single-end
									individualSequenceID2FilePairLs[individual_sequence.id].append([fileRecord])
		sys.stderr.write("%s individual sequences. Done.\n"%(len(individualSequenceID2FilePairLs)))
		return individualSequenceID2FilePairLs
	
	def getISQ_ID2LibrarySplitOrder2FileLs(self, individualSequenceIDList=None, dataDir=None, filtered=None, ignoreEmptyReadFile=True):
		"""
		2012.3.19
			add argument, ignoreEmptyReadFile
		2012.2.24
			filtered=None means "no filtering based on this field.".
			
		2012.2.10
			If for one (isq_id, librarySplitOrder), there is only one mate (single-end).
				The isq_id2LibrarySplitOrder2FileLs only stores one file object (FileLs is of length 1).
			Length of FileLs is commesurate with the number of ends.
		"""
		sys.stderr.write("Getting isq_id2LibrarySplitOrder2FileLs for %s isq entries ..."%(len(individualSequenceIDList)))
		isq_id2LibrarySplitOrder2FileLs = {}
		if not dataDir:
			dataDir = self.data_dir
		counter = 0
		from mapper.CountFastqReadBaseCount import CountFastqReadBaseCount
		
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = IndividualSequence.get(individualSequenceID)
			if not individual_sequence:	#not present in db, ignore
				continue
			for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
				path = os.path.join(dataDir, individual_sequence_file.path)
				if filtered is not None and individual_sequence_file.filtered!=filtered:	#skip entries that don't matched the filtered argument
					continue
				if ignoreEmptyReadFile:	#2012.3.19	ignore empty read files.
					if individual_sequence_file.read_count is None:	#calculate it on the fly
						baseCountData = CountFastqReadBaseCount.getReadBaseCount(path, onlyForEmptyCheck=True)
						read_count = baseCountData.read_count
					else:
						read_count = individual_sequence_file.read_count
					if read_count==0:
						continue
				if individualSequenceID not in isq_id2LibrarySplitOrder2FileLs:
					isq_id2LibrarySplitOrder2FileLs[individualSequenceID] = {}
				counter += 1
				LibrarySplitOrder2FileLs = isq_id2LibrarySplitOrder2FileLs[individualSequenceID]
				library = individual_sequence_file.library
				split_order = individual_sequence_file.split_order
				mate_id = individual_sequence_file.mate_id
				if mate_id is None:
					mate_id = 1
				key = (library, split_order)
				if key not in LibrarySplitOrder2FileLs:
					LibrarySplitOrder2FileLs[key] = []
				if len(LibrarySplitOrder2FileLs[key])<mate_id:
					for i in xrange(mate_id-len(LibrarySplitOrder2FileLs[key])):	#expand the list to match the number of mates
						LibrarySplitOrder2FileLs[key].append(None)
				isq_file_obj = PassingData(db_entry=individual_sequence_file, path=path)
				LibrarySplitOrder2FileLs[key][mate_id-1] = isq_file_obj
		
		sys.stderr.write("%s individual sequence files from %s isq entries.\n"%(counter, len(isq_id2LibrarySplitOrder2FileLs)))
		return isq_id2LibrarySplitOrder2FileLs
	
	def getISQDBEntryLsForAlignment(self, individualSequenceIDList=None, dataDir=None, filtered=None, ignoreEmptyReadFile=True):
		"""
		2012.9.19 similar to getISQ_ID2LibrarySplitOrder2FileLs()
			isqLs
				library2Data
					isqFileRawID2Index
					isqFileRawDBEntryLs
					splitOrder2Index
					fileObjectPairLs
		"""
		sys.stderr.write("Getting isqLs for %s isq entries ..."%(len(individualSequenceIDList)))
		isqLs = []
		if not dataDir:
			dataDir = self.data_dir
		counter = 0
		from mapper.CountFastqReadBaseCount import CountFastqReadBaseCount
		
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = IndividualSequence.get(individualSequenceID)
			if not individual_sequence:	#not present in db, ignore
				continue
			library2Data = {}	#2012.9.19
			for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
				path = os.path.join(dataDir, individual_sequence_file.path)
				if filtered is not None and individual_sequence_file.filtered!=filtered:	#skip entries that don't matched the filtered argument
					continue
				if ignoreEmptyReadFile:	#2012.3.19	ignore empty read files.
					if individual_sequence_file.read_count is None:	#calculate it on the fly
						baseCountData = CountFastqReadBaseCount.getReadBaseCount(path, onlyForEmptyCheck=True)
						read_count = baseCountData.read_count
					else:
						read_count = individual_sequence_file.read_count
					if read_count==0:
						continue
				isqFileRawID = individual_sequence_file.individual_sequence_file_raw_id
				
				
				counter += 1
				library = individual_sequence_file.library
				splitOrder = individual_sequence_file.split_order
				mate_id = individual_sequence_file.mate_id
				if library not in library2Data:
					library2Data[library] = PassingData(isqFileRawID2Index = {}, isqFileRawDBEntryLs=[], \
											splitOrder2Index={}, fileObjectPairLs=[])
				#for convenience
				isqFileRawID2Index = library2Data[library].isqFileRawID2Index
				isqFileRawDBEntryLs = library2Data[library].isqFileRawDBEntryLs
				splitOrder2Index = library2Data[library].splitOrder2Index
				fileObjectPairLs = library2Data[library].fileObjectPairLs
				if isqFileRawID not in isqFileRawID2Index:
					isqFileRawID2Index[isqFileRawID] = len(isqFileRawID2Index)
					isqFileRawDBEntryLs.append(individual_sequence_file.individual_sequence_file_raw)
				
				if splitOrder not in splitOrder2Index:
					splitOrder2Index[splitOrder] = len(splitOrder2Index)
					fileObjectPairLs.append([])
				fileObjectPairIndex = splitOrder2Index[splitOrder]
				if mate_id is None:
					mate_id = 1
				noOfFileObjects = len(fileObjectPairLs[fileObjectPairIndex])
				if noOfFileObjects<mate_id:
					for i in xrange(mate_id-noOfFileObjects):	#expand the list to match the number of mates
						fileObjectPairLs[fileObjectPairIndex].append(None)
				isq_file_obj = PassingData(db_entry=individual_sequence_file, path=path)
				fileObjectPairLs[fileObjectPairIndex][mate_id-1] = isq_file_obj
			#add it to the individual_sequence db entry
			individual_sequence.library2Data = library2Data
			isqLs.append(individual_sequence)
		sys.stderr.write("%s individual sequence files from %s isq entries.\n"%(counter, len(isqLs)))
		return isqLs
	
	@classmethod
	def get_data_dir(cls, ):
		"""
		2012.7.13
			classmethod. assuming the db connection is in place. usually used within methods of the table classes
			
			Example:
				data_dir = VervetDB.get_data_dir()
		"""
		dataDirEntry = README.query.filter_by(title='data_dir').first()
		if not dataDirEntry or not dataDirEntry.description or not os.path.isdir(dataDirEntry.description):
			# todo: need to test dataDirEntry.description is writable to the user
			sys.stderr.write("data_dir not available in db or not accessible on the harddisk. Raise exception.\n")
			raise BaseException
			return None
		else:
			return dataDirEntry.description
	
	def constructPedgreeGraphOutOfAlignments(self, alignmentLs):
		"""
		2012.4.20
			Warning if one individual appears in the graph >1 times. 
		2011-12-15
			copied from AlignmentToTrioCallPipeline.py
			
			construct a directed graph (edge: from parent to child) of which nodes are all from alignmentLs.
		"""
		sys.stderr.write("Construct pedigree out of %s alignments... "%(len(alignmentLs)))
		import networkx as nx
		DG=nx.DiGraph()
		
		individual_id2alignmentLs = {}
		for alignment in alignmentLs:
			individual_id = alignment.individual_sequence.individual_id
			if individual_id not in individual_id2alignmentLs:
				individual_id2alignmentLs[individual_id] = []
			else:
				sys.stderr.write("Warning: individual_id %s appears in >%s alignments.\n"%(individual_id, \
																	len(individual_id2alignmentLs[individual_id])))
			individual_id2alignmentLs[individual_id].append(alignment)
			if DG.has_node(individual_id):
				sys.stderr.write("Warning: individual_id %s already in the graph.\n"%(individual_id))
			else:
				DG.add_node(individual_id)
		
		
		for row in Ind2Ind.query:
			if row.individual1_id in individual_id2alignmentLs and row.individual2_id in individual_id2alignmentLs:
				DG.add_edge(row.individual1_id, row.individual2_id)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
															nx.number_connected_components(DG.to_undirected())))
		return PassingData(DG=DG, individual_id2alignmentLs=individual_id2alignmentLs)
	
	def constructPedgree(self, directionType=1):
		"""
		2012.9.6 add argument directionType
			1: parent -> child as edge
			2: child -> parent as edge
			3: undirected
		2012.1.23
		"""
		sys.stderr.write("Constructing pedigree from db directionType=%s..."%(directionType))
		import networkx as nx
		if directionType==3:
			DG = nx.Graph()
		else:
			DG = nx.DiGraph()
		
		for row in Ind2Ind.query:
			if directionType==2:
				DG.add_edge(row.individual2_id, row.individual1_id)
			else:
				DG.add_edge(row.individual1_id, row.individual2_id)
		
		"""
		#2012.6.19 initialization, but found unnecessary. increasePredecessorNoOfDescendantCount() can do it itself.
		for node in DG:
			DG.node[node]["noOfDescendant"] = 0	#initialize this
			DG.node[node]["noOfChildren"] = len(DG.successors(node))	#
		"""
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
															nx.number_connected_components(DG.to_undirected())))
		return DG
	
	def accumulatePredecessorAttributeRecursively(self, DG=None, source=None, attributeName='noOfDescendant'):
		"""
		2012.11.30
			bugfix. use copy.deepcopy around attributeInitValue or attributeIncrementValue everytime.
		2012.6.19
			attributeName=
				noOfDescendant: wrong, don't use it
				descendantList: treated as a set
				descendantSet: same as descendantList
			a recursive function to accumulate the attribute of the predecessors of the source node
				by absorbing the attribute of the source node.
			The algorithm starts from the bottom of the pedigree. called in calculateCumulativeAttributeForEachNodeInPedigree() 
		"""
		import copy
		predecessorList = DG.predecessors(source)
		if attributeName =='noOfDescendant':
			attributeInitValue = 0
			attributeIncrementValue = 1
		elif attributeName =='descendantList' or attributeName=='descendantSet':
			attributeInitValue = set()
			attributeIncrementValue = set([source])
		else:	#same as noOfDescendant
			attributeInitValue = 0
			attributeIncrementValue = 1
		if attributeName not in DG.node[source]:	#never accessed
			DG.node[source][attributeName] = copy.deepcopy(attributeInitValue)	#avoid "by-reference"
		for node in predecessorList:
			if attributeName not in DG.node[node]:
				DG.node[node][attributeName] = copy.deepcopy(attributeInitValue)
			if attributeName=='descendantList' or attributeName=='descendantSet':
				DG.node[node][attributeName] |= copy.deepcopy(DG.node[source][attributeName])
				DG.node[node][attributeName] |= copy.deepcopy(attributeIncrementValue)	#source itself
			else:
				DG.node[node][attributeName] += copy.deepcopy(DG.node[source][attributeName])
				DG.node[node][attributeName] += copy.deepcopy(attributeIncrementValue)	#source itself
			self.accumulatePredecessorAttributeRecursively(DG, source=node, attributeName=attributeName)
	
	def calculateCumulativeAttributeForEachNodeInPedigree(self, DG=None, attributeName='noOfDescendant'):
		"""
		2012.6.19
			the algorithm starts from the bottom (leaves) of the graph and do a reverse BFS
		"""
		leafNodeList = []
		for node in DG:
			if len(DG.successors(node))==0:
				leafNodeList.append(node)
		sys.stderr.write("%s leaf nodes.\n"%(len(leafNodeList)))
		for leafNode in leafNodeList:
			self.accumulatePredecessorAttributeRecursively(DG, source=leafNode, attributeName=attributeName)
		
	
	def getNoOfDescendant(self, individual_id=None, attributeName = 'noOfDescendant'):
		"""
		2012.6.19
			this function is wrong. don't use it. use getDescendantSet() to get a unique set of descendants and get its length.
			It's wrong because in BFS upwards traversal, if one parent has >1 children, its descendant count would be added to
				that of grandparent for >1 times, which is bad.
		"""
		if not getattr(self, 'pedigreeDG', None):
			self.pedigreeDG = self.constructPedgree()
		
		if individual_id not in self.pedigreeDG.node:
			return None
		
		firstNode = self.pedigreeDG.nodes()[0]
		nodeProperty = self.pedigreeDG.node[firstNode]
		if attributeName not in nodeProperty:
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeDG, attributeName=attributeName)
		return self.pedigreeDG.node[individual_id][attributeName]
	
	def getDescendantSet(self, individual_id=None, attributeName = 'descendantList'):
		"""
		2012.6.19
			descendantList is treated as descendantSet.
			similar to getNoOfDescendant() (which is wrong) but return a set of descendant nodes
		"""
		if not getattr(self, 'pedigreeDG', None):
			self.pedigreeDG = self.constructPedgree()
		
		if individual_id not in self.pedigreeDG.node:
			return None
		
		firstNode = self.pedigreeDG.nodes()[0]
		nodeProperty = self.pedigreeDG.node[firstNode]
		if attributeName not in nodeProperty:	#not run yet
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeDG, attributeName=attributeName)
		return self.pedigreeDG.node[individual_id][attributeName]
	
	def getAncestorSet(self, individual_id=None, attributeName = 'descendantList'):
		"""
		2012.11.30
			use a reverse pedigree DG
		2012.6.19
			descendantList is treated as descendantSet.
			similar to getNoOfDescendant() (which is wrong) but return a set of descendant nodes
		"""
		if not getattr(self, 'pedigreeReverseDG', None):
			self.pedigreeReverseDG = self.constructPedgree(directionType=2)
		
		if individual_id not in self.pedigreeReverseDG.node:
			return None
		
		firstNode = self.pedigreeReverseDG.nodes()[0]
		nodeProperty = self.pedigreeReverseDG.node[firstNode]
		if attributeName not in nodeProperty:	#not run yet
			self.calculateCumulativeAttributeForEachNodeInPedigree(self.pedigreeReverseDG, attributeName=attributeName)
		return self.pedigreeReverseDG.node[individual_id][attributeName]
	
	
	def findFamilyFromPedigreeGivenSize(self, DG=None, familySize=3, removeFamilyFromGraph=True, noOfOutgoingEdges=None):
		"""
		2012.4.2
			add noOfOutgoingEdges to find true singletons.
		2011.12.15
			copied from AlignmentToTrioCallPipeline.py
		2011-12.4
			DG is a directed graph.
			It views the graph from the offspring point of view because it's designed to find max number of trios (for trioCaller).
			Nuclear families (size>3) could not be handled by trioCaller and not considered in this algorithm.
			
			1. sort all nodes by their out degree ascendingly
				nodes with higher out degree would be checked later in the process to minimize the disruption of dependent trios.
				(if this offspring happens to be parent in other trios).
			2. if a node's incoming degree is same as familySize-1, then a family of given size is found.
		"""
		sys.stderr.write("Finding families of size %s .."%(familySize))
		familyLs = []	#each element of this list is either [singleton] or [parent, child] or [father, mother, child]
		allNodes = DG.nodes()
		out_degree_node_ls = []
		for node in allNodes:
			out_degree_node_ls.append([DG.out_degree(node), node])
		out_degree_node_ls.sort()	#nodes with low out degree would be considered first.
		
		for out_degree, node in out_degree_node_ls:
			if DG.has_node(node):	#it could have been removed as trios/duos were found
				#edges = DG.in_edges(node)
				noOfIncomingEdges = DG.in_degree(node)
				if noOfOutgoingEdges is not None and out_degree!=noOfOutgoingEdges:
					continue
				#noOfIncomingEdges = len(edges)
				if noOfIncomingEdges==familySize-1:	# a trio is found
					parents = DG.predecessors(node)
					family = parents + [node]
					familyLs.append(family)
					if removeFamilyFromGraph:
						DG.remove_nodes_from(family)
				elif noOfIncomingEdges>familySize-1:
					sys.stderr.write("Warning: noOfIncomingEdges for node %s is %s (family size =%s).\n"%(node, \
																noOfIncomingEdges, noOfIncomingEdges+1))
		sys.stderr.write(" found %s.\n"%(len(familyLs)))
		return familyLs
	
	@classmethod
	def parseAlignmentCompositeID(cls, compositeID):
		"""
		2012.5.4
			composite ID is used to designate members in trio, similar to read groups used in bam files but different.
			It's generated via IndividualAlignment.getCompositeID()
		
				compositeID = '%s_%s_%s_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.id, \
								self.individual_sequence.individual.code)
		"""
		split_id_ls = compositeID.split('_')
		individual_alignment_id = int(split_id_ls[0])
		ind_seq_id = int(split_id_ls[1])
		individual_id = int(split_id_ls[2])
		individual_code = split_id_ls[3]
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code)
	
	@classmethod
	def parseAlignmentReadGroupWithoutDB(cls, read_group):
		"""
		2012.5.11
			read-group is used to designate samples in alignment files and also in the ensuing multi-sample VCF files,
				similar to compositeID used in trio-representations but different.
			It's generated via IndividualAlignment.getReadGroup()
				read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
								
		"""
		individual_alignment_id = None
		ind_seq_id = None
		individual_code = None
		sequencer = None
		ref_ind_seq_id = None
		individual_id = None
		try:
			split_id_ls = read_group.split('_')
			individual_alignment_id = int(split_id_ls[0])
			ind_seq_id = int(split_id_ls[1])
			individual_code = split_id_ls[2]
			sequencer = split_id_ls[3]
			ref_ind_seq_id = int(split_id_ls[5])	#index 4 is "vs"
			individual_id = None
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			
			
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer, individualAlignment=None)
	
	def parseAlignmentReadGroup(self, read_group):
		"""
		2012.10.8 no more try except
		2012.5.10
			read-group is used to designate samples in alignment files and also in the ensuing multi-sample VCF files,
				similar to compositeID used in trio-representations but different.
			It's generated via IndividualAlignment.getReadGroup()
				read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.individual_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
								
		"""
		split_id_ls = read_group.split('_')
		try:
			individual_alignment_id = int(split_id_ls[0])
		except:	#only catch the string error, not catching the db error,
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			individualAlignment = None
			individual_alignment_id = None
			ind_seq_id = None
			individual_code = None
			sequencer = None
			ref_ind_seq_id = None
			individual_id = None
			return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer, individualAlignment=individualAlignment)
		
		individualAlignment = IndividualAlignment.get(individual_alignment_id)
		ind_seq_id = individualAlignment.individual_sequence.id
		individual_code = individualAlignment.individual_sequence.individual.code
		sequencer = individualAlignment.individual_sequence.sequencer
		ref_ind_seq_id = individualAlignment.ref_sequence.id
		individual_id = individualAlignment.individual_sequence.individual.id
		"""
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			individualAlignment = None
			individual_alignment_id = None
			ind_seq_id = None
			individual_code = None
			sequencer = None
			ref_ind_seq_id = None
			individual_id = None
		"""
		return PassingData(individual_alignment_id=individual_alignment_id, individual_sequence_id=ind_seq_id,\
						individual_id=individual_id, individual_code=individual_code, ref_ind_seq_id=ref_ind_seq_id,\
						sequencer=sequencer, individualAlignment=individualAlignment)
	
	def parseDBEntryGivenDBAffiliatedFilename(self, filename=None, TableClass=None):
		"""
		2012.8.15
			almost all db-affliated filenames are beginned with the db_entry's primary ID, separated by _.
			this is a generic function, to be called by others with different TableClass
		"""
		baseFilename = os.path.basename(filename)
		try:
			dbID = int(baseFilename.split('_')[0])
			db_entry = TableClass.get(dbID)
		except:
			sys.stderr.write("Warning: could not parse dbID from %s.\n"%(baseFilename))
			db_entry = None
		return db_entry
	
	def parseGenotypeFileGivenDBAffiliatedFilename(self, filename=None):
		"""
		2012.8.15
			call parseDBEntryGivenDBAffiliatedFilename
		"""
		return self.parseDBEntryGivenDBAffiliatedFilename(filename=filename, TableClass=GenotypeFile)
	
	def parseIndividualAlignmentGivenDBAffiliatedFilename(self, filename=None):
		"""
		2012.8.15
			call parseDBEntryGivenDBAffiliatedFilename
		"""
		return self.parseDBEntryGivenDBAffiliatedFilename(filename=filename, TableClass=IndividualAlignment)
	
	def findISQFoldersNotInDB(self, dataDir=None, subFolder='individual_sequence',):
		"""
		2012.7.11 scan the individual_sequence folder to see which subfolder is not in individual_sequence.path
			i.e.
				orphanPathList = db_vervet.findISQFoldersNotInDB(dataDir=None, subFolder='individual_sequence',)
				orphanPathList = db_vervet.findISQFoldersNotInDB(dataDir='/u/home/eeskin2/polyacti/NetworkData/vervet/db/', \
														subFolder='individual_sequence',)
		"""
		sys.stderr.write("Finding individual_sequence folders that are not in db (orphaned)... ")
		dbPath2dbEntry = {}
		for row in IndividualSequence.query:
			if row.path and row.path not in dbPath2dbEntry:
				dbPath2dbEntry[row.path] = row
		sys.stderr.write("%s items in db.\n"%(len(dbPath2dbEntry)))
		
		if not dataDir:
			dataDir = self.data_dir
		
		topFolder = os.path.join(dataDir, subFolder)
		orphanPathList = []
		for item in os.listdir(topFolder):
			itemAbsPath = os.path.join(topFolder, item)
			if os.path.isdir(itemAbsPath):
				isqPath = os.path.join(subFolder, item)
				if isqPath not in dbPath2dbEntry:
					orphanPathList.append(itemAbsPath)
		sys.stderr.write(" %s found.\n"%(len(orphanPathList)))
		return orphanPathList
	
	def findAlignmentFileNotInDB(self, dataDir=None, subFolder='individual_alignment',):
		"""
		2012.7.11 scan the individual_alignment folder to see which file is not in individual_alignment.path
			i.e.
				orphanPathList = db_vervet.findAlignmentFileNotInDB(dataDir=None, subFolder='individual_alignment',)
		"""
		sys.stderr.write("Finding individual_alignment folders that are not in db (orphaned)... ")
		dbPath2dbEntry = {}
		for row in IndividualAlignment.query:
			if row.path and row.path not in dbPath2dbEntry:
				dbPath2dbEntry[row.path] = row
		
		sys.stderr.write("%s items in db.\n"%(len(dbPath2dbEntry)))
		
		if not dataDir:
			dataDir = self.data_dir
		
		topFolder = os.path.join(dataDir, subFolder)
		orphanPathList = []
		from pymodule import utils
		for item in os.listdir(topFolder):
			itemAbsPath = os.path.join(topFolder, item)
			
			if os.path.isfile(itemAbsPath) and utils.getRealPrefixSuffixOfFilenameWithVariableSuffix(itemAbsPath)[1]=='.bam':
				#make sure the suffix is bam. loads of .bai files are not recorded in db.
				isqPath = os.path.join(subFolder, item)
				if isqPath not in dbPath2dbEntry:
					orphanPathList.append(itemAbsPath)
		sys.stderr.write(" %s found.\n"%(len(orphanPathList)))
		return orphanPathList
		
	def updateGenotypeMethodNoOfLoci(self, db_entry=None, format=None, no_of_chromosomes=None):
		"""
		2012.8.30
		2012.8.15
			add argument format and aggregate the number of loci according to format
		2012.7.17
		"""
		no_of_loci = 0
		formatAndNoOfChromosomes2NoOfLoci = {}
		query = GenotypeFile.query.filter_by(genotype_method_id=db_entry.id)
		if format:
			query = query.filter_by(format=format)
		if no_of_chromosomes:	#2012.8.30
			query = query.filter_by(no_of_chromosomes = no_of_chromosomes)
		
		for genotype_file in query:
			key = (genotype_file.format, genotype_file.no_of_chromosomes)
			if key not in formatAndNoOfChromosomes2NoOfLoci:
				formatAndNoOfChromosomes2NoOfLoci[key] = 0
			formatAndNoOfChromosomes2NoOfLoci[key] += genotype_file.no_of_loci
		noOfLociList = formatAndNoOfChromosomes2NoOfLoci.values()
		#do some checking here to make sure every number in noOfLociList agrees with each other
		db_entry.no_of_loci = max(noOfLociList)	#take the maximum
		self.session.add(db_entry)
		self.session.flush()
	
if __name__ == '__main__':
	main_class = VervetDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.setup(create_tables=True)
	import pdb
	pdb.set_trace()
