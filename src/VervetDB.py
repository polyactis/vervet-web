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

bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:	   #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
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
import os
import hashlib


__session__ = scoped_session(sessionmaker(autoflush=False, autocommit=True))
#__metadata__ = ThreadLocalMetaData()	#2008-10 not good for pylon

__metadata__ = MetaData()

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
	latitude = Field(Float)	#2011-3-1
	longitude = Field(Float)	#2011-3-1
	altitude = Field(Float)	#2011-4-28
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
		2011.12.4
			represent male as 1. represent female as 2.
		"""
		if self.sex[0]=='M':
			return 1
		else:
			return 2
	
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
	2011-11-28
		add pass_qc_read_base_count which counts the number of bases in all pass-QC reads (base quality>=20, read mapping quality >=30).
	2011-9-19
		add mean_depth, read_group_added
	2011-8-3
		add column median_depth, mode_depth
	2011-3-3
	"""
	ind_sequence = ManyToOne('IndividualSequence', colname='ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_sequence = ManyToOne('IndividualSequence', colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	aln_method = ManyToOne('AlignmentMethod', colname='aln_method_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype_method_ls = ManyToMany("GenotypeMethod",tablename='genotype_method2individual_alignment', local_colname='individual_alignment_id')
	path = Field(Text)
	format = Field(String(512))
	median_depth = Field(Float)	#2011-8-2
	mode_depth = Field(Float)	#2011-8-2
	mean_depth = Field(Float)	#2011-9-12
	pass_qc_read_base_count = Field(BigInteger)	#2011-11-28	QC = (base quality>=20, read mapping quality >=30). 
	read_group_added = Field(Integer, default=0)	# 2011-9-15 0=No, 1=Yes
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_alignment', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('ind_seq_id', 'ref_ind_seq_id', 'aln_method_id'))
	
	def getReadGroup(self):
		"""
		2011-11-02
			read group is this alignment's unique identifier in a sam/bam file.
		"""
		sequencer = self.ind_sequence.sequencer
		read_group = '%s_%s_%s_%s_vs_%s'%(self.id, self.ind_seq_id, self.ind_sequence.individual.code, \
								sequencer, self.ref_ind_seq_id)
		return read_group

	def getCompositeID(self):
		"""
		2012.1.25
			ID to be used to identify members of trios. almost same as getReadGroup()
		"""
		read_group = '%s_%s_%s_%s'%(self.id, self.ind_seq_id, self.ind_sequence.individual.id, \
								self.ind_sequence.individual.code)
		return read_group
	
	def constructRelativePath(self, subFolder='individual_alignment'):
		"""
		2012.2.24
			fix a bug
		2012.2.10
			moved from VervetDB.constructRelativePathForIndividualAlignment
		2011-8-29
			called by getAlignment() and other programs
		"""
		#'/' must not be put in front of the relative path.
		# otherwise, os.path.join(self.data_dir, dst_relative_path) will only take the path of dst_relative_path.
		dst_relative_path = '%s/%s_%s_vs_%s_by_%s.%s'%(subFolder, self.id, self.ind_seq_id,\
												self.ref_ind_seq_id, self.aln_method_id, self.format)
		
		return dst_relative_path
	
class IndividualSequence(Entity, TableClass):
	"""
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

class IndividualSequenceFile(Entity, TableClass):
	"""
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
	base_count = Field(BigInteger)	#2011-8-2
	path = Field(Text, unique=True)	#path to the actual file
	format = Field(String(512))	#fasta, fastq
	quality_score_format = Field(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	filtered = Field(Integer, default=0)	#0 means not. 1 means yes.
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


class IndividualSequenceFileRaw(Entity, TableClass):
	"""
	2012.1.26
		this table is used to store the bam files from WUSTL. The bam files hold either single-end or paired-end reads
			in one file.
		Technically, this table is parent of IndividualSequenceFile.
	"""
	individual_sequence = ManyToOne('IndividualSequence', colname='individual_sequence_id', ondelete='CASCADE', onupdate='CASCADE')
	individual_sequence_file_ls = OneToMany("IndividualSequenceFile")
	library = Field(Text)	#id for the preparation library
	base_count = Field(BigInteger)
	path = Field(Text)	#path to the actual file
	md5sum = Field(Text, unique=True)	#used to identify each raw file
	quality_score_format = Field(String(512), default='Standard')
		#Standard=Phred+33 (=Sanger), Illumina=Phred+64 (roughly, check pymodule/utils for exact formula)
		# Illumina1.8+ (after 2011-02) is Standard.
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='individual_sequence_file_raw', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('library', 'md5sum'))


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
	2011-4-5
		add table locus_method
	2011-2-3
		add ref_allele, alt_allele
	"""
	chromosome = Field(String(512))
	start = Field(Integer)
	stop = Field(Integer)
	ref_seq = ManyToOne('Sequence', colname='ref_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	alt_seq = ManyToOne('Sequence', colname='alt_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	ref_sequence = ManyToOne('IndividualSequence', colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	#locus_method = ManyToOne('LocusMethod', colname='locus_method_id', ondelete='CASCADE', onupdate='CASCADE')
	locus_method_ls = ManyToMany('LocusMethod',tablename='locus2locus_method', local_colname='locus_id', \
								remote_colname='locus_method_id')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'ref_ind_seq_id'))

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


class LocusMethod(Entity, TableClass):
	"""
	2011-4-5
		to mark different sets of loci
	"""
	short_name = Field(String(256), unique=True)
	description = Field(Text)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='locus_method', metadata=__metadata__, session=__session__)
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
	2011-2-4
		file format:
						locus1	locus2
			individual1	allele1/allele2
			individual2
	"""
	short_name = Field(String(256), unique=True)
	genotype_file_dir = Field(Text)
	bam_filename = Field(Text)
	vcf_filename = Field(Text)
	filename = Field(Text)
	description = Field(Text)
	individual_alignment_ls = ManyToMany("IndividualAlignment",tablename='genotype_method2individual_alignment', \
										local_colname='genotype_method_id')
	ref_sequence = ManyToOne('IndividualSequence', colname='ref_ind_seq_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genotype_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('short_name', 'ref_ind_seq_id'))

class GenotypeFile(Entity, TableClass):
	"""
	2011-2-4
		file format:
			locus.id	allele_order	allele_type	seq.id	score	target_locus
	"""
	individual = ManyToOne('Individual', colname='individual_id', ondelete='CASCADE', onupdate='CASCADE')
	filename = Field(Text)
	genotype_method = ManyToOne('GenotypeMethod', colname='genotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(String(4096))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='genotype_file', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('individual_id', 'filename', 'genotype_method_id'))

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
		if getattr(self, 'schema', None):	#for postgres
			for entity in entities:
				if entity.__module__==self.__module__:	#entity in the same module
					using_table_options_handler(entity, schema=self.schema)
		
		#2008-10-05 MySQL typically close connections after 8 hours resulting in a "MySQL server has gone away" error.
		__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle, echo=self.sql_echo)
		
		self.metadata = __metadata__
		self.session = __session__
	
	def setup(self, create_tables=True):
		"""
		2008-09-07
			expose option create_tables, default=True. assign it to False if no new table is to be created.
		"""
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()
		#2008-08-26 setup_all() would setup other databases as well if they also appear in the program. Seperate this to be envoked after initialization
		# to ensure the metadata of other databases is setup properly.
	
	def findGenotypeMethodGivenName(self, genotypeMethodName):
		"""
		2011-2-11
			create one if not-existent
		"""
		genotypeMethod = GenotypeMethod.query.filter_by(short_name=genotypeMethodName).first()
		if not genotypeMethod:
			genotypeMethod = GenotypeMethod(short_name=genotypeMethodName)
			self.session.add(genotypeMethod)
			self.session.flush()
		return genotypeMethod
	
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
	
	def getAlignment(self, individual_code=None, individual_sequence_id=None, path_to_original_alignment=None, sequencer='GA', \
					sequence_type='SR', sequence_format='fastq', \
					ref_individual_sequence_id=10, \
					alignment_method_name='bwa-short-read', alignment_format='bam', subFolder='individual_alignment', \
					createSymbolicLink=False, individual_sequence_filtered=0, read_group_added=None, dataDir=None):
		"""
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
		individual = self.getIndividual(individual_code)
		if individual_sequence_id:
			individual_sequence = IndividualSequence.get(individual_sequence_id)
		else:
			individual_sequence = self.getIndividualSequence(individual.id, sequencer=sequencer, sequence_type=sequence_type,\
								sequence_format=sequence_format, filtered=individual_sequence_filtered)
		alignment_method = self.getAlignmentMethod(alignment_method_name=alignment_method_name)
		query = IndividualAlignment.query.filter_by(ind_seq_id=individual_sequence.id).\
				filter_by(ref_ind_seq_id=ref_individual_sequence_id).\
				filter_by(aln_method_id=alignment_method.id)
		if alignment_format:
			query = query.filter_by(format=alignment_format)
		db_entry = query.first()
		if not db_entry:
			db_entry = IndividualAlignment(ind_seq_id=individual_sequence.id, ref_ind_seq_id=ref_individual_sequence_id,\
								aln_method_id=alignment_method.id, format=alignment_format, read_group_added=read_group_added)
			self.session.add(db_entry)
			self.session.flush()
			#copy the file over
			
			#'/' must not be put in front of the relative path.
			# otherwise, os.path.join(dataDir, dst_relative_path) will only take the path of dst_relative_path.
			dst_relative_path = db_entry.constructRelativePath(subFolder=subFolder)
			
			#update its path in db to the relative path
			db_entry.path = dst_relative_path
			
			dst_pathname = os.path.join(dataDir, dst_relative_path)
			from pymodule.utils import runLocalCommand
			if path_to_original_alignment and (os.path.isfile(path_to_original_alignment) or os.path.isdir(path_to_original_alignment)):
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
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	
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
					altitude=None, ucla_id=None, site=None, collection_date=None, collector=None, \
					approx_age_group_at_collection=None, tax_id=None, birthdate=None, vrc_founder=None, comment=None,\
					microchip_id=None):
		"""
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
		if ucla_id:	#ucla_id should be unique
			db_entry = Individual.query.filter_by(ucla_id=ucla_id).first()
		else:
			query = Individual.query.filter_by(code=code)
			if tax_id:
				query = query.filter_by(tax_id=tax_id)
			db_entry = query.first()
		if not db_entry:
			if sex is None or sex=='?' or sex=='':	#2011-4-29
				sex = None
			elif len(sex)>=1:	#2011-5-5 take the first letter
				sex = sex[0].upper()
			db_entry = Individual(code=code, sex=sex, age=age, age_cas=age_cas, latitude=latitude,\
						longitude=longitude, altitude=altitude, ucla_id=ucla_id, site=site, \
						collection_date=collection_date, collector=collector,\
						approx_age_group_at_collection=approx_age_group_at_collection, tax_id=tax_id, birthdate=birthdate,\
						vrc_founder=vrc_founder, comment=comment, microchip_id=microchip_id)
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequence(self, individual_id=None, sequencer=None, sequence_type=None,\
						sequence_format=None, path_to_original_sequence=None, tissue_name=None, coverage=None,\
						subFolder='individual_sequence', quality_score_format="Standard", filtered=0,\
						parent_individual_sequence_id=None, dataDir=None):
		"""
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
		if sequencer:
			query = query.filter_by(sequencer=sequencer)
		if sequence_type:
			query = query.filter_by(sequence_type=sequence_type)
		if sequence_format:
			query = query.filter_by(format=sequence_format)
		if tissue_name:
			tissue = self.getTissue(short_name=tissue_name)
			query = query.filter_by(tissue_id=tissue.id)
		if parent_individual_sequence_id:
			query = query.filter_by(parent_individual_sequence_id=parent_individual_sequence_id)
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
	
	def copyParentIndividualSequence(self, parent_individual_sequence=None, parent_individual_sequence_id=None,\
									subFolder='individual_sequence', quality_score_format='Standard', filtered=1):
		"""
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
		return db_entry
	
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
	
	def getIndividualSequenceFileRaw(self, individual_sequence_id, library=None, md5sum=None, path=None):
		"""
		2012.2.14
		"""
		#query first
		query = IndividualSequenceFileRaw.query.filter_by(individual_sequence_id=individual_sequence_id)
		if library:
			query = query.filter_by(library=library)
		if md5sum:
			query = query.filter_by(md5sum=md5sum)
		if path:
			query = query.filter_by(path=os.path.realpath(path))
		db_entry = query.first()
		if not db_entry:
			db_entry = IndividualSequenceFileRaw(individual_sequence_id=individual_sequence_id, library=library, md5sum=md5sum, \
										path=os.path.realpath(path))
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getIndividualSequenceFile(self, individual_sequence_id, library=None, mate_id=None, split_order=None, format=None,\
							filtered=0, parent_individual_sequence_file_id=None, individual_sequence_file_raw_id=None,\
							quality_score_format='Standard'):
		"""
		2012.2.14
		"""
		#query first
		query = IndividualSequenceFile.query.filter_by(individual_sequence_id=individual_sequence_id)
		if library:
			query = query.filter_by(library=library)
		if mate_id:
			query = query.filter_by(mate_id=mate_id)
		if split_order:
			query = query.filter_by(split_order=split_order)
		if format:
			query = query.filter_by(format=format)
		if filtered:
			query = query.filter_by(filtered=filtered)
		if parent_individual_sequence_file_id:
			query = query.filter_by(parent_individual_sequence_file_id=parent_individual_sequence_file_id)
		if individual_sequence_file_raw_id:
			query = query.filter_by(individual_sequence_file_raw_id=individual_sequence_file_raw_id)
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
	
	def getGenotypeFile(self, individual, genotypeMethod):
		"""
		2011-2-11
		"""
		if individual.id and genotypeMethod.id:
			db_entry = GenotypeFile.query.filter_by(individual_id=individual.id).filter_by(genotype_method_id=genotypeMethod.id).first()
		else:
			db_entry = None
		if not db_entry:
			db_entry = GenotypeFile()
			db_entry.individual = individual
			db_entry.genotype_method = genotypeMethod
			db_entry.filename = os.path.join(genotypeMethod.genotype_file_dir, '%s_%s.tsv'%(genotypeMethod.id, individual.id))
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getLocus(self, chr=None, start=None, stop=None, ref_seq=None, alt_seq=None):
		"""
		2011-2-11
		
		"""
		db_entry = Locus.query.filter_by(chromosome=chr).filter_by(start=start).filter_by(stop=stop).first()
		if not db_entry:
			db_entry = Locus(chromosome=chr, start=start, stop=stop)
			db_entry.ref_seq = ref_seq
			db_entry.alt_seq = alt_seq
			self.session.add(db_entry)
			self.session.flush()
		return db_entry
	
	def getSite(self, description=None, city=None, stateprovince=None, country_name=None, latitude=None, longitude=None,\
			altitude=None):
		"""
		2011-4-28
		"""
		country = self.getCountry(country_name=country_name)
		if description:
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
	
	def getIndividualSequenceID2FilePairLs(self, individualSequenceIDList, dataDir=None, needPair=True, checkOldPath=False):
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
	
	def getISQ_ID2LibrarySplitOrder2FileLs(self, individualSequenceIDList, dataDir=None, filtered=None):
		"""
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
		for individualSequenceID in individualSequenceIDList:
			individual_sequence = IndividualSequence.get(individualSequenceID)
			if not individual_sequence:	#not present in db, ignore
				continue
			for individual_sequence_file in individual_sequence.individual_sequence_file_ls:
				if filtered is not None and individual_sequence_file.filtered!=filtered:	#skip entries that don't matched the filtered argument
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
				path = os.path.join(dataDir, individual_sequence_file.path)
				isq_file_obj = PassingData(db_entry=individual_sequence_file, path=path)
				LibrarySplitOrder2FileLs[key][mate_id-1] = isq_file_obj
		
		sys.stderr.write("%s individual sequence files from %s isq entries.\n"%(counter, len(isq_id2LibrarySplitOrder2FileLs)))
		return isq_id2LibrarySplitOrder2FileLs
	
	@property
	def data_dir(self, ):
		"""
		2011-2-11
			get the master directory in which all files attached to this db are stored.
		"""
		dataDirEntry = README.query.filter_by(title='data_dir').first()
		if not dataDirEntry or not dataDirEntry.description or not os.path.isdir(dataDirEntry.description):
			# todo: need to test dataDirEntry.description is writable to the user
			sys.stderr.write("data_dir not available in db or not accessible on the harddisk. Raise exception.\n")
			raise
			return None
		else:
			return dataDirEntry.description
	
	def constructPedgreeGraphOutOfAlignments(self, alignmentLs):
		"""
		2011-12-15
			copied from AlignmentToTrioCallPipeline.py
			
			construct a directed graph (edge: from parent to child) of which nodes are all from alignmentLs.
		"""
		sys.stderr.write("Construct pedigree out of %s alignments... "%(len(alignmentLs)))
		import networkx as nx
		DG=nx.DiGraph()
		
		individual_id2alignmentLs = {}
		for alignment in alignmentLs:
			individual_id = alignment.ind_sequence.individual_id
			if individual_id not in individual_id2alignmentLs:
				individual_id2alignmentLs[individual_id] = []
			else:
				sys.stderr.write("Warning: individual_id %s appears >1 alignments.\n"%(individual_id))
			individual_id2alignmentLs[individual_id].append(alignment)
			DG.add_node(individual_id)
		
		
		for row in Ind2Ind.query:
			if row.individual1_id in individual_id2alignmentLs and row.individual2_id in individual_id2alignmentLs:
				DG.add_edge(row.individual1_id, row.individual2_id)
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
															nx.number_connected_components(DG.to_undirected())))
		return PassingData(DG=DG, individual_id2alignmentLs=individual_id2alignmentLs)
	
	def constructPedgree(self):
		"""
		2012.1.23
		"""
		sys.stderr.write("Constructing pedigree from db ...")
		import networkx as nx
		DG=nx.DiGraph()
		
		for row in Ind2Ind.query:
			DG.add_edge(row.individual1_id, row.individual2_id)
		
		sys.stderr.write("%s nodes. %s edges. %s connected components.\n"%(DG.number_of_nodes(), DG.number_of_edges(), \
															nx.number_connected_components(DG.to_undirected())))
		return DG
	
	def findFamilyFromPedigreeGivenSize(self, DG, familySize=3, removeFamilyFromGraph=True):
		"""
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
	
if __name__ == '__main__':
	main_class = VervetDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.setup(create_tables=True)
	import pdb
	pdb.set_trace()
