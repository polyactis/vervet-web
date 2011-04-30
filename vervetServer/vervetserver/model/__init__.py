import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
#import sqlalchemy
from vervet.src import VervetDB
from pymodule import GenomeDB


db_vervet = None
genome_db = None
gene_annotation = None