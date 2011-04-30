"""Pylons environment configuration"""
import os

from mako.lookup import TemplateLookup
from pylons.configuration import PylonsConfig
from pylons.error import handle_mako_error

import vervetserver.lib.app_globals as app_globals
import vervetserver.lib.helpers
from vervetserver.config.routing import make_map

import vervetserver.model as model


def load_environment(global_conf, app_conf):
	"""Configure the Pylons environment via the ``pylons.config``
	object
	"""
	config = PylonsConfig()
	
	# Pylons paths
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	paths = dict(root=root,
				 controllers=os.path.join(root, 'controllers'),
				 static_files=os.path.join(root, 'public'),
				 templates=[os.path.join(root, 'templates')])

	# Initialize config with the basic options
	config.init_app(global_conf, app_conf, package='vervetserver', paths=paths)

	config['routes.map'] = make_map(config)
	config['pylons.app_globals'] = app_globals.Globals(config)
	config['pylons.h'] = vervetserver.lib.helpers
	
	# Setup cache object as early as possible
	import pylons
	pylons.cache._push_object(config['pylons.app_globals'].cache)
	

	# Create the Mako TemplateLookup, with the default auto-escaping
	config['pylons.app_globals'].mako_lookup = TemplateLookup(
		directories=paths['templates'],
		error_handler=handle_mako_error,
		module_directory=os.path.join(app_conf['cache_dir'], 'templates'),
		input_encoding='utf-8', default_filters=['escape'],
		imports=['from webhelpers.html import escape'])

	# CONFIGURATION OPTIONS HERE (note: all config options will override
	# any Pylons config options)
	
	#2008-10-05 setup the database connection
	drivername = config['app_conf']['drivername']
	hostname = config['app_conf']['hostname']
	dbname = config['app_conf']['dbname']
	schema = config['app_conf']['schema']
	db_user = config['app_conf']['db_user']
	db_passwd = config['app_conf']['db_passwd']
	pool_recycle = int(config['app_conf']['pool_recycle'])
	sql_echo = False
	if ('sql_echo' in config['app_conf'] and config['app_conf']['sql_echo'] == 'True'):
		sql_echo = True
	if config['app_conf'].get('echo_pool', False)=='True':	#watch:  bool('False')= True
		echo_pool = True	#2010-9-20 to enable monitoring of the pool
	else:
		echo_pool = False
	
	model.db_vervet = model.VervetDB.VervetDB(drivername=drivername, username=db_user, password=db_passwd, \
						hostname=hostname, database=dbname, schema=schema, pool_recycle=pool_recycle, \
						sql_echo=sql_echo, echo_pool= echo_pool)
	
	model.db_vervet.setup(create_tables=False)

	model.genome_db = model.GenomeDB.GenomeDatabase(drivername=drivername, username=db_user, password=db_passwd, \
				hostname=hostname, database='genome', schema=schema, pool_recycle=pool_recycle)
	model.genome_db.setup(create_tables=False)
	
	return config
