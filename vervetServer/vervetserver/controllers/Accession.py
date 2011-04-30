import logging

from pylons import request, response, session, tmpl_context as c, url
from pylons.controllers.util import abort, redirect

from vervetserver.lib.base import BaseController, render

log = logging.getLogger(__name__)


from vervetserver import model
import simplejson
import vervetserver.lib.helpers as h
import gviz_api, datetime, re
from pymodule import PassingData
from pymodule.utils import processRegexpString
from pylons.decorators import jsonify


class AccessionController(BaseController):
	individual_central_view = 'view_individual'	#could also be ecotype_info.
	
	def index(self):
		"""
		2011-4-29
		"""
		c.find250kAccessionsURL = url(controller="Accession", action="findAllAccessions", id=None)
		c.findAccessionsByCountryURL = url(controller='Accession', action="findAccessionsByCountry", id=None)
		c.AccessionCountrySuggestOracleURL = url(controller='Accession', action="countryNameAutoComplete", id=None)
		c.AccessionNameSuggestOracleURL = url(controller='Accession', action="autoComplete", id=None)
		c.AccessionByNameURL = url(controller='Accession', action="findAccessionsByName", id=None)
		return render('/Accession.html')

	def findAccessionsByName(self):
		"""
		2011-4-29
		"""
		name = request.params.get('name', 'Algutsrum')	#default is 'Algutsrum'
		name = processRegexpString(name).p_str_sql
		condition = "code ~* '%s' or ucla_id ~* '%s' or name ~* '%s'"%\
						(name, name, name)
		return self.findAccessions(condition)
	
	def findAccessionsByID(self):
		"""
		2011-4-29
		"""
		ecotype_id = request.params.get('id', '100')	#default is 100
		condition = "id=%s"%(ecotype_id)
		return self.findAccessions(condition)
	
	def findAccessionsByCountry(self):
		"""
		2011-4-29
		"""
		country = request.params.get('country', 'USA')	#default is 100
		condition = "country='%s'"%(country)
		return self.findAccessions(condition)
	
	@classmethod
	def findAccessions(cls, condition=None, extra_tables=None):
		"""
		2011-4-29
		"""
		if condition:
			condition = 'where %s'%condition
		else:
			condition = ""
		table_str = '%s v'%cls.individual_central_view
		if extra_tables:
			table_str += ', %s'%extra_tables
		rows = model.db_vervet.metadata.bind.execute("select v.* from %s %s"%\
													(table_str, condition))
		
		#3rd finally construct the full data and turn it into json
		column_name_type_ls = [("tg_ecotypeid", ("number", "Individual ID")), ("nativename", ("string", "UCLA ID")), \
							("tax_id",("number", "Tax ID")), \
							("sex",("string", "sex")), \
							("age",("number", "age")), \
							("age_cas",("number", "age_cas")), \
							("approx_age_group_at_collection",("string", "Approx. Age")), \
							("latitude",("number", "Latitude")), \
							("longitude",("number", "Longitude")), \
							("city",("string", "City")),\
							("province", ("string", "Province")), ("country", ("string", "Country")),\
							("collector", ("string", "Collector")), ("collection_date", ("date", "Collection Date"))]
		
		description = dict(column_name_type_ls)
		return_ls = []
		for row in rows:
			entry = dict()
			for column_name_type in column_name_type_ls:
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type =='number':
					default_value = -1
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				
				#2011-4-29 fake some ids so that it can work with MapWithPhenotype from the client end
				if column_name=='tg_ecotypeid':
					column_value = getattr(row, 'id', default_value)
				elif column_name=='nativename':
					column_value = getattr(row, 'ucla_id', default_value)
				else:
					column_value = getattr(row, column_name, default_value)
				entry[column_name] = column_value
			return_ls.append(entry)
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		response.charset = 'utf-8'
		return json_result
	
	def findAccessionsNameLike(self, namelike):
		"""
		2011-4-29
			find all accessions with name beginned with 'namelike'
		"""
		name_processed_ob = processRegexpString(namelike)
		namelike = name_processed_ob.p_str_sql
		rows = model.db_vervet.metadata.bind.execute("select * from %s where code ~* '^%s' or ucla_id ~* '^%s'"%\
											(self.individual_central_view, namelike, namelike))
		name_set = set()
		name_p = re.compile(r'^%s'%name_processed_ob.p_str, re.IGNORECASE)
		columns_to_be_added = ['code', 'ucla_id', ]
		for row in rows:
			for col in columns_to_be_added:
				col_value = getattr(row, col)
				if col_value and name_p.search(col_value):
					name_set.add(col_value)
		return list(name_set)
	
	def findCountryLike(self, namelike):
		"""
		2011-4-30
			find all countries with name beginned with 'namelike'
		"""
		name_processed_ob = processRegexpString(namelike)
		namelike = name_processed_ob.p_str_sql
		rows = model.db_vervet.metadata.bind.execute("select * from %s where country ~* '^%s'"%\
											(self.individual_central_view, namelike))
		name_set = set()
		name_p = re.compile(r'^%s'%name_processed_ob.p_str, re.IGNORECASE)
		columns_to_be_added = ['country']
		for row in rows:
			for col in columns_to_be_added:
				col_value = getattr(row, col)
				if col_value and name_p.search(col_value):
					name_set.add(col_value)
		return list(name_set)
	
	#@jsonify	#2009-4-3, comment it out and using simplejson.dumps() directly because the default encoding 'utf-8' doesn't work due to some swedish letters.
	def autoComplete(self):
		"""
		2011-4-29
			autoComplete server end for AccessionByName.java
		"""
		namelike = request.params.get('namelike')
		name_ls = self.findAccessionsNameLike(namelike)
		name_ls.sort()
		if len(name_ls)>100:
			name_ls = name_ls[:100]
		response.headers['Content-Type'] = 'application/json'
		return simplejson.dumps(dict(result=name_ls), encoding='latin1')
	
	@jsonify
	def countryNameAutoComplete(self):
		"""
		2011-4-29
			auto complete server end for AccessionByCountry.java
		"""
		query = request.params.get('country')
		name_ls = self.findCountryLike(query)
		name_ls.sort()
		if len(name_ls)>100:
			name_ls = name_ls[:100]
		return dict(result=name_ls)
	
	def findAllAccessions(self):
		"""
		2011-4-29
		"""
		condition = "latitude is not null and longitude is not null"
		response.headers['Content-Type'] = 'application/json'
		return self.findAccessions(condition)
	