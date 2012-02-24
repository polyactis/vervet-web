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
		c.accessionAttributeNameURL = url(controller='Accession', action="getAccessionAttributeNameLs", id=None)
		c.accessionAttributeDataURL = url(controller='Accession', action='getAccessionAttributeValue', id=None)
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
		2011-9-13
			split country by "|" first
		2011-4-29
		"""
		country_ls = request.params.get('country', 'USA')	#default is USA
		country_ls = country_ls.split('|')
		condition_ls = []
		for country in country_ls:
			if country.strip():	#not empty after stripping
				condition_ls.append("country='%s'"%(country.strip()))
		
		condition = " or ".join(condition_ls)
		return self.findAccessions(condition)
	
	@classmethod
	def findAccessions(cls, condition=None, extra_tables=None):
		"""
		2011-5-3
			return "code" info as well
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
		column_name_type_ls = [("accession_id", ("number", "Individual ID")), \
							("nativename", ("string", "Code")),\
							("ucla_id", ("string", "UCLA ID")), \
							("tax_id",("number", "Tax ID")), \
							("sex",("string", "sex")), \
							("age",("number", "age")), \
							("age_cas",("number", "age_cas")), \
							("approx_age_group_at_collection",("string", "Approx. Age")), \
							("latitude",("number", "Latitude")), \
							("longitude",("number", "Longitude")), \
							("site", ("string", "Site")), ("city",("string", "City")),\
							("country", ("string", "Country")),\
							("collector", ("string", "Collector")), ("collection_date", ("date", "Collection Date")),\
							("seqCoverage", ("number", "Seq Coverage")), \
							("seqCoverageFilter", ("number", "Filtered Seq Coverage"))]
		
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
				if column_name=='accession_id':
					column_value = getattr(row, 'id', default_value)
				elif column_name=='nativename':
					column_value = getattr(row, 'code', default_value)
				elif column_name=='seqCoverage':
					isq = model.VervetDB.IndividualSequence.query.filter_by(individual_id=row.id).filter_by(filtered=0).first()
					if isq:
						column_value = isq.coverage
					else:
						column_value = default_value
				elif column_name=='seqCoverageFilter':
					isq = model.VervetDB.IndividualSequence.query.filter_by(individual_id=row.id).filter_by(filtered=1).first()
					if isq:
						column_value = isq.coverage
					else:
						column_value = default_value
				elif column_name=='site':
					column_value = getattr(row, "site_name", default_value)
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
	
	@jsonify
	def getAccessionAttributeNameLs(self):
		"""
		2011-5-15
		"""
		accessionAttributeNameLs = [[-3, 'Gender'], [-2, 'Age'], [-1, "Size"], [0, "Same"]]
		for row in model.VervetDB.PhenotypeMethod.query:
			label = '%s_%s'%(row.id, row.short_name)
			accessionAttributeNameLs.append([row.id, label])
		return accessionAttributeNameLs
	
	def getAccessionPhenotypeValueStructure(self, phenotype_method_id=None):
		"""
		2011-5-15
		"""
		accession_id2phenotype_value = {}
		min_value = None
		max_value = None
		query = model.VervetDB.Phenotype.query.filter_by(phenotype_method_id=phenotype_method_id)
		for row in query:
			if row.value is not None:
				phenotype_value = row.value
				accession_id2phenotype_value[row.individual_id] = phenotype_value
				
				if min_value == None or phenotype_value<min_value:
					min_value = phenotype_value
				if max_value == None or phenotype_value>max_value:
					max_value = phenotype_value
		dc = dict(min_value=min_value, max_value=max_value, accession_id2phenotype_value=accession_id2phenotype_value)
		return dc
	
	def getAccessionAttributeValue(self, returnJson=True, transform_func = None):
		"""
		2011-5-15
		"""
		attribute_id = int(request.params.get('attribute_id', 0))
		
		if attribute_id>0:
			dc = self.getAccessionPhenotypeValueStructure(phenotype_method_id=attribute_id)
			#dc["accession_id2attribute_value"] = dc["accession_id2phenotype_value"]
			dc = dict(min_value=dc["min_value"], max_value=dc["max_value"], accession_id2attribute_value=dc["accession_id2phenotype_value"])
		elif attribute_id<=0:	#same size for each cluster of accessions. so each accession attribute
			accession_id2attribute_value = {}
			min_value = None
			max_value = None
			query = model.VervetDB.Individual.query
			for row in query:
				if attribute_id==0:	#each cluster has the same size.
					attribute_value = 0
				elif attribute_id ==-1:	#size of cluster corresponds to the number of individuals included.
					attribute_value = 1
				elif attribute_id == -2:	#age
					attribute_value = row.age
				elif attribute_id == -3:	#gender
					attribute_value = None
					if row.sex:
						if row.sex[0]=='M':#male is coded as 0, 
							attribute_value = 0
						elif row.sex[0]=='F':
							attribute_value = 1
				if attribute_value is not None:
					if min_value == None or attribute_value<min_value:
						min_value = attribute_value
					if max_value == None or attribute_value>max_value:
						max_value = attribute_value
					accession_id2attribute_value[row.id] = attribute_value
			dc = dict(min_value=min_value, max_value=max_value, accession_id2attribute_value=accession_id2attribute_value)
			
		if transform_func is not None:
			dc["accession_id2attribute_value"] = transform_func(dc["accession_id2attribute_value"], dc['min_value'],\
													dc['max_value']) 
		if returnJson:
			response.headers['Content-Type'] = 'application/json'
			return simplejson.dumps(dc)
		else:
			return dc