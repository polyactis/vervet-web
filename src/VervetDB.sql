-- 2011-4-29

drop view view_individual;
create or replace view view_individual as select i.id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
	i.approx_age_group_at_collection, u.realname as collector, i.collection_date, i.latitude, 
	i.longitude, s.short_name as site_name,
	s.city, s.stateprovince as province, c.name as country from individual i, site s, country c,
	acl_user u where i.site_id=s.id and s.country_id=c.id and i.collector_id=u.id;