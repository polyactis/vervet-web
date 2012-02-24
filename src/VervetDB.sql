-- 2011-4-29

drop view view_individual;
create or replace view view_individual as select i.id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
	i.approx_age_group_at_collection, u.realname as collector, i.collection_date, i.latitude, 
	i.longitude, s.short_name as site_name,
	s.city, s.stateprovince as province, c.name as country from individual i, site s, country c,
	acl_user u where i.site_id=s.id and s.country_id=c.id and i.collector_id=u.id;
	
drop view view_alignment;
create or replace view view_alignment as select i.id, i.code, i.ucla_id, i.tax_id, i.sex, i.age, 
	i.collection_date, i.latitude, i.longitude, isq.id as isq_id, isq.sequencer, isq.sequence_type, isq.base_count, 
	isq.coverage as raw_coverage, ia.id as alignment_id, ia.ref_ind_seq_id, ia.aln_method_id,
	ia.median_depth, ia.mean_depth, ia.mode_depth
	from individual i, individual_alignment ia, individual_sequence isq where isq.individual_id=i.id and isq.id=ind_seq_id
	order by id, isq.sequencer, isq.sequence_type, ref_ind_seq_id;
	
drop view view_individual_sequence;
create or replace view view_individual_sequence as select isq.id as isq_id, isq.sequencer, isq.sequence_type, isq.format, isq.coverage, isq.path, 
	i.id as individual_id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
	i.approx_age_group_at_collection, i.collection_date, i.latitude, 
	i.longitude, s.short_name as site_name, s.city, s.stateprovince as province
	from individual i, individual_sequence isq, site s
	where i.site_id=s.id and i.id=isq.individual_id order by code;