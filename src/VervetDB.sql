-- 2011-4-29

drop view view_individual;
create or replace view view_individual as select i.id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
    i.approx_age_group_at_collection, u.realname as collector, i.collection_date, i.latitude, 
    i.longitude,i.site_id, s.short_name as site_name,
    s.city, s.stateprovince as province, c.name as country from individual i, site s, country c,
    acl_user u where i.site_id=s.id and s.country_id=c.id and i.collector_id=u.id;
    
drop view view_alignment;
create or replace view view_alignment as select i.id as individual_id, i.code, i.ucla_id, i.tax_id, i.sex, i.age, 
    i.site_id, i.collection_date, i.latitude, i.longitude, isq.id as isq_id, 
    isq.filtered, isq.sequencer, isq.sequence_type, isq.tissue_id, isq.base_count, 
    isq.coverage as raw_coverage, ia.id as alignment_id, ia.ref_ind_seq_id, ia.aln_method_id,
    ia.median_depth, ia.mean_depth, ia.mode_depth, ia.outdated_index
    from individual i, individual_alignment ia, individual_sequence isq where isq.individual_id=i.id and isq.id=ind_seq_id
    order by individual_id, isq.sequencer, isq.sequence_type, ref_ind_seq_id, alignment_id;
    
drop view view_individual_sequence;
create or replace view view_individual_sequence as select newt.*, t.short_name as tissue 
    from ( select isq.id as individual_sequence_id, isq.sequencer, 
    isq.sequence_type, isq.format, i.target_coverage, isq.coverage, isq.path, isq.filtered, isq.quality_score_format, isq.tissue_id, 
    i.id as individual_id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
    i.approx_age_group_at_collection, i.collection_date, i.latitude, 
    i.longitude, s.id as site_id, s.short_name as site_name, s.city, c.name as country, isq.date_created, isq.date_updated
    from individual i, individual_sequence isq, site s, country c where i.site_id=s.id and i.id=isq.individual_id and s.country_id=c.id) 
    as newt left join tissue t on t.id=newt.tissue_id order by code;

drop view view_individual_sequence_with_raw;
create or replace view view_individual_sequence_with_raw as select v.*, isqr.id as isqr_id, isqr.library, 
    isqr.path as raw_path, isqr.file_size, isqr.date_created as raw_date_created 
    from view_individual_sequence v, individual_sequence_file_raw isqr where v.individual_sequence_id=isqr.individual_sequence_id
    order by individual_sequence_id, filtered, raw_date_created;

-- 2012.7.12 view isq_file
drop view view_individual_sequence_file;
create or replace view view_individual_sequence_file as select v.*, isqf.id as isqf_id, isqf.library, isqf.split_order, isqf.mate_id, 
    isqf.base_count, isqf.date_created as file_date_created 
    from view_individual_sequence v, individual_sequence_file isqf 
    where v.individual_sequence_id=isqf.individual_sequence_id 
    order by individual_sequence_id, filtered, split_order, mate_id;

-- 2012.7.12 view isq_file with raw
drop view view_individual_sequence_file_with_raw;
create or replace view view_individual_sequence_file_with_raw as select isqr.id as isqr_id, isqr.library, 
    isqr.path as raw_path, isqr.file_size, isqr.date_created as raw_date_created, 
    isqf.id as isqf_id, isqf.individual_sequence_id, isqf.filtered, isqf.library as isqf_library, 
    isqf.split_order, isqf.mate_id, isqf.base_count, isqf.read_count
    from individual_sequence_file_raw isqr, individual_sequence_file isqf 
    where isqf.individual_sequence_file_raw_id=isqr.id
    order by individual_sequence_file_raw_id, filtered, raw_date_created;

-- 2012.7.12 view isqr summary
drop view view_isqr_summary;
create or replace view view_isqr_summary as 
    select isqr_id, library, individual_sequence_id, filtered, file_size, sum(base_count) as sum_base_count, sum(read_count) as sum_read_count
    from view_individual_sequence_file_with_raw 
    group by isqr_id, library, individual_sequence_id, filtered, file_size
    order by individual_sequence_id, filtered, library;

-- 2012.6.3 this view is not right.
drop view view_compare_seq_coverage;
create or replace view view_compare_seq_coverage as select i1.individual_id, i1.path, i1.coverage as coverage_after_filter, 
    i1.tissue_id,   
    i2.coverage, i1.coverage/i2.coverage as ratio from individual_sequence i1 right join individual_sequence i2 
    on i1.parent_individual_sequence_id=i2.id order by ratio;

drop view view_sequence_before_after_filter;
create or replace view view_sequence_before_after_filter as select i.code, i.ucla_id, t1.*
    from (select s2.individual_id, s1.id as isq_id, s1.coverage, s1.tissue_id, s1.date_created, s2.id as unfilter_isq_id, s2.coverage as unfilter_coverage, 
    float8(s1.base_count)/s2.base_count as ratio, s2.date_created as unfilter_date_created from (select * from individual_sequence where filtered =1) as s1 right join 
    (select * from individual_sequence where filtered =0) as s2 on s1.parent_individual_sequence_id=s2.id) as t1,
    individual i where i.id=t1.individual_id;

drop view viewInd2Ind;
create or replace view viewInd2Ind as select i1.code as i1_code, i2i.individual1_id, i2.code as i2_code , i2i.individual2_id, i2i.relationship_type_id
    from ind2ind i2i, individual i1, individual i2 where  i1.id=i2i.individual1_id and i2.id=i2i.individual2_id;

-- 2012.6.26 view on all the alignments
drop view view_seq_comparison_with_alignment;
create or replace view view_seq_comparison_with_alignment as select v.*, va1.alignment_id, va1.aln_method_id,
    va1.ref_ind_seq_id, va1.median_depth, va1.outdated_index
    from view_sequence_before_after_filter v left join view_alignment va1 on v.isq_id=va1.isq_id;

--  2012.7.3 check whether each mate in individual_sequence_file has its corresponding entry;
drop view check_isq_file_parity;
create or replace view check_isq_file_parity as select individual_sequence_id, library, split_order, count(mate_id) 
 from (select distinct individual_sequence_id, library, split_order, mate_id from individual_sequence_file  )
 as newt group by individual_sequence_id, library, split_order 
 order by count, individual_sequence_id, library, split_order;
