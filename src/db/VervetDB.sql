 
-- 2013.2.18 good whole-genome alignment
-- added "individual_sequence_file_raw_id is null" because there are some  sub-alignments (alignment using one lane/library's reads out of all)
drop view view_good_alignment cascade;
create or replace view view_good_alignment as select * from view_alignment 
    where filtered=1 and outdated_index=0 and alignment_method_id=2 and ref_ind_seq_id=3280
    and individual_sequence_file_raw_id is null;


-- 2011-4-29

drop view view_individual cascade;
create or replace view view_individual as select i.id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
    i.approx_age_group_at_collection, u.realname as collector, i.collection_date, s.latitude, 
    s.longitude, i.site_id, i.target_coverage, s.short_name as site_name,
    s.city, s.stateprovince as province, s.country_id, c.name as country from individual i, site s, country c,
    acl_user u where i.site_id=s.id and s.country_id=c.id and i.collector_id=u.id;
    
drop view view_alignment cascade;
create or replace view view_alignment as select i.id as individual_id, i.code, i.ucla_id, i.tax_id, i.sex, i.age, 
    i.site_id, i.collection_date, i.latitude, i.longitude, isq.id as isq_id, 
    isq.filtered, isq.sequencer_id, isq.sequence_type_id, isq.tissue_id, isq.base_count, 
    isq.coverage as raw_coverage, isq.is_contaminated, isq.outdated_index as isq_outdated_index,
    ia.id as alignment_id, ia.ref_ind_seq_id, ia.alignment_method_id,
    ia.median_depth, ia.mean_depth, ia.mode_depth, ia.outdated_index, ia.individual_sequence_file_raw_id,
    ia.file_size, ia.total_no_of_reads, ia.parent_individual_alignment_id, ia.mask_genotype_method_id, ia.local_realigned,
    ia.date_created as date_alignment_created, ia.date_updated as date_alignment_updated
    from view_individual i, individual_alignment ia, individual_sequence isq where isq.individual_id=i.id and isq.id=ind_seq_id
    order by individual_id, isq.sequencer_id, isq.sequence_type_id, ref_ind_seq_id, alignment_id;

    
drop view view_alignment_with_country cascade;
create or replace view view_alignment_with_country as select i.id as individual_id, i.code, i.ucla_id, i.tax_id, i.sex, i.age, 
    i.site_id, i.collection_date, i.latitude, i.longitude, i.country_id, i.country, isq.id as isq_id, 
    isq.filtered, isq.sequencer_id, isq.sequence_type_id, isq.tissue_id, isq.base_count, 
    isq.coverage as raw_coverage, ia.id as alignment_id, ia.ref_ind_seq_id, ia.alignment_method_id,
    ia.median_depth, ia.mean_depth, ia.mode_depth, ia.outdated_index, ia.individual_sequence_file_raw_id,
    ia.total_no_of_reads, 
    ia.date_created, ia.date_updated
    from view_individual i, individual_alignment ia, individual_sequence isq where isq.individual_id=i.id and isq.id=ind_seq_id
    order by individual_id, isq.sequencer_id, isq.sequence_type_id, ref_ind_seq_id, alignment_id;

-- 2012.10.1
drop view view_genotype_method_alignment;
create or replace view view_genotype_method_alignment as select va.*, g2a.genotype_method_id, g.short_name,
    g.no_of_individuals, g.no_of_loci
    from view_alignment va, genotype_method2individual_alignment g2a, genotype_method g
    where va.alignment_id=g2a.individual_alignment_id and g2a.genotype_method_id=g.id;

drop view view_individual_sequence cascade;
create or replace view view_individual_sequence as select ntt.*, t2.batch_id_list, t2.coverage_list
    from (select newt.*, t.short_name as tissue  from ( select isq.id as individual_sequence_id, isq.sequencer_id, 
    isq.sequence_type_id, isq.format, i.target_coverage, isq.coverage, isq.path, isq.filtered, isq.quality_score_format, isq.tissue_id,
    isq.read_count, 
    i.id as individual_id, i.code, i.name, i.ucla_id, i.tax_id, i.sex, i.age, i.age_cas, 
    i.approx_age_group_at_collection, i.collection_date, s.latitude, 
    s.longitude, s.id as site_id, s.short_name as site_name, s.city, s.country_id, c.name as country, isq.date_created, isq.date_updated
    from individual i, individual_sequence isq, site s, country c 
    where i.site_id=s.id and i.id=isq.individual_id and s.country_id=c.id) 
    as newt left join tissue t on t.id=newt.tissue_id) as ntt, 
    (select i.id as individual_id, string_agg(cast(b.sequence_batch_id as text),',' order by sequence_batch_id) as batch_id_list,
       string_agg(cast(b.coverage as text),',' order by coverage) as coverage_list 
       from individual i left join 
         (select ib.individual_id, ib.sequence_batch_id, b.coverage from individual2batch ib, sequence_batch b where ib.sequence_batch_id=b.id) as b
         on i.id=b.individual_id group by i.id) as t2
    where ntt.individual_id=t2.individual_id
    order by code ;

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
create or replace view view_sequence_before_after_filter as select t2.*, visq.batch_id_list, visq.coverage_list from (select i.code, i.ucla_id, t1.*
    from (select s2.individual_id, s1.id as isq_id, s1.coverage, s1.tissue_id, s1.date_created, s2.id as unfilter_isq_id, s2.coverage as unfilter_coverage, 
    float8(s1.base_count)/s2.base_count as ratio, s2.date_created as unfilter_date_created from (select * from individual_sequence where filtered =1) as s1 right join 
    (select * from individual_sequence where filtered =0) as s2 on s1.parent_individual_sequence_id=s2.id) as t1,
    individual i where i.id=t1.individual_id) as t2, view_individual_sequence visq where t2.isq_id=visq.individual_sequence_id;

-- 2013.2.1 updated with sex and relationship info
drop view view_ind2ind;
create or replace view view_ind2ind as select i2i.id as i2i_id, i2i.individual1_id, i1.code as i1_code, i1.sex as i1_sex, 
    i2i.individual2_id, i2.code as i2_code , i2.sex as i2_sex, 
    i2i.relationship_type_id, r.short_name as relationship
    from ind2ind i2i, individual i1, individual i2, relationship_type r where  i1.id=i2i.individual1_id and i2.id=i2i.individual2_id
    and r.id=i2i.relationship_type_id;

-- 2012.6.26 view on all the alignments
drop view view_seq_comparison_with_alignment;
create or replace view view_seq_comparison_with_alignment as select v.*, va1.alignment_id, va1.alignment_method_id,
    va1.ref_ind_seq_id, va1.median_depth, va1.outdated_index
    from view_sequence_before_after_filter v left join view_alignment va1 on v.isq_id=va1.isq_id;

--  2012.7.3 check whether each mate in individual_sequence_file has its corresponding entry;
drop view check_isq_file_parity;
create or replace view check_isq_file_parity as select individual_sequence_id, library, split_order, count(mate_id) 
 from (select distinct individual_sequence_id, library, split_order, mate_id from individual_sequence_file  )
 as newt group by individual_sequence_id, library, split_order 
 order by count, individual_sequence_id, library, split_order;

-- 2012.8.6 compare no_of_loci, file_size of two genotype methods side by side
select t1.genotype_method_id, t1.chromosome, t1.no_of_loci, t1.file_size, 
    t2.genotype_method_id, t2.chromosome, t2.no_of_loci, t2.file_size
    from (select * from genotype_file where genotype_method_id =10) t1 full join 
        (select * from genotype_file where genotype_method_id =12) t2 on t1.chromosome=t2.chromosome
    order by t1.no_of_loci desc;
