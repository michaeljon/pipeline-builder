create database variant_store;

drop table if exists variant_reference;
create unlogged table variant_reference (
    sequence_id character varying(20),
    position integer,
    ref text,
    alt text,
    hgvs_root character varying(256),
    hgvs character varying(8192),
    variant_type character varying(100),
    common boolean,
    gene_info character varying(8192),
    rs integer
);

copy variant_reference from '/Users/michaeljon/src/ovation/variant-store/variant_reference.tsv' delimiter E'\t' csv header;

create index ndx_variant_ref_seq on variant_reference (sequence_id);
create index ndx_variant_ref_seq_pos on variant_reference (sequence_id, position);
create index ndx_variant_ref_hgvs on variant_reference (hgvs);
create index ndx_variant_ref_hgvs_root on variant_reference (hgvs_root);
create index ndx_variant_ref_rs on variant_reference (rs);
create index ndx_variant_ref_gene on variant_reference (gene_info);
create index ndx_variant_ref_seq_pos_alt on variant_reference (sequence_id, position, alt);
create index ndx_variant_ref_gene_vartype on variant_reference (gene_info, variant_type);

drop table if exists genome_features;
create unlogged table genome_features (
    sequence_id character varying(20),
    feature_type character varying(100),
    start_position int,
    stop_position int,
    id character varying(200),
    name character varying(200),
    description character varying(200),
    gene character varying(200),
    start_key bigint,
    stop_key bigint,
    local_region int4range,
    global_region int8range
);

copy genome_features from '/Users/michaeljon/src/ovation/variant-store/genome_features.tsv' delimiter E'\t' csv header;

update genome_features set local_region = int4range(start_position, stop_position, '[]');
update genome_features set global_region = int8range(start_key, stop_key, '[]');

create index ndx_feature_seq on genome_features (sequence_id);
create index ndx_feature_seq_start_stop on genome_features (sequence_id, start_position, stop_position);
create index ndx_feature_id on genome_features (id);

create index ndx_feature_gene on genome_features (gene);
create index ndx_feature_feature_type on genome_features (feature_type);

create index ndx_feature_local on genome_features using gist (local_region);
create index ndx_feature_global on genome_features using gist (global_region);

with 
s as (select length(hgvs) / 10 bin from variant_reference),
t as (
    select s.bin bin, count(*) instances from s group by 1
    union all
    select n, 0 from generate_series(0, 75) n
)
select bin * 10 bin, max(instances)
from t
group by 1
order by 1;


select sequence_id, 
       feature_type, 
       id, 
       name, 
       gene,
       start_position, 
       stop_position 
from genome_features 
where sequence_id = 'NC_000001.11' 
  and 43440063 <@ local_region 
order by start_position, stop_position desc;

select sequence_id, 
       feature_type, 
       id, 
       name, 
       gene,
       start_position, 
       stop_position 
from genome_features 
where sequence_id = 'NC_000001.11' 
  and start_position <= 30380 
  and 30380 <= stop_position 
order by start_position, stop_position desc;


select sequence_id, 
       feature_type, 
       id, 
       name, 
       gene,
       start_position, 
       stop_position 
from genome_features 
where sequence_id = 'NC_000001.11' 
  and '[30380,30480]'::int4range && local_region
order by start_position, stop_position;


select sequence_id, 
       feature_type, 
       id, 
       name, 
       gene,
       start_position, 
       stop_position 
from genome_features 
where sequence_id = 'NC_000001.11' 
  and (30380 <= stop_position and start_position <= 30480)
order by start_position, stop_position;

select g.sequence_id, 
       g.feature_type, 
       g.id, 
       c.position,
       g.start_position, 
       g.stop_position,
       g.gene,
       c.ref,
       c.alt,
       c.variant_type,
       c.hgvs,
       c.rs
from genome_features g, variant_reference c
where g.sequence_id = 'NC_000001.11'
  and g.sequence_id = c.sequence_id 
  and g.local_region @> c.position
  and c.rs = 1570700741
limit 100;

-- by rs #
select g.sequence_id, 
       g.feature_type, 
       g.id, 
       c.position,
       g.start_position, 
       g.stop_position,
       g.gene,
       c.hgvs,
       c.ref,
       c.alt,
       c.rs
from genome_features g, variant_reference c
where c.sequence_id = 'NC_000001.11'
  and g.sequence_id = c.sequence_id 
  and g.local_region @> c.position
  and c.rs = 1570700741
limit 100;

-- by hgvs
select g.sequence_id, 
       g.feature_type, 
       g.id, 
       c.position,
       g.start_position, 
       g.stop_position,
       g.gene,
       c.hgvs_root,
       c.ref,
       c.alt,
       c.rs
from genome_features g, variant_reference c
where c.sequence_id = 'NC_000001.11'
  and g.sequence_id = c.sequence_id 
  and g.local_region @> c.position
  and c.hgvs_root = 'NC_000001.11:g.100187648'
order by g.start_position, g.stop_position
limit 100;

select sequence_id, 
       feature_type, 
       id, 
       name, 
       gene,
       start_position, 
       stop_position 
from genome_features 
where sequence_id = 'NC_000001.11' 
  and '[43438699,43438817]'::int4range && local_region
order by start_position, stop_position;
