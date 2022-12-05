create database variant_store;

-- copy variant_reference from '$LOC/GRCh37.p13.variants.tsv' delimiter E'\t' csv header;
-- copy variant_reference from '$LOC/GRCh38.p14.variants.tsv' delimiter E'\t' csv header;

-- copy genome_features from '$LOC/GRCh37.p13.features.tsv' delimiter E'\t' csv header;
-- copy genome_features from '$LOC/GRCh38.p14.features.tsv' delimiter E'\t' csv header;

drop table if exists sample_variant;
create table sample_variant (
    sample character varying(20),
    reference character varying(20),
    sequence_id character varying(20),
    position integer,
    rs text,
    ref text,
    alt text,
    locus character varying(50)
);

drop table if exists variant_reference;
create unlogged table variant_reference (
    reference character varying(20),
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

drop table if exists genome_features;
create unlogged table genome_features (
    reference character varying(20),
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

--
-- sample
--
create index ndx_sample on sample_variant (sample);
create index ndx_sample_position on sample_variant (sample, reference, sequence_id, position);
create index ndx_position on sample_variant (reference, sequence_id, position);
create index ndx_locus on sample_variant (sample, reference, locus);

--
-- variant reference
--
create index ndx_variant_ref_ref on variant_reference (reference);
create index ndx_variant_ref_refseq on variant_reference (reference, sequence_id);
create index ndx_variant_ref_refseq_pos on variant_reference (reference, sequence_id, position);
create index ndx_variant_ref_refhgvs on variant_reference (reference, hgvs);
create index ndx_variant_refref_hgvs_root on variant_reference (reference, hgvs_root);
create index ndx_variant_ref_rs on variant_reference (rs);
create index ndx_variant_ref_gene on variant_reference (gene_info);
create index ndx_variant_ref_refseq_pos_alt on variant_reference (reference, sequence_id, position, alt);
create index ndx_variant_ref_gene_vartype on variant_reference (gene_info, variant_type);

--
-- genomic structure
--
update genome_features set local_region = int4range(start_position, stop_position, '[]');
update genome_features set global_region = int8range(start_key, stop_key, '[]');

create index ndx_feature_refseq on genome_features (reference, sequence_id);
create index ndx_feature_refseq_start_stop on genome_features (reference, sequence_id, start_position, stop_position);
create index ndx_feature_refid on genome_features (reference, id);

create index ndx_feature_refgene on genome_features (reference, gene);
create index ndx_feature_refgene_type on genome_features (reference, gene, feature_type);
create index ndx_feature_reffeature_type on genome_features (reference, feature_type);

create index ndx_feature_local on genome_features using gist (local_region);
create index ndx_feature_global on genome_features using gist (global_region);
