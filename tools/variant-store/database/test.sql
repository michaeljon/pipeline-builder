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


with A as (
    select sample, sequence_id, position, rs, locus, ref, alt
    from variant
    where sample = 'platinum'
      and sequence_id = 'chr5'
      and (position >= 10000 and position < 99999)
),
P as (
    select sample, sequence_id, position, rs, ref, alt
    from A
    where locus not in (
        select locus
        from variant
        where sample = 'sampled'
          and sequence_id = 'chr5'
          and (position >= 10000 and position < 99999)
    )
),
B as (
    select sample, sequence_id, position, rs, locus, ref, alt
    from variant
    where sample = 'sampled'
      and sequence_id = 'chr5'
      and (position >= 10000 and position < 99999)
),
S as (
    select sample, sequence_id, position, rs, ref, alt
    from B
    where locus not in (
        select locus
        from variant
        where sample = 'platinum'
          and sequence_id = 'chr5'
          and (position >= 10000 and position < 99999)
    )
),
R as (
    select * from P
    union all 
    select * from S
)
select *
from R
order by sequence_id, position


with C as
(
    select sample, sequence_id, position, rs, locus, ref, alt
    from variant
    where sample = 'sampled'
      and sequence_id = 'chr22'
      and (position >= 1 and position < 100000)
)
select sample, sequence_id, position, rs, ref, alt
from C
where locus not in (
    select locus
    from variant
    where sample = 'platinum'
)

select sample, sequence_id, 1 + (position / 25000000 )::int as intv, count(*) 
from variant 
group by 1, 2, 3
order by 2, 3, 1;

with a as (
    select sequence_id, 1 + (position / 25000000 )::int intv, count(*) c
    from variant 
    group by 1, 2
)
select sequence_id, intv, avg(c)::int
from a
group by 1, 2
order by 1, 2

with a as (
    select sequence_id, 1 + (position / 25000000 )::int intv, count(*) calls, min(position) minimum, max(position) maximum
    from variant 
    group by 1, 2
)
select sequence_id, intv, calls, minimum, maximum
from a
order by 1, 2

copy (with a as (
    select sequence_id, 1 + (pos / 25000000 )::int intv, count(*) calls, min(position) minimum, max(position) maximum
    from variant 
    group by 1, 2
)
select sequence_id, intv, calls, minimum, maximum, (intv - 1) * 25000000 + 1 intv_low, (intv) * 25000000 intv_high
from a
order by 1, 2) to '/Users/michaeljon/vcf/spread.csv' csv header;

with x as (
    select sample, sequence_id, count(*) cnt
    from variant 
    group by 1, 2
),
y as (
    select p.sample p_sample, q.sample q_sample, p.sequence_id, 
           p.cnt platinum, 
           q.cnt sampled,
           abs(p.cnt - q.cnt) diff, 
           case substring(p.sequence_id from 4) 
            when 'X' then 23 
            when 'Y' then 24 
            else substring(p.sequence_id from 4)::int 
           end chr
    from x p right join x q on p.sequence_id = q.sequence_id
)
select sequence_id, 
       chr,
       case when platinum > sampled then 'platinum' else 'sampled' end larger, 
       platinum platinum_variants,
       sampled sampled_variants,
       diff delta
from y
where diff != 0
  and p_sample = 'platinum'
order by chr



with a as (
    select sample,            
           case substring(sequence_id from 4) 
            when 'X' then 23 
            when 'Y' then 24 
            else substring(sequence_id from 4)::int 
           end chr, 
           count(*) calls, min(pos) minimum, max(pos) maximum
    from variant 
    group by 1, 2
),
b as (
    select sample, chr, calls p_calls, minimum p_min, maximum p_max, 0 s_calls, 0 s_min, 0 s_max
    from a
    where sample = 'platinum'

    union all

    select sample, chr, 0 p_calls, 0 p_min, 0 p_max, calls s_calls, minimum s_min, maximum s_max
    from a
    where sample = 'sampled'
),
c as (
    select chr, 
        sum(p_calls) p_calls, sum(s_calls) s_calls, 
        sum(p_min) p_min_loc, sum(s_min) s_min_loc, 
        sum(p_max) p_max_loc, sum(s_max) s_max_loc
    from b
    group by 1
)
select chr, 
       p_calls, s_calls,
       abs(p_calls - s_calls) diff,
       p_min_loc, s_min_loc,
       p_max_loc, s_max_loc
from c





select p.sample, p.sequence_id, p.position, p.rs, p.locus, p.ref, p.alt, s.ref, s.alt
from variant p join variant s on p.sample = 'platinum' 
                             and s.sample = 'sampled'
                             and p.locus = s.locus
where p.sequence_id = 'chr5'
  and (p.position >= 1 and p.position < 100000)





with c as (
    select sequence_id, position, count(*) cnt 
    from variant 
    group by sequence_id, position
    ) 
select cnt, count(*) 
from c 
group by cnt 
order by cnt asc;

select *
from crosstab(
    'select sequence_id, sample, count(*) as variant_count
     from variant 
     group by sequence_id, sample
     order by sequence_id, sample',
    'select distinct sample from variant order by 1'
) as (
    sequence_id varchar(20),
    "zr8912_1_S69" varchar(20),
    "zr8912_2_S70" varchar(20),
    "zr8912_3_S71" varchar(20),
    "zr8912_4_S72" varchar(20),
    "zr8937_1_S73" varchar(20),
    "zr8937_2_S74" varchar(20),
    "zr8937_3_S75" varchar(20),
    "zr8937_4_S76" varchar(20)
)


-- contributions
select *
from crosstab(
    $$
        with singles as (
            select v.sequence_id, v.position --, v.alt
            from variant v
            group by v.sequence_id, v.position --, v.alt 
            having count(*) = 7
        ),
        samples as (
            select v.sequence_id, v.sample, count(*) as variant_count
            from variant v join singles s on v.sequence_id = s.sequence_id and v.position = s.position -- and s.alt = any(string_to_array(v.alt, ','))
            group by v.sequence_id, v.sample
            order by v.sequence_id, v.sample
        )
        select *
        from samples
    $$,
    $$
        select samples 
        from ( 
            values
                ('zr8912_1_S69'), ('zr8912_2_S70'), ('zr8912_3_S71'), ('zr8912_4_S72'),
                ('zr8937_1_S73'), ('zr8937_2_S74'), ('zr8937_3_S75'), ('zr8937_4_S76')
            ) s (samples)
    $$ 
) as (
    sequence_id varchar(20),
    zr8912_1_S69 int,
    zr8912_2_S70 int,
    zr8912_3_S71 int,
    zr8912_4_S72 int,
    zr8937_1_S73 int,
    zr8937_2_S74 int,
    zr8937_3_S75 int,
    zr8937_4_S76 int
)


select *
from crosstab(
    'with singles as (
        select v.sequence_id, v.position
        from variant v
        group by v.sequence_id, v.position 
        having count(*) = 8
    ),
    samples as (
        select v.sequence_id, v.sample, count(*) as variant_count
        from variant v join singles s on v.sequence_id = s.sequence_id and v.position = s.position and v.alt = any(string_to_array(r.alt, ','))
        group by v.sequence_id, v.sample
        order by v.sequence_id, v.sample
    )
    select *
    from samples',
    'select distinct sample from variant order by 1'
) as (
    sequence_id varchar(20),
    "zr8912_1_S69" int,
    "zr8912_2_S70" int,
    "zr8912_3_S71" int,
    "zr8912_4_S72" int,
    "zr8937_1_S73" int,
    "zr8937_2_S74" int,
    "zr8937_3_S75" int,
    "zr8937_4_S76" int
)



select *
from crosstab(
    'select v.sequence_id, v.sample, count(*) as variant_count
     from variant v
     group by v.sequence_id, v.sample
     order by v.sequence_id, v.sample',
    'select distinct sample from variant order by 1'
) as (
    sequence_id varchar(20),
    "zr8912_1_S69" int,
    "zr8912_2_S70" int,
    "zr8912_3_S71" int,
    "zr8912_4_S72" int,
    "zr8937_1_S73" int,
    "zr8937_2_S74" int,
    "zr8937_3_S75" int,
    "zr8937_4_S76" int
)



with singles as (
    select v.sequence_id, v.position
    from variant v
    group by v.sequence_id, v.position 
    having count(*) = 1
),
samples as (
    select v.sequence_id, v.sample, v.position
    from variant v join singles s on v.sequence_id = s.sequence_id and v.position = s.position
),
genes as (
    select distinct 
        s.sample,
        g.sequence_id, 
        g.gene,
        g.start_position, 
        g.stop_position
    from genome_features g join samples s on s.sequence_id = g.sequence_id and s.position <@ g.local_region 
    where g.feature_type = 'gene'
    order by g.start_position, g.stop_position desc
)
select sample, sequence_id, gene
from genes;



with samples as (
    select v.sequence_id, v.sample, v.position, v.alt
    from variant v
    where v.sample = 'zr8937_3_S75'
),
genes as (
    select distinct 
        s.sample,
        g.sequence_id, 
        g.gene,
        g.start_position, 
        g.stop_position
    from genome_features g join samples s on s.sequence_id = g.sequence_id and s.position <@ g.local_region 
    where g.feature_type = 'gene'
    order by g.start_position, g.stop_position desc
)
select sample, sequence_id, gene
from genes;


select v.sample, v.sequence_id, v.position, r.ref, v.alt, g.gene, r.rs
from variant v join genome_features g on v.sequence_id = g.sequence_id and v.position <@ g.local_region
               left join variant_reference r on v.sequence_id = r.sequence_id and v.position = r.position and v.alt = r.alt
where v.sample = 'zr8937_3_S75'
  and g.feature_type = 'gene'
limit 1000;







select *
from crosstab(
$$
select (v.sample || ':' || g.gene)::varchar(100), coalesce(r.variant_type, 'other'), coalesce(count(*), 0)::decimal
from variant v join genome_features g on v.sequence_id = g.sequence_id and v.position <@ g.local_region
               left join variant_reference r on v.sequence_id = r.sequence_id and v.position = r.position
where g.gene in (
'ATG16L1', 'CARD15', 'CDKAL1', 'CNR1', 'CRP', 'CX3CR1', 'CXCL9', 'DLG5', 'FHIT', 'FcRL3', 'HLA', 'HLA-B*35',
'HLA-G', 'HLAB*27HLA-B*44', 'HLADRB*103', 'HLADRB1', 'HSP70-2', 'IBD5', 'IL-10', 'IL-12B', 'IL-6',
'IL23R', 'IRGM', 'LRRK2', 'LTA', 'MDR1', 'MHC', 'MIF', 'NAT',  'NOD2', 'OCTN', 'PAI-1', 'STAT3', 'TCF-4(TCF7L2)',
'TLR', 'TLR1', 'TLR1-2', 'TLR1-6', 'TLR4', 'TNF', 'TNF-1031C', 'TNFRSF6B', 'TNFα-308A', 'TPMT', 'ZNF365')
  and g.feature_type = 'gene'
group by v.sample || ':' || g.gene, r.variant_type
$$
) as (
    sample varchar(100),
    del decimal,
    indel decimal,
    ins decimal,
    snv decimal,
    other decimal
)


select *
from crosstab(
$$
    select v.sample::varchar(100), g.gene, coalesce(count(*), 0)::decimal
    from variant v join genome_features g on v.sequence_id = g.sequence_id and v.position <@ g.local_region
                left join variant_reference r on v.sequence_id = r.sequence_id and v.position = r.position and v.alt = any(string_to_array(r.alt, ','))
    where g.gene in (
    'CDKAL1', 'CNR1', 'CRP', 
    'CX3CR1', 'CXCL9', 'DLG5', 'FHIT', 'HLA-G', 'IL23R', 'IRGM', 
    'LRRK2', 'NOD2', 'STAT3', 'TLR1', 'TLR4', 'TNFRSF6B', 'TPMT', 'ZNF365')
    and g.feature_type = 'gene'
    and r.variant_type = 'indel'
    group by v.sample, g.gene
    order by v.sample, g.gene
$$,
$$
    select genes 
    from ( 
        values ('CDKAL1'), ('CNR1'), ('CRP'),
               ('CX3CR1'), ('CXCL9'), ('DLG5'), ('FHIT'), ('HLA-G'), ('IL23R'), ('IRGM'),
               ('LRRK2'), ('NOD2'), ('STAT3'), ('TLR1'), ('TLR4'), ('TNFRSF6B'), ('TPMT'), ('ZNF365')
        ) g(genes)
$$     
) as (
    sample varchar(100),
    "CDKAL1" decimal, "CNR1" decimal, "CRP" decimal, 
    "CX3CR1" decimal, "CXCL9" decimal, "DLG5" decimal, "FHIT" decimal, "HLA-G" decimal, "IL23R" decimal, "IRGM" decimal, 
    "LRRK2" decimal, "NOD2" decimal, "STAT3" decimal, "TLR1" decimal, "TLR4" decimal, "TNFRSF6B" decimal, "TPMT" decimal, "ZNF365" decimal
)


with genes as (
     select g.reference, g.sequence_id, g.gene, g.local_region
     from genome_features g
     where g.reference = 'GRCh37.p13'
       and g.gene in ( 'BRCA1', 'BRCA2' )
       and g.feature_type = 'gene'
 )
 select g.gene, coalesce(count(*), 0)::decimal
 from sample_variant s join genes g on s.reference = g.reference and s.sequence_id = g.sequence_id and s.position <@ g.local_region
                       left join variant_reference r on g.reference = r.reference and g.sequence_id = r.sequence_id and s.position = r.position and s.alt = any(string_to_array(r.alt, ','))
 where s.sample = 'mj'
 group by g.gene
 order by g.gene


select g.reference, g.sequence_id, g.gene, s.position, g.local_region, r.ref, s.ref, s.alt
from sample_variant s join genome_features g on s.reference = g.reference and s.sequence_id = g.sequence_id and s.position <@ g.local_region
                      left join variant_reference r on g.reference = r.reference and s.rs = r.rs || ''
where g.reference = 'GRCh37.p13'
  and g.gene in ( 'BRCA1', 'BRCA2' )
  and g.feature_type = 'gene'
  and s.sample = 'mj'

