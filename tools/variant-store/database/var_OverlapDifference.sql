drop function if exists var_OverlapDifference;

create function var_OverlapDifference(_primary varchar(10), _secondary varchar(10), _chrom varchar(20), _startpos bigint, _endpos bigint)
returns table (
    sample varchar(10),
    seq_id varchar(20),
    pos bigint,
    rs text,
    ref text,
    alt text
) as
$$
declare
begin
    return query (
        with A as (
            select v.sample, v.seq_id, v.pos, v.rs, v.locus, v.ref, v.alt
            from variant v
            where v.sample = _primary
              and v.seq_id = _chrom
              and (v.pos >= _startpos and v.pos < _endpos)
        ),
        P as (
            select A.sample, A.seq_id, A.pos, A.rs, A.ref, A.alt
            from A
            where A.locus not in (
                select i.locus
                from variant i
                where i.sample = _secondary
                  and i.seq_id = _chrom
                  and (i.pos >= _startpos and i.pos < _endpos)
            )
        ),
        B as (
            select v.sample, v.seq_id, v.pos, v.rs, v.locus, v.ref, v.alt
            from variant v
            where v.sample = _secondary
              and v.seq_id = _chrom
              and (v.pos >= _startpos and v.pos < _endpos)
        ),
        S as (
            select B.sample, B.seq_id, B.pos, B.rs, B.ref, B.alt
            from B
            where B.locus not in (
                select i.locus
                from variant i
                where i.sample = _primary
                  and i.seq_id = _chrom
                  and (i.pos >= _startpos and i.pos < _endpos)
            )
        ),
        R as (
            select * from P
            union all 
            select * from S
        )
        select *
        from R
        order by R.seq_id, R.pos
    );
end
$$ LANGUAGE plpgsql;