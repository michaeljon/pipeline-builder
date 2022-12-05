drop function if exists var_Overlap;

create function var_Overlap(_primary varchar(10), _secondary varchar(10), _chrom varchar(20), _startpos bigint, _endpos bigint)
returns table (
    seq_id varchar(20),
    pos bigint,
    rs text,
    ref text,
    p_alt text,
    s_alt text
) as
$$
declare
begin
    return query (
        select p.seq_id, p.pos, p.rs, p.ref, p.alt, s.alt
        from variant p join variant s on p.sample = _primary 
                                    and s.sample = _secondary
                                    and p.locus = s.locus
        where p.seq_id = _chrom
          and (p.pos >= _startpos and p.pos < _endpos)
          and p.alt != s.alt
        order by p.seq_id, p.pos
    );
end
$$ LANGUAGE plpgsql;