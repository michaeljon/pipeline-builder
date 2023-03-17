BEGIN { 
    OFS="\t"
    CONVFMT="%.17g"
    PREC=100
}

/^##contig=/ { 
    match($0, /ID=([A-Z_0-9.]+),length=([0-9]+)>/, m)
    ss[m[1]] = m[2]
}

! /^#/ { 
    counts[$1] += 1;

    if (length($4) == 1 && length($5) == 1) {
        if ($3 ~ /\./) n_snv[$1] += 1; else snv[$1] += 1;
    } else if (length($4) == length($5)) {
        if ($3 ~ /\./) n_mnv[$1] += 1; else mnv[$1] += 1;
    } else if (length($4) < length($5)) {
        if ($3 ~ /\./) n_ins[$1] += 1; else ins[$1] += 1;
    } else {
        if ($3 ~ /\./) n_del[$1] += 1; else del[$1] += 1;
    }
} 

! /^#/ && ! ($3 ~ /\./) { 
    rs[$1] += 1;
} 

END { 
    print "sequence", "sample", "sequence", "seq-size", "variant-count", "rs-assigned", "rs-unassigned", "snv", "mnv", "ins", "del", "unk-snv", "unk-mnv", "unk-ins", "unk-del", "base/var", "var/base", "normal"

    for (a in counts) {
        ns[a] = ss[a]/counts[a]
        # ns[a] = counts[a]/ss[a]
    }

    mn = min(ns)
    mx = max(ns)

    for (a in counts) 
        print SEQUENCE, SAMPLE, a, ss[a], counts[a], rs[a], counts[a] - rs[a], snv[a] + 0, mnv[a] + 0, ins[a] + 0, del[a] + 0, n_snv[a] + 0, n_mnv[a] + 0, n_ins[a] + 0, n_del[a] + 0, ss[a]/counts[a], counts[a]*1.0/ss[a]*1.0, normal(ns, a, mn, mx) * 100.0
}

function min(ar,  m) {
    m = 9999999999999
    for (a in ar) {
        if (ar[a] < m) {
            m = ar[a];
        }
    }
    return m
}

function max(ar,  m) {
    m = -1 
    for (a in ar) {
        if (ar[a] > m) m = ar[a];
    }
    return m
}

function normal(ar, a, mn, mx) {
    return (ar[a] - mn) / (mx - mn)
}