# Stuff...

We had to monkey-patch vep (for now) `~/bin/ensembl-vep/modules/Bio/EnsEMBL/VEP/Stats.pm` at lines 912 and 925.

```perl
join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} sort {
    my $aterm = $a =~ s/chr//r;
    my $bterm = $b =~ s/chr//r;
    $aterm = ord($aterm) unless $aterm =~ /^\d+$/;
    $bterm = ord($bterm) unless $bterm =~ /^\d+$/;
    return $aterm <=> $bterm;
} @keys),
```