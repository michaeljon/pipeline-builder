---
title: Handling non-SNV variants at interval boundaries
margin-left: 0.75in
margin-top: 0.5in
---

Our current informatics pipeline preprocesses the sequence data (FASTQ files), then aligns and sorts the resulting BAM files. Those BAMs are then partitoned across intervals for quality adjustment and variant calling. Our partitioning process uses a brute-force method to define the intervals at a fixed partition point. The exception is the _last_ interval in a chromosome where an expansion factor is used to include the remaining base pairs in the interval. (For example, a final interval that is less than the factor size would be folded into the previously considered final interval. This interval will now be larger than the prior intervals. This was done to remove too-small final intervals for given interval sizes.)

The above interval construction approach introduces hard artificial boundaries into the variant calling process because each interval is called independent of other intervals. For example, `NC_000001.11` (chr1) might be partitioned as follows (50mpb intervals with a 25% combiner factor).

    NC_000001.11:1-50000000
    NC_000001.11:50000001-100000000
    NC_000001.11:100000001-150000000
    NC_000001.11:150000001-200000000
    NC_000001.11:200000001-248956422

While the alignment and sorting operations currently operate on the whole genome and input FASTQ, the variant calling process will be assigned to one of the above intervals. At the "edges" of those intervals (50, 100, 150, and 200mbp) several cases arise where incorrect variants may be called. This discussion does not take into account chromosome boundaries or at known break points (centromeres) as those regions are already under-characterized.

There are several cases to consider. Let's call the intervals here $I_1$ and $I_2$.

- **Case 1**: A pair of SNVs called, one at the end of $I_1$ and one at the beginning of $I_2$ where these two SNVs comprise an MNV instead.

- **Case 2**: An SNV called at the end of $I_1$, but a structural variant called at the start of $I_2$, whre the two variants should be combined.

- **Case 3**: A structural variant at the end of $I_1$ with an SNV called at the start of $I_2$, where the two variants should be combined.

- **Case 4**: A structural variant called at the end of $I_1$ with a second structural variant called at the beginning of $I_2$, where the two variants comprise a single variant.

- **Case 5**: A deletion overlapping the boundary that's missed as a variant because only the beginning of it is called in $I_1$.

- **Case 6**: An insertion overlapping the boundary that's missed as a variant because only the beginning of it is called in $I_1$.

Note, in cases 5 and 6 variants might be called on either side of the boundary where those variants do not match a well-known or previously identified variant.

Each of the above cases results in a "missed" or otherwise incorrect variant.

One proposal is to construct small overlapping intervals $I_O$ (of concern here is the length of that overlap, it needs to be small enough to be processed quickly but large enough to catch large structural variants in the boundary region).

The following shows one option that introduces a ±50kbp overlap (a 100kbp interval centered on the boundary) marked with `->`.

    NC_000001.11:1-50000000
    -> NC_000001.11:49950000-50050001
    NC_000001.11:50000001-100000000
    -> NC_000001.11:99950000-100050001
    NC_000001.11:100000001-150000000
    -> NC_000001.11:149950000-150050001
    NC_000001.11:150000001-200000000
    -> NC_000001.11:199950000-200050001
    NC_000001.11:200000001-248956422

We need an algorithm to reconcile conflicting variants that now arise when recombining the data. At some point looking "backward" and "forward" we should find a place where the called variants start to "match". Variants beyond this point can be pulled from either $I_1$, $I_2$ or the overlap $I_O$. Where we find differences between the numbered intervals and $I_O$ we might (and here's where the tricky part comes in) choose

- The larger of the variants
- The smaller of the variants
- Only variants in $I_O$
- To vary the selection based on the variant type and source interval

The way that we recombine our called variants causes the above problem. Variants called on $I_1$ are distinct from those called on $I_2$. When merged into the complete VCF for the sample variants that meet one of the above cases abut 

We want some of the functionality available in [bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm) to handle overlapping or adjoining variants. However this still assumes that a variant would be called at boundary positions so that `bcftools` could "normalize" it.


<!--
  pandoc -o boundary-variants.pdf -f markdown+inline_notes+yaml_metadata_block+fancy_lists --standalone -t latex boundary-variants.md
-->