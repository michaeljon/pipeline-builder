---
title: Handling non-SNV variants at interval boundaries
author: MJ Miller (mj@ovation.io)
margin-left: 1.5in
abstract: |
  Partitioning variant calling to decrease costs can miss or incorrectly call variants that intersect or otherwise overlap with an interval boundary. This paper describes the approaches Ovation has taken to calling variants that introduce the boundary problem and proposes a solution to that problem which maintains our low compute cost / low time to process.
  \newpage
---

## The situation

There are two coarse-grained partitioning steps in the Ovation variant calling pipeline. The first is an optimization for performing read alignment to the reference genome. Briefly this optimization selects batches of reads and distributes those reads to many individual alignment processes. Each process aligns the reads to the reference and returns the alignment details. Those aligned reads are gathered after all reads have been processed^[We say these reads are query-sorted, but it's probably safer to think of them as unordered. A non-partitioning alignment would result in an output that is essentially one-to-one with the reads found in the source FASTQ data.], sorted, and a duplicate detection / marking process is applied. The sorting process results in a what is known as coordinate-sorted reads. That is the ordering of the reads follows a monotonically-increasing ordering by [ chromosome, position ]. The output file is in [BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).

It is this data that we'll turn to for the rest of this paper.^[Note, this paper mainly talks about WGS informatics. While many of the issues apply to RNA-seq we won't cover any RNA-seq topics.]

Before talking about the current approach, issues with that approach, and possible solutions, a quick review of other partitioning models that Ovation has used might be helpful.

### Naive model - call variants on entire BAM

Our first alignment processes used no partitioning at all. We provided the entire BAM to our alignment tooling [GATK](https://gatk.broadinstitute.org/hc/en-us) and ran that pipeline on a single EC2 instance. This model should be considered the "most correct" in that it doesn't suffer from any of the other partitioning problems that we'll talk about in the rest of the paper. However, GATK is essentially a single-threaded application and there are three sub-steps that must be taken to go from BAM to called variants. Each of those steps (determination of base quality, application and normalization of base quality, variant calling against the reference and a set of known variants) must run sequentially and completely over the BAM. In this case we were effectively "wasting" a large EC2 instance and spending upwards of 24 hours to call variants _per sample_.

### Partition by chromosome

With that in mind the first partitioning approach was to extract all reads for each chromosome into individual batches. We took those 23 files and ran the same GATK pipeline against each, but in parallel on a large EC2 instance. While this reduced the time to process the sample there was still a large "waste" of CPU resources because a few large chromosomes would dominate the processing time. While the smaller ones completed early the overall variant pipeline (gathering all variants into a single output file and annotating those variants) completion would be delayed and most of the CPUs would be idle. In this case we might still consider this partitioning scheme to result in a "correct" variant calling outcome because we don't expect to call variants that cross chromosome "boundaries".

### Even the workload - "smart" partitioning

The next approach was to combine chromosome batches to create most evenly-sized input files.^[There are two ways to think of "size" here. The first is the number of base pairs in an input file which is easily determined from the reference data. The second is the number of reads to process and to use as inputs into the variant calling process. This number is generally non-deterministic but is related to the chromosome size in that a larger chromosome _should have_ a corresponding larger number of reads.] Performing the batching in a deterministic way meant basing the input sizes on the base pair size of each chromosome and then doing a one-time bin packing process to find the best and most even distribution. However, this model had the side effect of causing the variant calling process to still take as long as the longest input file, but would idle more CPUs than the per-chromosome model because there were fewer input files.

### The current model - brute force partitioning

Our current production partitioning process uses a brute-force method to define the intervals at a fixed partition point. The exception is the _last_ interval in a chromosome where an expansion factor is used to include the remaining base pairs in the interval. For example, a final interval that is less than the factor size would be folded into the previously considered final interval. This interval will now be larger than the prior intervals. This was done to remove too-small final intervals for given interval sizes.

Our default partitioning uses a 25mbp partiting size with a 0.25 load factor for the final partition. This results in 133 intervals over the 23 chromosomes. We then split the GATK pipeline into three distinct steps and parallelize those individual steps. In this way we're able to keep a 72 core EC2 running at nearly 100% capacity and we reduce the overall variant. For a low read depth FASTQ we might see the following variant calling times

| step      | start time | stop time  | seconds |
| --------- | ---------- | ---------- | ------- |
| calibrate | 1668890681 | 1668891052 | 371     |
| bqsr      | 1668889870 | 1668890591 | 721     |
| call      | 1668891580 | 1668895267 | 3687    |

With this approach we're able to quickly and cheaply call variants. In the example above we can call variants on a 12x whole genome in around 80 real minutes.

```sh
for f in *.calibrate.log *.bqsr.log *.call.log; do
	min=$(awk '/^[0-9]/{printf("%f\n", $3+$4)}' $f | head -n 1 | sed -r 's/\.[0-9]+$//')
	max=$(awk '/^[0-9]/{printf("%f\n", $3+$4)}' $f | tail -n 1 | sed -r 's/\.[0-9]+$//')
	seconds=$(expr $max - $min)

	echo \| $f\| $min\| $max\| $seconds \|
done
```

### The problem

The above interval construction approach introduces hard artificial boundaries into the variant calling process because each interval is called independent of other intervals. For example, `NC_000001.11` (chr1) might be partitioned as follows (50mpb intervals with a 25% combiner factor).

    NC_000001.11:1-50000000
    NC_000001.11:50000001-100000000
    NC_000001.11:100000001-150000000
    NC_000001.11:150000001-200000000
    NC_000001.11:200000001-248956422

While the alignment and sorting operations currently operate on the whole genome and input FASTQ, the variant calling process will be assigned to one of the above intervals. At the "edges" of those intervals (50, 100, 150, and 200mbp) several cases arise where incorrect variants may be called. This discussion does not take into account chromosome boundaries or at known break points (centromeres) as those regions are already under-characterized.

There are several cases to consider. Let's call the intervals here $I_1$ and $I_2$.

- **Case 1**: A pair of SNVs called, one at the end of $I_1$ and one at the beginning of $I_2$ where these two SNVs comprise an MNV instead.

  ![Case 1](../boundary-images/case-1)

- **Case 2**: An SNV called at the end of $I_1$, but a structural variant called at the start of $I_2$, whre the two variants should be combined.

  ![Case 2](../boundary-images/case-2)

- **Case 3**: A structural variant at the end of $I_1$ with an SNV called at the start of $I_2$, where the two variants should be combined.

  ![Case 3](../boundary-images/case-3)

- **Case 4**: A structural variant called at the end of $I_1$ with a second structural variant called at the beginning of $I_2$, where the two variants comprise a single variant.

  ![Case 4](../boundary-images/case-4)

- **Case 5**: A deletion overlapping the boundary that's missed as a variant because only the beginning of it is called in $I_1$.

  ![Case 5](../boundary-images/case-5)

- **Case 6**: An insertion overlapping the boundary that's missed as a variant because only the beginning of it is called in $I_1$.

  ![Case 6](../boundary-images/case-6)

Note, in cases 5 and 6 variants might be called on either side of the boundary where those variants do not match a well-known or previously identified variant.

Each of the above cases results in a "missed" or otherwise incorrect variant.

### What we're going to do

One proposal is to construct small overlapping intervals $I_O$ (of concern here is the length of that overlap, it needs to be small enough to be processed quickly but large enough to catch large structural variants in the boundary region).

![Overlapped interval](../boundary-images/overlapped)

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
  pandoc -o boundary-variants.pdf -f markdown+inline_notes+yaml_metadata_block+fancy_lists --lua-filter=./pagebreak.lua --standalone -t latex --default-image-extension=pdf --resource-path=.:boundary-images boundary-variants.md
-->
