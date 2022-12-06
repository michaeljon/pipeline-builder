# Outline for boundary variants presentation

- introduction: reminder on how we process sequences

  - two-pass process
  - fragment the sequence file(s) / FASTQ by read and scatter them
  - recombine the aligned reads, sort them
  - fragment the resulting alignments
  - perform quality recalibration
  - call variants
  - merge everything back together and annotate

- where the costs are

  - fragmenting the FASTQ (the files are big)
  - aligning the reads (the more alignment processes the better)
  - the scatter / gather process
  - calling variants
  - mapping variants to "known sites" (somewhere around 1 billion sites)

- how we've called variants so far

  - naive model - call variants across the intact BAM
  - partitioning by chromosome
  - segue - batch smaller chromosomes to even the load
  - the current model - brute force
    - pre-generate partitions
    - generate partitions based on alignments

- some numbers, for some fun

  - 25mbp partition sizes
    | step | start time | stop time | seconds |
    | --------- | ---------- | ---------- | ------- |
    | calibrate | 1668890681 | 1668891052 | 371 |
    | bqsr | 1668889870 | 1668890591 | 721 |
    | call | 1668891580 | 1668895267 | 3687 |

  - 50 mbp partition sizes
    | step | start time | stop time | seconds |
    | --------- | ---------- | ---------- | ------- |
    | calibrate | 1669176745 | 1669177122 | 377 |
    | bqsr | 1669175857 | 1669176663 | 806 |
    | call | 1669178065 | 1669182661 | 4596 |

- where do those numbers come from

  - GRCh37.p13
    - 297 contigs / sequences / regions -> interval
    - 1,074,344,881 documented variants
  - GRCh38.p14
    - 705 contigs / sequences / regions -> interval
    - 1,112,554,591 documented variants
  - T2T - CHM13v2.0
    - 24 contigs / sequences / regions -> interval
    - 0(!) documented variants

- opportunities for making this faster

  - let's create more partitions
  - is there a lower bound / maximum number of partitions?
  - but calibrate and bqsr is effectively fixed
  - there might be an opportunity to reduce those (partitioning known sites)

- so, what could possibly go wrong?

  - brute force partitioning breaks things
  - it stomps on interesting "features"

    - GRCh37.p13

      - Affected genes are: AGBL4 NR5A2 EIF5B LYPD6B RBM6 TBC1D23 LINC01214 ADH5 PARP8 SYNPO CCNC LATS1 ZPBP ZCWPW1 WDFY4 R3HCC1L CNTN5 ANKS1B CAB39L UBAC2 CA10 DCC CCNB3
      - Affected genomic feature types: mRNA D_loop region type lnc_RNA gene match pseudogene transcript

    - GRCh38.p14
      - Affected genes are: AGBL4 LOC124904230 SLC35A3 OTUD7B NRXN1 AFF3 LOC124906112 RBM6 CMSS1 FILIP1L H2AZ1-DT TIGD6 MCHR2-AS1 ZPBP LOC105375423 SNTG1 RGS22 ERP44 TIMM23B-AGAP6 AGAP6 DNMBP CNTN5 RACGAP1 DLEU2 TRIM13 LINC01588 ATP8B4 MYO5B CLCN5
      - Affected genomic feature types: mRNA D_loop exon region type lnc_RNA gene centromere match pseudogene transcript

- we can't reliably call variants in those regions

  - case 1a - two snvs split at the cut point and are really an mnv
  - cases 2a and 2b - an snv and an mnv straddling the cut point
  - cases 3a, 3b, and 3c - structural variants (ins, del, indel / delins)
  - variant / structure have no match in the known sites database

- what we're going to do - experiment

  - brute force with overlaps
    - how big should the overlap be
    - where does it start and stop
    - how do we "merge" the called variants
  - smart partitioning shifting the cut point based on features
    - what type of feature is "safe" to cut
    - where do we make the cuts (± some reasonable number of bases)
  - keep in mind that the cuts are based on the reference and not the sample
    - the reference is fixed and the features "work"
    - the sample is still an unknown, so we need to be smart about shifting features

- next steps

  - carve the known sites database into fragments
    - this is a very expensive process, so we need to make it work
    - we might have a different interval set here than the fragments we create
    - there are 100s of intervals, so we need to combine some anyway
  - identify a set of safe-to-cut features
  - pre-calculate safe brute-force cut points
    - this is a smarter version of our current algorithm
  - start working on variant merge
    - identify a set of experiments to determine overlap size
    - identify a means of merging called variants
  - identify a model for collapsing smaller contigs
    - 700 intervals are too many to deal with
    - average size is around 300 kilobases, expensive to ship
    - many contigs are alternate defintions of known regions, merge them with parent
