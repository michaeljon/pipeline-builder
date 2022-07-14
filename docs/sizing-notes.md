# Some sizing notes

The basis on which we're working

- Each individual genomic dataset is approximately $3\times10^9$ base pairs
- Each individual shares approximately 99% of their genome leaving approximately $3\times10^7$ variations
- We base our variation discovery off a _reference genome_, in our case GRCh38.p13

Assuming short read (151bp), paired sequence data, each read captured in the FASTQ comprises

- One line capturing the sequence identifier and metadata of approximately 70 characteres
- One line capturing the individual bases within the read (the ACTGs) at 151 characters
- One line of sequence quality data captured per read at 151 characters
- Some additional trivial housekeeping data

At 30x read depth / coverage, we assume an average of 30 reads per base pair, which gives us approximately $100\times10^9$ bases recorded across the FASTQ file pair. Note that 30x is considered the low-end of usefulness and we are likely looking at 50x or 100x coverage for our production runs. Our example FASTQs, at 30x depth, have $1.6\times10^9$ lines per file and we have two FASTQ files per sequence meaning we have $3.2\times10^9$ lines and $8\times10^8$ "reads" to work with. The compressed size on disk is approximately 30Gb per file; uncompressed they are around 154Gb for a nearly 5:1 compression ratio.

Some additional back-of-the-napkin math tells us that we'll see around 1Gb per "read depth" for a 30x sequence and that we're only covering about $1/3$ of the entire genome.

(For the following I'm ignoring structural variants, they do not significantly impact file size, but they do drive what the assembled sequence looks like and how it compares to the reference.)

Ignoring the variant call file's headers and metadata (noise compared to the overall file size and accounting for about 3500 lines of text) and assuming an unannotated file, we might see approximately $4.5\times10^6$ variant lines. Of these a variant without annotations has an average line length of approximately 256 characters. This results in an uncompressed file a little over 900Mb. Compressing this (using gzip) results in a 150Mb file, which is a 6:1 compression ratio.

Annotating variants leaves the variant entry count the same as it's simply adding to each "row" the known information about that variation. However, each annotation can inject around 4kb which is approximately 16 times more information. That gives us an uncompressed size of around 14Gb.