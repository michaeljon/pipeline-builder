# How to build index files

If / when we pull a new reference, for example GRCh38.p14 or CHM13, we'll need to construct a number of indexes for the tools. For the most part this process applies to all references. In some cases the pipelines will construct additional indexes during their bootstrap phase (these are very lightweight and fast, so pre-building them doesn't make the overall process faster).

For now we "install" our references in `~/reference/${ORGANISM}`. The current list includes two human genomes (GRCh38.p13 and .p14), SARS-CoV-2, and our four endemic coronaviruses (HCoV-229E, HCoV-HKU1, HCoV-NL63, and HCoV-OC43). The `--reference` parameter passed to the builder should point at this "top-level" directory.

### BWA indexes

First, we create our main index. This one feeds `bwa-mem2` which is our primary aligner.

```bash
bwa-mem2 index GCF_000001405.40_GRCh38.p14_genomic.fna
```

That command will build the BWT indexes from the reference. The filename for the reference will be used to construct the index file names: `.0123`, `.amb`, `.ann`, `.pac`, and `.bwt.2bit.64`. This last file is based on the processor installed on the EC2 instance and its support for various Intel extensions.

This particular step can take upwards of an hour depending on the EC2 instance. This is why we build the files and burn them into the AMI. We can afford the download / launch costs for the instance far easier than we can the build time (plus, we can reuse these files for each alignment). The next step, creating the index for `bwa` can also take about an hour.

We also create a set of index files for our "old" aligner `bwa` for comparison.

```bash
bwa index GCF_000001405.40_GRCh38.p14_genomic.fna
```

In this case some of the same files are constructed, but there are two new ones: `.sa` and `.bwt`.

### hisat and minimap indexes (for virus data)

For virome analysis we have a number of alternate aligners. BWA has already been covered above, but we need to build index files for HISAT and MiniMap. This process, like the others for virus data, is extremely fast. If you blink you might miss it. These two processes generate files named `*.[1-8].ht2` and `*.mmi` to match the reference name

```bash
hisat2-build -p 16 reference.fasta reference
minimap2 -d reference.mmi reference.fasta
```


### samtools indexes

```bash
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna
```

This command creates a single index file with a `.fai` extension.

### The GATK sequence dictionary

This index can be created by the pipeline, but if running many pipelines in parallel on the same EC2 you'll want to construct this ahead of time to avoid race conditions.

```bash
gatk CreateSequenceDictionary \
    -R GCF_000001405.40_GRCh38.p14_genomic.fna \
    -O GCF_000001405.40_GRCh38.p14_genomic.dict
```

The output filename is used directly by the pipeline tooling so it needs to follow the above convention directly.

### The GATK interval list

This file is built per-organism right now. For humans the build is based off the chromosome list that's stored in the associated JSON file. For other, simpler, organisms we're hard-coding the build step in the pipeline itself. For those organisms that have a single reference in their associated FASTA (the coronaviruses for example), the build looks something like this:

```bash
    egrep '(MN908947.3)\\s' ${REFERENCE}/GCA_009858895.3_ASM985889v3_genomic.fna.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat ${REFERENCE}/GCA_009858895.3_ASM985889v3_genomic.dict - >${REFERENCE}/ref_genome_autosomal.interval_list
```

But for the human FASTA which contain all of the individual chromosome references the process involves parsing the chromosome "names" from the JSON to build a larger regular expression.

First, in python

```python
    chromosomes = ["chr" + str(c) for c in loadChromosomeList(options["chromosomeSizes"])]
    regex = "|".join(chromosomes)
```

Then, later inside the pipeline itself

```bash
    egrep '({REGEX})\\s' {REFERENCE}/{ASSEMBLY}.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/{ASSEMBLY}.dict - >{REFERENCE}/{ASSEMBLY}_autosomal.interval_list
```

### The GATK feature file index

HaplotypeCaller can make use of an index on the "known sites" reference. We create this index using GATK. This process requires that the reference dictionary was created in an early step (see the `CreateSequenceDictionary` step above) because it's the file that's effectively "merged" into the feature VCF.

```bash
# cache the reference copy so we can reuse the name consistently
mv GCF_000001405.39.gz GCF_000001405.39.non-indexed.vcf.gz

# update the sequence dictionary in the vcf
gatk UpdateVcfSequenceDictionary \
    -I GCF_000001405.39.non-indexed.vcf.gz \
    -SD GCF_000001405.40_GRCh38.p14_genomic.dict \
    -O GCF_000001405.39.gz

# then index the feature file (still in vcf format)
gatk IndexFeatureFile -I GCF_000001405.39.gz
```

## Additional reference data


```bash
nextclade dataset get --name='sars-cov-2' --output-dir=$HOME/reference/sars-cov-2/nextclade-data/sars-cov-2
```