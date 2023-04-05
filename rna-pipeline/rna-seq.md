# RNA-seq processing for Ovation IBD data

This document describes how individual FASTQ pairs undergo the "front end" of the RNA-seq pipeline. In particular, the processes of identifying, running initial QC, adapter (and other) trimming, splice- and junction-aware alignment, and finally feature counting.

The following `bash` scripts assume a "base" file name, which maps to the initial Zymo sample or run identifier (typically starting with 'zr'), an assigned sequence identifier for later use in the Ovation Sequence Store, and a sequence run identifier. These are `$SAMPLE`, `$R1`, `$R2`, `$R1_PATH`, `$R2_PATH`, `$SEQUENCEID`, and `$SEQUENCERUNID`. The driver script will, for this exercise, reuse those values. For now the Ovation `BBID` will be unused but can be assigned at the end of the process.

It's also assumed for this process that all indexes are present and the base directory for them is in the `$REFERENCE` environment variable.

Source FASTQ files will be identified using pattern matching from the `ovation-prod-sequence-data` S3 bucket and the resulting output files will be stored in `ovation-prod-rna-sequence-results` using the following structure:

```
ovation-prod-rna-sequence-results
  - zymo
    - $SEQUENCEID
      - $SEQUENCERUNID
        - stats
        - data
        - counts
        - multiqc
```

In the `stats` folder will be the summary output from `featureCounts`, the raw and JSON outputs from `cutadapt`, the HTML and compressed JSON from `fastqc`, and the two log files from `STAR`. In the `data` directory will be the resulting adapter-trimmed friles from `cutadapt` and the aligned BAM file from `STAR`. The `counts` directory will hold the raw `featureCounts` output with all seven columns (Geneid, Chr, Start, End, Strand, Length, and $SEQUENCEID). The last bit of information will be the per-run MultiQC HTML and "data" files.

The above format should let use easily identify subsets of the RNA files for reporting such that we can run deSEQ2 across a collection of files and generate the resulting heatmap.

```bash
#!/usr/bin/env bash

export SEQUENCE_S3=s3://ovation-prod-sequence-data/zymo
export RESULTS_S3=s3://ovation-prod-rna-sequence-results/zymo
```

### Gathering the files

This script can be used to retrieve the files. Each row is a 5-tuple of { sample, r1, r2, r1_path, and r2_path }. This lends itself nicely to running using GNU parallel. Each row will be parsed into the various environment variables. If a given FASTQ pair has not yet been assigned a `$SEQUENCEID` or `$SEQUENCERUNID` that will be done at this time and the values will be written to the `sequence_map.tsv` and `sequence_run_map.tsv` files.

```bash
aws s3 ls --recursive ${SEQUENCE_S3} |
    grep --color=never '/RiboFree/' |
    grep --color=never '.fastq.gz$' |
    sort -t ' ' -k 3n,3 |
    awk -vREMOVE="zymo/" '
        BEGIN {
            OFS = "\t"
            print "sample", "r1", "r2", "s3_path"
        }
        /_R1_/ {
            s3_path = $4;
            nps = split(s3_path, ps, "/");
            r1 = r2 = ps[nps];

            sub(/_R1_/, "_R2_", r2);
            sub(r1, "", s3_path);
            sub(REMOVE, "", s3_path)

            nsc = split(r1, sc, "_");
            sample = sc[1];
            for (n = 2; n <= nsc; n++) {
                if (sc[n] == "R1") {
                    break;
                }
                sample = sample "_" sc[n]
            }
            print sample, r1, r2, s3_path
    }'
```

Each file pair in the above list will be downloaded to the local directory and the following steps will be executed in order. The `$R1` and `$R2` environment variables will be set based on the `$SEQUENCEID` as:

```bash
R1=${SEQUENCEID}_R1.fastq.gz
R2=${SEQUENCEID}_R2.fastq.gz

R1_TRIMMED=data/${SEQUENCEID}_trimmed_R1.fastq.gz
R2_TRIMMED=data/${SEQUENCEID}_trimmed_R2.fastq.gz
```

While it seems like we're going to a lot of trouble to renumber and rename things, this is an opportunity for us to a) normalize the naming conventions, but more importantly to b) add another layer of de-identification to the patient assets.

### Folder setup for the pipeline

```bash
for output in stats data counts multiqc; do
    rm -rf ${output}
    mkdir -p ${output}
do
```

### Initial QC

Our RNA initial QC process is a standard `fastqc` run. We take the original FASTQ pair in `$R1` and `$R2` and run a full `fastqc` pass over them.

```bash
fastqc --threads 2 ${R1} ${R2} -o stats
```

### Adapter trimming

We're using `cutadapt` for RNA-seq adapter trimming and hard-clipping. This matches well with the output results we've received from Zymo and allows us to pass in the specific trimming parameters that we need.

```bash
# the output file names are deliberately named in this weird way
cutadapt \
    --json stats/${SEQUENCERUNID}-cutadapt.json \
    --report full \
    -j 0 \
    -a NNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
    -U 10 \
    -m 20 \
    -o ${R1_TRIMMED} \
    -p ${R2_TRIMMED} \
    ${R1} \
    ${R2} >stats/${SEQUENCERUNID}.cutadapt.txt
```

### Post-trimming QC

Our RNA post-trimming QC process is also a standard `fastqc` run except that we're going to use the "trimmed" files from `cutadapt` to verify that our trimming process isn't being overly aggressive.

```bash
fastqc --threads 2 ${R1_TRIMMED} ${R2_TRIMMED} -o stats
```

### Alignment

```bash
STAR \
    --genomeDir ${REFERENCE} \
    --readFilesIn ${R1_TRIMMED} ${R2_TRIMMED} \
    --runThreadN $(nproc) \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 64324509440 \
    --readFilesCommand zcat \
    --outSAMmultNmax 1 \
    --outFileNamePrefix ${SEQUENCERUNID} \
    --peOverlapNbasesMin 10 \
    --peOverlapMMp 0.01

mv ${SEQUENCERUNID}Aligned.sortedByCoord.out.bam data/${SEQUENCERUNID}.bam

mv ${SEQUENCERUNID}Log.final.out stats
mv ${SEQUENCERUNID}Log.out stats
mv ${SEQUENCERUNID}Log.progress.out stats
```

### Feature counting

```bash
featureCounts \
    -T 24 \
    -F GTF \
    -a ${REFERENCE}/GCF_000001405.40_GRCh38.p14_genomic.gtf \
    -p \
    --countReadPairs \
    --primary \
    -o ${SEQUENCERUNID} \
    -s 2 \
    -t exon \
    -g gene_id \
    data/${SEQUENCERUNID}.bam

mv ${SEQUENCERUNID} counts
mv ${SEQUENCERUNID}.summary stats
```

## Per-sample quality control

- STAR ERCC
- RSeQC
- Preseq
- Markdup / dupRadar
- Qualimap
- gProfiler (uses featureCounts output)

### Overall QC

Finally we'll run MultiQC over the single RNA sequence

```bash
cd stats

multiqc \
    --verbose \
    --force \
    --cl-config "use_filename_as_sample_name: true" \
    --cl-config "fn_clean_sample_names: false" \
    ../multiqc
```

### Cleanup, packaging

When we've completed running all QC we'll remove the local source FASTQ pair and use AWS again to push the files to our output S3 bucket.

```bash
rm -f ${R1}
rm -f ${R2}

aws s3 sync . ${RESULTS_S3}/${SEQUENCEID}/${SEQUENCERUNID}/

echo "sample: ${SAMPLE}" > manifest.txt
echo "sequence: ${SEQUENCEID}" >> manifest.txt
echo "rundate: $(date --utc --iso-8601=seconds)" >> manifest.txt

aws s3 cp manifest.txt ${RESULTS_S3}/${SEQUENCEID}/${SEQUENCERUNID}/
```
