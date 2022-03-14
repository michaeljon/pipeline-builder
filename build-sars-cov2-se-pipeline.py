#!/usr/bin/env python

from io import TextIOWrapper
import argparse
import math

from argparse import Namespace
from typing import Dict, Tuple, Any
from datetime import datetime
from os.path import exists, expandvars, basename
from os import cpu_count, system

# This script generates an Illumina paired-read VCF annotation pipeline.
# There are some short-circuit options generated as part of this script,
# but they will not normally be used in production. They are meant to
# prevent a developer from re-running various steps in the pipeline in
# the case that the output files are already present.

# a good example of how to run this might be:
# python ./build-pipeline.py \
#   --sample S119813_no_sonication \
#   --work-dir $WORKING \
#   --temp-dir $WORKING/stats/temp \
#   --script ../pipeline-runner \
#   --watchdog 240 \
#   --skip-fastp \
#   --align-only
#
# which sets the working and temp dir (important for large sequences)
# specifies the target script
# sets the watchdog timer
# skips fastp processing
# only runs the alignment steps

OptionsDict = Dict[str, Any]
FastqSet = Tuple[str, str, str, str]


def getFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_L001_R1_001.fastq.gz".format(
            FASTQ_DIR=fastq_dir, SAMPLE=sample
        ),
        "{FASTQ_DIR}/{SAMPLE}_L002_R1_001.fastq.gz".format(
            FASTQ_DIR=fastq_dir, SAMPLE=sample
        ),
        "{FASTQ_DIR}/{SAMPLE}_L003_R1_001.fastq.gz".format(
            FASTQ_DIR=fastq_dir, SAMPLE=sample
        ),
        "{FASTQ_DIR}/{SAMPLE}_L004_R1_001.fastq.gz".format(
            FASTQ_DIR=fastq_dir, SAMPLE=sample
        ),
    )


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    bin = options["bin"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    script.write(
        """
if [[ ! -f {REFERENCE}/covid_reference.dict ]]; then
    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/covid_reference.fasta \\
        -O {REFERENCE}/covid_reference.dict
else
    echo "Reference dictionary {REFERENCE}/covid_reference.dict ${{green}}already present${{reset}}"
fi

if [[ ! -f {REFERENCE}/sars_cov2_ref_genome_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    egrep '(NC_045512.2)\\s' {REFERENCE}/covid_reference.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/covid_reference.dict - >{REFERENCE}/sars_cov2_ref_genome_autosomal.interval_list
else
    echo "Interval list {REFERENCE}/sars_cov2_ref_genome_autosomal.interval_list ${{green}}already present${{reset}}"
fi

""".format(
            REFERENCE=reference, BIN=bin
        )
    )


def combineLaneData(script: TextIOWrapper, files: FastqSet, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz ]]; then
    echo "Combining lane data files"
    
    zcat {L1} \\
         {L2} \\
         {L3} \\
         {L4} | pigz >{PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz
else
    echo "{PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz found, ${{green}}not combining${{reset}}"
fi
""".format(
            SAMPLE=sample,
            PIPELINE=pipeline,
            L1=files[0],
            L2=files[1],
            L3=files[2],
            L4=files[3],
        )
    )


def alignWithFastpAndBWA(
    script: TextIOWrapper,
    options: OptionsDict,
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    nonRepeatable = options["non-repeatable"]
    readLimit = int(options["read-limit"])

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --report_title "fastp report for sample {SAMPLE}" \\
            --in1 {PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz \\
            --verbose {LIMITREADS} \\
            --adapter_sequence CTGTCTCTTATACACATCT \\
            --stdout \\
            --thread 8 \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html | \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/covid_reference.fasta \\
            - | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            DASHK=""
            if nonRepeatable == True
            else "-K " + str((10_000_000 * int(threads))),
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def alignWithFastpAndMinimap(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    readLimit = int(options["read-limit"])

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --report_title "fastp report for sample {SAMPLE}" \\
            --in1 {PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz \\
            --verbose {LIMITREADS} \\
            --adapter_sequence CTGTCTCTTATACACATCT \\
            --stdout \\
            --thread 8 \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html | \\
        minimap2 -a -t {THREADS} {REFERENCE}/covid_reference.mmi  \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            - | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def alignWithoutFastpUsingBWA(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    nonRepeatable = options["non-repeatable"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/covid_reference.fasta \\
            {PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz \\
            | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            DASHK=""
            if nonRepeatable == True
            else "-K " + str((10_000_000 * int(threads))),
            BIN=bin,
        )
    )


def alignWithoutFastpUsingMinimap(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
         minimap2 -a -t {THREADS} {REFERENCE}/covid_reference.mmi  \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {PIPELINE}/{SAMPLE}.combined_lanes.fastq.gz | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            BIN=bin,
        )
    )


def sortAlignedAndMappedData(script: TextIOWrapper, options: OptionsDict, output: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    alignOnly = options["alignOnly"]

    # this is one of the more time-consuming operations, so, if we've already
    # done the alignment operation for this sample, we'll skip this
    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        unpigz --stdout {PIPELINE}/{SAMPLE}.aligned.sam.gz |
        bamsormadup \\
            SO=coordinate \\
            threads={THREADS} \\
            level=6 \\
            tmpfile={TEMP}/{SAMPLE} \\
            inputformat=sam \\
            indexfilename={SORTED}.bai \\
            M={STATS}/{SAMPLE}.duplication_metrics >{SORTED}'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed sort / dup process errno = $status"
        rm -f {SORTED}
        rm -f {SORTED}.bai
        exit $status
    fi

    # force the index to look "newer" than its source
    touch {SORTED}.bai
else
    echo "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi
{EXIT_IF_ALIGN_ONLY}
""".format(
            SAMPLE=sample,
            THREADS=threads,
            SORTED=output,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
            EXIT_IF_ALIGN_ONLY="exit" if alignOnly == True else "",
        )
    )


def alignAndSort(
    script: TextIOWrapper, files: FastqSet, options: OptionsDict, output: str
):
    skipFastp = options["skipPreprocess"]
    useAlternateAligner = options["alternateAligner"]

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    combineLaneData(script, files, options)

    if skipFastp == False:
        if useAlternateAligner == False:
            alignWithFastpAndBWA(script, options)
        else:
            alignWithFastpAndMinimap(script, options)
    else:
        if useAlternateAligner == False:
            alignWithoutFastpUsingBWA(script, options)
        else:
            alignWithoutFastpUsingMinimap(script, options)

    sortAlignedAndMappedData(script, options, output)


def callVariants(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    consensus: str,
    bam: str,
    vcf: str,
):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        gatk HaplotypeCaller --java-options '-Xmx8g' \\
            -R {REFERENCE}/covid_reference.fasta \\
            -I {BAM} \\
            -O {VCF} \\
            --verbosity ERROR \\
            --pairHMM FASTEST_AVAILABLE \\
            --native-pair-hmm-threads 4

        bcftools view --output-type z <{VCF} >{VCF}.gz
        bcftools index {VCF}.gz
        bcftools consensus \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
            {VCF}.gz \\
        | sed '/>/ s/$/ | {SAMPLE}/' >{CONSENSUS}
    else
        echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
    fi
""".format(
            REFERENCE=reference, BAM=bam, VCF=vcf, CONSENSUS=consensus, SAMPLE=sample
        )
    )


def callVariants2(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    consensus: str,
    bam: str,
    vcf: str,
    ploidy: str,
):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        echo Starting variant calling for {BAM}

        if [[ ! -f {PLOIDY} ]]; then
            echo '* * * * 1' >{PLOIDY}
        fi

        bcftools mpileup \\
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
            --max-idepth 1000000 \\
            --max-depth 1000000 \\
            --count-orphans \\
            --full-BAQ \\
            --redo-BAQ \\
            --threads 4 \\
            --output-type u \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
            {BAM} 2>/dev/null | \\
        bcftools call \\
            --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
            --variants-only \\
            --keep-alts \\
            --multiallelic-caller \\
            --ploidy-file {PLOIDY} \\
            --threads 4 \\
            --output-type v  \\
            --output {VCF} 2>/dev/null

        bcftools view --output-type z <{VCF} >{VCF}.gz
        bcftools index {VCF}.gz
        bcftools consensus \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
            {VCF}.gz \\
        | sed '/>/ s/$/ | {SAMPLE}/' >{CONSENSUS}

        echo Completed variant calling for {BAM}
    else
        echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
    fi
""".format(
            REFERENCE=reference,
            BAM=bam,
            VCF=vcf,
            PLOIDY=ploidy,
            CONSENSUS=consensus,
            SAMPLE=sample,
        )
    )


def producePileup(
    script: TextIOWrapper,
    reference: str,
    bam: str,
    pileup: str,
    ploidy: str,
):
    script.write(
        """
    # create pileup
    if [[ ! -f {PILEUP} ]]; then
        echo Generate pileup for {BAM}

        if [[ ! -f {PLOIDY} ]]; then
            echo '* * * * 1' >{PLOIDY}
        fi

        samtools mpileup \\
            --max-depth 1000000 \\
            --count-orphans \\
            --redo-BAQ \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
            {BAM} | gzip >{PILEUP} 2>/dev/null

        echo Pileup completed for {BAM}
    else
        echo "Pileup already finished for {BAM}, ${{green}}skipping${{reset}}"
    fi
""".format(
            REFERENCE=reference,
            BAM=bam,
            PILEUP=pileup,
            PLOIDY=ploidy,
        )
    )


def annotate(
    script: TextIOWrapper,
    options: OptionsDict,
    vcf: str,
    annotated: str
):
    sample = options["sample"]
    pipeline = options["pipeline"]
    bin = options["bin"]

    script.write(
        """
    # annotate
    if [[ ! -f {PIPELINE}/{SAMPLE}.nirvana.json.gz || ! -f {PIPELINE}/{SAMPLE}.nirvana.json.gz.jsi ]]; then
        echo Start nirvana annotation for {VCF}

        dotnet {BIN}/nirvana/Nirvana.dll \\
            -c {BIN}/nirvana/Data/Cache/SARS-CoV-2/SARS-CoV-2 \\
            -sd {BIN}/nirvana/Data/SupplementaryAnnotation/SARS-CoV-2 \\
            --enable-dq \\
            -r {BIN}/nirvana/Data/References/SARS-CoV-2.ASM985889v3.dat \\
            -i {VCF} \\
            -o {SAMPLE}.nirvana

        mv {SAMPLE}.nirvana.json.gz {PIPELINE}
        mv {SAMPLE}.nirvana.json.gz.jsi {PIPELINE}

        echo Annotion complete for {VCF}
    else
        echo "Annotations already for {VCF}, ${{green}}skipping${{reset}}"
    fi

    if [[ ! -f {PIPELINE}/{SAMPLE}.annotated.vcf ]]; then
        echo Start snpEff annotation for {VCF}

        java -jar ~/bin/snpEff/snpEff.jar NC_045512.2 {VCF} >{ANNOTATED}

        echo snpEff annotion complete for {VCF}
    else
        echo "snpEff annotations already for {VCF}, ${{green}}skipping${{reset}}"
    fi
""".format(
            VCF=vcf, SAMPLE=sample, BIN=bin, PIPELINE=pipeline, ANNOTATED=annotated
        )
    )


def runPipeline(script: TextIOWrapper, options: OptionsDict, prefix: str):
    useAlternateCaller = options["alternateCaller"]
    sample = options["sample"]

    bam = """{PREFIX}.sorted.bam""".format(PREFIX=prefix)
    vcf = """{PREFIX}.vcf""".format(PREFIX=prefix)
    annotated = """{PREFIX}.annotated.vcf""".format(PREFIX=prefix)
    pileup = """{PREFIX}.pileup.gz""".format(PREFIX=prefix)
    consensus = """{PREFIX}.consensus.fa""".format(PREFIX=prefix)
    ploidy = """{PREFIX}.ploidy""".format(PREFIX=prefix)

    script.write("\n")
    script.write("#\n")
    script.write("# Running pipeline for {BAM}\n".format(BAM=bam))
    script.write("#\n")
    script.write("(")

    if useAlternateCaller == False:
        callVariants(script, options["reference"], sample, consensus, bam, vcf)
    else:
        callVariants2(script, options["reference"], sample, consensus, bam, vcf, pileup)

    producePileup(script, options["reference"], bam, pileup, ploidy)
    annotated = """{PREFIX}.annotated.vcf""".format(PREFIX=prefix)
    annotate(script, options, vcf, annotated)

    script.write(") &\n")
    script.write("\n")


def doVariantQC(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    script.write(
        """
#
# RUN Variant QC process
# 
echo "Starting Variant QC processes"

# if [[ ! -f {STATS}/{SAMPLE}.variant_calling_detail_metrics ]]; then
#    gatk CollectVariantCallingMetrics \\
#        --VERBOSITY ERROR \\
#        -I {PIPELINE}/{SAMPLE}.vcf \\
#        -O {STATS}/{SAMPLE} &
# else
#     echo "Variant metrics already run, ${{green}}skipping${{reset}}"
# fi

#
# we need to quiet vcftools here because it's stupid chatty and doesn't have an option to quiet
#
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --freq2 --out {STATS}/{SAMPLE} --max-alleles 2 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --depth --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --site-mean-depth --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --site-quality --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --missing-indv --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --missing-site --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --het --out {STATS}/{SAMPLE} 2>/dev/null &
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )


def runAlignmentQC(
    script: TextIOWrapper, options: OptionsDict, sorted: str, aligned: str
):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]

    script.write(
        """
#
# RUN QC process
# 
echo "Starting QC processes"

if [[ ! -f {STATS}/{SAMPLE}.flagstat.txt ]]; then
    samtools flagstat \\
        {SORTED} >{STATS}/{SAMPLE}.flagstat.txt &
else
    echo "samtools flagstat already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.alignment_metrics.txt ]]; then
    gatk CollectAlignmentSummaryMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        --ADAPTER_SEQUENCE CTGTCTCTTATACACATCT \\
        -R {REFERENCE}/covid_reference.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    echo "Alignment metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    gatk CollectGcBiasMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/covid_reference.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    echo "GC bias metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then
    samtools stats -@ 8 \\
        -r reference/covid_reference.fasta \\
        {SORTED} >{STATS}/{SAMPLE}.samstats &
else
    echo "samtools stats already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samidx ]]; then
    samtools idxstats -@ 8 \\
        {SORTED} >{STATS}/{SAMPLE}.samidx &
else
    echo "samtools idxstats already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.coverage ]]; then
    bedtools genomecov -d -ibam \\
        {SORTED} >{STATS}/{SAMPLE}.coverage &
else
    echo "bedtools genomecov already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.sorted_fastqc.zip || ! -f {STATS}/{SAMPLE}.sorted_fastqc.html ]]; then
    fastqc \\
        --outdir {STATS} \\
        --noextract \\
        {SORTED} &
else
    echo "FASTQC already run, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
            SORTED=sorted,
            ALIGNED=aligned,
        )
    )


def runMultiQC(script: TextIOWrapper, options: OptionsDict):
    stats = options["stats"]
    sample = options["sample"]

    script.write(
        """
#
# Run MultiQC across everything
# 

## Clean up any old stats
cd {STATS}
rm -rf *multiqc*

# Setup working space
rm -rf {STATS}/qc
mkdir -p {STATS}/qc
cd {STATS}/qc

# Run multiqc
multiqc --tag RNA -f {STATS}

# Save the output
mv {STATS}/qc/multiqc_data {STATS}/{SAMPLE}_multiqc_data
mv {STATS}/qc/multiqc_report.html {STATS}/{SAMPLE}_multiqc_report.html
rm -rf {STATS}/qc

""".format(
            STATS=stats, SAMPLE=sample
        )
    )


def cleanup(script: TextIOWrapper, prefix: str, options: OptionsDict):
    pipeline = options["pipeline"]

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all intermediate files\n")
    script.write("#\n")

    script.write("# TODO: Delete all the thingz\n")

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all remaining intermediate files\n")
    script.write("#\n")
    script.write("\n")
    for type in ["snps", "indels"]:
        script.write(
            "rm -f {PIPELINE}/merge.{TYPE}.list\n".format(PIPELINE=pipeline, TYPE=type)
        )

    script.write("\n")


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastqSet):
    script.write("#\n")
    script.write(
        "# generated at {TIME}\n".format(
            TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
    )
    script.write("#\n")
    script.write("# Parameters\n")
    for opt in options.keys():
        script.write("#   {OPTION} = {VALUE}\n".format(OPTION=opt, VALUE=options[opt]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Input files\n")
    script.write("#   L001 = {F}\n".format(F=filenames[0]))
    script.write("#   L002 = {F}\n".format(F=filenames[1]))
    script.write("#   L003 = {F}\n".format(F=filenames[2]))
    script.write("#   L004 = {F}\n".format(F=filenames[3]))
    script.write("#\n")


def writeVersions(script: TextIOWrapper):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...\n")
    script.write("#\n")


def writeEnvironment(script: TextIOWrapper, options: OptionsDict):
    script.write("#\n")
    script.write(
        """
# color shortcuts
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
reset=$(tput sgr0)

# deal with really large fastq files
ulimit -n 8192

# perl stuff
export PATH=/home/ubuntu/perl5/bin:$PATH
export PERL5LIB=/home/ubuntu/perl5/lib/perl5:$PERL5LIB
export PERL_LOCAL_LIB_ROOT=/home/ubuntu/perl5:$PERL_LOCAL_LIB_ROOT

# shared library stuff
export LD_LIBRARY_PATH={WORKING}/bin:/usr/lib64:/usr/local/lib/:$LB_LIBRARY_PATH
export LD_PRELOAD={WORKING}/bin/libz.so.1.2.11.zlib-ng

# dotnet stuff
export PATH=/home/ubuntu/.dotnet


# handy path
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.3.0:{WORKING}/bin:$PATH\n""".format(
            WORKING=options["working"]
        )
    )
    script.write("\n")


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.set_defaults(doQC=False, cleanIntermediateFiles=True)
    parser.add_argument(
        "-s",
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>__L00[1-4]_R1_001.fastq.gz",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        required=True,
        action="store",
        metavar="WORKING_DIR",
        dest="working",
        help="Working directory, e.g. base for $WORKING/pipeline, $WORKING/stats",
    )
    parser.add_argument(
        "-q",
        "--skip-qc",
        action="store_false",
        dest="doQC",
        default=True,
        help="Skip QC process on input and output files",
    )
    parser.add_argument(
        "-P",
        "--skip-fastp",
        action="store_true",
        dest="skipPreprocess",
        default=False,
        help="Skip running fastp on input file(s)",
    )
    parser.add_argument(
        "-c",
        "--clean",
        action="store_true",
        dest="cleanIntermediateFiles",
        default=False,
        help="Clean up the mess we make",
    )

    parser.add_argument(
        "-a",
        "--align-only",
        action="store_true",
        dest="alignOnly",
        default=False,
        help="Only run alignment and sorting processes",
    )

    parser.add_argument(
        "-A",
        "--alternate-caller",
        action="store_true",
        dest="alternateCaller",
        default=False,
        help="Use alternate caller bcftools",
    )

    parser.add_argument(
        "-C",
        "--alternate-aligner",
        action="store_true",
        dest="alternateAligner",
        default=False,
        help="Use minimap2 as the alignment tool",
    )

    parser.add_argument(
        "-r",
        "--reference-dir",
        action="store",
        metavar="REFERENCE_DIR",
        dest="reference",
        help="Location of the reference genome files",
    )
    parser.add_argument(
        "-p",
        "--pipeline-dir",
        action="store",
        metavar="PIPELINE_DIR",
        dest="pipeline",
        help="Output directory for processing",
    )
    parser.add_argument(
        "-F",
        "--fastq-dir",
        action="store",
        metavar="FASTQ_DIR",
        dest="fastq_dir",
        help="Location of L00[1234] files",
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        action="store",
        metavar="TEMP_DIR",
        dest="temp",
        help="Temporary storage for alignment",
    )
    parser.add_argument(
        "-o",
        "--stats-dir",
        action="store",
        metavar="STATS_DIR",
        dest="stats",
        help="Destination for statistics and QC files",
    )
    parser.add_argument(
        "-b",
        "--bin-dir",
        action="store",
        metavar="BIN_DIR",
        dest="bin",
        default="$HOME/bin",
        help="Install location of all tooling",
    )

    parser.add_argument(
        "--script",
        action="store",
        metavar="SHELL_SCRIPT",
        dest="script",
        default="pipeline-runner",
        help="Filename of bash shell to create",
    )

    parser.add_argument(
        "-z",
        "--cores",
        action="store",
        dest="cores",
        metavar="CPU_COUNT",
        default=cpu_count(),
        help="Specify the number of available CPU",
    )

    parser.add_argument(
        "-N",
        "--non-repeatable",
        action="store_true",
        dest="non-repeatable",
        default=False,
        help="Force the -K option to BWA, will use 10,000,000bp per thread",
    )

    parser.add_argument(
        "--read-limit",
        action="store",
        dest="read-limit",
        metavar="READ_LIMIT",
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    parser.add_argument(
        "-d",
        "--watchdog",
        action="store",
        metavar="WATCHDOG_TIMEOUT",
        dest="watchdog",
        default=150,
        help="Specify a watchdog timeout for the alignment process. Value is in minutes",
    )

    return parser.parse_args()


def fixupPathOptions(opts: Namespace) -> OptionsDict:
    options = vars(opts)

    if options["working"] == None:
        options["working"] = "$HOME"
    if options["reference"] == None:
        options["reference"] = "{WORKING}/reference".format(WORKING=options["working"])
    if options["pipeline"] == None:
        options["pipeline"] = "{WORKING}/pipeline".format(WORKING=options["working"])
    if options["fastq_dir"] == None:
        options["fastq_dir"] = "{WORKING}/pipeline".format(WORKING=options["working"])
    if options["stats"] == None:
        options["stats"] = "{WORKING}/stats".format(WORKING=options["working"])
    if options["temp"] == None:
        options["temp"] = "{WORKING}/pipeline".format(WORKING=options["working"])

    if options["sample"] == None:
        print("--sample is a required option")
        quit(1)

    for opt in [
        "working",
        "reference",
        "pipeline",
        "stats",
        "temp",
        "bin",
        "script",
    ]:
        options[opt] = expandvars(options[opt])

    return options


def verifyOptions(options: OptionsDict):
    filenames = getFileNames(options)

    if (
        exists(expandvars(filenames[0])) == False
        or exists(expandvars(filenames[1])) == False
        or exists(expandvars(filenames[2])) == False
        or exists(expandvars(filenames[3])) == False
    ):
        print(
            "Unable to locate the FASTQ files at {L1}, {L2}, {L3}, or {L4}".format(
                L1=filenames[0],
                L2=filenames[1],
                L3=filenames[2],
                L4=filenames[3],
            )
        )
        print("Check your --sample and --work-dir parameters")
        quit(1)

    if exists(options["bin"]) == False:
        print(
            "Unable to find your --bin-dir directory at {PATH}".format(
                PATH=options["bin"]
            )
        )
        quit(1)

    if exists(options["working"]) == False:
        print(
            "Unable to find your --work-dir directory at {PATH}".format(
                PATH=options["working"]
            )
        )
        quit(1)

    if exists(options["temp"]) == False:
        print(
            "Unable to find your --temp-dir directory at {PATH}".format(
                PATH=options["temp"]
            )
        )
        quit(1)

    if exists(options["reference"]) == False:
        print(
            "Unable to find your --reference-dir directory at {PATH}".format(
                PATH=options["reference"]
            )
        )
        quit(1)

    if exists(options["pipeline"]) == False:
        print(
            "Unable to find your --pipeline-dir directory at {PATH}".format(
                PATH=options["pipeline"]
            )
        )
        quit(1)

    if exists(options["stats"]) == False:
        print(
            "Unable to find your --stats-dir directory at {PATH}".format(
                PATH=options["stats"]
            )
        )
        quit(1)


def main():
    opts = defineArguments()
    options = fixupPathOptions(opts)
    verifyOptions(options)

    filenames = getFileNames(options)
    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )
    prefix = "{PIPELINE}/{SAMPLE}".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )
    aligned = "{PIPELINE}/{SAMPLE}.aligned.sam.gz".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"]
    )

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write(
            "touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"])
        )
        script.write("\n")

        updateDictionary(script, options)
        filenames = getFileNames(options)

        alignAndSort(
            script,
            filenames,
            options,
            sorted,
        )

        if options["doQC"]:
            runAlignmentQC(script, options, sorted, aligned)

        runPipeline(script, options, prefix)

        if options["doQC"]:
            script.write(
                """
echo ${yellow}Waiting for variant calling pipeline to complete before starting variant qc and multiqc${reset}
wait
echo ${green}Pipeline processed${reset}
                """
            )

            doVariantQC(script, options)
            runMultiQC(script, options)

        script.write(
            """
\necho ${{yellow}}Waiting for any outstanding processes to complete, this might return immediately, it might not.${{reset}}
wait
\necho -e "${{green}}Done processing${{reset}} {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"\n""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                PIPELINE=options["pipeline"],
            )
        )

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, prefix, options)

        script.write("\n")
        script.write(
            "touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"])
        )
        script.write("\n")

    system("chmod +x " + options["script"])


if __name__ == "__main__":
    main()
