#!/usr/bin/env python

from io import TextIOWrapper
import json
import argparse
import math

from argparse import Namespace
from typing import Dict, Tuple, Sequence, List, Any
from datetime import datetime
from math import ceil
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

ChromosomeSizeList = Dict[str, str]
ChromosomeNames = List[str]
ChromosomeList = List[Tuple[str, str]]
OptionsDict = Dict[str, Any]
FastaPair = Tuple[str, str]


def loadIntervals(chromosomeSizes: str) -> ChromosomeSizeList:
    with open(chromosomeSizes, "r") as file:
        chromosomeSizesDict = json.load(file)
        return chromosomeSizesDict


def loadChromosomeList(chromosomeSizes: str) -> ChromosomeNames:
    with open(chromosomeSizes, "r") as file:
        chromosomeSizesDict = json.load(file)
        return chromosomeSizesDict.keys()


def computeIntervals(options: OptionsDict) -> ChromosomeList:
    intervals = []

    chromosomeSizes = options["chromosomeSizes"]
    segmentSize = options["segmentSize"]
    lastBlockMax = math.floor(options["segmentSize"] * options["factor"])

    with open(chromosomeSizes, "r") as file:
        chromosomeSizes = json.load(file)

        for c in chromosomeSizes:
            remainder = chromosomeSizes[c] - segmentSize
            segments = ceil(chromosomeSizes[c] / segmentSize)
            segment = 0

            while remainder > lastBlockMax:
                lower = segment * segmentSize + 1
                upper = (segment + 1) * segmentSize
                intervals.append("chr{c}:{lower}-{upper}".format(c=c, lower=lower, upper=upper))

                segment += 1
                remainder -= segmentSize

            if remainder > 0:
                lower = (segments - 2) * segmentSize
                if lower % 10 == 0:
                    lower += 1

                intervals.append(
                    "chr{c}:{lower}-{upper}".format(
                        c=c,
                        lower=lower,
                        upper=chromosomeSizes[c],
                    )
                )
            else:
                lower = (segments - 1) * segmentSize
                if lower % 10 == 0:
                    lower += 1

                intervals.append(
                    "chr{c}:{lower}-{upper}".format(
                        c=c,
                        lower=lower,
                        upper=chromosomeSizes[c],
                    )
                )

    return [(interval, interval.replace(":", "_").replace("-", "_")) for interval in intervals]


def getFileNames(options: OptionsDict) -> FastaPair:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_R1.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    bin = options["bin"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    chromosomes = ["chr" + str(c) for c in loadChromosomeList(options["chromosomeSizes"])]
    regex = "|".join(chromosomes)

    script.write(
        """
if [[ ! -f {REFERENCE}/Homo_sapiens_assembly38.dict ]]; then
    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -O {REFERENCE}/Homo_sapiens_assembly38.dict
else
    logthis "Reference dictionary {REFERENCE}/Homo_sapiens_assembly38.dict ${{green}}already completed${{reset}}"
fi

if [[ ! -f {REFERENCE}/hg38_ref_genome_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    logthis "Building {REFERENCE}/hg38_ref_genome_autosomal.interval_list"

    egrep '({REGEX})\\s' {REFERENCE}/Homo_sapiens_assembly38.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/Homo_sapiens_assembly38.dict - >{REFERENCE}/hg38_ref_genome_autosomal.interval_list

    logthis "Building {REFERENCE}/hg38_ref_genome_autosomal.interval_list finished"
else
    logthis "Interval list {REFERENCE}/hg38_ref_genome_autosomal.interval_list ${{green}}already completed${{reset}}"
fi

""".format(
            REFERENCE=reference, BIN=bin, REGEX=regex
        )
    )


def alignAndSort(script: TextIOWrapper, r1: str, r2: str, options: OptionsDict, output: str):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]
    nonRepeatable = options["non-repeatable"]
    readLimit = int(options["read-limit"])
    skipPreprocess = options["skipPreprocess"]
    alignOnly = options["alignOnly"]

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    if skipPreprocess == False:
        script.write(
            """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    logthis "Running fastp and bwa-mem on {SAMPLE}"

    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --in1 {R1} \\
            --in2 {R2} \\
            --verbose {LIMITREADS} \\
            --stdout \\
            --thread 8 \\
            --detect_adapter_for_pe \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html | \\
        bwa-mem2 mem -t {THREADS} \\
            -Y \\
            -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            - | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        logthis "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi

    logthis "fastp and bwa-mem run complete on {SAMPLE}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}already completed${{reset}}"
fi
""".format(
                R1=r1,
                R2=r2,
                REFERENCE=reference,
                SAMPLE=sample,
                THREADS=threads,
                PIPELINE=pipeline,
                TIMEOUT=timeout,
                STATS=stats,
                BIN=bin,
                TEMP=temp,
                DASHK="" if nonRepeatable == True else "-K " + str((10_000_000 * int(threads))),
                LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
            )
        )
    else:
        script.write(
            """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.sam.gz ]]; then
    logthis "Running bwa-mem on {SAMPLE}"

    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            {R1} \\
            {R2} \\
            | pigz >{PIPELINE}/{SAMPLE}.aligned.sam.gz'

    status=$?
    if [ $status -ne 0 ]; then
        logthis "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.sam.gz
        exit $status
    fi

    logthis "bwa-mem run complete on {SAMPLE}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.sam.gz, aligned temp file found, ${{green}}already completed${{reset}}"
fi

""".format(
                R1=r1,
                R2=r2,
                REFERENCE=reference,
                SAMPLE=sample,
                THREADS=threads,
                PIPELINE=pipeline,
                TIMEOUT=timeout,
                BIN=bin,
                DASHK="" if nonRepeatable == True else "-K " + str((10_000_000 * int(threads))),
            )
        )

    # this is one of the more time-consuming operations, so, if we've already
    # done the alignment operation for this sample, we'll skip this
    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} ]]; then
    logthis "Sorting and marking duplicates for {SAMPLE}"

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
        logthis "Watchdog timer killed sort / dup process errno = $status"
        rm -f {SORTED}
        rm -f {SORTED}.bai
        exit $status
    fi

    logthis "Sorting and marking duplicates complete for {SAMPLE}"
else
    logthis "Sorting and duplicate marking done for {SORTED}, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {SORTED}.bai ]]; then
    logthis "Indexing {SORTED}"
    samtools index -@ {THREADS} {SORTED}
else
    logthis "Indexing already completed for {SORTED}, ${{green}}already completed${{reset}}"
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


def splinter(script: TextIOWrapper, bam: str, sorted: str, interval: str):
    script.write(
        """
\n# Scatter interval {INTERVAL}        
if [[ ! -f {BAM} || ! -f {BAM}.bai ]]; then
    logthis "Creating interval {INTERVAL}"
    (samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM} && samtools index -@ 4 {BAM})&
    logthis "Interval {INTERVAL} created"
else
    logthis "Splinter for {INTERVAL} has been computed, ${{green}}already completed${{reset}}"
fi""".format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )


def genBQSR(script: TextIOWrapper, reference: str, interval: str, bam: str, bqsr: str):
    script.write(
        """
    # run base quality score recalibration - build the bqsr table
    if [[ ! -f {BQSR}.table ]]; then
        logthis "Generating {BQSR}.table"

        gatk BaseRecalibrator --java-options '-Xmx8g' \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BAM} \\
            -O {BQSR}.table \\
            --verbosity ERROR \\
            --preserve-qscores-less-than 6 \\
            --known-sites {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            --known-sites {REFERENCE}/Homo_sapiens_assembly38.known_indels.vcf \\
            --known-sites {REFERENCE}/Mills_and_1000G_gold_standard.indels.hg38.vcf \\
            -L {INTERVAL}

        logthis "{BQSR}.table completed"
    else
        logthis "BQSR table generation for {INTERVAL} ${{green}}already completed${{reset}}"
    fi

    # run base quality score recalibration - apply the bqsr table
    if [[ ! -f {BQSR} || ! -f {BQSR}.bai ]]; then
        logthis "Applying calibration for {BQSR}"

        gatk ApplyBQSR --java-options '-Xmx8g' \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BAM} \\
            -O {BQSR} \\
            --verbosity ERROR \\
            --emit-original-quals true \\
            --preserve-qscores-less-than 6 \\
            --static-quantized-quals 10 \\
            --static-quantized-quals 20 \\
            --static-quantized-quals 30 \\
            --bqsr-recal-file {BQSR}.table \\
            -L {INTERVAL}

        # index that file, this is our target for getting vcf
        samtools index -@ 4 {BQSR}

        logthis "Calibation and indexing completed for for {BQSR}"
    else
        logthis "BQSR application for {INTERVAL} ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr
        )
    )


def callVariants(script: TextIOWrapper, reference: str, interval: str, bqsr: str, vcf: str):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        logthis "Starting variant calling for {INTERVAL}"

        gatk HaplotypeCaller --java-options '-Xmx8g' \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BQSR} \\
            -O {VCF} \\
            --verbosity ERROR \\
            --dbsnp {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            --pairHMM FASTEST_AVAILABLE \\
            --native-pair-hmm-threads 4 \\
            -L {INTERVAL}

        logthis "Completed variant calling for {INTERVAL}"
    else
        logthis "Variants already called for {INTERVAL}, ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def callVariants2(script: TextIOWrapper, reference: str, interval: str, bqsr: str, vcf: str):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        logthis "Starting variant calling for {INTERVAL}"

        bcftools mpileup \\
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
            --max-depth 250 \\
            --no-BAQ \\
            --threads 4 \\
            --output-type u \\
            --regions {INTERVAL} \\
            --fasta-ref {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            {BQSR} 2>/dev/null | \\
        bcftools call \\
            --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
            --variants-only \\
            --multiallelic-caller \\
            --ploidy GRCh38 \\
            --threads 4 \\
            --output-type v  \\
            --output {VCF} 2>/dev/null

        logthis "Completed variant calling for {INTERVAL}"
    else
        logthis "Variants already called for {INTERVAL}, ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def annotate(
    script: TextIOWrapper,
    options: OptionsDict,
    vep: str,
    interval: str,
    input: str,
    output: str,
    summary: str,
):
    reference = options["reference"]
    stats = options["stats"]
    chromosome = interval.split(":")[0]

    script.write(
        """
    if [[ ! -f {OUTPUT} || ! -f {STATS}/{SUMMARY} ]]; then
        logthis "Starting annotation for {INTERVAL}"
        
        vep --dir {VEP} \\
            --cache \\
            --format vcf \\
            --vcf \\
            --merged \\
            --fork {FORKS} \\
            --offline \\
            --use_given_ref \\
            --verbose \\
            --force_overwrite {CHROMOSOME} \\
            --exclude_null_alleles \\
            --symbol \\
            --coding_only \\
            --fasta {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            --input_file {INPUT} \\
            --output_file {OUTPUT} \\
            --stats_file {STATS}/{SUMMARY}

        logthis "Completed annotation for {INTERVAL}"
    else
        logthis "Annotations for {INTERVAL} already completed, ${{green}}already completed${{reset}}"
    fi
""".format(
            VEP=vep,
            REFERENCE=reference,
            INPUT=input,
            OUTPUT=output,
            SUMMARY=summary,
            INTERVAL=interval,
            CHROMOSOME="--chr " + chromosome if interval != "FINAL" else "",
            STATS=stats,
            FORKS=8 if interval != "FINAL" else 64,
        )
    )


def scatter(script: TextIOWrapper, options: OptionsDict, prefix: str, sorted: str):
    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        splinter(script, bam, sorted, interval[0])

    script.write(
        """
\n
logthis "${yellow}Waiting for scattering to complete${reset}"
wait
logthis "${green}Scattering completed${reset}"
\n
        """
    )


def runIntervals(script: TextIOWrapper, options: OptionsDict, prefix: str):
    #
    # for now commenting out individual region annotation as there's
    # really no downstream consumer of that data. instead we'll
    # just keep the individual region's calls for later merging and
    # bulk annotation (which provides a complete summary annotation)
    # working = options["working"]
    useAlternateCaller = options["alternateCaller"]

    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        bqsr = """{PREFIX}.{INTERVAL}_bqsr.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        vcf = """{PREFIX}.{INTERVAL}.vcf""".format(PREFIX=prefix, INTERVAL=interval[1])

        script.write("\n")
        script.write("#\n")
        script.write("# Run interval {INTERVAL}\n".format(INTERVAL=interval[0]))
        script.write("#\n")
        script.write("(")

        genBQSR(script, options["reference"], interval[0], bam, bqsr)

        if useAlternateCaller == False:
            callVariants(script, options["reference"], interval[0], bqsr, vcf)
        else:
            callVariants2(script, options["reference"], interval[0], bqsr, vcf)

        # annotate(
        #     script,
        #     options,
        #     "{WORKING}/vep_data".format(WORKING=working),
        #     interval[0],
        #     vcf,
        #     vcf.replace(".vcf", ".annotated.vcf"),
        #     basename(vcf.replace(".vcf", ".annotated.vcf_summary.html")),
        # )

        script.write(") &\n")
        script.write("\n")

    script.write(
        """
logthis "${yellow}Waiting for intervals to complete${reset}"
wait
logthis "${green}Intervals processed${reset}"
        """
    )


def generateConsensus(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
if [[ ! -f {PIPELINE}/{SAMPLE}.consensus.fasta ]]; then
    logthis "${{yellow}}Building consensus {PIPELINE}/{SAMPLE}.unannotated.vcf.gz${{reset}}"
    bcftools consensus \\
        --fasta-ref {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
    | sed '/>/ s/$/ | {SAMPLE}/' >{PIPELINE}/{SAMPLE}.consensus.fasta
else
    logthis "Consensus fasta already generated for {PIPELINE}/{SAMPLE}.consensus.fasta, ${{green}}already completed${{reset}}"
fi
        """.format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample
        )
    )

    pass


def gather(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# Gather interval data and recombine(s)
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "Building merge list for {SAMPLE}"
    /bin/ls -1 {PIPELINE}/{SAMPLE}.chr[0-9MXY]*.vcf | sort -k1,1V >{PIPELINE}/{SAMPLE}.merge.list

    logthis "Concatenating intermediate VCFs into final {PIPELINE}/{SAMPLE}.unannotated.vcf.gz"

    vcf-concat --files {PIPELINE}/{SAMPLE}.merge.list \\
        | bgzip >{PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "VCF concatenation complete for {PIPELINE}/{SAMPLE}.unannotated.vcf.gz"
else
    logthis "VCFs already merged, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz.tbi ]]; then
    logthis "Indexing {PIPELINE}/{SAMPLE}.unannotated.vcf.gz"

    tabix -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "Indexing {PIPELINE}/{SAMPLE}.unannotated.vcf.gz complete"
else
    logthis "VCFs index already created, ${{green}}already completed${{reset}}"
fi

logthis "${{yellow}}Waiting for VCF merge to complete${{reset}}"
wait
logthis "${{green}}Final VCF merge completed${{reset}}"
    """.format(
            PIPELINE=pipeline, SAMPLE=sample
        )
    )


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
logthis "Starting Variant QC processes"

(
    if [[ ! -f {STATS}/{SAMPLE}.variant_calling_detail_metrics ]]; then
        logthis "Collection variant calling metrics"

        gatk CollectVariantCallingMetrics \\
            --VERBOSITY ERROR \\
            --DBSNP {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            -I {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
            -O {STATS}/{SAMPLE}

        logthis "Variant calling metrics collected"
    else
        logthis "Variant metrics already run, ${{green}}already completed${{reset}}"
    fi
) &

#
# we need to quiet vcftools here because it's stupid chatty and doesn't have an option to quiet
#
logthis "Running vcftools statistics"

if [[ ! -f {STATS}/{SAMPLE}.frq ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --freq2 --out {STATS}/{SAMPLE} --max-alleles 2 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.idepth ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --depth --out {STATS}/{SAMPLE} 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.ldepth.mean ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --site-mean-depth --out {STATS}/{SAMPLE} 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.lqual ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --site-quality --out {STATS}/{SAMPLE} 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.imiss ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --missing-indv --out {STATS}/{SAMPLE} 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.lmiss ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --missing-site --out {STATS}/{SAMPLE} 2>/dev/null &
fi

if [[ ! -f {STATS}/{SAMPLE}.het ]]; then
    vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --het --out {STATS}/{SAMPLE} 2>/dev/null &
fi

logthis "${{yellow}}Waiting for variant QC metrics to complete${{reset}}"
job
wait
logthis "${{green}}Variant QC metrics completed${{reset}}"
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )


def runAlignmentQC(script: TextIOWrapper, options: OptionsDict, sorted: str):
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
logthis "Starting QC processes"

if [[ ! -f {STATS}/{SAMPLE}.flagstat.txt ]]; then
    logthis "Starting samtools flagstat on {SAMPLE}"

    samtools flagstat \\
        {SORTED} >{STATS}/{SAMPLE}.flagstat.txt &
else
    logthis "samtools flagstat already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.alignment_metrics.txt ]]; then
    logthis "Starting alignment summary metrics on {SAMPLE}"

    gatk CollectAlignmentSummaryMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    logthis "Alignment metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    logthis "Starting GC bais metrics on {SAMPLE}"

    gatk CollectGcBiasMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    logthis "GC bias metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then
    logthis "Starting WGS metrics on {SAMPLE}"

    gatk CollectWgsMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.wgs_metrics.txt \\
        --READ_LENGTH 151 \\
        -INTERVALS {REFERENCE}/hg38_ref_genome_autosomal.interval_list \\
        --USE_FAST_ALGORITHM \\
        --INCLUDE_BQ_HISTOGRAM &
else
    logthis "WGS metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then
    logthis "Starting samtools stats on {SAMPLE}"

    samtools stats -@ 8 \\
        -r {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        {SORTED} >{STATS}/{SAMPLE}.samstats &
else
    logthis "samtools stats already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samidx ]]; then
    logthis "Starting samtools idxstats on {SAMPLE}"

    samtools idxstats -@ 8 \\
        {SORTED} >{STATS}/{SAMPLE}.samidx &
else
    logthis "samtools idxstats already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.sorted_fastqc.zip || ! -f {STATS}/{SAMPLE}.sorted_fastqc.html ]]; then
    logthis "Starting fastqc for {SAMPLE}"

    fastqc \\
        --outdir {STATS} \\
        --noextract \\
        {SORTED} &
else
    logthis "FASTQC already run, ${{green}}already completed${{reset}}"
fi
""".format(
            REFERENCE=reference,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
            SORTED=sorted,
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
logthis "Starting MultiQC on {SAMPLE}"

## Clean up any old stats
cd {STATS}
rm -rf *multiqc*

# Setup working space
rm -rf {STATS}/qc
mkdir -p {STATS}/qc
cd {STATS}/qc

# Run multiqc
multiqc --tag DNA --verbose -f {STATS}

# Save the output
mv {STATS}/qc/multiqc_data {STATS}/{SAMPLE}_multiqc_data
mv {STATS}/qc/multiqc_report.html {STATS}/{SAMPLE}_multiqc_report.html
rm -rf {STATS}/qc

logthis "MultiQC for {SAMPLE} is complete"
""".format(
            STATS=stats, SAMPLE=sample
        )
    )


def cleanup(script: TextIOWrapper, prefix: str, options: OptionsDict):
    sample = options["sample"]

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all intermediate interval files\n")
    script.write("#\n")

    script.write("rm -f {PREFIX}/{SAMPLE}.chr*\n".format(PREFIX=prefix, SAMPLE=sample))


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastaPair):
    script.write("#\n")
    script.write("# generated at {TIME}\n".format(TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    script.write("#\n")
    script.write("# Parameters\n")
    for opt in options.keys():
        script.write("#   {OPTION} = {VALUE}\n".format(OPTION=opt, VALUE=options[opt]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Input files\n")
    script.write("#   R1 = {F}\n".format(F=filenames[0]))
    script.write("#   R2 = {F}\n".format(F=filenames[1]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Assumed chromosome sizes (from {SIZES})\n".format(SIZES=options["chromosomeSizes"]))
    intervals = loadIntervals(options["chromosomeSizes"])
    for interval in intervals.keys():
        script.write("#   {CHROME} = {SIZE}\n".format(CHROME=interval, SIZE=intervals[interval]))
    script.write("#\n")

    script.write("#\n")
    script.write("# Intervals\n")
    intervals = computeIntervals(options)
    for i, f in intervals:
        script.write("#   {INTERVAL} -> {FILE}\n".format(INTERVAL=i, FILE=f))
    script.write("#\n")

    script.write("#\n")
    script.write("# Split parameters\n")
    script.write("#   segmentSize = {P}\n".format(P=options["segmentSize"]))
    script.write("#   factor = {P}\n".format(P=options["factor"]))
    script.write("#   last block max size = {P}\n".format(P=math.floor(options["segmentSize"] * options["factor"])))


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

function logthis() {{
  NOW=$(date "+%Y-%m-%d %H:%M:%S")
  echo -e "[${{NOW}}] ${{1}}"
}}

# deal with really large fastq files
ulimit -n 8192

# perl stuff
export PATH=/home/ubuntu/perl5/bin:$PATH
export PERL5LIB=/home/ubuntu/perl5/lib/perl5:$PERL5LIB
export PERL_LOCAL_LIB_ROOT=/home/ubuntu/perl5:$PERL_LOCAL_LIB_ROOT

# shared library stuff
export LD_LIBRARY_PATH={WORKING}/bin:/usr/lib64:/usr/local/lib/:$LB_LIBRARY_PATH
export LD_PRELOAD={WORKING}/bin/libz.so.1.2.11.zlib-ng

# handy path
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.3.0:{WORKING}/bin:$PATH\n""".format(
            WORKING=options["working"]
        )
    )
    script.write("\n")


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>_R[12].fastq.gz",
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
        "-F",
        "--fastq-dir",
        action="store",
        metavar="FASTQ_DIR",
        dest="fastq_dir",
        help="Location of R1 and R2 files",
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
        help="Location of R1 and R2 files",
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

    parser.add_argument(
        "--sizes",
        action="store",
        metavar="CHROME_SIZES",
        dest="chromosomeSizes",
        default="chromosomeSizes.json",
        help="Name of JSON file containing chromosome sizes (ADVANCED)",
    )

    parser.add_argument(
        "--segment",
        action="store",
        metavar="SEGMENT_SIZE",
        dest="segmentSize",
        default=50_000_000,
        help="Size of interval partition (ADVANCED)",
    )
    parser.add_argument(
        "--factor",
        action="store",
        metavar="FACTOR",
        dest="factor",
        default=0.25,
        help="Interval remainder buffer (between 0.10 and 0.50) (ADVANCED)",
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
        "chromosomeSizes",
    ]:
        options[opt] = expandvars(options[opt])

    return options


def verifyOptions(options: OptionsDict):
    filenames = getFileNames(options)

    if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
        print(
            "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                R1=filenames[0],
                R2=filenames[1],
            )
        )
        print("Check your --sample and --work-dir parameters")
        quit(1)

    if exists(options["bin"]) == False:
        print("Unable to find your --bin-dir directory at {PATH}".format(PATH=options["bin"]))
        quit(1)

    if exists(options["working"]) == False:
        print("Unable to find your --work-dir directory at {PATH}".format(PATH=options["working"]))
        quit(1)

    if exists(options["temp"]) == False:
        print("Unable to find your --temp-dir directory at {PATH}".format(PATH=options["temp"]))
        quit(1)

    if exists(options["reference"]) == False:
        print("Unable to find your --reference-dir directory at {PATH}".format(PATH=options["reference"]))
        quit(1)

    if exists(options["pipeline"]) == False:
        print("Unable to find your --pipeline-dir directory at {PATH}".format(PATH=options["pipeline"]))
        quit(1)

    if exists(options["stats"]) == False:
        print("Unable to find your --stats-dir directory at {PATH}".format(PATH=options["stats"]))
        quit(1)


def main():
    opts = defineArguments()
    options = fixupPathOptions(opts)

    verifyOptions(options)

    filenames = getFileNames(options)
    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])
    prefix = "{PIPELINE}/{SAMPLE}".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])
    cleantarget = options["pipeline"]

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write("touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

        updateDictionary(script, options)
        filenames = getFileNames(options)

        alignAndSort(
            script,
            filenames[0],
            filenames[1],
            options,
            sorted,
        )

        # the order of the next two operations is important in this script
        # because scatter does its job but waits on all the background
        # operations to complete. this doesn't take long, so waiting is
        # ok.
        #
        # however, if we were to start the alignment qc process, which runs
        # in the background, we would start them and then block until they
        # completed while we're waiting to scatter intervals around
        #
        # this makes scatter a standalone process for now
        scatter(
            script,
            options,
            prefix,
            sorted,
        )

        # if we've been asked to run qc at all we can actually start
        # the alignment qc processes very early (right after we've
        # completed alignment and sorting). so, we start them here. there's
        # an overall impact on the interval runs while these qc processes
        # are running, but we're will to deal with that since the overall
        # qc process takes a long time
        #
        # in particular it's the fastqc process itself which chews up
        # so much time. it's a single-threaded limited-memory (250mb)
        # process so we're not too worried about it's impact on the rest of
        # the variant calling processes
        if options["doQC"]:
            runAlignmentQC(script, options, sorted)

        runIntervals(script, options, prefix)
        gather(script, options)

        annotate(
            script,
            options,
            "{WORKING}/vep_data".format(WORKING=options["working"]),
            "FINAL",
            "{PIPELINE}/{SAMPLE}.unannotated.vcf.gz".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"]),
            "{PIPELINE}/{SAMPLE}.annotated.vcf.gz".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"]),
            "{SAMPLE}.annotated.vcf_summary.html".format(SAMPLE=options["sample"]),
        )

        # generateConsensus(script, options)

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, cleantarget, options)

        if options["doQC"]:
            doVariantQC(script, options)
            runMultiQC(script, options)

        script.write(
            """
logthis "${{yellow}}Waiting for any outstanding processes to complete, this might return immediately, it might not.${{reset}}"
wait
logthis "${{green}}Done processing${{reset}} {SAMPLE}"
""".format(
                SAMPLE=options["sample"],
            )
        )

        script.write("\n")
        script.write("touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

    system("chmod +x {SCRIPT}".format(SCRIPT=options["script"]))


if __name__ == "__main__":
    main()
