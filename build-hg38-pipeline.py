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
from os import cpu_count

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
                intervals.append(
                    "chr{c}:{lower}-{upper}".format(c=c, lower=lower, upper=upper)
                )

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

    return [
        (interval, interval.replace(":", "_").replace("-", "_"))
        for interval in intervals
    ]


def getFileNames(options: OptionsDict) -> FastaPair:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    )


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    bin = options["bin"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    chromosomes = [
        "chr" + str(c) for c in loadChromosomeList(options["chromosomeSizes"])
    ]
    regex = "|".join(chromosomes)

    script.write(
        """
if [[ ! -f {REFERENCE}/Homo_sapiens_assembly38.dict ]]; then
    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -O {REFERENCE}/Homo_sapiens_assembly38.dict
else
    echo "Reference dictionary {REFERENCE}/Homo_sapiens_assembly38.dict ${{green}}already present${{reset}}"
fi

if [[ ! -f {REFERENCE}/ref_genome_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    egrep '({REGEX})\\s' {REFERENCE}/Homo_sapiens_assembly38.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/Homo_sapiens_assembly38.dict - >{REFERENCE}/ref_genome_autosomal.interval_list
else
    echo "Interval list {REFERENCE}/ref_genome_autosomal.interval_list ${{green}}already present${{reset}}"
fi

""".format(
            REFERENCE=reference, BIN=bin, REGEX=regex
        )
    )


def alignAndSort(
    script: TextIOWrapper, r1: str, r2: str, options: OptionsDict, output: str
):
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
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            -i {R1} \\
            -I {R2} \\
            --verbose {LIMITREADS} \\
            --stdout \\
            --thread 8 \\
            --detect_adapter_for_pe \\
            -j {STATS}/{SAMPLE}-fastp.json \\
            -h {STATS}/{SAMPLE}-fastp.html | \\
        bwa-mem2 mem -t {THREADS} \\
            -Y -M {DASHK} \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/Homo_sapiens_assembly38.fasta \\
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
                DASHK=""
                if nonRepeatable == True
                else "-K " + str((10_000_000 * threads)),
                LIMITREADS="--reads_to_process " + str(readLimit)
                if readLimit > 0
                else "",
            )
        )
    else:
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
            {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            {R1} \\
            {R2} \\
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
                R1=r1,
                R2=r2,
                REFERENCE=reference,
                SAMPLE=sample,
                THREADS=threads,
                PIPELINE=pipeline,
                TIMEOUT=timeout,
                BIN=bin,
                DASHK=""
                if nonRepeatable == True
                else "-K " + str((10_000_000 * threads)),
            )
        )

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


def splinter(script: TextIOWrapper, bam: str, sorted: str, interval: str):
    script.write(
        """
\n# Scatter interval {INTERVAL}        
if [[ ! -f {BAM} || ! -f {BAM}.bai ]]; then
    echo "Creating interval {INTERVAL}"
    (samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM} && samtools index -@ 4 {BAM})&
else
    echo "Splinter for {INTERVAL} has been computed, ${{green}}skipping${{reset}}"
fi""".format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )


def genBQSR(script: TextIOWrapper, reference: str, interval: str, bam: str, bqsr: str):
    script.write(
        """
    # run base quality score recalibration - build the bqsr table
    if [[ ! -f {BQSR}.table ]]; then
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
    else
        echo "BQSR table generation for {INTERVAL} ${{green}}skipped${{reset}}"
    fi

    # run base quality score recalibration - apply the bqsr table
    if [[ ! -f {BQSR} || ! -f {BQSR}.bai ]]; then
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
    else
        echo "BQSR application for {INTERVAL} ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr
        )
    )


def callVariants(
    script: TextIOWrapper, reference: str, interval: str, bqsr: str, vcf: str
):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        gatk HaplotypeCaller --java-options '-Xmx8g' \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BQSR} \\
            -O {VCF} \\
            --verbosity ERROR \\
            --dbsnp {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            --pairHMM FASTEST_AVAILABLE \\
            --native-pair-hmm-threads 4 \\
            -L {INTERVAL}
    else
        echo "Variants already called for {INTERVAL}, ${{green}}skipping${{reset}}"
    fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def callVariants2(
    script: TextIOWrapper, reference: str, interval: str, bqsr: str, vcf: str
):
    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        echo Starting variant calling for {INTERVAL}

        bcftools mpileup \\
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
            --max-depth 250 \\
            --no-BAQ \\
            --threads 4 \\
            --output-type u \\
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

        echo Completed variant calling for {INTERVAL}
    else
        echo "Variants already called for {INTERVAL}, ${{green}}skipping${{reset}}"
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
        echo Starting annotation for {INTERVAL}
        
        vep --dir {VEP} \\
            --cache \\
            --format vcf \\
            --vcf \\
            --merged \\
            --fork {FORKS} \\
            --offline \\
            --use_given_ref \\
            --verbose \\
            --force_overwrite \\
            --chr {CHROMOSOME} \\
            --exclude_null_alleles \\
            --symbol \\
            --coding_only \\
            --fasta {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            --input_file {INPUT} \\
            --output_file {OUTPUT} \\
            --stats_file {STATS}/{SUMMARY}

        echo Completed annotation for {INTERVAL}
    else
        echo "Annotations for {INTERVAL} already completed, ${{green}}skipping${{reset}}"
    fi
""".format(
            VEP=vep,
            REFERENCE=reference,
            INPUT=input,
            OUTPUT=output,
            SUMMARY=summary,
            INTERVAL=interval,
            CHROMOSOME=chromosome,
            STATS=stats,
            FORKS=8,
        )
    )


def scatter(script: TextIOWrapper, options: OptionsDict, prefix: str, sorted: str):
    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        splinter(script, bam, sorted, interval[0])

    script.write(
        """
\necho ${yellow}Waiting for scattering to complete${reset}
wait
echo ${green}Scattering completed${reset}\n
        """
    )


def runIntervals(script: TextIOWrapper, options: OptionsDict, prefix: str):
    working = options["working"]
    useAlternateCaller = options["alternateCaller"]

    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        bqsr = """{PREFIX}.{INTERVAL}_bqsr.bam""".format(
            PREFIX=prefix, INTERVAL=interval[1]
        )
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

        annotate(
            script,
            options,
            "{WORKING}/vep_data".format(WORKING=working),
            interval[0],
            vcf,
            vcf.replace(".vcf", ".annotated.vcf"),
            basename(vcf.replace(".vcf", ".annotated.vcf_summary.html")),
        )

        script.write(") &\n")
        script.write("\n")

    script.write(
        """
echo ${yellow}Waiting for intervals to complete${reset}
wait
echo ${green}Intervals processed${reset}
        """
    )


def gather(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# Gather interval data and recombine(s)
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.final.vcf ]]; then
    /bin/ls -1 {PIPELINE}/*.annotated.vcf | sort -k1,1V >{PIPELINE}/merge.list

    gatk MergeVcfs \\
        --VERBOSITY ERROR \\
        -I {PIPELINE}/merge.list \\
        -O {PIPELINE}/{SAMPLE}.final.vcf

    gatk IndexFeatureFile \\
        --verbosity ERROR \\
        -I {PIPELINE}/{SAMPLE}.final.vcf
else
    echo "VCFs already merged and index created, ${{green}}skipping${{reset}}"
fi

echo ${{yellow}}Waiting for VCF merge to complete${{reset}}
wait
echo ${{green}}Final VCF merge completed${{reset}}
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
echo "Starting Variant QC processes"

if [[ ! -f {STATS}/{SAMPLE}.variant_calling_detail_metrics ]]; then
    (
        gatk CollectVariantCallingMetrics \\
            --VERBOSITY ERROR \\
            --DBSNP {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            -I {PIPELINE}/{SAMPLE}.final.vcf \\
            -O {STATS}/{SAMPLE}
     ) &
else
    echo "Variant metrics already run, ${{green}}skipping${{reset}}"
fi

#
# we need to quiet vcftools here because it's stupid chatty and doesn't have an option to quiet
#
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --freq2 --out {STATS}/{SAMPLE} --max-alleles 2 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --depth --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --site-mean-depth --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --site-quality --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --missing-indv --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --missing-site --out {STATS}/{SAMPLE} 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.final.vcf --het --out {STATS}/{SAMPLE} 2>/dev/null &

echo ${{yellow}}Waiting for variant QC metrics to complete${{reset}}
wait
echo ${{green}}Variant QC metrics completed${{reset}}
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
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    echo "Alignment metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    gatk CollectGcBiasMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    echo "GC bias metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then
    gatk CollectWgsMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.wgs_metrics.txt \\
        --READ_LENGTH 151 \\
        -INTERVALS {REFERENCE}/ref_genome_autosomal.interval_list \\
        --USE_FAST_ALGORITHM \\
        --INCLUDE_BQ_HISTOGRAM &
else
    echo "WGS metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then
    samtools stats -@ 8 \\
        -r reference/Homo_sapiens_assembly38.fasta \\
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
multiqc --tag DNA --verbose -f {STATS}

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
    sample = options["sample"]

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all intermediate interval files\n")
    script.write("#\n")

    bam = """{PREFIX}.chr*.bam""".format(PREFIX=prefix)
    bqsr = """{PREFIX}.chr*_bqsr.bam""".format(PREFIX=prefix)
    vcf = """{PREFIX}.chr*.vcf""".format(PREFIX=prefix)

    script.write("rm -f {BAM}\n".format(BAM=bam))
    script.write("rm -f {BAM}.bai\n".format(BAM=bam))

    script.write("\n")
    script.write("rm -f {BQSR}\n".format(BQSR=bqsr))
    script.write("rm -f {BQSR}.bai\n".format(BQSR=bqsr.replace(".bam", "")))
    script.write("rm -f {BQSR}.table\n".format(BQSR=bqsr))

    script.write("\n")
    script.write("rm -f {VCF}\n".format(VCF=vcf))
    script.write("rm -f {SAMPLE}.chr*.vcf.idx".format(SAMPLE=sample))
    script.write("rm -f {VCF}.idx\n".format(VCF=vcf))
    script.write("rm -f {VCF}\n".format(VCF=vcf).replace(".vcf", ".html"))

    for type in ["snps", "indels"]:
        script.write(
            """
rm -f {SAMPLE}.chr*.{TYPE}.vcf.idx
            """.format(
                PIPELINE=pipeline, TYPE=type, SAMPLE=sample
            )
        )

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all remaining intermediate files\n")
    script.write("#\n")
    script.write("\n")
    for type in ["snps", "indels"]:
        script.write(
            "rm -f {PIPELINE}/merge.{TYPE}.list\n".format(PIPELINE=pipeline, TYPE=type)
        )

    for type in ["snps", "indels"]:
        script.write(
            """
rm -f {PIPELINE}/{SAMPLE}.{TYPE}.final.vcf
rm -f {PIPELINE}/{SAMPLE}.{TYPE}.final.vcf.idx
""".format(
                PIPELINE=pipeline, SAMPLE=sample, TYPE=type
            )
        )

    script.write("\n")


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastaPair):
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
    script.write("#   R1 = {F}\n".format(F=filenames[0]))
    script.write("#   R2 = {F}\n".format(F=filenames[1]))
    script.write("#\n")

    script.write("#\n")
    script.write(
        "# Assumed chromosome sizes (from {SIZES})\n".format(
            SIZES=options["chromosomeSizes"]
        )
    )
    intervals = loadIntervals(options["chromosomeSizes"])
    for interval in intervals.keys():
        script.write(
            "#   {CHROME} = {SIZE}\n".format(CHROME=interval, SIZE=intervals[interval])
        )
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
    script.write(
        "#   last block max size = {P}\n".format(
            P=math.floor(options["segmentSize"] * options["factor"])
        )
    )


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

    if (
        exists(expandvars(filenames[0])) == False
        or exists(expandvars(filenames[1])) == False
    ):
        print(
            "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                R1=filenames[0],
                R2=filenames[1],
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

        if options["doQC"]:
            doVariantQC(script, options)
            runMultiQC(script, options)

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, prefix, options)

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

        script.write("\n")
        script.write(
            "touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"])
        )
        script.write("\n")


if __name__ == "__main__":
    main()
