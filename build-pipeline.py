#!/usr/bin/env python

import json
import argparse
import math
from datetime import datetime
from math import ceil
from os.path import exists, expandvars
from os import cpu_count

# This script generates an Illumina paired-read VCF annotation pipeline.
# There are some short-circuit options generated as part of this script,
# but they will not normally be used in production. They are meant to
# prevent a developer from re-running various steps in the pipeline in
# the case that the output files are already present.


def loadIntervals(chromosomeSizes):
    with open(chromosomeSizes, "r") as file:
        chromosomeSizes = json.load(file)
        return chromosomeSizes


def loadChromosomeList(chromosomeSizes):
    with open(chromosomeSizes, "r") as file:
        chromosomeSizes = json.load(file)
        return chromosomeSizes.keys()


def computeIntervals(options):
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


def getFileNames(options):
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    )


def updateDictionary(script, options):
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
    echo "Reference dictionary {REFERENCE}/Homo_sapiens_assembly38.dict already present"
fi

if [[ ! -f {REFERENCE}/ref_genome_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    egrep '({REGEX})\\s' {REFERENCE}/Homo_sapiens_assembly38.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/Homo_sapiens_assembly38.dict - >{REFERENCE}/ref_genome_autosomal.interval_list
else
    echo "Interval list {REFERENCE}/ref_genome_autosomal.interval_list already present"
fi

""".format(
            REFERENCE=reference, BIN=bin, REGEX=regex
        )
    )


def alignAndSort(script, r1, r2, options, output):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    timeout = options["watchdog"]
    stats = options["stats"]
    bin = options["bin"]
    nonRepeatable = options["non-repeatable"]
    readLimit = int(options["read-limit"])

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    # this is one of the more time-consuming operations, so, if we've already
    # done the alignment operation for this sample, we'll skip this
    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
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
            - >{PIPELINE}/{SAMPLE}.aligned.bam'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed alignment process errno = $status"
        rm -f {SAMPLE}.aligned.bam
        exit $status
    fi
else
    echo "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, skipping"
fi

if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    timeout {TIMEOUT}m bash -c \\
        'LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        bamsormadup \\
            SO=coordinate \\
            threads={THREADS} \\
            level=6 \\
            tmpfile={PIPELINE}/{SAMPLE} \\
            inputformat=sam \\
            indexfilename={SORTED}.bai \\
            M={STATS}/{SAMPLE}.duplication_metrics <{PIPELINE}/{SAMPLE}.aligned.bam >{SORTED}'

    status=$?
    if [ $status -ne 0 ]; then
        echo "Watchdog timer killed sort / dup process errno = $status"
        rm -f {SORTED}
        rm -f {SORTED}.bai
        exit $status
    fi
else
    echo "{SORTED}, index, and metrics found, not aligning"
fi

""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            SORT_THREADS=threads,
            IO_THREADS=threads // 2,
            SORTED=output,
            PIPELINE=pipeline,
            TIMEOUT=timeout,
            STATS=stats,
            BIN=bin,
            DASHK="" if nonRepeatable == True else "-K " + str((10_000_000 * threads)),
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def splinter(script, bam, sorted, interval):
    script.write(
        """
\n# Scatter interval {INTERVAL}        
if [[ ! -f {BAM} || ! -f {BAM}.bai ]]; then
    echo "Creating interval {INTERVAL}"
    (samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM} && samtools index -@ 4 {BAM})&
else
    echo "Splinter for {INTERVAL} has been computed, not re-splintering"
fi""".format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )


def genBQSR(script, reference, interval, bam, bqsr):
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
            echo "BQSR table generation for {INTERVAL} skipped"
        fi

        # run base quality score recalibration - apply the bqsr table
        if [[ ! -f {BQSR} || ! -f {BQSR}.bai ]]; then
            gatk ApplyBQSR --java-options '-Xmx8g' \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -I {BAM} \\
                -O {BQSR} \\
                --verbosity ERROR \\
                --preserve-qscores-less-than 6 \\
                --static-quantized-quals 10 \\
                --static-quantized-quals 20 \\
                --static-quantized-quals 30 \\
                --bqsr-recal-file {BQSR}.table \\
                -L {INTERVAL}

            # index that file, this is our target for getting vcf
            samtools index -@ 4 {BQSR}
        else
            echo "BQSR application for {INTERVAL} already completed"
        fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr
        )
    )


def callVariants(script, reference, interval, bqsr, vcf):
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
            echo "Variants already called for {INTERVAL}, skipping"
        fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def filterSNPs(script, reference, vcf, interval, snps, filtered):
    script.write(
        """
        # pull snps out of out called variants and annotate them
        if [[ ! -f {SNPS} ]]; then
            gatk SelectVariants --java-options '-Xmx8g' \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {VCF} \\
                --verbosity ERROR \\
                -select-type SNP \\
                -O {SNPS}
        else
            echo "SNPs already selected for {INTERVAL}, skipping"
        fi

        if [[ ! -f {FILTERED} ]]; then
            gatk VariantFiltration --java-options '-Xmx8g' \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {SNPS} \\
                --verbosity ERROR \\
                --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \\
                --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \\
                --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \\
                --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \\
                --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \\
                -O {FILTERED}
        else
            echo "SNPs already filtered for {INTERVAL}, skipping"
        fi
""".format(
            REFERENCE=reference,
            VCF=vcf,
            SNPS=snps,
            FILTERED=filtered,
            INTERVAL=interval,
        )
    )


def filterINDELs(script, reference, vcf, interval, indels, filtered):
    script.write(
        """
        # pull indels out of out called variants and annotate them
        if [[ ! -f {INDELS} ]]; then
            gatk SelectVariants --java-options '-Xmx8g' \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {VCF} \\
                --verbosity ERROR \\
                -select-type INDEL \\
                -O {INDELS}
        else
            echo "INDELs already selected for {INTERVAL}, skipping"
        fi

        if [[ ! -f {FILTERED} ]]; then
            gatk VariantFiltration --java-options '-Xmx8g' \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {INDELS} \\
                --verbosity ERROR \\
                --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \\
                --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \\
                --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS_lt_n20" \\
                -O {FILTERED}
        else
            echo "INDELs already filtered for {INTERVAL}, skipping"
        fi
""".format(
            REFERENCE=reference,
            VCF=vcf,
            INDELS=indels,
            FILTERED=filtered,
            INTERVAL=interval,
        )
    )


def filterVariants(script, reference, interval, vcf):
    snps = vcf.replace(".vcf", ".snps.vcf")
    filtered = vcf.replace(".vcf", ".snps.filtered.vcf")
    filterSNPs(script, reference, vcf, interval, snps, filtered)

    indels = vcf.replace(".vcf", ".indels.vcf")
    filtered = vcf.replace(".vcf", ".indels.filtered.vcf")
    filterINDELs(script, reference, vcf, interval, indels, filtered)


def annotate(script, options, vep, type, interval, input, output, summary):
    reference = options["reference"]
    chromosome = interval.split(":")[0]

    script.write(
        """
        if [[ ! -f {INPUT} || ! -f {SUMMARY} ]]; then
            echo Starting {TYPE} annotation for {INTERVAL}
            
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
                --fasta {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                --input_file {INPUT} \\
                --output_file {OUTPUT} \\
                --stats_file {SUMMARY}

            echo Completed {TYPE} annotation for {INTERVAL}
        else
            echo "{TYPE} annotations for {INTERVAL} already completed, skipping"
        fi
""".format(
            VEP=vep,
            REFERENCE=reference,
            TYPE=type,
            INPUT=input,
            OUTPUT=output,
            SUMMARY=summary,
            INTERVAL=interval,
            CHROMOSOME=chromosome,
            FORKS=8,
        )
    )


def annotateVariants(script, options, vcf, interval):
    working = options["working"]

    vep = "{WORKING}/vep_data".format(WORKING=working)

    filtered = vcf.replace(".vcf", ".indels.filtered.vcf")
    annotated = vcf.replace(".vcf", ".indels.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".indels.filtered.annotated.vcf_summary.html")
    annotate(script, options, vep, "INDEL", interval, filtered, annotated, summary)

    filtered = vcf.replace(".vcf", ".snps.filtered.vcf")
    annotated = vcf.replace(".vcf", ".snps.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".snps.filtered.annotated.vcf_summary.html")
    annotate(script, options, vep, "SNP", interval, filtered, annotated, summary)


def scatter(script, options, prefix, sorted):
    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        splinter(script, bam, sorted, interval[0])

    script.write("\n\necho Waiting for scattering to complete\n")
    script.write("wait\n")


def runIntervals(script, options, prefix):
    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        bqsr = """{PREFIX}.{INTERVAL}_bqsr.bam""".format(
            PREFIX=prefix, INTERVAL=interval[1]
        )
        vcf = """{PREFIX}.{INTERVAL}.vcf""".format(PREFIX=prefix, INTERVAL=interval[1])

        script.write("\n")
        script.write("    #\n")
        script.write("    # Run interval {INTERVAL}\n".format(INTERVAL=interval[0]))
        script.write("    #\n")
        script.write("    (\n")

        genBQSR(script, options["reference"], interval[0], bam, bqsr)
        callVariants(script, options["reference"], interval[0], bqsr, vcf)
        filterVariants(script, options["reference"], interval[0], vcf)
        annotateVariants(script, options, vcf, interval[0])

        script.write("    ) &\n")
        script.write("\n")

    script.write("echo Waiting for intervals to complete\n")
    script.write("wait\n")


def gather(script, options):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# Gather interval data and recombine(s)
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.snps.final.vcf ]]; then
    /bin/ls -1 {PIPELINE}/*.snps.filtered.annotated.vcf >{PIPELINE}/merge.snps.list

    gatk MergeVcfs \\
        --VERBOSITY ERROR \\
        -I {PIPELINE}/merge.snps.list \\
        -O {PIPELINE}/{SAMPLE}.snps.final.vcf &
else
    echo "SNPs already merged, skipping"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.indels.final.vcf ]]; then
    /bin/ls -1 {PIPELINE}/*.indels.filtered.annotated.vcf >{PIPELINE}/merge.indels.list

    gatk MergeVcfs \\
        --VERBOSITY ERROR \\
        -I {PIPELINE}/merge.indels.list \\
        -O {PIPELINE}/{SAMPLE}.indels.final.vcf &
else
    echo "INDELs already merged, skipping"
fi

echo Waiting for SNP and INDEL merge to complete
wait
""".format(
            PIPELINE=pipeline, SAMPLE=sample
        )
    )


def mergeFinal(script, options):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# Create final VCF file(s)
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.final.unfiltered.vcf ]]; then
    gatk MergeVcfs \\
        --VERBOSITY ERROR \\
        -I {PIPELINE}/{SAMPLE}.snps.final.vcf \\
        -I {PIPELINE}/{SAMPLE}.indels.final.vcf \\
        -O {PIPELINE}/{SAMPLE}.final.unfiltered.vcf
else
    echo "Unfiltered INDEL and SNP VCFs already merged, skipping"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.final.filtered.vcf ]]; then
    gatk SelectVariants \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -V {PIPELINE}/{SAMPLE}.final.unfiltered.vcf \\
        -O {PIPELINE}/{SAMPLE}.final.filtered.vcf \\
        --verbosity ERROR \\
        --exclude-filtered \\
        --exclude-non-variants \\
        --remove-unused-alternates
else
    echo "Filtered INDEL and SNP VCFs already merged, skipping"
fi
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def doVariantQC(script, options):
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

if [[ ! -f {STATS}/{SAMPLE}_unfilt.variant_calling_detail_metrics ]]; then
    (gatk CollectVariantCallingMetrics \\
        --VERBOSITY ERROR \\
        --DBSNP {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
        -I {PIPELINE}/{SAMPLE}.final.unfiltered.vcf \\
        -O {STATS}/{SAMPLE}_unfilt
    
     sed -i 's/^{SAMPLE}/{SAMPLE}_unfiltered/' {STATS}/{SAMPLE}_unfilt.variant_calling_detail_metrics
     sed -i 's/^{SAMPLE}/{SAMPLE}_unfiltered/' {STATS}/{SAMPLE}_unfilt.variant_calling_summary_metrics) &
else
    echo "Variant (unfiltered) metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}_filt.variant_calling_detail_metrics ]]; then
    (gatk CollectVariantCallingMetrics \\
         --VERBOSITY ERROR \\
         --DBSNP {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
         -I {PIPELINE}/{SAMPLE}.final.filtered.vcf \\
         -O {STATS}/{SAMPLE}_filt

     sed -i 's/^{SAMPLE}/{SAMPLE}_filtered/' {STATS}/{SAMPLE}_filt.variant_calling_detail_metrics
     sed -i 's/^{SAMPLE}/{SAMPLE}_filtered/' {STATS}/{SAMPLE}_filt.variant_calling_summary_metrics) &
else
    echo "Variant (filtered) metrics already run, skipping"
fi
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )


def runQC(script, options, sorted):
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
    echo "samtools flagstat already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.alignment_metrics.txt ]]; then
    gatk CollectAlignmentSummaryMetrics \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    echo "Alignment metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    gatk CollectGcBiasMetrics \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    echo "GC bias metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then
    gatk CollectWgsMetrics \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.wgs_metrics.txt \\
        --READ_LENGTH 151 \\
        -INTERVALS {REFERENCE}/ref_genome_autosomal.interval_list \\
        --USE_FAST_ALGORITHM \\
        --INCLUDE_BQ_HISTOGRAM &
else
    echo "WGS metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.sorted_fastqc.zip || ! -f {STATS}/{SAMPLE}.sorted_fastqc.html ]]; then
    fastqc \\
        --outdir {STATS} \\
        --noextract \\
        {SORTED} &
else
    echo "FASTQC already run, skipping"
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


def runMultiQC(script, options):
    stats = options["stats"]
    sample = options["sample"]

    script.write(
        """
#
# Run MultiQC across everything
# 
(
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
)    

""".format(
            STATS=stats, SAMPLE=sample
        )
    )


def cleanup(script, prefix, options):
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
rm -f {SAMPLE}.chr*.{TYPE}.filtered.vcf.idx
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


def writeHeader(script, options, filenames):
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


def writeVersions(script):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...\n")
    script.write("#\n")


def main():
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
        "-c",
        "--clean",
        action="store_true",
        dest="cleanIntermediateFiles",
        default=False,
        help="Clean up the mess we make",
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
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    parser.add_argument(
        "-d",
        "--watchdog",
        action="store",
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

    opts = parser.parse_args()

    options = vars(opts)

    if options["working"] == None:
        options["working"] = "$HOME"
    if options["reference"] == None:
        options["reference"] = "{WORKING}/reference".format(WORKING=options["working"])
    if options["pipeline"] == None:
        options["pipeline"] = "{WORKING}/pipeline".format(WORKING=options["working"])
    if options["stats"] == None:
        options["stats"] = "{WORKING}/stats".format(WORKING=options["working"])

    if options["sample"] == None:
        print("--sample is a required option")
        return

    for opt in [
        "working",
        "reference",
        "pipeline",
        "stats",
        "bin",
        "script",
        "chromosomeSizes",
    ]:
        options[opt] = expandvars(options[opt])

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
        return

    if exists(options["bin"]) == False:
        print(
            "Unable to find your --bin-dir directory at {PATH}".format(
                PATH=options["bin"]
            )
        )
        return

    if exists(options["working"]) == False:
        print(
            "Unable to find your --work-dir directory at {PATH}".format(
                PATH=options["working"]
            )
        )
        return

    if exists(options["reference"]) == False:
        print(
            "Unable to find your --reference-dir directory at {PATH}".format(
                PATH=options["reference"]
            )
        )
        return

    if exists(options["pipeline"]) == False:
        print(
            "Unable to find your --pipeline-dir directory at {PATH}".format(
                PATH=options["pipeline"]
            )
        )
        return

    if exists(options["stats"]) == False:
        print(
            "Unable to find your --stats-dir directory at {PATH}".format(
                PATH=options["stats"]
            )
        )
        return

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

        script.write("#\n")
        script.write(
            """
# perl stuff
export PATH=/home/ubuntu/perl5/bin:$PATH
export PERL5LIB=/home/ubuntu/perl5/lib/perl5:$PERL5LIB
export PERL_LOCAL_LIB_ROOT=/home/ubuntu/perl5:$PERL_LOCAL_LIB_ROOT

# shared library stuff
export LD_LIBRARY_PATH={WORKING}/bin:/usr/lib64:/usr/local/lib/:$LB_LIBRARY_PATH

# handy path
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.3.0:{WORKING}/bin:$PATH\n""".format(
                WORKING=options["working"]
            )
        )
        script.write("\n")

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

        scatter(
            script,
            options,
            prefix,
            sorted,
        )
        runIntervals(script, options, prefix)
        gather(script, options)
        mergeFinal(script, options)

        if options["doQC"]:
            doVariantQC(script, options)
            runQC(script, options, sorted)

            script.write("\necho Waiting for QC metrics to complete\n")
            script.write("wait\n")

            runMultiQC(script, options)

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, prefix, options)

        script.write(
            """\necho -e "Done processing {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"\n""".format(
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
