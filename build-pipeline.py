#!/usr/bin/env python

import json
import argparse
import math
from datetime import datetime
from math import ceil
from os.path import exists, expandvars
from os import cpu_count, uname, getpid

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
                intervals.append(
                    "chr{c}:{lower}-{upper}".format(
                        c=c,
                        lower=(segments - 2) * segmentSize,
                        upper=chromosomeSizes[c],
                    )
                )
            else:
                intervals.append(
                    "chr{c}:{lower}-{upper}".format(
                        c=c,
                        lower=(segments - 1) * segmentSize,
                        upper=chromosomeSizes[c],
                    )
                )

    return [
        (interval, interval.replace(":", "_").replace("-", "_"))
        for interval in intervals
    ]


def getFileNames(options, trimmed):
    sample = options["sample"]
    pipeline = options["pipeline"]

    if trimmed == False:
        return (
            "{PIPELINE}/{SAMPLE}_R1.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
            "{PIPELINE}/{SAMPLE}_R2.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        )
    else:
        return (
            "{PIPELINE}/{SAMPLE}_R1_trimmed.fq.gz".format(
                PIPELINE=pipeline, SAMPLE=sample
            ),
            "{PIPELINE}/{SAMPLE}_R2_trimmed.fq.gz".format(
                PIPELINE=pipeline, SAMPLE=sample
            ),
        )


def genTrimmer(script, r1, r2, options):
    pipeline = options["pipeline"]
    stats = options["stats"]
    threads = options["cores"]

    script.write("#\n")
    script.write("# Trim reads\n")
    script.write("#\n")

    # because we're in the trimmer we can assume that we want
    # trimmed output files, so we need to get those names
    # and check for a short-circuit if appropriate
    filenames = getFileNames(options, True)

    script.write(
        """
if [[ ! -f {TRIMMED_R1} || ! -f {TRIMMED_R2} ]]; then
    trim_galore \\
        --illumina \\
        --cores {THREADS} \\
        --output_dir {PIPELINE} \\
        {R1} \\
        {R2}

    mv {TRIMMED_R1}_trimming_report.txt {STATS}
    mv {TRIMMED_R2}_trimming_report.txt {STATS}
else
    echo "{TRIMMED_R1} and {TRIMMED_R2} found, not trimming"
fi        
""".format(
            R1=r1,
            R2=r2,
            PIPELINE=pipeline,
            THREADS=threads,
            TRIMMED_R1=filenames[0],
            TRIMMED_R2=filenames[1],
            STATS=stats,
        )
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
    echo "Reference dictionary {REFERENCE}/Homo_sapiens_assembly38.dict already present"
fi

""".format(
            REFERENCE=reference, BIN=bin, REGEX=regex
        )
    )


def genBWA(script, r1, r2, options, output):
    reference = options["reference"]
    sample = options["sample"]
    working = options["working"]
    stats = options["stats"]
    pipeline = options["pipeline"]
    threads = options["cores"]

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    # each side of the pipeline also wants a read / write thread to push data
    # through the UNIX pipe, we'll provide one here by giving up one
    # processing thread per side
    threads -= 1

    # this is one of the more time-consuming operations, so, if we've already
    # done the alignment operation for this sample, we'll skip this
    script.write(
        """
if [[ ! -f {SORTED} || ! -f {SORTED}.bai || ! -f {STATS}/{SAMPLE}.duplication_metrics ]]; then
    bwa-mem2 mem -t {THREADS} \\
        {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        {R1} \\
        {R2} \\
        -Y \\
        -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:MJS.SEQUENCER.7\\tLB:{SAMPLE}\\tSM:{SAMPLE}" |
    bamsormadup \\
        SO=coordinate \\
        threads={THREADS} \\
        level=0 \\
        tmpfile={TMP}/bamsormapdup_{NODENAME}_{PID} \\
        inputformat=sam \\
        indexfilename={SORTED}.bai \\
        M={STATS}/{SAMPLE}.duplication_metrics >{SORTED}
else
    echo "{SORTED}, index, and metrics found, not aligning"
fi
""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            WORKING=working,
            STATS=stats,
            THREADS=threads,
            SORTED=output,
            TMP=pipeline,
            NODENAME=uname().nodename,
            PID=getpid(),
        )
    )


def splinter(script, bam, sorted, interval):
    script.write(
        """
        # Splinter our interval off the main file
        if [[ ! -f {BAM} || ! -f {BAM}.bai ]]; then
            samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM}
            samtools index -@ 4 {BAM}
        else
            echo "Splinter for {INTERVAL} has been computed, not re-splintering"
        fi
    """.format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )


def genBQSR(script, reference, interval, bam, bqsr):
    script.write(
        """
        # run base quality score recalibration - build the bqsr table
        if [[ ! -f {BQSR}.table ]]; then
            gatk BaseRecalibrator \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -I {BAM} \\
                -O {BQSR}.table \\
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
            gatk ApplyBQSR \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -I {BAM} \\
                -O {BQSR} \\
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


def callVariants(script, reference, interval, bam, bqsr, vcf):
    script.write(
        """
        # call variants
        if [[ ! -f {VCF} ]]; then
            gatk HaplotypeCaller \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -I {BQSR} \\
                -O {VCF} \\
                --pairHMM AVX_LOGLESS_CACHING_OMP \\
                --native-pair-hmm-threads 8 \\
                -L {INTERVAL}
        else
            echo "Variants already called for {INTERVAL}, skipping"
        fi
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr, VCF=vcf
        )
    )


def filterSNPs(script, reference, vcf, interval, snps, filtered):
    script.write(
        """
        # pull snps out of out called variants and annotate them
        if [[ ! -f {SNPS} || ! -f {FILTERED} ]]; then
            gatk SelectVariants \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {VCF} \\
                -select-type SNP \\
                -O {SNPS}

            gatk VariantFiltration \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {SNPS} \\
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
        if [[ ! -f {INDELS} || ! -f {FILTERED} ]]; then
            gatk SelectVariants \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {VCF} \\
                -select-type INDEL \\
                -O {INDELS}

            gatk VariantFiltration \\
                -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
                -V {INDELS} \\
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


def annotate(script, options, vep, type, vcf, interval, filtered, annotated, summary):
    reference = options["reference"]
    pipeline = options["pipeline"]

    input = filtered.split("/")[-1]
    output = annotated.split("/")[-1]

    script.write(
        """
        # run vep in a docker container to annotate the vcf
        if [[ ! -f {ANNOTATED} || ! -f {SUMMARY} ]]; then
            sudo docker run \\
                -v {VEP}:/opt/vep/.vep:Z \\
                -v {PIPELINE}:/opt/vep/.vep/input:Z \\
                -v {PIPELINE}:/opt/vep/.vep/output:Z \\
                -v {REFERENCE}:/opt/vep/.vep/reference:Z \\
                ensemblorg/ensembl-vep \\
                    ./vep --cache --format vcf --merged --offline --use_given_ref --vcf --verbose \\
                    --fasta /opt/vep/.vep/reference/Homo_sapiens_assembly38.fasta \\
                    --input_file /opt/vep/.vep/input/{INPUT} \\
                    --output_file /opt/vep/.vep/output/{OUTPUT}
        else
            echo "{TYPE} annotations for {INTERVAL} already completed, skipping"
        fi
""".format(
            VEP=vep,
            PIPELINE=pipeline,
            REFERENCE=reference,
            TYPE=type,
            VCF=vcf,
            ANNOTATED=annotated,
            FILTERED=filtered,
            SUMMARY=summary,
            INPUT=input,
            OUTPUT=output,
            INTERVAL=interval,
        )
    )


def annotateVariants(script, options, vcf, interval):
    working = options["working"]

    vep = "{WORKING}/vep_data".format(WORKING=working)

    filtered = vcf.replace(".vcf", ".indels.filtered.vcf")
    annotated = vcf.replace(".vcf", ".indels.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".indels.filtered.annotated.vcf_summary.html")
    annotate(script, options, vep, "INDEL", vcf, interval, filtered, annotated, summary)

    filtered = vcf.replace(".vcf", ".snps.filtered.vcf")
    annotated = vcf.replace(".vcf", ".snps.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".snps.filtered.annotated.vcf_summary.html")
    annotate(script, options, vep, "SNP", vcf, interval, filtered, annotated, summary)


def runIntervals(script, options, prefix, sorted):
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

        splinter(script, bam, sorted, interval[0])
        genBQSR(script, options["reference"], interval[0], bam, bqsr)
        callVariants(script, options["reference"], interval[0], bam, bqsr, vcf)
        filterVariants(script, options["reference"], interval[0], vcf)
        annotateVariants(script, options, vcf, interval[0])

        script.write("    ) &\n")
        script.write("\n")

    script.write("echo Waiting for intervals to complete\n")
    script.write("wait\n")


def merge(script, options):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# Create final VCF file(s)
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.snps.final.vcf ]]; then
    /bin/ls -1 {PIPELINE}/*.snps.filtered.annotated.vcf >{PIPELINE}/merge.snps.list

    gatk MergeVcfs \\
        -I {PIPELINE}/merge.snps.list \\
        -O {PIPELINE}/{SAMPLE}.snps.final.vcf &
else
    echo "SNPs already merged, skipping"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.indels.final.vcf ]]; then
    /bin/ls -1 {PIPELINE}/*.indels.filtered.annotated.vcf >{PIPELINE}/merge.indels.list

    gatk MergeVcfs \\
        -I {PIPELINE}/merge.indels.list \\
        -O {PIPELINE}/{SAMPLE}.indels.final.vcf &
else
    echo "INDELs already merged, skipping"
fi

echo Waiting for SNP and INDEL merge to complete
wait

if [[ ! -f {PIPELINE}/{SAMPLE}.final.unfiltered.vcf ]]; then
    gatk MergeVcfs \\
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
        --exclude-filtered
else
    echo "Filtered INDEL and SNP VCFs already merged, skipping"
fi
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def runQC(script, options):
    reference = options["reference"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]

    filenames = getFileNames(options, False)

    script.write(
        """
#
# RUN QC process
# 
if [[ ! -f {PIPELINE}/{SAMPLE}.merged.bqsr.bam || ! -f {PIPELINE}/{SAMPLE}.merged.bqsr.bam.bai ]]; then
    samtools merge -@ {THREADS} -o {PIPELINE}/{SAMPLE}.merged.bqsr.bam {PIPELINE}/*bqsr*.bam
    samtools index -@ {THREADS} {PIPELINE}/{SAMPLE}.merged.bqsr.bam
else
    echo "Intermediate BQSR files already merged and indexed, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.bqsr.flagstat.txt ]]; then
    samtools flagstat \\
        {PIPELINE}/{SAMPLE}.merged.bqsr.bam >{STATS}/{SAMPLE}.bqsr.flagstat.txt &
else
    echo "samtools flagstat already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.bqsr.insert_metrics.txt || ! -f {STATS}/{SAMPLE}.bqsr.insert_metrics.pdf ]]; then
    gatk CollectInsertSizeMetrics \\
        -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
        -O {STATS}/{SAMPLE}.bqsr.insert_metrics.txt \\
        -H {STATS}/{SAMPLE}.bqsr.insert_metrics.pdf \\
        -M 0.5 &
else
    echo "Insert metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.bqsr.alignment_metrics.txt ]]; then
    gatk CollectAlignmentSummaryMetrics \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
        -O {STATS}/{SAMPLE}.bqsr.alignment_metrics.txt &
else
    echo "Alignment metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.bqsr.gc_bias_summary.txt ]]; then
    gatk CollectGcBiasMetrics \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
        -O {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.bqsr.gc_bias_summary.txt &
else
    echo "GC bias metrics already run, skipping"
fi

if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then
    gatk CollectWgsMetrics \\
        -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
        -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
        -O {STATS}/{SAMPLE}.wgs_metrics.txt \\
        --READ_LENGTH 151 \\
        -INTERVALS {REFERENCE}/ref_genome_autosomal.interval_list \\
        --USE_FAST_ALGORITHM \\
        --INCLUDE_BQ_HISTOGRAM &
else
    echo "WGS metrics already run, skipping"
fi
""".format(
            R1=filenames[0],
            R2=filenames[1],
            REFERENCE=reference,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
        )
    )

    if options["trim"] == True:
        trimmed = getFileNames(options, True)

        script.write(
            """
fastqc \\
    --threads=5 \\
    --outdir {STATS} \\
    --noextract \\
    {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    {R1} \\
    {R2} \\
    {TRIMMED_R1} \\
    {TRIMMED_R2} &
""".format(
                R1=filenames[0],
                R2=filenames[1],
                TRIMMED_R1=trimmed[0],
                TRIMMED_R2=trimmed[1],
                REFERENCE=reference,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
            )
        )
    else:
        script.write(
            """
fastqc \\
    --threads=3 \\
    --outdir {STATS} \\
    --noextract \\
    {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    {R1} \\
    {R2} &
""".format(
                R1=filenames[0],
                R2=filenames[1],
                REFERENCE=reference,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
            )
        )

    script.write("\necho Waiting for QC metrics to complete\n")
    script.write("wait\n")


def runMultiQC(script, options):
    stats = options["stats"]

    script.write(
        """
#
# Run MultiQC across everything
# 
(
    ## TODO - find out how to short-circuit this
    mkdir -p {STATS}/qc
    cd {STATS}/qc
    multiqc -f {STATS}
)    

""".format(
            STATS=stats
        )
    )


def cleanup(script, prefix, options):
    pipeline = options["pipeline"]
    sample = options["sample"]

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
    script.write("rm -f {BQSR}.bai\n".format(BQSR=bqsr))
    script.write("rm -f {BQSR}.table\n".format(BQSR=bqsr))

    script.write("\n")
    script.write("rm -f {VCF}\n".format(VCF=vcf))
    script.write("rm -f {VCF}\n".format(VCF=vcf).replace(".vcf", ".html"))

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
    for type in ["snps", "indels"]:
        script.write(
            "rm -f {PIPELINE}/{SAMPLE}.{TYPE}.final.vcf\n".format(
                PIPELINE=pipeline, SAMPLE=sample, TYPE=type
            )
        )

    script.write("\n")
    script.write(
        "rm -f {PIPELINE}/{SAMPLE}.merged.bqsr.bam\n".format(
            PIPELINE=pipeline, SAMPLE=sample
        )
    )
    script.write(
        "rm -f {PIPELINE}/{SAMPLE}.merged.bqsr.bam.bai\n".format(
            PIPELINE=pipeline, SAMPLE=sample
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
    parser.set_defaults(doQC=False, trim=True, cleanIntermediateFiles=True)
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
        "-t",
        "--trim",
        action="store_true",
        dest="trim",
        default=False,
        help="Run trim-galore on the input FASTQ. Use this with great caution because it can blow up the read pair matching.",
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
        "-z",
        "--cores",
        action="store",
        dest="cores",
        default=cpu_count(),
        help="Specify the number of available CPU",
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
        "--sizes",
        action="store",
        metavar="CHROME_SIZES",
        dest="chromosomeSizes",
        default="chromosomeSizes.json",
        help="Name of JSON file containing chromosome sizes",
    )

    parser.add_argument(
        "--segment",
        action="store",
        metavar="SEGMENT_SIZE",
        dest="segmentSize",
        default=50_000_000,
        help="Size of interval partition",
    )
    parser.add_argument(
        "--factor",
        action="store",
        metavar="FACTOR",
        dest="factor",
        default=0.25,
        help="Interval remainder buffer (between 0.10 and 0.50)",
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

    filenames = getFileNames(options, False)

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
            """export PATH={WORKING}/bin/FastQC:{WORKING}/bin/TrimGalore-0.6.7:{WORKING}/bin/gatk-4.2.3.0:{WORKING}/bin:$PATH\n""".format(
                WORKING=options["working"]
            )
        )
        script.write("\n")

        updateDictionary(script, options)

        # assume we're not trimming, this gets the original R1/R2
        filenames = getFileNames(options, False)

        if options["trim"] == True:
            genTrimmer(script, filenames[0], filenames[1], options)

            # trimming writes new files under different names, this will grab those
            filenames = getFileNames(options, options["trim"])

        genBWA(
            script,
            filenames[0],
            filenames[1],
            options,
            sorted,
        )
        runIntervals(
            script,
            options,
            prefix,
            sorted,
        )
        merge(script, options)

        if options["doQC"]:
            runQC(script, options)
            runMultiQC(script, options)

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, prefix, options)

        script.write(
            """\necho "Done processing {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"\n""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                PIPELINE=options["pipeline"],
            )
        )


if __name__ == "__main__":
    main()
