#!/usr/bin/env python

from io import TextIOWrapper
import json
import argparse
import math

from argparse import Namespace
from typing import Dict, Tuple, Sequence, List, Any
from datetime import datetime
from math import ceil
from os.path import exists, expandvars
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
FastqSet = Tuple[str, str]


def loadIntervals(rulesFile: str):
    with open(rulesFile, "r") as file:
        chromosomeSizes = json.load(file)
        return chromosomeSizes


def computeIntervals(options: OptionsDict):
    intervals = []

    rulesFile = options["chromosomeSizes"]
    segmentSize = options["segmentSize"]
    lastBlockMax = math.floor(options["segmentSize"] * options["factor"])

    with open(rulesFile, "r") as file:
        rules = json.load(file)

        for rule in rules:
            remainder = rule["length"] - segmentSize
            segments = ceil(rule["length"] / segmentSize)
            segment = 0

            while remainder > lastBlockMax:
                lower = segment * segmentSize + 1
                upper = (segment + 1) * segmentSize
                intervals.append(
                    "{accession}:{lower}-{upper}".format(accession=rule["accession"], lower=lower, upper=upper)
                )

                segment += 1
                remainder -= segmentSize

            if remainder > 0:
                lower = (segments - 2) * segmentSize
                if lower % 10 == 0:
                    lower += 1

                intervals.append(
                    "{accession}:{lower}-{upper}".format(
                        accession=rule["accession"],
                        lower=lower,
                        upper=rule["length"],
                    )
                )
            else:
                lower = (segments - 1) * segmentSize
                if lower % 10 == 0:
                    lower += 1

                intervals.append(
                    "{accession}:{lower}-{upper}".format(
                        accession=rule["accession"],
                        lower=lower,
                        upper=rule["length"],
                    )
                )

    return [(interval, interval.replace(":", "_").replace("-", "_")) for interval in intervals]


def writeMergeList(options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]
    intervals = computeIntervals(options)

    mergeList = "{PIPELINE}/{SAMPLE}.merge.list".format(PIPELINE=pipeline, SAMPLE=sample)
    with open(mergeList, "w") as ml:
        ml.truncate()
        for interval in intervals:
            ml.write(
                "{PIPELINE}/{SAMPLE}.{INTERVAL}.vcf\n".format(INTERVAL=interval[1], PIPELINE=pipeline, SAMPLE=sample)
            )


def getFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_R1.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    )


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    bin = options["bin"]
    assembly = options["referenceAssembly"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    chromosomes = [c["accession"] for c in loadIntervals(options["chromosomeSizes"])]
    regex = "|".join(chromosomes)

    script.write(
        """
if [[ ! -f {REFERENCE}/{ASSEMBLY}.dict ]]; then
    logthis "${{yellow}}Creating sequence dictionary${{reset}}"

    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -O {REFERENCE}/{ASSEMBLY}.dict
else
    logthis "Reference dictionary {REFERENCE}/{ASSEMBLY}.dict ${{green}}already completed${{reset}}"
fi

if [[ ! -f {REFERENCE}/{ASSEMBLY}_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    logthis "Building {REFERENCE}/{ASSEMBLY}.interval_list"

    egrep '^({REGEX})\\s' {REFERENCE}/{ASSEMBLY}.fna.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/{ASSEMBLY}.dict - >{REFERENCE}/{ASSEMBLY}_autosomal.interval_list

    logthis "Building {REFERENCE}/{ASSEMBLY}_autosomal.interval_list finished"
else
    logthis "Interval list {REFERENCE}/{ASSEMBLY}_autosomal.interval_list ${{green}}already completed${{reset}}"
fi

""".format(
            REFERENCE=reference, ASSEMBLY=assembly, BIN=bin, REGEX=regex
        )
    )


def runIdentityPreprocessor(script: TextIOWrapper, r1: str, r2: str, o1: str, o2: str):
    script.write(
        """
#
# run the identity preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running identity preprocessor${{reset}}"

    ln -s {R1} {O1}
    ln -s {R2} {O2}
else
    echo "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
        )
    )


def runFastpPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    reference = options["reference"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    readLimit = int(options["read-limit"])

    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running fastp preprocessor${{reset}}"

    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
    fastp \\
        --report_title "fastp report for sample {SAMPLE}" \\
        --in1 {R1} \\
        --in2 {R2} \\
        --out1 {O1} \\
        --out2 {O2} \\
        --detect_adapter_for_pe \\
        --verbose {LIMITREADS} \\
        --thread 8 \\
        -j {STATS}/{SAMPLE}-fastp.json \\
        -h {STATS}/{SAMPLE}-fastp.html
else
    echo "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            R2=r2,
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def preprocessFASTQ(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    preprocessor = options["preprocessor"]

    if preprocessor == "none":
        runIdentityPreprocessor(script, r1, r2, o1, o2)
    elif preprocessor == "fastp":
        runFastpPreprocessor(script, r1, r2, o1, o2, options)
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)


def runBwaAligner(
    script: TextIOWrapper,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    nonRepeatable = options["non-repeatable"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    bwa-mem2 mem -t {THREADS} \\
        -Y -M {DASHK} \\
        -v 1 \\
        -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        {REFERENCE}/{ASSEMBLY}.fna \\
        {O1} \\
        {O2} | 
    samtools view -Sb - >{PIPELINE}/{SAMPLE}.aligned.bam
else
    echo "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            DASHK="" if nonRepeatable == True else "-K " + str((10_000_000 * int(threads))),
        )
    )


def alignFASTQ(
    script: TextIOWrapper,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    aligner = options["aligner"]

    if aligner == "bwa":
        runBwaAligner(script, o1, o2, options)
    else:
        print("Unexpected value {ALIGNER} given for the --aligner option".format(ALIGNER=aligner))
        quit(1)

    pass


def sortWithBiobambam(script: TextIOWrapper, options: OptionsDict, output: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    logthis "${{yellow}}Sorting and marking duplicates${{reset}}"

    bamsormadup \\
        SO=coordinate \\
        threads={THREADS} \\
        level=6 \\
        tmpfile={TEMP}/{SAMPLE} \\
        inputformat=bam \\
        indexfilename={SORTED}.bai \\
        M={STATS}/{SAMPLE}.duplication_metrics <{PIPELINE}/{SAMPLE}.aligned.bam >{SORTED}

    # force the index to look "newer" than its source
    touch {SORTED}.bai
else
    echo "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            THREADS=threads,
            SORTED=output,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
        )
    )

    pass


def sortWithSamtools(script: TextIOWrapper, options: OptionsDict, output: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    unmarked = "{PIPELINE}/{SAMPLE}.unmarked.bam".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {UNMARKED} ]]; then
    logthis "${{yellow}}Sorting aligned file${{reset}}"

    samtools sort {PIPELINE}/{SAMPLE}.aligned.bam -o {UNMARKED}
else
    echo "{UNMARKED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    logthis "${{yellow}}Marking duplicates${{reset}}"

    java -Xmx8g -jar {BIN}/picard.jar MarkDuplicates \\
        --TAGGING_POLICY All \\
        --REFERENCE_SEQUENCE {REFERENCE}/hcov-oc43.fasta \\
        -I {UNMARKED} \\
        -O {SORTED} \\
        -M {STATS}/{SAMPLE}_marked_dup_metrics.txt    

    # generate an index on the result
    samtools index -b {SORTED} {SORTED}.bai
else
    echo "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi
    """.format(
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            UNMARKED=unmarked,
            SORTED=output,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
        )
    )


def sortAlignedAndMappedData(script: TextIOWrapper, options: OptionsDict, output: str):
    sorter = options["sorter"]

    if sorter == "biobambam":
        sortWithBiobambam(script, options, output)
    else:
        sortWithSamtools(script, options, output)


def extractUmappedReads(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# extract unmapped reads
#
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped_R1.fastq || ! -f {PIPELINE}/{SAMPLE}_unmapped_R2.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads into initial FASTQ${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}_unmapped_R1.fastq \\
        -2 {PIPELINE}/{SAMPLE}_unmapped_R2.fastq \\
        {PIPELINE}/{SAMPLE}.aligned.bam
else
    echo "Unmapped reads extracted to initial FASTQ, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            PIPELINE=pipeline,
        )
    )


def alignAndSort(script: TextIOWrapper, options: OptionsDict, output: str):
    processUnmapped = options["processUnmapped"]
    filenames = getFileNames(options)
    trimmedFilenames = getTrimmedFileNames(options)
    alignOnly = options["alignOnly"]

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filenames[0], filenames[1], trimmedFilenames[0], trimmedFilenames[1], options)
    alignFASTQ(script, trimmedFilenames[0], trimmedFilenames[1], options)

    sortAlignedAndMappedData(script, options, output)

    if processUnmapped == True:
        extractUmappedReads(script, options)

    if alignOnly == True:
        script.write(
            """

# align-only flag set
logthis "align-only set, exiting pipeline early"
exit
"""
        )


def genBQSR(script: TextIOWrapper, options: OptionsDict, interval: str, bam: str, bqsr: str):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    knownSites = options["knownSites"]

    script.write(
        """
    # run base quality score recalibration - build the bqsr table
    if [[ ! -f {BQSR}.table ]]; then
        logthis "Generating {BQSR}.table"

        gatk BaseRecalibrator --java-options '-Xmx4g' \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
            -I {BAM} \\
            -O {BQSR}.table \\
            --verbosity ERROR \\
            --preserve-qscores-less-than 6 \\
            --known-sites {REFERENCE}/{KNOWN_SITES} \\
            -L {INTERVAL}

        logthis "{BQSR}.table completed"
    else
        logthis "BQSR table generation for {INTERVAL} ${{green}}already completed${{reset}}"
    fi

    # run base quality score recalibration - apply the bqsr table
    if [[ ! -f {BQSR} ]]; then
        logthis "Applying calibration for {BQSR}"

        gatk ApplyBQSR --java-options '-Xmx4g' \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
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

        logthis "Calibation completed for for {BQSR}"
    else
        logthis "BQSR application for {INTERVAL} ${{green}}already completed${{reset}}"
    fi

    # run base quality score recalibration - apply the bqsr table
    if [[ ! -f {BQSR}.bai ]]; then
        logthis "Creating index for {BQSR}"

        samtools index -@ 4 {BQSR}

        logthis "Indexing completed for for {BQSR}"
    else
        logthis "BQSR index for {INTERVAL} ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, KNOWN_SITES=knownSites, INTERVAL=interval, BAM=bam, BQSR=bqsr
        )
    )


def callVariantsUsingGatk(script: TextIOWrapper, options: OptionsDict, interval: str, bqsr: str, vcf: str):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    knownSites = options["knownSites"]

    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        logthis "Starting variant calling for {INTERVAL}"

        gatk HaplotypeCaller --java-options '-Xmx4g' \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
            -I {BQSR} \\
            -O {VCF} \\
            --verbosity ERROR \\
            --dbsnp {REFERENCE}/{KNOWN_SITES} \\
            --pairHMM FASTEST_AVAILABLE \\
            --native-pair-hmm-threads 4 \\
            -L {INTERVAL}

        logthis "Completed variant calling for {INTERVAL}"
    else
        logthis "Variants already called for {INTERVAL}, ${{green}}already completed${{reset}}"
    fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, KNOWN_SITES=knownSites, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def callVariantsUsingBcftools(script: TextIOWrapper, options: OptionsDict, interval: str, bqsr: str, vcf: str):
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
    # call variants
    if [[ ! -f {VCF} ]]; then
        logthis "Starting variant calling for {INTERVAL}"

        bcftools mpileup \\
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
            --max-depth 500 \\
            --no-BAQ \\
            --threads 4 \\
            --output-type u \\
            --regions {INTERVAL} \\
            --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
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
            REFERENCE=reference, ASSEMBLY=assembly, INTERVAL=interval, BQSR=bqsr, VCF=vcf
        )
    )


def annotate(
    script: TextIOWrapper,
    options: OptionsDict,
    vep: str,
    input: str,
    output: str,
    summary: str,
):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    stats = options["stats"]

    script.write(
        """
if [[ ! -f {OUTPUT} || ! -f {STATS}/{SUMMARY} ]]; then
    logthis "Starting annotation"
    
    vep --dir {VEP} \\
        --cache \\
        --format vcf \\
        --vcf \\
        --compress_output gzip \\
        --merged \\
        --fork 64 \\
        --offline \\
        --use_given_ref \\
        --verbose \\
        --force_overwrite \\
        --symbol \\
        --fasta {REFERENCE}/{ASSEMBLY}.fna \\
        --input_file {INPUT} \\
        --output_file {OUTPUT} \\
        --stats_file {STATS}/{SUMMARY}

    logthis "Completed annotation"
else
    logthis "Annotations already completed, ${{green}}already completed${{reset}}"
fi
""".format(
            VEP=vep,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            INPUT=input,
            OUTPUT=output,
            SUMMARY=summary,
            STATS=stats,
        )
    )


def scatter(script: TextIOWrapper, options: OptionsDict, prefix: str, sorted: str):
    intervals = computeIntervals(options)

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])

        script.write(
            """
\n# Scatter interval {INTERVAL}        
if [[ ! -f {BAM} || ! -f {BAM}.bai ]]; then
    logthis "Creating interval {INTERVAL}"
    (
        samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM}
        samtools index -@ 4 {BAM}
    )&
else
    logthis "Splinter for {INTERVAL} has been computed, ${{green}}already completed${{reset}}"
fi

""".format(
                SORTED=sorted, INTERVAL=interval[0], BAM=bam
            )
        )

    script.write('logthis "${yellow}Waiting for scattering processes to complete${reset}"\n')
    script.write("wait\n")
    script.write('logthis "${green}Scattering completed${reset}"\n')


def runIntervals(script: TextIOWrapper, options: OptionsDict, prefix: str):
    #
    # for now commenting out individual region annotation as there's
    # really no downstream consumer of that data. instead we'll
    # just keep the individual region's calls for later merging and
    # bulk annotation (which provides a complete summary annotation)
    # working = options["working"]
    caller = options["caller"]
    skipBQSR = options["skipBQSR"]

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

        if skipBQSR == False:
            genBQSR(script, options, interval[0], bam, bqsr)
        else:
            bqsr = bam

        if caller == "gatk":
            callVariantsUsingGatk(script, options, interval[0], bqsr, vcf)
        elif caller == "bcftools":
            callVariantsUsingBcftools(script, options, interval[0], bqsr, vcf)
        else:
            print("Unexpected value {CALLER} given for the --caller option".format(CALLER=caller))
            quit(1)

        script.write(") &\n")
        script.write("interval_processing_pids+=($!)\n")
        script.write("\n")

    script.write('logthis "${yellow}Waiting on interval processing to complete${reset}"\n')
    script.write("wait\n")
    script.write('logthis "${green}Intervals processed${reset}"\n')


def generateConsensus(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
if [[ ! -f {PIPELINE}/{SAMPLE}.consensus.fasta ]]; then
    logthis "${{yellow}}Building consensus {PIPELINE}/{SAMPLE}.unannotated.vcf.gz${{reset}}"

    bcftools consensus \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
    | sed '/>/ s/$/ | {SAMPLE}/' >{PIPELINE}/{SAMPLE}.consensus.fasta
else
    logthis "Consensus fasta already generated for {PIPELINE}/{SAMPLE}.consensus.fasta, ${{green}}already completed${{reset}}"
fi
        """.format(
            REFERENCE=reference, ASSEMBLY=assembly, PIPELINE=pipeline, SAMPLE=sample
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

    # for now this is hard-coded, but we need to get this from the json rules
    /bin/ls -1 {PIPELINE}/{SAMPLE}.NC_*.vcf | sort -k1,1V >{PIPELINE}/{SAMPLE}.merge.list

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
    """.format(
            PIPELINE=pipeline, SAMPLE=sample
        )
    )


def doVariantQC(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    knownSites = options["knownSites"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    script.write(
        """
#
# RUN Variant QC process
# 
logthis "Starting Variant QC processes"

#
# this is commented out right now because CollectVariantCallingMetrics crashes while reading
# the variants after about 45 minutes. i'll ping my contact at the broad to check if it's a
# known issue and get it resolved
#
# if [[ ! -f {STATS}/{SAMPLE}.variant_calling_detail_metrics ]]; then
#     logthis "Collection variant calling metrics"

#     gatk CollectVariantCallingMetrics \\
#         --VERBOSITY ERROR \\
#         --REFERENCE_SEQUENCE {REFERENCE}/{ASSEMBLY}.fna \\
#         --SEQUENCE_DICTIONARY {REFERENCE}/{ASSEMBLY}.dict \\
#         --DBSNP {REFERENCE}/{KNOWN_SITES} \\
#         --INPUT {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
#         --OUTPUT {STATS}/{SAMPLE}.variant_calling_detail_metrics &
# else
#     logthis "Variant metrics already run, ${{green}}already completed${{reset}}"
# fi

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

echo "VCF QC processes complete"

#
# run some other stats on the vcf file
#
if [[ ! -d {STATS}/{SAMPLE}_bcfstats ]]; then
    bcftools stats --fasta-ref {REFERENCE}/{ASSEMBLY}.fna {PIPELINE}/{SAMPLE}.unannotated.vcf.gz > {STATS}/{SAMPLE}.chk
    plot-vcfstats --prefix {STATS}/{SAMPLE}_bcfstats {STATS}/{SAMPLE}.chk
else
    echo "bcftools stats and plots already run, ${{green}}skipping${{reset}}"
fi

""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            KNOWN_SITES=knownSites,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
        )
    )


def startAlignmentQC(script: TextIOWrapper, options: OptionsDict, sorted: str):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
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

    gatk CollectAlignmentSummaryMetrics --java-options '-Xmx4g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    logthis "Alignment metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    logthis "Starting GC bais metrics on {SAMPLE}"

    gatk CollectGcBiasMetrics --java-options '-Xmx4g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    logthis "GC bias metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then
    logthis "Starting WGS metrics on {SAMPLE}"

    gatk CollectWgsMetrics --java-options '-Xmx4g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.wgs_metrics.txt \\
        --MINIMUM_BASE_QUALITY 20 \\
        --MINIMUM_MAPPING_QUALITY 20 \\
        --COVERAGE_CAP 250 \\
        --READ_LENGTH 151 \\
        --INTERVALS {REFERENCE}/{ASSEMBLY}_autosomal.interval_list \\
        --USE_FAST_ALGORITHM \\
        --INCLUDE_BQ_HISTOGRAM &
else
    logthis "WGS metrics already run, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then
    logthis "Starting samtools stats on {SAMPLE}"

    ( 
        samtools stats -@ 8 \\
            -r {REFERENCE}/{ASSEMBLY}.fna \\
                {SORTED} >{STATS}/{SAMPLE}.samstats

        plot-bamstats --prefix {SAMPLE}_samstats/ {STATS}/{SAMPLE}.samstats
    ) &
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

# 
# we'll bring this back in when it makes sense. for human sequences this
# generates a file that's enormous and likely can't be processed by anything
# anyway, so..
#
# if [[ ! -f {STATS}/{SAMPLE}.bedtools.coverage.gz ]]; then
#     bedtools genomecov -d -ibam {SORTED} \\
#         | gzip >{STATS}/{SAMPLE}.bedtools.coverage.gz &
# else
#     echo "bedtools genomecov already run, ${{green}}skipping${{reset}}"
# fi

if [[ ! -f {STATS}/{SAMPLE}.samtools.coverage ]]; then
    samtools coverage -d 0 \\
        --reference {REFERENCE}/{ASSEMBLY}.fna \\
        {SORTED} >{STATS}/{SAMPLE}.samtools.coverage &
else
    echo "samtools coverage already run, ${{green}}skipping${{reset}}"
fi

# if [[ ! -f {STATS}/{SAMPLE}.sorted_fastqc.zip || ! -f {STATS}/{SAMPLE}.sorted_fastqc.html ]]; then
#     logthis "Starting fastqc for {SAMPLE}"

#     fastqc \\
#         --outdir {STATS} \\
#         --noextract \\
#         {SORTED} &
# else
#     logthis "FASTQC already run, ${{green}}already completed${{reset}}"
# fi
""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
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
    pipeline = options["pipeline"]
    sample = options["sample"]
    intervals = computeIntervals(options)

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all intermediate interval files\n")
    script.write("#\n")

    for interval in intervals:
        script.write(
            """
rm -f {PIPELINE}/{SAMPLE}.{INTERVAL}*
rm -f {PIPELINE}/{SAMPLE}.merge.list
""".format(INTERVAL=interval[1], PIPELINE=pipeline, SAMPLE=sample)
        )

    script.write("\n")


def writeHeader(script: TextIOWrapper, options: OptionsDict, filenames: FastqSet):
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
    for interval in intervals:
        script.write(
            "#   {CHROME} => {ACCESSION} {SIZE}\n".format(
                CHROME=interval["chromosome"], ACCESSION=interval["accession"], SIZE=interval["length"]
            )
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
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.6.1:{WORKING}/bin:$PATH\n""".format(
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
        "-Q",
        "--skip-multi-qc",
        action="store_false",
        dest="doMultiQc",
        default=True,
        help="Skip running multiqc QC process on input and output files",
    )
    parser.add_argument(
        "--skip-variant-qc",
        action="store_false",
        dest="doVariantQc",
        default=True,
        help="Skip running variant QC process on input and output files",
    )
    parser.add_argument(
        "--skip-alignment-qc",
        action="store_false",
        dest="doAlignmentQc",
        default=True,
        help="Skip alignment QC process on input and output files",
    )
    parser.add_argument(
        "--skip-bqsr",
        action="store_true",
        dest="skipBQSR",
        default=False,
        help="Skip running BQSR processing on input file(s)",
    )
    parser.add_argument(
        "--skip-interval-processing",
        action="store_true",
        dest="skipIntervalProcessing",
        default=False,
        help="Skip all interval processing (scatter, call, gather)",
    )

    parser.add_argument(
        "-P",
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["trimmomatic", "fastp", "none"],
        help="Optionally run a FASTQ preprocessor",
    )

    parser.add_argument(
        "-A",
        "--aligner",
        action="store",
        dest="aligner",
        default="bwa",
        choices=["bwa", "hisat2"],
        help="Use 'bwa' or 'hisat2' as the aligner.",
    )

    parser.add_argument(
        "-S",
        "--sorter",
        action="store",
        dest="sorter",
        default="biobambam",
        choices=["biobambam", "samtools"],
        help="Use 'biobambam' or 'samtools' as the sorter.",
    )

    parser.add_argument(
        "-V",
        "--caller",
        action="store",
        dest="caller",
        default="bcftools",
        choices=["bcftools", "gatk"],
        help="Use `bcftools` or `gatk` as the variant caller",
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
        "-u",
        "--unmapped",
        action="store_true",
        dest="processUnmapped",
        default=False,
        help="Extract unmapped reads into secondary _R1 and _R2 FASTQ files",
    )

    parser.add_argument(
        "-X",
        "--align-only",
        action="store_true",
        dest="alignOnly",
        default=False,
        help="Only run alignment and sorting processes",
    )

    parser.add_argument(
        "-R",
        "--reference-assembly",
        required=True,
        action="store",
        metavar="REFERENCE_ASSEMBLY",
        dest="referenceAssembly",
        help="Base name of the reference assembly",
    )

    parser.add_argument(
        "-k",
        "--known-sites",
        required=True,
        action="store",
        metavar="KNOWN_SITES",
        dest="knownSites",
        help="Name of the 'known sites' VCF",
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
        default="hg38-chromosomeSizes.json",
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

    skipIntervalProcessing = options["skipIntervalProcessing"]

    if skipIntervalProcessing == False:
        writeMergeList(options)

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

        alignAndSort(script, options, sorted)

        # the order of the next two operations is important in this script
        # because scatter does its job but waits on all the background
        # operations to complete. this doesn't take long, so waiting is
        # ok.
        #
        # however, if we were to start the alignment qc process, which runs
        # in the background, we would start them and then block until they
        # completed while we're waiting to scatter intervals around
        #

        if skipIntervalProcessing == False:
            # this makes scatter a standalone process for now
            scatter(script, options, prefix, sorted)

        # if we've been asked to run qc at all we can actually start
        # the alignment qc processes very early (right after we've
        # completed alignment and sorting). so, we start them here. there's
        # an overall impact on the interval runs while these qc processes
        # are running, but we're willing to deal with that since the overall
        # qc process takes a very long time
        #
        # in particular it's the fastqc process itself which chews up
        # so much time. it's a single-threaded limited-memory (250mb)
        # process so we're not too worried about it's impact on the rest of
        # the variant calling processes
        if options["doAlignmentQc"] == True:
            startAlignmentQC(script, options, sorted)

        if skipIntervalProcessing == False:
            # process the intervals
            runIntervals(script, options, prefix)
            gather(script, options)

            # and cleanup the interval files so we have some space
            if options["cleanIntermediateFiles"] == True:
                cleanup(script, cleantarget, options)


        # we can start variant qc here because we're going to run against
        # the unannotated vcf (we'll get the same metrics either way), and
        # parts of this process take a long time
        if options["doVariantQc"] == True:
            doVariantQC(script, options)

        annotate(
            script,
            options,
            "{WORKING}/vep_data".format(WORKING=options["working"]),
            "{PIPELINE}/{SAMPLE}.unannotated.vcf.gz".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"]),
            "{PIPELINE}/{SAMPLE}.annotated.vcf.gz".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"]),
            "{SAMPLE}.annotated.vcf_summary.html".format(SAMPLE=options["sample"]),
        )

        generateConsensus(script, options)

        # we'll wait here to make sure all the background stuff is done before we
        # run multiqc and cleanup
        script.write('logthis "${yellow}Waiting for background processes to complete.${reset}"\n')
        script.write("wait\n")
        script.write('logthis "${green}Done processing${reset}"\n')

        if options["doMultiQc"] == True:
            runMultiQC(script, options)

        script.write('logthis "${yellow}Waiting for any outstanding processes to complete.${reset}"\n')
        script.write("wait\n")
        script.write('logthis "${green}Done processing${reset}"\n')

        script.write("\n")
        script.write("touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

    system("chmod +x {SCRIPT}".format(SCRIPT=options["script"]))


if __name__ == "__main__":
    main()
