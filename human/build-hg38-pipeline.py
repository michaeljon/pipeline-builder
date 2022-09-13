#!/usr/bin/env python

from io import TextIOWrapper
import json
import argparse
import math

from argparse import Namespace
from optparse import Option
from typing import Dict, Tuple, Sequence, List, Any
from datetime import datetime
from math import ceil
from os.path import exists, expandvars
from os import cpu_count, system

OptionsDict = Dict[str, Any]
FastqSet = Tuple[str, str]


def writeFragmentList(options: OptionsDict):
    fragmentCount = int(options["fragmentCount"])

    if fragmentCount == 1:
        return

    pipeline = options["pipeline"]
    sample = options["sample"]

    fragmentList = "{PIPELINE}/{SAMPLE}.fragments.tsv".format(PIPELINE=pipeline, SAMPLE=sample)
    with open(fragmentList, "w") as fl:
        fl.truncate()
        fl.write("fragment\tr1\tr2\taligned\tunmarked\tsorted\n")

        for fragment in range(1, fragmentCount + 1):
            fl.write(
                "{FRAGMENT:03n}\t{SAMPLE}_R1.trimmed_{FRAGMENT:03n}.fastq.gz\t{SAMPLE}_R2.trimmed_{FRAGMENT:03n}.fastq.gz\t{SAMPLE}_001_{FRAGMENT:03n}.aligned.bam\t{SAMPLE}_001_{FRAGMENT:03n}.unmarked.bam\t{SAMPLE}_001_{FRAGMENT:03n}.sorted.bam\n".format(
                    FRAGMENT=fragment, SAMPLE=sample
                )
            )


fallback_warning_shown = False


def getFileNames(options: OptionsDict) -> FastqSet:
    global fallback_warning_shown

    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    # assume we have the _001 pattern first
    filenames = (
        "{FASTQ_DIR}/{SAMPLE}_R1_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )

    if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
        if fallback_warning_shown == False:
            print("Falling back to shortened fastq file names")
            fallback_warning_shown = True

        # if that didn't work, try for the redacted names
        filenames = (
            "{FASTQ_DIR}/{SAMPLE}_R1.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
            "{FASTQ_DIR}/{SAMPLE}_R2.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        )

        if exists(expandvars(filenames[0])) == False or exists(expandvars(filenames[1])) == False:
            print(
                "Unable to locate the R1 or R2 files at {R1} and {R2}".format(
                    R1=filenames[0],
                    R2=filenames[1],
                )
            )
            print("Check your --sample and --fastq-dir parameters")
            quit(1)

    return filenames


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
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

    logthis "${{yellow}}Identity preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
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
    readLimit = options["read-limit"]

    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running FASTP preprocessor${{reset}}"

    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
    fastp \\
        --report_title "fastp report for sample {SAMPLE}" \\
        --in1 {R1} \\
        --in2 {R2} \\
        --out1 {O1} \\
        --out2 {O2} \\
        --detect_adapter_for_pe \\
        --verbose {LIMITREADS} \\
        --thread 16 \\
        -j {STATS}/{SAMPLE}-fastp.json \\
        -h {STATS}/{SAMPLE}-fastp.html

    logthis "${{yellow}}FASTP preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
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
    aligner = options["aligner"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    {ALIGNER} mem -t {THREADS} \\
        -Y -M \\
        -v 1 \\
        -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        {REFERENCE}/{ASSEMBLY}.fna \\
        {O1} \\
        {O2} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            ALIGNER=aligner,
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
        )
    )


def alignFASTQ(
    script: TextIOWrapper,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    aligner = options["aligner"]

    if aligner == "bwa" or aligner == "bwa-mem2":
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

    logthis "${{yellow}}Sorting and marking duplicates completed${{reset}}"
else
    logthis "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
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


def sortWithSamtools(script: TextIOWrapper, options: OptionsDict, output: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
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

    samtools sort -@ {THREADS} {PIPELINE}/{SAMPLE}.aligned.bam -o {UNMARKED} --verbosity 5

    logthis "${{yellow}}Sorting aligned file completed${{reset}}"
else
    logthis "{UNMARKED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    logthis "${{yellow}}Marking duplicates${{reset}}"

    java -Xmx8g -jar {BIN}/picard.jar MarkDuplicates \\
        --TAGGING_POLICY All \\
        --REFERENCE_SEQUENCE {REFERENCE}/{ASSEMBLY}.fna \\
        -I {UNMARKED} \\
        -O {SORTED} \\
        -M {STATS}/{SAMPLE}_marked_dup_metrics.txt    

    # generate an index on the result
    samtools index -b {SORTED} {SORTED}.bai

    logthis "${{yellow}}Marking duplicates completed${{reset}}"
else
    logthis "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi
    """.format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
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


def fragment(script: TextIOWrapper, r1: str, r2: str, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    fragmentCount = options["fragmentCount"]

    script.write(
        """
logthis "${{yellow}}Fragmenting trimmed FASTQ${{reset}}"

Ovation.Pipeline.FastqProcessor split \\
    --format FastqGz \\
    --splits {FRAGMENTS} \\
    --output {PIPELINE} \\
    --in1 {R1} \\
    --in2 {R2}

logthis "${{yellow}}Fragmenting trimmed FASTQ completed${{reset}}"
""".format(
            PIPELINE=pipeline, R1=r1, R2=r2, SAMPLE=sample, FRAGMENTS=fragmentCount
        )
    )


def alignFragments(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]

    cores = options["cores"]
    threads = 4
    jobs = int(cores / threads)

    fragmentList = "{PIPELINE}/{SAMPLE}.fragments.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
#
# align the input files
#
logthis "${{yellow}}Running aligners${{reset}}"
parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.alignment.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
        bwa-mem2 mem -t {THREADS} \\
            -Y -M \\
            -v 1 \\
            -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
            {REFERENCE}/{ASSEMBLY}.fna \\
            {PIPELINE}/{{r1}} \\
            {PIPELINE}/{{r2}} | \\
        samtools view -Sb -@ 4 - >{PIPELINE}/{{aligned}}
    fi' :::: {FRAGMENT_LIST}
logthis "${{yellow}}Fragment alignment complete${{reset}}"
""".format(
            JOBS=jobs,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            FRAGMENT_LIST=fragmentList,
        )
    )


def sortParallelWithBiobambam(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    cores = options["cores"]
    threads = 4
    jobs = int(cores / threads)

    fragmentList = "{PIPELINE}/{SAMPLE}.fragments.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
#
# sort and mark duplicates
#
logthis "${{yellow}}Fragment sorting and marking duplicates${{reset}}"
parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.sort.log --header --colsep $'\\t' \\
    'if [[ ! -f {{sorted}} ]]; then
        bamsormadup \\
            SO=coordinate \\
            threads={THREADS} \\
            level=6 \\
            tmpfile={TEMP}/{SAMPLE} \\
            inputformat=bam \\
            indexfilename={{sorted}}.bai \\
            M={STATS}/{SAMPLE}_{{fragment}}.duplication_metrics <{PIPELINE}/{{aligned}} >{PIPELINE}/{{sorted}}
    fi' :::: {FRAGMENT_LIST}

# force the index to look "newer" than its source
parallel --header --colsep $'\\t' \\
    touch {{sorted}}.bai

logthis "${{yellow}}Fragment sorting and marking duplicates completed${{reset}}"
""".format(
            JOBS=jobs,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
            FRAGMENT_LIST=fragmentList,
        )
    )


def sortParallelWithSamtools(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    fragmentList = "{PIPELINE}/{SAMPLE}.fragments.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
#
# sort and mark duplicates
#
logthis "${{yellow}}Sorting aligned file${{reset}}"
parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.sort.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{{unmarked}} ]]; then
        samtools sort {PIPELINE}/{{aligned}} -o {PIPELINE}/{{unmarked}}
     fi' :::: {FRAGMENT_LIST}

logthis "${{yellow}}Marking duplicates${{reset}}"
parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.mark.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{{sorted}} ]]; then
        java -Xmx8g -jar {BIN}/picard.jar MarkDuplicates \\
            --TAGGING_POLICY All \\
            --REFERENCE_SEQUENCE {REFERENCE}/{ASSEMBLY}.fna \\
            -I {PIPELINE}/{{unmarked}} \\
            -O {PIPELINE}/{{sorted}} \\
            -M {STATS}/{SAMPLE}_{{fragment}}_marked_dup_metrics.txt    
     fi' :::: {FRAGMENT_LIST}

# generate an index on the result
parallel --header --colsep $'\\t' \\
    samtools index -b {PIPELINE}/{{sorted}} {PIPELINE}/{{sorted}}.bai

logthis "${{yellow}}Sorting and marking duplicates compeleted${{reset}}"
    """.format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
            FRAGMENT_LIST=fragmentList,
        )
    )


def sortFragments(script: TextIOWrapper, options: OptionsDict):
    sorter = options["sorter"]

    if sorter == "biobambam":
        sortParallelWithBiobambam(script, options)
    else:
        sortParallelWithSamtools(script, options)


def combine(script: TextIOWrapper, options: OptionsDict, output: str):
    # samtools merge -f -o {OUTPUT} -c -p <*.bam>
    pass


def extractUmappedReads(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# extract unmapped reads
#
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped_R1.fastq || ! -f {PIPELINE}/{SAMPLE}_unmapped_R2.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}_unmapped_R1.fastq \\
        -2 {PIPELINE}/{SAMPLE}_unmapped_R2.fastq \\
        {PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Unmapped read extraction completed${{reset}}"
else
    logthis "Unmapped reads extracted to initial FASTQ, ${{green}}skipping${{reset}}"
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
    fragmentCount = options["fragmentCount"]

    bin = options["bin"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    segmentSize = options["segmentSize"]
    factor = options["factor"]

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filenames[0], filenames[1], trimmedFilenames[0], trimmedFilenames[1], options)

    if fragmentCount == 1:
        alignFASTQ(script, trimmedFilenames[0], trimmedFilenames[1], options)
        sortAlignedAndMappedData(script, options, output)
    else:
        fragment(script, trimmedFilenames[0], trimmedFilenames[1], options)
        alignFragments(script, options)
        sortFragments(script, options)
        combine(script, options, output)

    script.write(
        """
# Write intervals file computed from aligned and sorted BAM
samtools view --header-only {BAM} | \\
    grep --color=never '^@SQ' | \\
    python {BIN}/interval-builder.py --segment {SEGMENT} --factor {FACTOR} --process intervalList > {PIPELINE}/{SAMPLE}.intervals.tsv

# Write VCF merge list file computed from aligned and sorted BAM
samtools view --header-only {BAM} | \\
    grep --color=never '^@SQ' | \\
    python {BIN}/interval-builder.py --segment {SEGMENT} --factor {FACTOR} --process mergeList --root {PIPELINE}/{SAMPLE} > {PIPELINE}/{SAMPLE}.merge.list
""".format(
            BIN=bin, BAM=output, PIPELINE=pipeline, SAMPLE=sample, SEGMENT=segmentSize, FACTOR=factor
        )
    )

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


def genBQSRTables(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    reference = options["reference"]
    assembly = options["referenceAssembly"]
    knownSites = options["knownSites"]

    intervalList = "{PIPELINE}/{SAMPLE}.intervals.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
logthis "${{yellow}}Generating BQSR tables${{reset}}"
parallel --joblog {PIPELINE}/{SAMPLE}.bqsr.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}_bqsr.table ]]; then
        gatk BaseRecalibrator --java-options -Xmx4g \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
            -I {PIPELINE}/{SAMPLE}.{{root}}.bam \\
            -O {PIPELINE}/{SAMPLE}.{{root}}_bqsr.table \\
            --verbosity ERROR \\
            --preserve-qscores-less-than 6 \\
            --known-sites {REFERENCE}/{KNOWN_SITES} \\
            -L {{interval}}
     fi' :::: {INTERVAL_LIST}
logthis "${{green}}BQSR table generation complete${{reset}}"
""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            INTERVAL_LIST=intervalList,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            KNOWN_SITES=knownSites,
        )
    )


def applyBQSRTables(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    reference = options["reference"]
    assembly = options["referenceAssembly"]

    intervalList = "{PIPELINE}/{SAMPLE}.intervals.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
logthis "${{yellow}}Applying BQSR calibration${{reset}}"
parallel --joblog {PIPELINE}/{SAMPLE}.calibrate.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}_bqsr.bam ]]; then
        gatk ApplyBQSR --java-options -Xmx4g \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
            -I {PIPELINE}/{SAMPLE}.{{root}}.bam \\
            -O {PIPELINE}/{SAMPLE}.{{root}}_bqsr.bam \\
            --verbosity ERROR \\
            --emit-original-quals true \\
            --preserve-qscores-less-than 6 \\
            --static-quantized-quals 10 \\
            --static-quantized-quals 20 \\
            --static-quantized-quals 30 \\
            --bqsr-recal-file {PIPELINE}/{SAMPLE}.{{root}}_bqsr.table \\
            -L {{interval}}
     fi' :::: {INTERVAL_LIST}
logthis "${{green}}BQSR calibration completed${{reset}}"
""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            INTERVAL_LIST=intervalList,
            REFERENCE=reference,
            ASSEMBLY=assembly,
        )
    )


def genBQSR(script: TextIOWrapper, options: OptionsDict):
    genBQSRTables(script, options)
    applyBQSRTables(script, options)


def callVariantsUsingGatk(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]
    intervalList = "{PIPELINE}/{SAMPLE}.intervals.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    # haplotypecaller _can_ use threads, but doesn't do so in a really
    # good way. after a bunch of experimentation we found that limiting
    # it to a _single_ thread and running as many jobs as possible leads
    # to a much better load average and overall runtime

    # read limit 100,000,000
    # interval size 50,000,000
    # 1 hmm threads / 72 jobs
    # [2022-07-12 18:31:50] Calling variants using GATK
    # [2022-07-12 19:06:45] GATK variant calling completed
    # 2 hmm threads / 36 jobs
    # [2022-07-12 16:59:59] Calling variants using GATK
    # [2022-07-12 17:45:17] GATK variant calling completed
    # 4 hmm threads / 18 jobs
    # [2022-07-12 02:29:11] Calling variants using GATK
    # [2022-07-12 03:36:31] GATK variant calling completed
    # 6 hmm threads / 12 jobs
    # [2022-07-12 15:31:14] Calling variants using GATK
    # [2022-07-12 16:57:06] GATK variant calling completed
    # 8 hmm threads / 9 jobs
    # [2022-07-11 23:04:40] Calling variants using GATK
    # [2022-07-12 00:56:43] GATK variant calling completed
    cores = options["cores"]
    hmmThreads = 1
    jobs = int(cores / hmmThreads)

    reference = options["reference"]
    assembly = options["referenceAssembly"]
    knownSites = options["knownSites"]

    script.write(
        """
logthis "${{yellow}}Calling variants using GATK${{reset}}"

parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.call.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}.vcf ]]; then
        gatk HaplotypeCaller --java-options -Xmx4g \\
            -R {REFERENCE}/{ASSEMBLY}.fna \\
            -I {PIPELINE}/{SAMPLE}.{{root}}_bqsr.bam \\
            -O {PIPELINE}/{SAMPLE}.{{root}}.vcf \\
            --verbosity ERROR \\
            --dbsnp {REFERENCE}/{KNOWN_SITES} \\
            --pairHMM FASTEST_AVAILABLE \\
            --native-pair-hmm-threads {THREADS} \\
            -L {{interval}}
    fi' :::: {INTERVAL_LIST}

logthis "${{green}}GATK variant calling completed${{reset}}"
""".format(
            JOBS=jobs,
            THREADS=hmmThreads,
            PIPELINE=pipeline,
            SAMPLE=sample,
            INTERVAL_LIST=intervalList,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            KNOWN_SITES=knownSites,
        )
    )


def callVariantsUsingBcftools(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]
    intervalList = "{PIPELINE}/{SAMPLE}.intervals.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    reference = options["reference"]
    assembly = options["referenceAssembly"]

    # bcftools uses threads, so we limit the number of jobs, in this case
    # we'll assign the number of threads equally to each part of the pipeline
    # and limit the number of jobs to account for 2x those threads being in use
    cores = options["cores"]
    threads = 4
    jobs = int(cores / (threads * 2))

    script.write(
        """
logthis "Calling variants using BCFTOOLS"

parallel -j {JOBS} --joblog {PIPELINE}/{SAMPLE}.call.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}.vcf ]]; then
        bcftools mpileup \\
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
            --max-depth 500 \\
            --no-BAQ \\
            --threads {THREADS} \\
            --output-type u \\
            --regions {{interval}} \\
            --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
            {PIPELINE}/{SAMPLE}.{{root}}_bqsr.bam 2>/dev/null | \\
        bcftools call \\
            --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
            --variants-only \\
            --multiallelic-caller \\
            --ploidy GRCh38 \\
            --threads {THREADS} \\
            --output-type v  \\
            --output {PIPELINE}/{SAMPLE}.{{root}}.vcf 2>/dev/null
    fi' :::: {INTERVAL_LIST}

logthis "BCFTOOLS variant calling completed"
""".format(
            JOBS=jobs,
            THREADS=threads,
            PIPELINE=pipeline,
            SAMPLE=sample,
            INTERVAL_LIST=intervalList,
            REFERENCE=reference,
            ASSEMBLY=assembly,
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

    # vep is an interesting beast when it comes to threading. it spends a ton
    # of time forking and waiting on processes that are then waiting on I/O
    # so, we grab the number of cores and double it, then assign that to vep
    cores = options["cores"]
    forks = int(cores * 2)

    script.write(
        """
if [[ ! -f {OUTPUT} || ! -f {STATS}/{SUMMARY} ]]; then
    logthis "Starting annotation"
    
    vep --dir {VEP} \\
        --cache \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --merged \\
        --fork {FORKS} \\
        --offline \\
        --use_given_ref \\
        --verbose \\
        --force_overwrite \\
        --symbol \\
        --hgvs \\
        --protein \\
        --af \\
        --fasta {REFERENCE}/{ASSEMBLY}.fna \\
        --input_file {INPUT} \\
        --output_file {OUTPUT} \\
        --stats_file {STATS}/{SUMMARY}

    tabix -p vcf {OUTPUT}

    logthis "Completed annotation"
else
    logthis "Annotations already completed, ${{green}}already completed${{reset}}"
fi
""".format(
            FORKS=forks,
            VEP=vep,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            INPUT=input,
            OUTPUT=output,
            SUMMARY=summary,
            STATS=stats,
        )
    )


def scatter(script: TextIOWrapper, options: OptionsDict, sorted: str):
    pipeline = options["pipeline"]
    sample = options["sample"]
    intervalList = "{PIPELINE}/{SAMPLE}.intervals.tsv".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
logthis "${{yellow}}Waiting for scattering processes to complete${{reset}}"

parallel --joblog {PIPELINE}/{SAMPLE}.scatter.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}.bam ]]; then
        samtools view -bh {SORTED} {{interval}} --output {PIPELINE}/{SAMPLE}.{{root}}.bam
     fi' :::: {INTERVAL_LIST}

parallel --joblog {PIPELINE}/{SAMPLE}.scatindex.log --header --colsep $'\\t' \\
    'if [[ ! -f {PIPELINE}/{SAMPLE}.{{root}}.bam.bai ]]; then
        samtools index {PIPELINE}/{SAMPLE}.{{root}}.bam
     fi' :::: {INTERVAL_LIST}

logthis "${{green}}Scattering completed${{reset}}"
""".format(
            PIPELINE=pipeline, SAMPLE=sample, SORTED=sorted, INTERVAL_LIST=intervalList
        )
    )


def runIntervals(script: TextIOWrapper, options: OptionsDict, prefix: str):
    caller = options["caller"]

    script.write('logthis "${yellow}Processing intervals${reset}"\n')

    genBQSR(script, options)

    if caller == "gatk":
        callVariantsUsingGatk(script, options)
    elif caller == "bcftools":
        callVariantsUsingBcftools(script, options)
    else:
        print("Unexpected value {CALLER} given for the --caller option".format(CALLER=caller))
        quit(1)

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
        --sample {SAMPLE} \\
    | sed '/>/ s/$/ | {SAMPLE}/' >{PIPELINE}/{SAMPLE}.consensus.fasta

    logthis "${{yellow}}Consensus completed${{reset}}"
else
    logthis "Consensus fasta already generated for {PIPELINE}/{SAMPLE}.consensus.fasta, ${{green}}already completed${{reset}}"
fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, PIPELINE=pipeline, SAMPLE=sample
        )
    )


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

    checks = []

    # this is commented out right now because CollectVariantCallingMetrics crashes while reading
    # the variants after about 45 minutes. i'll ping my contact at The Broad to check if it's a
    # known issue and get it resolved
    # checks.append(
    #     """'if [[ ! -f {STATS}/{SAMPLE}.variant_calling_detail_metrics ]]; then gatk CollectVariantCallingMetrics --VERBOSITY ERROR --REFERENCE_SEQUENCE {REFERENCE}/{ASSEMBLY}.fna --SEQUENCE_DICTIONARY {REFERENCE}/{ASSEMBLY}.dict --DBSNP {REFERENCE}/{KNOWN_SITES} --INPUT {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --OUTPUT {STATS}/{SAMPLE}.variant_calling_detail_metrics; fi' \\\n""".format(
    #         REFERENCE=reference,
    #         ASSEMBLY=assembly,
    #         KNOWN_SITES=knownSites,
    #         PIPELINE=pipeline,
    #         SAMPLE=sample,
    #         STATS=stats,
    #     )
    # )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.frq ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --freq2 --out {STATS}/{SAMPLE} --max-alleles 2 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.idepth ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --depth --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.ldepth.mean ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --site-mean-depth --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.lqual ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --site-quality --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.imiss ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --missing-indv --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.lmiss ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --missing-site --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.het ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz --het --out {STATS}/{SAMPLE} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )

    checks.append(
        """'if [[ ! -d {STATS}/{SAMPLE}_bcfstats.stats.txt ]]; then bcftools stats --fasta-ref {REFERENCE}/{ASSEMBLY}.fna {PIPELINE}/{SAMPLE}.unannotated.vcf.gz > {STATS}/{SAMPLE}_bcfstats.stats.txt; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            KNOWN_SITES=knownSites,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
        )
    )

    return checks


def doAlignmentQC(script: TextIOWrapper, options: OptionsDict, sorted: str):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]
    skipFastQc = options["skipFastQc"]

    filenames = getTrimmedFileNames(options)

    checks = []

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.flagstat.txt ]]; then samtools flagstat -@ 8 {SORTED} >{STATS}/{SAMPLE}.flagstat.txt; fi' \\\n""".format(
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.alignment_metrics.txt ]]; then gatk CollectAlignmentSummaryMetrics --java-options -Xmx4g --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.alignment_metrics.txt; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then gatk CollectGcBiasMetrics --java-options -Xmx4g --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.gc_bias_metrics.txt -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf -S {STATS}/{SAMPLE}.gc_bias_summary.txt; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then gatk CollectWgsMetrics --java-options -Xmx4g --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.wgs_metrics.txt --MINIMUM_BASE_QUALITY 20 --MINIMUM_MAPPING_QUALITY 20 --COVERAGE_CAP 250 --READ_LENGTH 151 --INTERVALS {REFERENCE}/{ASSEMBLY}_autosomal.interval_list --USE_FAST_ALGORITHM --INCLUDE_BQ_HISTOGRAM; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.stats.txt ]]; then samtools stats -@ 8 -r {REFERENCE}/{ASSEMBLY}.fna {SORTED} >{STATS}/{SAMPLE}.stats.txt; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.idxstats.txt ]]; then samtools idxstats -@ 8 {SORTED} >{STATS}/{SAMPLE}.idxstats.txt; fi' \\\n""".format(
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.samtools.coverage ]]; then samtools coverage -d 0 --reference {REFERENCE}/{ASSEMBLY}.fna {SORTED} >{STATS}/{SAMPLE}.samtools.coverage; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    if skipFastQc == False:
        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.html || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.html ]]; then fastqc --threads 2 --outdir {STATS} --noextract {R1} {R2}; fi' \\\n""".format(
                SAMPLE=sample,
                STATS=stats,
                R1=filenames[0],
                R2=filenames[1],
            )
        )

    return checks


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
multiqc \\
    --verbose \\
    --force \\
    --cl_config 'custom_logo: "{STATS}/ovationlogo.png"' \\
    --cl_config 'custom_logo_url: "https://www.ovation.io"' \\
    --cl_config 'custom_logo_title: "Ovation"' \\
    {STATS}

# Save the output
mv {STATS}/qc/multiqc_data {STATS}/{SAMPLE}_multiqc_data
mv {STATS}/qc/multiqc_report.html {STATS}/{SAMPLE}_multiqc_report.html
rm -rf {STATS}/qc

logthis "MultiQC for {SAMPLE} is complete"
""".format(
            STATS=stats, SAMPLE=sample
        )
    )


def doQualityControl(script: TextIOWrapper, options: OptionsDict, sorted: str):
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    plot_vcfstats = """plot-vcfstats --prefix {STATS}/{SAMPLE}_bcfstats {STATS}/{SAMPLE}_bcfstats.stats.txt""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    # this doesn't have a test, it's fast enough that we can afford to run it
    plot_bamstats = """plot-bamstats --prefix {STATS}/{SAMPLE}_samstats/ {STATS}/{SAMPLE}.stats.txt""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    alignment_checks = doAlignmentQC(script, options, sorted)
    variant_checks = doVariantQC(script, options)

    cmd = ""
    if options["doAlignmentQc"] == True or options["doVariantQc"] == True:
        cmd = "parallel --joblog {PIPELINE}/{SAMPLE}.qc.log ::: \\\n".format(PIPELINE=pipeline, SAMPLE=sample)

        if options["doAlignmentQc"] == True:
            for qc in alignment_checks:
                cmd += "    " + qc

        if options["doVariantQc"] == True:
            for qc in variant_checks:
                cmd += "    " + qc

        cmd += "    'echo End of stats'\n"

        script.write(
            """
#
# RUN QC processes, if there are any
# 
logthis "Starting QC processes"
{PARALLEL}

# can't run plot-bamstats until the statistics files are run above
{PLOT_BAMSTATS}

# can't run plot-vcfstats until the statistics files are run above
{PLOT_VCFSTATS}
""".format(
                PARALLEL=cmd,
                PLOT_BAMSTATS=plot_bamstats if options["doAlignmentQc"] == True else "",
                PLOT_VCFSTATS=plot_vcfstats if options["doVariantQc"] == True else "",
            )
        )

    if options["doMultiQc"] == True:
        runMultiQC(script, options)


def cleanup(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write("\n")
    script.write("#\n")
    script.write("# Clean up all intermediate interval files\n")
    script.write("#\n")

    script.write(
        """
for i in $(cut -d $'\\t' -f 2 zr5654_9_S5.intervals.tsv); do        
    rm -f {PIPELINE}/{SAMPLE}.${{i}}*
done
""".format(
            PIPELINE=pipeline, SAMPLE=sample
        )
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

    segmentSize = options["segmentSize"]
    factor = options["factor"]

    script.write("#\n")
    script.write("# Split parameters\n")
    script.write("#   segmentSize = {P}\n".format(P=segmentSize))
    script.write("#   factor = {P}\n".format(P=factor))
    script.write("#   last block max size = {P}\n".format(P=math.floor(segmentSize * factor)))

    fragmentCount = options["fragmentCount"]

    script.write("#\n")
    script.write("# Fragments\n")
    script.write("#   fragmentCount = {P}\n".format(P=fragmentCount))


def writeVersions(script: TextIOWrapper):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...\n")
    script.write("#\n")


def writeEnvironment(script: TextIOWrapper, options: OptionsDict):
    noColor = options["noColor"]

    script.write("#\n")

    if noColor == True:
        script.write(
            """
# color shortcuts disabled
red=$(echo -n "")
green=$(echo -n "")
yellow=$(echo -n "")
reset=$(echo -n "")
"""
        )
    else:
        script.write(
            """
# color shortcuts enabled
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
reset=$(tput sgr0)
"""
        )

    script.write(
        """
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
export LD_LIBRARY_PATH={WORKING}/lib:{WORKING}/bin:/usr/lib64:/usr/local/lib/:$LB_LIBRARY_PATH
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
        "--skip-fastqc",
        action="store_true",
        dest="skipFastQc",
        default=False,
        help="Skip running FASTQC statistics",
    )

    parser.add_argument(
        "--skip-interval-processing",
        action="store_true",
        dest="skipIntervalProcessing",
        default=False,
        help="Skip all interval processing (scatter, call, gather)",
    )

    parser.add_argument(
        "--skip-annotation",
        action="store_true",
        dest="skipAnnotation",
        default=False,
        help="Skip all VCF annotation processes",
    )

    parser.add_argument(
        "--generate-consensus",
        action="store_true",
        dest="generateConsensus",
        default=False,
        help="Generate a consensus FASTA file",
    )

    parser.add_argument(
        "--no-color",
        action="store_true",
        dest="noColor",
        default=False,
        help="Turn off colorized log output",
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
        default="bwa-mem2",
        choices=["bwa", "bwa-mem2"],
        help="Use 'bwa' or 'bwa-mem2' as the aligner.",
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
        "--vep-dir",
        action="store",
        metavar="VEP_DATA",
        dest="vep_data_dir",
        help="Location of the current VEP data",
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
        type=int,
        action="store",
        dest="cores",
        metavar="CPU_COUNT",
        default=cpu_count(),
        help="Specify the number of available CPU",
    )

    parser.add_argument(
        "--read-limit",
        type=int,
        action="store",
        dest="read-limit",
        metavar="READ_LIMIT",
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    parser.add_argument(
        "-d",
        "--watchdog",
        type=int,
        action="store",
        metavar="WATCHDOG_TIMEOUT",
        dest="watchdog",
        default=150,
        help="Specify a watchdog timeout for the alignment process. Value is in minutes",
    )

    parser.add_argument(
        "--fragments",
        type=int,
        action="store",
        metavar="FRAGMENT_COUNT",
        dest="fragmentCount",
        default=1,
        help="Number of FASTQ fragments to align",
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
        type=int,
        action="store",
        metavar="SEGMENT_SIZE",
        dest="segmentSize",
        default=50_000_000,
        help="Size of interval partition (ADVANCED)",
    )
    parser.add_argument(
        "--factor",
        type=float,
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

    writeFragmentList(options)

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        script.write("set -ep\n")
        writeHeader(script, options, filenames)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write("touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

        alignAndSort(script, options, sorted)

        if skipIntervalProcessing == False:
            # this makes scatter a standalone process for now
            scatter(script, options, sorted)
            runIntervals(script, options, prefix)
            gather(script, options)

            # and cleanup the interval files so we have some space
            if options["cleanIntermediateFiles"] == True:
                cleanup(script, options)

        if options["skipAnnotation"] == False:
            annotate(
                script,
                options,
                "{VEP_DATA}".format(VEP_DATA=options["vep_data_dir"]),
                "{PIPELINE}/{SAMPLE}.unannotated.vcf.gz".format(
                    PIPELINE=options["pipeline"], SAMPLE=options["sample"]
                ),
                "{PIPELINE}/{SAMPLE}.annotated.vcf.gz".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"]),
                "{SAMPLE}.annotated.vcf_summary.html".format(SAMPLE=options["sample"]),
            )

        if options["generateConsensus"] == True:
            generateConsensus(script, options)

        # we'll wait here to make sure all the background stuff is done before we
        # run multiqc and cleanup
        script.write('logthis "${green}Done with front-end processing${reset}"\n')

        doQualityControl(script, options, sorted)

        script.write('logthis "${green}Done with back-end processing${reset}"\n')

        script.write("\n")
        script.write("touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

    system("chmod +x {SCRIPT}".format(SCRIPT=options["script"]))


if __name__ == "__main__":
    main()
