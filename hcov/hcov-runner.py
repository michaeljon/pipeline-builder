#!/usr/bin/env python

from io import TextIOWrapper
import argparse
from argparse import Namespace
from typing import Tuple, Dict, Any
from datetime import datetime
from os.path import exists, expandvars
from os import cpu_count, system

FastqSet = Tuple[str, str]
OptionsDict = Dict[str, Any]


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    bin = options["bin"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    chromosome = options["referenceName"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    script.write(
        """
if [[ ! -f {REFERENCE}/{ASSEMBLY}.dict ]]; then
    logthis "${{yellow}}Creating sequence dictionary${{reset}}"

    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -O {REFERENCE}/{ASSEMBLY}.dict

    logthis "Sequence dictionary completed"
else
    logthis "Reference dictionary {REFERENCE}/{ASSEMBLY}.dict ${{green}}already completed${{reset}}"
fi

if [[ ! -f {REFERENCE}/{ASSEMBLY}_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    logthis "${{yellow}}Building {REFERENCE}/{ASSEMBLY}.interval_list${{reset}}"

    egrep '({CHROMOSOME})\\s' {REFERENCE}/{ASSEMBLY}.fna.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/{ASSEMBLY}.dict - >{REFERENCE}/{ASSEMBLY}_autosomal.interval_list

    logthis "Building {REFERENCE}/{ASSEMBLY}_autosomal.interval_list finished"
else
    logthis "Interval list {REFERENCE}/{ASSEMBLY}_autosomal.interval_list ${{green}}already completed${{reset}}"
fi

""".format(
            REFERENCE=reference, ASSEMBLY=assembly, BIN=bin, CHROMOSOME=chromosome
        )
    )


def sortWithBiobambam(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])

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
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
            SORTED=sorted,
        )
    )

    pass


def sortWithSamtools(script: TextIOWrapper, options: OptionsDict):
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

    samtools sort {PIPELINE}/{SAMPLE}.aligned.bam -o {UNMARKED}

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
        -O {PIPELINE}/{SAMPLE}.sorted.bam \\
        -M {STATS}/{SAMPLE}_marked_dup_metrics.txt    

    # generate an index on the result
    samtools index -b {PIPELINE}/{SAMPLE}.sorted.bam {PIPELINE}/{SAMPLE}.sorted.bam.bai

    logthis "${{yellow}}Marking duplicates completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.sorted.bam, index, and metrics found, ${{green}}skipping${{reset}}"
fi
    """.format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            UNMARKED=unmarked,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
        )
    )


def sortAlignedAndMappedData(script: TextIOWrapper, options: OptionsDict):
    sorter = options["sorter"]

    if sorter == "biobambam":
        sortWithBiobambam(script, options)
    else:
        sortWithSamtools(script, options)


def callVariantsUsingGatk(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using GATK${{reset}}"

    gatk HaplotypeCaller --java-options '-Xmx8g' \\
        -R {REFERENCE}/{ASSEMBLY}.fna \\
        -I {PIPELINE}/{SAMPLE}.sorted.bam \\
        -O {PIPELINE}/{SAMPLE}.unannotated.vcf \\
        --verbosity ERROR \\
        --pairHMM FASTEST_AVAILABLE \\
        --native-pair-hmm-threads 4

    bgzip {PIPELINE}/{SAMPLE}.unannotated.vcf
    tabix -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "${{green}}GATK variant calling completed${{reset}}"
else
    logthis "Variants already called for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, SAMPLE=sample
        )
    )


def callVariantsUsingLofreq(
    script: TextIOWrapper,
    options: OptionsDict,
):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using lofreq${{reset}}"

    lofreq indelqual \\
        --dindel \\
        --ref {REFERENCE}/{ASSEMBLY}.fna \\
        --verbose \\
        --out - \\
        {PIPELINE}/{SAMPLE}.sorted.bam |
    lofreq call \\
        --ref {REFERENCE}/{ASSEMBLY}.fna \\
        --call-indels \\
        --no-default-filter \\
        --max-depth 1000000 \\
        --force-overwrite \\
        --verbose \\
        --out {PIPELINE}/{SAMPLE}.unannotated.vcf \\
        -

    bgzip {PIPELINE}/{SAMPLE}.unannotated.vcf
    tabix -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "${{green}}lofreq variant calling completed${{reset}}"
else
    logthis "Variants already called for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, SAMPLE=sample, PIPELINE=pipeline
        )
    )


def callVariantsUsingBcftools(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using bcftools${{reset}}"

    if [[ ! -f {PIPELINE}/{SAMPLE}.ploidy ]]; then
        echo '* * * * 1' >{PIPELINE}/{SAMPLE}.ploidy
    fi

    bcftools mpileup \\
        --annotate FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SCR,FORMAT/SP,INFO/AD,INFO/SCR \\
        --max-depth 1000000 \\
        --max-idepth 1000000 \\
        --threads 4 \\
        --output-type u \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.sorted.bam 2>/dev/null | \\
    bcftools call \\
        --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
        --variants-only \\
        --keep-alts \\
        --multiallelic-caller \\
        --ploidy-file {PIPELINE}/{SAMPLE}.ploidy \\
        --prior 0.05 \\
        --threads 4 \\
        --output-type v  \\
        --output {PIPELINE}/{SAMPLE}.unannotated.vcf.tmp 2>/dev/null

    bcftools +fill-tags \\
        {PIPELINE}/{SAMPLE}.unannotated.vcf.tmp \\
        --output-type v  \\
        --output {PIPELINE}/{SAMPLE}.unannotated.vcf \\
        -- --tags AC,AN,AF,VAF,MAF 2>/dev/null

    rm -f {PIPELINE}/{SAMPLE}.unannotated.vcf.tmp

    bgzip {PIPELINE}/{SAMPLE}.unannotated.vcf
    tabix -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "${{green}}bcftools variant calling completed${{reset}}"
else
    logthis "Variants already called for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def callVariants(script: TextIOWrapper, options: OptionsDict):
    caller = options["caller"]

    if caller == "gatk":
        callVariantsUsingGatk(script, options)
    elif caller == "bcftools":
        callVariantsUsingBcftools(script, options)
    elif caller == "lofreq":
        callVariantsUsingLofreq(script, options)
    else:
        print("Unexpected value {CALLER} given for the --caller option".format(CALLER=caller))
        quit(1)


def produceConsensusUsingBcftools(
    script: TextIOWrapper,
    options: OptionsDict,
):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# produce consensus using bcftools
if [[ ! -f {PIPELINE}/{SAMPLE}.consensus.fa ]]; then
    logthis "${{yellow}}Building consensus {PIPELINE}/{SAMPLE}.unannotated.vcf.gz${{reset}}"

    bcftools index {PIPELINE}/{SAMPLE}.unannotated.vcf.gz
    bcftools consensus \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
    | sed '/>/ s/$/ | {SAMPLE}/' >{PIPELINE}/{SAMPLE}.consensus.fa

    logthis "${{yellow}}Consensus completed${{reset}}"
else
    logthis "Consensus generation already complete for {PIPELINE}/{SAMPLE}.unannotated.vcf, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
        )
    )


def produceConsensusUsingIvar(
    script: TextIOWrapper,
    options: OptionsDict,
):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# produce consensus using ivar
if [[ ! -f {PIPELINE}/{SAMPLE}.consensus.fa ]]; then
    logthis "${{yellow}}Building consensus {PIPELINE}/{SAMPLE}.sorted.bam${{reset}}"

    bcftools mpileup \\
        -aa \\
        -A \\
        --max-depth 0 \\
        --max-idepth 0 \\
        --count-orphans \\
        --min-BQ 0 \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.sorted.bam | \\
    ivar consensus -t 0 -m 3 -p {PIPELINE}/{SAMPLE}.consensus.fa

    sed -i 's/Consensus_{SAMPLE}.consensus_threshold_0_quality_20/{SAMPLE}/g' {PIPELINE}/{SAMPLE}.consensus.fa

    logthis "${{yellow}}Consensus completed${{reset}}"
else
    logthis "Consensus generation already complete for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
        )
    )


def producePileup(
    script: TextIOWrapper,
    options: OptionsDict,
):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    sample = options["sample"]
    pipeline = options["pipeline"]

    script.write(
        """
# create pileup
if [[ ! -f {PIPELINE}/{SAMPLE}.pileup.gz ]]; then
    logthis "Generate pileup for {PIPELINE}/{SAMPLE}.sorted.bam"

    bcftools mpileup \\
        --max-depth 0 \\
        --max-idepth 0 \\
        --count-orphans \\
        --min-BQ 0 \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.sorted.bam | gzip >{PIPELINE}/{SAMPLE}.pileup.gz 2>/dev/null

    logthis "Pileup completed for {PIPELINE}/{SAMPLE}.sorted.bam"
else
    logthis "Pileup already finished for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, ASSEMBLY=assembly, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def annotate(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    bin = options["bin"]

    referenceName = options["referenceName"]

    # we only run this for sars-cov-2 right now
    if options["__canAnnotateVariants"] == False:
        print(
            "Reference "
            + referenceName
            + " was requested. Variant annotation is only possible for SARS-CoV-2 right now."
        )
        return

    script.write(
        """
# annotate
if [[ ! -f {PIPELINE}/{SAMPLE}.nirvana.json.gz || ! -f {PIPELINE}/{SAMPLE}.nirvana.json.gz.jsi ]]; then
    logthis "Starting nirvana annotation"

    dotnet {BIN}/nirvana/Nirvana.dll \\
        -c {BIN}/nirvana/Data/Cache/SARS-CoV-2/SARS-CoV-2 \\
        -sd {BIN}/nirvana/Data/SupplementaryAnnotation/SARS-CoV-2 \\
        --enable-dq \\
        -r {BIN}/nirvana/Data/References/SARS-CoV-2.ASM985889v3.dat \\
        -i {PIPELINE}/{SAMPLE}.unannotated.vcf.gz \\
        -o {SAMPLE}.nirvana

    mv {SAMPLE}.nirvana.json.gz {PIPELINE}
    mv {SAMPLE}.nirvana.json.gz.jsi {PIPELINE}

    logthis "Completed nirvana annotation"
else
    logthis "nirvana annotations already completed, ${{green}}already completed${{reset}}"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.annotated.vcf.gz ]]; then
    logthis "Starting snpeff annotation"

    java -jar {BIN}/snpEff/snpEff.jar ann \\
        -htmlStats {PIPELINE}/{SAMPLE}.snpeff.html \\
        -noLog \\
        -verbose \\
        {REFERENCE_NAME} {PIPELINE}/{SAMPLE}.unannotated.vcf.gz | \\
    bgzip >{PIPELINE}/{SAMPLE}.annotated.vcf.gz

    tabix -p vcf {PIPELINE}/{SAMPLE}.annotated.vcf.gz

    logthis "Completed snpeff annotation"
else
    logthis "snpeff annotations already completed, ${{green}}already completed${{reset}}"
fi
""".format(
            SAMPLE=sample, BIN=bin, PIPELINE=pipeline, REFERENCE_NAME=referenceName
        )
    )


def assignClade(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    referenceName = options["referenceName"]

    # we only run this for sars-cov-2 right now
    if options["__canAssignClades"] == False:
        print(
            "Reference "
            + referenceName
            + " was requested. Clade assignment is only possible for SARS-CoV-2 right now."
        )
        return

    script.write(
        """
# assign clade
if [[ ! -f {PIPELINE}/{SAMPLE}.nextclade.tsv ]]; then
    logthis "Starting nextclade characterization for {SAMPLE}"

    nextclade run \\
        --input-dataset {REFERENCE}/nextclade-data \\
        --output-all {PIPELINE}/ \\
        --output-basename {SAMPLE}.nextclade \\
        -- \\
        {PIPELINE}/{SAMPLE}.consensus.fa

    logthis "Clade assignment complete for {SAMPLE}"
else
    logthis "Clade assignment already complete for {SAMPLE}, ${{green}}skipping${{reset}}"
fi
    """.format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def runVariantPipeline(script: TextIOWrapper, options: OptionsDict):
    referenceName = options["referenceName"]

    script.write("\n")
    script.write("#\n")
    script.write("# Running pipeline\n")
    script.write("#\n")

    callVariants(script, options)
    producePileup(script, options)
    annotate(script, options)

    script.write("\n")


def generateConsensus(script: TextIOWrapper, options: OptionsDict):
    consensusGenerator = options["consensusGenerator"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    if consensusGenerator == "ivar":
        produceConsensusUsingIvar(script, options)
    elif consensusGenerator == "bcftools":
        produceConsensusUsingBcftools(script, options)

    script.write(
        """
if [[ ! -f {REFERENCE}/{ASSEMBLY}.flat.fna ]]; then
    tail -n +2 {REFERENCE}/{ASSEMBLY}.fna | 
        grep -v '^>' | 
        tr -d '\\n' | 
        sed 's/\\(.\\)/\\1 /g' | 
        tr ' ' '\\n' > {REFERENCE}/{ASSEMBLY}.flat.fna
else
    logthis "Flat reference already generated, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.diff ]]; then
    logthis "Generating consensus difference for {SAMPLE}"

    tail -n +2 {PIPELINE}/{SAMPLE}.consensus.fa | 
        grep -v '^>' | 
        tr -d '\\n' | 
        sed 's/\\(.\\)/\\1 /g' | 
        tr ' ' '\\n' >{PIPELINE}/{SAMPLE}.consensus.flat.fna
    dwdiff -L8 -s -3 \\
        {REFERENCE}/{ASSEMBLY}.flat.fna \\
        {PIPELINE}/{SAMPLE}.consensus.flat.fna >{PIPELINE}/{SAMPLE}.diff

    logthis "Consensus difference already completed for {SAMPLE}, ${{green}}skipping${{reset}}"
else
    logthis "Consensus difference already generated, ${{green}}skipping${{reset}}"
fi
        """.format(
            REFERENCE=reference, ASSEMBLY=assembly, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def doVariantQC(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    checks = []

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
        """'if [[ ! -d {STATS}/{SAMPLE}_bcfstats ]]; then bcftools stats --fasta-ref {REFERENCE}/{ASSEMBLY}.fna {PIPELINE}/{SAMPLE}.unannotated.vcf.gz > {STATS}/{SAMPLE}.chk; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
        )
    )

    return checks


def doAlignmentQC(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]
    skipFastQc = options["skipFastQc"]

    filenames = getTrimmedFileNames(options)

    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])

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
        """'if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then samtools stats -@ 8 -r {REFERENCE}/{ASSEMBLY}.fna {SORTED} >{STATS}/{SAMPLE}.samstats; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.samidx ]]; then samtools idxstats -@ 8 {SORTED} >{STATS}/{SAMPLE}.samidx; fi' \\\n""".format(
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


def doQualityControl(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])

    plot_vcfstats = """plot-vcfstats --prefix {STATS}/{SAMPLE}_bcfstats {STATS}/{SAMPLE}.chk""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    # this doesn't have a test, it's fast enough that we can afford to run it
    plot_bamstats = """plot-bamstats --prefix {STATS}/{SAMPLE}_samstats/ {STATS}/{SAMPLE}.samstats""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    alignment_checks = doAlignmentQC(script, options)
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


def writeVersions(script: TextIOWrapper):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...\n")
    script.write("#\n")


def writeEnvironment(script: TextIOWrapper, options: OptionsDict):
    noColor = options["noColor"]

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

    script.write("#\n")
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
export LD_LIBRARY_PATH={WORKING}/bin:/usr/lib64:/usr/local/lib/:$LB_LIBRARY_PATH
export LD_PRELOAD={WORKING}/bin/libz.so.1.2.11.zlib-ng

# bcftools
export BCFTOOLS_PLUGINS={WORKING}/bin/plugins

# handy path
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.6.1:{WORKING}/bin:$PATH\n""".format(
            WORKING=options["working"]
        )
    )
    script.write("\n")


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

    for opt in ["working", "reference", "pipeline", "stats", "temp", "bin", "script", "adapters"]:
        options[opt] = expandvars(options[opt])

    return options


def verifyOptions(options: OptionsDict):
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

    ## this is a bit of a hack, but necessary since the tools and databases
    ## to do clade assignment and variant annotation don't exist right now
    referenceName = options["referenceName"]
    options["__canAssignClades"] = referenceName == "MN908947.3"
    options["__canAnnotateVariants"] = referenceName == "MN908947.3"


def getFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    fastq_dir = options["fastq_dir"]

    return (
        "{FASTQ_DIR}/{SAMPLE}_R1_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
        "{FASTQ_DIR}/{SAMPLE}_R2_001.fastq.gz".format(FASTQ_DIR=fastq_dir, SAMPLE=sample),
    )


def getTrimmedFileNames(options: OptionsDict) -> FastqSet:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return (
        "{PIPELINE}/{SAMPLE}_R1.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
        "{PIPELINE}/{SAMPLE}_R2.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample),
    )


def runIdentityPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
):
    script.write(
        """
#
# run the fastp preprocessor
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
    adapters = options["adapters"]
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
    logthis "${{yellow}}Running FASTP preprocessor${{reset}}"

    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        fastp \\
            --report_title "fastp report for sample {SAMPLE}" \\
            --in1 {R1} \\
            --in2 {R2} \\
            --out1 {O1} \\
            --out2 {O2} \\
            {ADAPTERS} \\
            --verbose {LIMITREADS} \\
            --thread 8 \\
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
            ADAPTERS="--adapter_fasta " + adapters if adapters != "" else "--detect_adapter_for_pe",
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            LIMITREADS="--reads_to_process " + str(readLimit) if readLimit > 0 else "",
        )
    )


def runTrimmomaticPreprocessor(
    script: TextIOWrapper,
    r1: str,
    r2: str,
    o1: str,
    o2: str,
    options: OptionsDict,
):
    bin = options["bin"]

    script.write(
        """
#
# run the trimmomatic preprocessor
#
if [[ ! -f {O1} || ! -f {O2} ]]; then
    logthis "${{yellow}}Running trimmomatic preprocessor${{reset}}"

    LD_PRELOAD={BIN}/libz.so.1.2.11.zlib-ng \\
        java -jar {BIN}/trimmomatic-0.39.jar PE \\
            {R1} \\
            {R2} \\
            {O1} /dev/null \\
            {O2} /dev/null \\
            ILLUMINACLIP:{BIN}/adapters/NexteraPE-PE.fa:2:30:10 \\
            LEADING:5 \\
            TRAILING:5 \\
            SLIDINGWINDOW:4:20 \\
            MINLEN:30

    logthis "${{yellow}}trimmomatic preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1, R2=r2, O1=o1, O2=o2, BIN=bin
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
    elif preprocessor == "trimmomatic":
        runTrimmomaticPreprocessor(script, r1, r2, o1, o2, options)
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)

    pass


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
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    bwa-mem2 mem -t {THREADS} \\
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
            O1=o1,
            O2=o2,
            ASSEMBLY=assembly,
            REFERENCE=reference,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            BIN=bin,
        )
    )


def runHisatAligner(
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
    bin = options["bin"]

    script.write(
        """
#
# align the input files
#
if [[ ! -f {PIPELINE}/{SAMPLE}.aligned.bam ]]; then
    logthis "${{yellow}}Running aligner${{reset}}"

    hisat2 \\
            -x {REFERENCE}/{ASSEMBLY} \\
            -1 {O1} \\
            -2 {O2} \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:ILLUMINA" \\
            --rg "PU:unspecified" \\
            --rg "LB:{SAMPLE}" \\
            --rg "SM:{SAMPLE}" \\
            --no-spliced-alignment \\
            --no-unal \\
            --threads {THREADS} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi
""".format(
            O1=o1,
            O2=o2,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            THREADS=threads,
            PIPELINE=pipeline,
            BIN=bin,
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
    elif aligner == "hisat2":
        runHisatAligner(script, o1, o2, options)
    else:
        print("Unexpected value {ALIGNER} given for the --aligner option".format(ALIGNER=aligner))
        quit(1)

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


def alignAndSort(script: TextIOWrapper, options: OptionsDict):
    processUnmapped = options["processUnmapped"]
    alignOnly = options["alignOnly"]

    filenames = getFileNames(options)
    trimmedFilenames = getTrimmedFileNames(options)

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filenames[0], filenames[1], trimmedFilenames[0], trimmedFilenames[1], options)
    alignFASTQ(script, trimmedFilenames[0], trimmedFilenames[1], options)
    sortAlignedAndMappedData(script, options)

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


def defineArguments() -> Namespace:
    parser = argparse.ArgumentParser()
    parser.set_defaults(doQC=False, cleanIntermediateFiles=True)
    parser.add_argument(
        "--sample",
        required=True,
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>_R[12].fastq.gz",
    )
    parser.add_argument(
        "--work-dir",
        required=True,
        action="store",
        metavar="WORKING_DIR",
        dest="working",
        help="Working directory, e.g. base for $WORKING/pipeline, $WORKING/stats",
    )

    parser.add_argument(
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
        "--preprocessor",
        action="store",
        dest="preprocessor",
        default="none",
        choices=["trimmomatic", "fastp", "none"],
        help="Optionally run a FASTQ preprocessor",
    )

    parser.add_argument(
        "--unmapped",
        action="store_true",
        dest="processUnmapped",
        default=False,
        help="Extract unmapped reads into secondary _R1 and _R2 FASTQ files",
    )

    parser.add_argument(
        "--align-only",
        action="store_true",
        dest="alignOnly",
        default=False,
        help="Only run alignment and sorting processes",
    )

    parser.add_argument(
        "--aligner",
        action="store",
        dest="aligner",
        default="bwa",
        choices=["bwa", "hisat2"],
        help="Use 'bwa' or 'hisat2' as the aligner.",
    )

    parser.add_argument(
        "--sorter",
        action="store",
        dest="sorter",
        default="biobambam",
        choices=["biobambam", "samtools"],
        help="Use 'biobambam' or 'samtools' as the sorter.",
    )

    parser.add_argument(
        "--caller",
        action="store",
        dest="caller",
        default="bcftools",
        choices=["bcftools", "gatk", "lofreq"],
        help="Use `bcftools`, `gatk`, or `lofreq` as the variant caller",
    )

    parser.add_argument(
        "--skip-annotation",
        action="store_true",
        dest="skipAnnotation",
        default=False,
        help="Skip all VCF annotation processes",
    )

    parser.add_argument(
        "--no-color",
        action="store_true",
        dest="noColor",
        default=False,
        help="Turn off colorized log output",
    )

    parser.add_argument(
        "--consensus-generator",
        action="store",
        dest="consensusGenerator",
        default="bcftools",
        choices=["ivar", "bcftools"],
        help="Choice of consensus FASTA generator",
    )

    parser.add_argument(
        "--reference-assembly",
        required=True,
        action="store",
        metavar="REFERENCE_ASSEMBLY",
        dest="referenceAssembly",
        help="Base name of the reference assembly",
    )

    help = (
        "'Chromosome' name from reference assembly "
        + "(MN908947.3, sars-cov-2), "
        + "(AF304460.1, hcov-229e), "
        + "(AY597011.2, hcov-hku1), "
        + "(AY567487.2, hcov-nl63), "
        + "(AY585228.1, hcov-oc43)"
    )

    parser.add_argument(
        "--reference-name",
        required=True,
        action="store",
        metavar="REFERENCE_NAME",
        dest="referenceName",
        choices=["MN908947.3", "AF304460.1", "AY597011.2", "AY567487.2", "AY585228.1"],
        help=help,
    )

    parser.add_argument(
        "--reference-dir",
        action="store",
        metavar="REFERENCE_DIR",
        dest="reference",
        help="Location of the reference genome files",
    )
    parser.add_argument(
        "--pipeline-dir",
        action="store",
        metavar="PIPELINE_DIR",
        dest="pipeline",
        help="Output directory for processing",
    )
    parser.add_argument(
        "--fastq-dir",
        action="store",
        metavar="FASTQ_DIR",
        dest="fastq_dir",
        help="Location of R1 and R2 files",
    )
    parser.add_argument(
        "--temp-dir",
        action="store",
        metavar="TEMP_DIR",
        dest="temp",
        help="Temporary storage for alignment",
    )
    parser.add_argument(
        "--stats-dir",
        action="store",
        metavar="STATS_DIR",
        dest="stats",
        help="Destination for statistics and QC files",
    )
    parser.add_argument(
        "--bin-dir",
        action="store",
        metavar="BIN_DIR",
        dest="bin",
        default="$HOME/bin",
        help="Install location of all tooling",
    )

    parser.add_argument(
        "--adapter-fasta",
        action="store",
        metavar="ADAPTER_FASTA",
        dest="adapters",
        default="",
        help="Optional FASTA file containing list of adapters to pass to fastp",
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
        "--cores",
        action="store",
        dest="cores",
        metavar="CPU_COUNT",
        default=cpu_count(),
        help="Specify the number of available CPU",
    )

    parser.add_argument(
        "--read-limit",
        action="store",
        dest="read-limit",
        metavar="READ_LIMIT",
        default=0,
        help="Limit the number of reads that fastp processes, for debugging",
    )

    return parser.parse_args()


def verifyFileNames(options: OptionsDict):
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


def main():
    opts = defineArguments()
    options = fixupPathOptions(opts)
    verifyFileNames(options)
    verifyOptions(options)

    filenames = getFileNames(options)

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

        alignAndSort(script, options)
        runVariantPipeline(script, options)
        generateConsensus(script, options)
        assignClade(script, options)

        # we'll wait here to make sure all the background stuff is done before we
        # run multiqc and cleanup
        script.write('logthis "${green}Done with front-end processing${reset}"\n')

        doQualityControl(script, options)

        script.write('logthis "${green}Done with back-end processing${reset}"\n')

        script.write("\n")
        script.write("touch {PIPELINE}/01-completed\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

    system("chmod +x " + options["script"])


if __name__ == "__main__":
    main()
