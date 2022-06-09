from io import TextIOWrapper
from typing import Dict, Any
from os.path import exists, expandvars
from argparse import Namespace

OptionsDict = Dict[str, Any]


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    bin = options["bin"]

    script.write("#\n")
    script.write("# Build the reference dictionary and interval list\n")
    script.write("#\n")

    script.write(
        """
if [[ ! -f {REFERENCE}/hcov-nl63.dict ]]; then
    java -jar {BIN}/picard.jar CreateSequenceDictionary \\
        -R {REFERENCE}/hcov-nl63.fasta \\
        -O {REFERENCE}/hcov-nl63.dict
else
    echo "Reference dictionary {REFERENCE}/hcov-nl63.dict ${{green}}already present${{reset}}"
fi

if [[ ! -f {REFERENCE}/ref_genome_autosomal.interval_list ]]; then
    # build the interval list, this is only done in the case where we're
    # processing a partial set of chromosomes. in the typical case this would
    # be a WGS collection.

    egrep '(NC_005831.2)\\s' {REFERENCE}/hcov-nl63.fasta.fai |
        awk '{{print $1"\\t1\\t"$2"\\t+\\t"$1}}' |
        cat {REFERENCE}/hcov-nl63.dict - >{REFERENCE}/ref_genome_autosomal.interval_list
else
    echo "Interval list {REFERENCE}/ref_genome_autosomal.interval_list ${{green}}already present${{reset}}"
fi

""".format(
            REFERENCE=reference, BIN=bin
        )
    )


def sortWithBiobambam(script: TextIOWrapper, options: OptionsDict, output: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    alignOnly = options["alignOnly"]

    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
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
    {EXIT_IF_ALIGN_ONLY}
    """.format(
            SAMPLE=sample,
            THREADS=threads,
            SORTED=output,
            PIPELINE=pipeline,
            STATS=stats,
            BIN=bin,
            TEMP=temp,
            EXIT_IF_ALIGN_ONLY="exit" if alignOnly == True else "",
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

    alignOnly = options["alignOnly"]

    unmarked = "{PIPELINE}/{SAMPLE}.unmarked.bam".format(PIPELINE=pipeline, SAMPLE=sample)

    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    samtools sort {PIPELINE}/{SAMPLE}.aligned.bam -o {UNMARKED}

    java -Xmx8g -jar {BIN}/picard.jar MarkDuplicates \\
        --TAGGING_POLICY All \\
        --REFERENCE_SEQUENCE {REFERENCE}/hcov-nl63.fasta \\
        -I {UNMARKED} \\
        -O {SORTED} \\
        -M {STATS}/{SAMPLE}_marked_dup_metrics.txt    

    # generate an index on the result
    samtools index -b {SORTED} {SORTED}.bai
else
    echo "{SORTED}, index, and metrics found, ${{green}}skipping${{reset}}"
fi
    {EXIT_IF_ALIGN_ONLY}
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
            EXIT_IF_ALIGN_ONLY="exit" if alignOnly == True else "",
        )
    )


def sortAlignedAndMappedData(script: TextIOWrapper, options: OptionsDict, output: str):
    sorter = options["sorter"]

    if sorter == "biobambam":
        sortWithBiobambam(script, options, output)
    else:
        sortWithSamtools(script, options, output)


def callVariantsUsingGatk(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    bam: str,
    vcf: str,
):
    script.write(
        """
# call variants
if [[ ! -f {VCF} ]]; then
    gatk HaplotypeCaller --java-options '-Xmx8g' \\
        -R {REFERENCE}/hcov-nl63.fasta \\
        -I {BAM} \\
        -O {VCF} \\
        --verbosity ERROR \\
        --pairHMM FASTEST_AVAILABLE \\
        --native-pair-hmm-threads 4

    echo Completed variant calling for {BAM}
else
    echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, BAM=bam, VCF=vcf, SAMPLE=sample
        )
    )


def callVariantsUsingLofreq(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    bam: str,
    vcf: str,
):
    script.write(
        """
# call variants
if [[ ! -f {VCF} ]]; then
    lofreq indelqual \\
        --dindel \\
        --ref {REFERENCE}/hcov-nl63.fasta \\
        --verbose \\
        --out - \\
        {BAM} |
    lofreq call \\
        --ref {REFERENCE}/hcov-nl63.fasta \\
        --call-indels \\
        --no-default-filter \\
        --max-depth 1000000 \\
        --force-overwrite \\
        --verbose \\
        --out {VCF} \\
        -

    echo Completed variant calling for {BAM}
else
    echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, BAM=bam, VCF=vcf, SAMPLE=sample
        )
    )


def callVariantsUsingBcftools(script: TextIOWrapper, reference: str, bam: str, vcf: str, gvcf: str, ploidy: str):
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
        --max-depth 1000000 \\
        --max-idepth 1000000 \\
        --threads 4 \\
        --output-type u \\
        --fasta-ref {REFERENCE}/hcov-nl63.fasta \\
        {BAM} 2>/dev/null | \\
    bcftools call \\
        --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
        --variants-only \\
        --keep-alts \\
        --multiallelic-caller \\
        --ploidy-file {PLOIDY} \\
        --prior 0.05 \\
        --threads 4 \\
        --output-type v  \\
        --output {VCF}.tmp 2>/dev/null

    bcftools +fill-tags \\
        {VCF}.tmp \\
        --output-type v  \\
        --output {VCF} \\
        -- --tags AC,AN,AF,VAF,MAF 2>/dev/null

    rm -f {VCF}.tmp
else
    echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference, BAM=bam, VCF=vcf, GVCF=gvcf, PLOIDY=ploidy
        )
    )


def produceConsensusUsingBcftools(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    consensus: str,
    vcf: str,
):
    script.write(
        """
# call variants
if [[ ! -f {CONSENSUS} ]]; then
    echo Creating consensus for {VCF}

    bcftools view --output-type z <{VCF} >{VCF}.gz
    bcftools index {VCF}.gz
    bcftools consensus \\
        --fasta-ref {REFERENCE}/hcov-nl63.fasta \\
        {VCF}.gz \\
    | sed '/>/ s/$/ | {SAMPLE}/' >{CONSENSUS}

    echo Completed consensus generation for {VCF}
else
    echo "Consensus generation already complete for {VCF}, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            VCF=vcf,
            CONSENSUS=consensus,
            SAMPLE=sample,
        )
    )


def produceConsensusUsingIvar(
    script: TextIOWrapper,
    reference: str,
    sample: str,
    consensus: str,
    bam: str,
):
    script.write(
        """
# call variants
if [[ ! -f {CONSENSUS} ]]; then
    echo Creating consensus for {BAM}

    samtools mpileup \\
        -aa \\
        -A \\
        --max-depth 0 \\
        --max-idepth 0 \\
        --count-orphans \\
        --min-BQ 0 \\
        --fasta-ref {REFERENCE}/hcov-nl63.fasta \\
        {BAM} | \\
    ivar consensus -t 0 -m 3 -p {CONSENSUS}

    sed -i 's/Consensus_{SAMPLE}.consensus_threshold_0_quality_20/{SAMPLE}/g' {CONSENSUS}

    echo Completed consensus generation for {BAM}
else
    echo "Consensus generation already complete for {BAM}, ${{green}}skipping${{reset}}"
fi
""".format(
            REFERENCE=reference,
            BAM=bam,
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
        -aa \\
        -A \\
        --max-depth 0 \\
        --max-idepth 0 \\
        --count-orphans \\
        --min-BQ 0 \\
        --fasta-ref {REFERENCE}/hcov-nl63.fasta \\
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


def annotate(script: TextIOWrapper, options: OptionsDict, vcf: str, annotated: str):
    sample = options["sample"]
    pipeline = options["pipeline"]
    bin = options["bin"]

    script.write(
        """
# annotate
if [[ ! -f {PIPELINE}/{SAMPLE}.annotated.vcf ]]; then
    echo Start snpEff annotation for {VCF}

    java -jar ~/bin/snpEff/snpEff.jar -htmlStats {PIPELINE}/{SAMPLE}.snpeff.html NC_005831.2 {VCF} >{ANNOTATED}

    echo snpEff annotion complete for {VCF}
else
    echo "snpEff annotations already for {VCF}, ${{green}}skipping${{reset}}"
fi
""".format(
            VCF=vcf, SAMPLE=sample, BIN=bin, PIPELINE=pipeline, ANNOTATED=annotated
        )
    )


def runPipeline(script: TextIOWrapper, options: OptionsDict, prefix: str):
    caller = options["caller"]
    useAlternateConsensus = options["alternateConsensus"]
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]

    bam = """{PREFIX}.sorted.bam""".format(PREFIX=prefix)
    vcf = """{PREFIX}.vcf""".format(PREFIX=prefix)
    gvcf = """{PREFIX}.gvcf""".format(PREFIX=prefix)
    annotated = """{PREFIX}.annotated.vcf""".format(PREFIX=prefix)
    pileup = """{PREFIX}.pileup.gz""".format(PREFIX=prefix)
    consensus = """{PREFIX}.consensus.fa""".format(PREFIX=prefix)
    ploidy = """{PREFIX}.ploidy""".format(PREFIX=prefix)

    script.write("\n")
    script.write("#\n")
    script.write("# Running pipeline for {BAM}\n".format(BAM=bam))
    script.write("#\n")

    if caller == "gatk":
        callVariantsUsingGatk(script, options["reference"], sample, bam, vcf)
    elif caller == "bcftools":
        callVariantsUsingBcftools(script, options["reference"], bam, vcf, gvcf, ploidy)
    elif caller == "lofreq":
        callVariantsUsingLofreq(script, options["reference"], sample, bam, vcf)
    else:
        print("Unexpected value {CALLER} given for the --caller option".format(CALLER=caller))
        quit(1)

    producePileup(script, options["reference"], bam, pileup, ploidy)
    annotate(script, options, vcf, annotated)

    if useAlternateConsensus == False:
        produceConsensusUsingBcftools(script, options["reference"], sample, consensus, vcf)
    else:
        produceConsensusUsingIvar(script, options["reference"], sample, consensus, bam)

    script.write(
        """
if [[ ! -f {REFERENCE}/hcov-nl63.flat.fasta ]]; then
    tail -n +2 {REFERENCE}/hcov-nl63.fasta | 
        grep -v '^>' | 
        tr -d '\\n' | 
        sed 's/\\(.\\)/\\1 /g' | 
        tr ' ' '\\n' > {REFERENCE}/hcov-nl63.flat.fasta
else
    echo "Flat reference already generated, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {PIPELINE}/{SAMPLE}.diff ]]; then
    echo Generating consensus difference for {SAMPLE}

    tail -n +2 {PIPELINE}/{SAMPLE}.consensus.fa | 
        grep -v '^>' | 
        tr -d '\\n' | 
        sed 's/\\(.\\)/\\1 /g' | 
        tr ' ' '\\n' >{PIPELINE}/{SAMPLE}.consensus.flat.fasta
    dwdiff -L8 -s -3 \\
        {REFERENCE}/hcov-nl63.flat.fasta \\
        {PIPELINE}/{SAMPLE}.consensus.flat.fasta >{PIPELINE}/{SAMPLE}.diff

    echo Consensus difference already completed for {SAMPLE}, ${{green}}skipping${{reset}}
else
    echo "Consensus difference already generated, ${{green}}skipping${{reset}}"
fi
        """.format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, CONSENSUS=consensus
        )
    )

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
#        -O {STATS}/{SAMPLE}
# else
#     echo "Variant metrics already run, ${{green}}skipping${{reset}}"
# fi

#
# we need to quiet vcftools here because it's stupid chatty and doesn't have an option to quiet
#
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --freq2 --out {STATS}/{SAMPLE}_vcfstats --max-alleles 2 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --depth --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --site-mean-depth --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --site-quality --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --missing-indv --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --missing-site --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &
vcftools --vcf {PIPELINE}/{SAMPLE}.vcf --het --out {STATS}/{SAMPLE}_vcfstats 2>/dev/null &

echo "Waiting for VCF QC processed to complete"
wait
echo "VCF QC processes complete"

#
# run some other stats on the vcf file
#
if [[ ! -d {STATS}/{SAMPLE}_bcfstats ]]; then
    bcftools stats --fasta-ref {REFERENCE}/hcov-nl63.fasta {PIPELINE}/{SAMPLE}.vcf > {STATS}/{SAMPLE}.chk
    plot-vcfstats --prefix {STATS}/{SAMPLE}_bcfstats {STATS}/{SAMPLE}.chk
else
    echo "bcftools stats and plots already run, ${{green}}skipping${{reset}}"
fi

""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, STATS=stats
        )
    )


def runAlignmentQC(script: TextIOWrapper, options: OptionsDict, sorted: str, aligned: str):
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
        -R {REFERENCE}/hcov-nl63.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.alignment_metrics.txt &
else
    echo "Alignment metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then
    gatk CollectGcBiasMetrics --java-options '-Xmx8g' \\
        --VERBOSITY ERROR \\
        -R {REFERENCE}/hcov-nl63.fasta \\
        -I {SORTED} \\
        -O {STATS}/{SAMPLE}.gc_bias_metrics.txt \\
        -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf \\
        -S {STATS}/{SAMPLE}.gc_bias_summary.txt &
else
    echo "GC bias metrics already run, ${{green}}skipping${{reset}}"
fi

if [[ ! -f {STATS}/{SAMPLE}.samstats ]]; then
    (samtools stats -@ 8 \\
        -r {REFERENCE}/hcov-nl63.fasta \\
        {SORTED} >{STATS}/{SAMPLE}.samstats

    plot-bamstats --prefix {STATS}/{SAMPLE}_samstats {STATS}/{SAMPLE}.samstats) &
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

echo "Waiting for QC processed to complete"
wait
echo "QC processes complete"

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

if [[ ! -d {STATS}/{SAMPLE}_multiqc_data ]]; then
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
fi
""".format(
            STATS=stats, SAMPLE=sample
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

# bcftools
export BCFTOOLS_PLUGINS={WORKING}/bin/plugins

# handy path
export PATH={WORKING}/bin/ensembl-vep:{WORKING}/bin/FastQC:{WORKING}/bin/gatk-4.2.3.0:{WORKING}/bin:$PATH\n""".format(
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
