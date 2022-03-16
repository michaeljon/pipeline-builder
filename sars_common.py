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


def callVariants(
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
            -R {REFERENCE}/covid_reference.fasta \\
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


def callVariants2(
    script: TextIOWrapper,
    reference: str,
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
            --max-depth 10000 \\
            --threads 4 \\
            --output-type u \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
            {BAM} 2>/dev/null | \\
        bcftools call \\
            --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \\
            --variants-only \\
            --keep-alts \\
            --multiallelic-caller \\
            --prior 0.05 \\
            --threads 4 \\
            --output-type v  \\
            --output {VCF} 2>/dev/null

        echo Completed variant calling for {BAM}
    else
        echo "Variants already called for {BAM}, ${{green}}skipping${{reset}}"
    fi
""".format(
            REFERENCE=reference,
            BAM=bam,
            VCF=vcf,
            PLOIDY=ploidy,
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
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
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
            --count-orphans \\
            --min-BQ 0 \\
            --fasta-ref {REFERENCE}/covid_reference.fasta \\
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
            --count-orphans \\
            --min-BQ 0 \\
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


def annotate(script: TextIOWrapper, options: OptionsDict, vcf: str, annotated: str):
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

        java -jar ~/bin/snpEff/snpEff.jar -htmlStats {PIPELINE}/{SAMPLE}.snpeff.html NC_045512.2 {VCF} >{ANNOTATED}

        echo snpEff annotion complete for {VCF}
    else
        echo "snpEff annotations already for {VCF}, ${{green}}skipping${{reset}}"
    fi
""".format(
            VCF=vcf, SAMPLE=sample, BIN=bin, PIPELINE=pipeline, ANNOTATED=annotated
        )
    )


def assignClade(
    script: TextIOWrapper, options: OptionsDict, reference: str, consensus: str
):
    sample = options["sample"]
    pipeline = options["pipeline"]

    script.write(
        """
    # assign clade
    if [[ ! -f {PIPELINE}/{SAMPLE}.nextclade.tsv ]]; then
        nextclade \\
            --in-order \\
            --input-fasta {CONSENSUS} \\
            --input-dataset {REFERENCE}/nextclade-data/sars-cov-2 \\
            --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \\
            --output-json {PIPELINE}/{SAMPLE}.nextclade.json \\
            --output-tsv {PIPELINE}/{SAMPLE}.nextclade.tsv \\
            --output-tree {PIPELINE}/{SAMPLE}.nextclade.auspice.json \\
            --output-dir {PIPELINE}/ \\
            --output-basename {SAMPLE}.nextclade

        echo "Clade assignment complete for {SAMPLE}"
    else
        echo "Clade assignment already complete for {SAMPLE}, ${{green}}skipping${{reset}}"
    fi
    """.format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample, CONSENSUS=consensus
        )
    )


def runPipeline(script: TextIOWrapper, options: OptionsDict, prefix: str):
    useAlternateCaller = options["alternateCaller"]
    useAlternateConsensus = options["alternateConsensus"]
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
        callVariants(script, options["reference"], sample, bam, vcf)
    else:
        callVariants2(script, options["reference"], bam, vcf, ploidy)

    producePileup(script, options["reference"], bam, pileup, ploidy)
    annotate(script, options, vcf, annotated)

    if useAlternateConsensus == False:
        produceConsensusUsingBcftools(
            script, options["reference"], sample, consensus, vcf
        )
    else:
        produceConsensusUsingIvar(script, options["reference"], sample, consensus, bam)

    assignClade(script, options, options["reference"], consensus)

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
    (samtools stats -@ 8 \\
        -r {REFERENCE}/covid_reference.fasta \\
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
