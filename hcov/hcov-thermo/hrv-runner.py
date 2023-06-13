#!/usr/bin/env python

from io import TextIOWrapper
import argparse
from argparse import Namespace
from typing import Tuple, Dict, Any
from datetime import datetime
from os.path import exists, expandvars
from os import cpu_count, system
from sys import platform

OptionsDict = Dict[str, Any]

panel_choices = ["NC_038311.1", "NC_038312.1", "NC_038878.1", "panel"]
panel_choice_help = (
    "'Chromosome' name from reference assembly "
    + "(NC_038311.1, hrv-a), "
    + "(NC_038312.1, hrv-b), "
    + "(NC_038878.1, hrv-c), "
    + "(panel, combined panel of all organisms)"
)


def getLibraryPath(options: OptionsDict) -> str:
    working = options["working"]

    if platform == "linux" or platform == "linux2":
        return "LD_LIBRARY_PATH={WORKING}/lib".format(WORKING=working)

    elif platform == "darwin":
        return "DYLD_LIBRARY_PATH={WORKING}/lib".format(WORKING=working)

    elif platform == "win32":
        print("What are you doing trying to run this on a Windows device? Only MacOS and Linux right now.")
        quit(1)

    # not reached
    return ""


def updateDictionary(script: TextIOWrapper, options: OptionsDict):
    bin = options["bin"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    chromosome = options["referenceName"]

    if chromosome == "panel":
        chromosome = "|".join([chr for chr in panel_choices if chr != "panel"])

    script.write("#\n")
    script.write("# Build the reference dictionary\n")
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

    {LIBRARY} {BIN}/bamsormadup \\
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
            LIBRARY=getLibraryPath(options),
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


def callVariants(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using bcftools${{reset}}"

    bcftools mpileup \\
        --annotate FORMAT/DP,FORMAT/SCR \\
        --max-depth 1000000 \\
        --max-idepth 1000000 \\
        --threads 8 \\
        --output-type u \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.sorted.bam 2>/dev/null | \\
    bcftools call \\
        --variants-only \\
        --multiallelic-caller \\
        --ploidy 1 \\
        --prior 0.05 \\
        --threads 8 \\
        --output-type u  \\
        --output {PIPELINE}/{SAMPLE}.bcf.tmp 2>/dev/null

    echo "Normalizing and filtering"
    bcftools norm -f {REFERENCE}/{ASSEMBLY}.fna \\
        {PIPELINE}/{SAMPLE}.bcf.tmp \\
        -Ou \\
        -o {PIPELINE}/{SAMPLE}.norm.bcf.tmp >/dev/null
    bcftools filter --IndelGap 5 \\
        {PIPELINE}/{SAMPLE}.norm.bcf.tmp \\
        -Ou \\
        -o {PIPELINE}/{SAMPLE}.unannotated.bcf

    echo "Indexing temporary BCF and calling consensus"
    bcftools index --force {PIPELINE}/{SAMPLE}.unannotated.bcf

    bcftools consensus \\
        --sample {SAMPLE} \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        -o {PIPELINE}/{SAMPLE}.consensus.fa \\
        {PIPELINE}/{SAMPLE}.unannotated.bcf

    echo "Converting temporary BCF to VCF"
    bcftools convert \\
        -Ov {PIPELINE}/{SAMPLE}.unannotated.bcf \\
        -o {PIPELINE}/{SAMPLE}.unannotated.vcf

    bgzip --force --keep {PIPELINE}/{SAMPLE}.unannotated.vcf
    tabix --force -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "${{green}}bcftools variant calling completed${{reset}}"
else
    logthis "Variants already called for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
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
            "Reference '"
            + referenceName
            + "' was requested. Variant annotation is only possible for SARS-CoV-2 right now."
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

    tabix --force -p vcf {PIPELINE}/{SAMPLE}.annotated.vcf.gz

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
            "Reference '"
            + referenceName
            + "' was requested. Clade assignment is only possible for SARS-CoV-2 right now."
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
    script.write("\n")
    script.write("#\n")
    script.write("# Running pipeline\n")
    script.write("#\n")

    callVariants(script, options)
    annotate(script, options)

    script.write("\n")


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
    doPicardQc = options["doPicardQc"]

    filename = getTrimmedFileName(options)

    sorted = "{PIPELINE}/{SAMPLE}.sorted.bam".format(PIPELINE=options["pipeline"], SAMPLE=options["sample"])

    checks = []

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.flagstat.txt ]]; then samtools flagstat -@ 8 {SORTED} >{STATS}/{SAMPLE}.flagstat.txt; fi' \\\n""".format(
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    if doPicardQc == True:
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
            """'if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then gatk CollectWgsMetrics --java-options -Xmx4g --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.wgs_metrics.txt --MINIMUM_BASE_QUALITY 10 --MINIMUM_MAPPING_QUALITY 10 --COVERAGE_CAP 2500 --USE_FAST_ALGORITHM --INCLUDE_BQ_HISTOGRAM; fi' \\\n""".format(
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

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.bedtools.coverage ]]; then samtools view -bq 15 -F 1284 {SORTED} | bedtools genomecov -d -ibam stdin | awk "\\$2 % 100 == 0 {{print \\$1,\\$2,\\$3}}" | sed "s/NC_038311.1/hrv_a/;s/NC_038312.1/hrv_b/;s/NC_038878.1/hrv_c/" >{STATS}/{SAMPLE}.bedtools.coverage; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    if skipFastQc == False:
        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}.fastqc.zip ]]; then fastqc --svg --threads 1 --outdir {STATS} --noextract {FILENAME}; fi' \\\n""".format(
                SAMPLE=sample,
                STATS=stats,
                FILENAME=filename,
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
    --cl-config 'custom_logo: "{STATS}/ovationlogo.png"' \\
    --cl-config 'custom_logo_url: "https://www.ovation.io"' \\
    --cl-config 'custom_logo_title: "Ovation"' \\
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
export PATH={WORKING}/perl5/bin:$PATH
export PERL5LIB={WORKING}/perl5/lib/perl5:$PERL5LIB
export PERL_LOCAL_LIB_ROOT={WORKING}/perl5:$PERL_LOCAL_LIB_ROOT

# shared library stuff (linux)
export LD_LIBRARY_PATH={WORKING}/lib:{WORKING}/bin:/usr/lib64:/usr/local/lib/:$LD_LIBRARY_PATH

# shared library stuff (darwin)
export DYLD_FALLBACK_LIBRARY_PATH={WORKING}/lib:{WORKING}/bin

# bcftools
export BCFTOOLS_PLUGINS={WORKING}/libexec/bcftools

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
    if options["bin"] == None:
        options["bin"] = "{WORKING}/bin".format(WORKING=options["working"])
    if options["lib"] == None:
        options["lib"] = "{WORKING}/lib".format(WORKING=options["working"])
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

    if exists(options["lib"]) == False:
        print("Unable to find your --lib-dir directory at {PATH}".format(PATH=options["lib"]))
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


fallback_warning_shown = False


def getFileName(options: OptionsDict) -> str:
    """
    Returns the input / source filename from the options. First, attempts the --fastq
    option as a full path. If that doesn't exist, construct a path from --fastq-dir
    as well and check it. If neither are present, this function exits the process.
    """
    global fallback_warning_shown

    fastq = options["fastq"]
    fastq_dir = options["fastq_dir"]

    if exists(expandvars(fastq)) == True:
        return fastq

    filename = "{FASTQ_DIR}/{FASTQ}".format(FASTQ_DIR=fastq_dir, FASTQ=fastq)

    if exists(expandvars(filename)) == True:
        return filename

    print("Check your --fastq and --fastq-dir parameters")
    quit(1)


def getTrimmedFileName(options: OptionsDict) -> str:
    sample = options["sample"]
    pipeline = options["pipeline"]

    return "{PIPELINE}/{SAMPLE}.trimmed.fastq.gz".format(PIPELINE=pipeline, SAMPLE=sample)


def runIdentityPreprocessor(
    script: TextIOWrapper,
    r1: str,
    o1: str,
):
    script.write(
        """
#
# run the fastp preprocessor
#
if [[ ! -f {O1} ]]; then
    logthis "${{yellow}}Running identity preprocessor${{reset}}"

    ln -s {R1} {O1}

    logthis "${{yellow}}Identity preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1,
            O1=o1,
        )
    )


def runCutadaptPreprocessor(
    script: TextIOWrapper,
    r1: str,
    o1: str,
    options: OptionsDict,
):
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]

    script.write(
        """
#
# run the cutadapt preprocessor
#
if [[ ! -f {O1} ]]; then
    logthis "${{yellow}}Running cutadapt preprocessor${{reset}}"

    cutadapt \\
        --cores {THREADS} \\
        --length 300 \\
        --cut 15 \\
        --minimum-length 30 \\
        --quality-cutoff 18,20 \\
        --trim-n \\
        --adapter ATCACCGACTGCCCATAGAGAGGCTGAGAC --times 3 \\
        --adapter AAAAAAAAAA$ --times 3 \\
        --adapter GGGGGGGGGG$ --times 3 \\
        --output {O1} \\
        {R1} >{STATS}/{SAMPLE}.cutadapt.log

    logthis "${{yellow}}cutadapt preprocessor completed${{reset}}"
else
    logthis "Preprocessor already run, ${{green}}skipping${{reset}}"
fi
""".format(
            R1=r1, O1=o1, SAMPLE=sample, STATS=stats, THREADS=threads
        )
    )


def preprocessFASTQ(
    script: TextIOWrapper,
    r1: str,
    o1: str,
    options: OptionsDict,
):
    preprocessor = options["preprocessor"]

    if preprocessor == "none":
        runIdentityPreprocessor(script, r1, o1)
    elif preprocessor == "cutadapt":
        runCutadaptPreprocessor(script, r1, o1, options)
    else:
        print("Unexpected value {PREPROCESSOR} given for the --preprocessor option".format(PREPROCESSOR=preprocessor))
        quit(1)


def runBwaAligner(
    script: TextIOWrapper,
    r1: str,
    options: OptionsDict,
):
    aligner = options["aligner"]
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

    {ALIGNER} mem \\
        -R "@RG\\tID:{SAMPLE}\\tPL:IONTORRENT\\tPU:unspecified\\tLB:{SAMPLE}\\tSM:{SAMPLE}" \\
        -t {THREADS} \\
        -Y \\
        -M \\
        -v 1 \\
        -S \\
        -P \\
        {REFERENCE}/{ASSEMBLY}.fna \\
        {R1} | 
    samtools view -Sb -@ 4 - >{PIPELINE}/{SAMPLE}.aligned.bam

    logthis "${{yellow}}Alignment completed${{reset}}"
else
    logthis "{PIPELINE}/{SAMPLE}.aligned.bam, aligned temp file found, ${{green}}skipping${{reset}}"
fi

""".format(
            ALIGNER=aligner,
            R1=r1,
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
    r1: str,
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
            -1 {R1} \\
            --rg-id "ID:{SAMPLE}" \\
            --rg "PL:IONTORRENT" \\
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
            R1=r1,
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
    r1: str,
    options: OptionsDict,
):
    aligner = options["aligner"]

    if aligner == "bwa" or aligner == "bwa-mem2":
        runBwaAligner(script, r1, options)
    elif aligner == "hisat2":
        runHisatAligner(script, r1, options)
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
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped_R1.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}_unmapped.fastq \\
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


def generateDepth(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# calculate depth by position
#
if [[ ! -f {PIPELINE}/{SAMPLE}.depth.gz ]]; then
    logthis "${{yellow}}Calculating depth by position${{reset}}"

    samtools depth {PIPELINE}/{SAMPLE}.sorted.bam | gzip >{PIPELINE}/{SAMPLE}.depth.gz

    logthis "${{yellow}}Depth calculation complete${{reset}}"
else
    logthis "Depth calculation already complete, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            PIPELINE=pipeline,
        )
    )


def alignAndSort(script: TextIOWrapper, options: OptionsDict):
    processUnmapped = options["processUnmapped"]
    alignOnly = options["alignOnly"]

    filename = getFileName(options)
    trimmedFilename = getTrimmedFileName(options)

    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    preprocessFASTQ(script, filename, trimmedFilename, options)
    alignFASTQ(script, trimmedFilename, options)
    sortAlignedAndMappedData(script, options)
    generateDepth(script, options)

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


def writeHeader(script: TextIOWrapper, options: OptionsDict, filename: str):
    script.write("#\n")
    script.write("# generated at {TIME}\n".format(TIME=datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    script.write("#\n")
    script.write("# Parameters\n")
    for opt in options.keys():
        script.write("#   {OPTION} = {VALUE}\n".format(OPTION=opt, VALUE=options[opt]))
    script.write("#\n")

    script.write("#\n")
    script.write("# FASTQ = {F}\n".format(F=filename))
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
        help="Short name of sample, all files will use this as an alias",
    )
    parser.add_argument(
        "--fastq",
        required=True,
        action="store",
        metavar="FASTQ",
        dest="fastq",
        help="Path name of FASTQ relative to --fastq-dir",
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
        "--skip-picard-qc",
        action="store_false",
        dest="doPicardQc",
        default=True,
        help="Skip any picard QC process on input and output files",
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
        choices=["cutadapt", "none"],
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
        choices=["bwa", "bwa-mem2", "hisat2"],
        help="Use 'bwa', 'bwa-mem2', or 'hisat2' as the aligner.",
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
        "--reference-assembly",
        required=True,
        action="store",
        metavar="REFERENCE_ASSEMBLY",
        dest="referenceAssembly",
        help="Base name of the reference assembly",
    )

    parser.add_argument(
        "--reference-name",
        required=True,
        action="store",
        metavar="REFERENCE_NAME",
        dest="referenceName",
        choices=panel_choices,
        help=panel_choice_help,
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
        help="Install location of all tooling",
    )
    parser.add_argument(
        "--lib-dir",
        action="store",
        metavar="LIB_DIR",
        dest="lib",
        help="Install location of all libraries",
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
    fastq = options["fastq"]
    fastq_dir = options["fastq_dir"]

    if exists(expandvars(fastq)) == True:
        return fastq

    filename = "{FASTQ_DIR}/{FASTQ}".format(FASTQ_DIR=fastq_dir, FASTQ=fastq)

    if exists(expandvars(filename)) == True:
        return filename

    print("Check your --fastq and --fastq-dir parameters")
    quit(1)


def main():
    opts = defineArguments()
    options = fixupPathOptions(opts)
    verifyFileNames(options)
    verifyOptions(options)

    filename = getFileName(options)

    with open(options["script"], "w+") as script:
        script.truncate(0)

        script.write("#!/usr/bin/env bash\n")
        writeHeader(script, options, filename)
        writeVersions(script)
        writeEnvironment(script, options)

        script.write("\n")
        script.write("touch {PIPELINE}/00-started\n".format(PIPELINE=options["pipeline"]))
        script.write("\n")

        updateDictionary(script, options)

        alignAndSort(script, options)
        runVariantPipeline(script, options)
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
