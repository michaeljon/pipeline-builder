from io import TextIOWrapper

from bio_types import *


def sortAlignedAndMappedData(script: TextIOWrapper, options: OptionsDict, reference):
    sample = options["sample"]
    pipeline = options["pipeline"]
    threads = options["cores"]
    stats = options["stats"]
    bin = options["bin"]
    temp = options["temp"]

    aligned = "{PIPELINE}/{SAMPLE}.{ORGANISM}.aligned.bam".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"], ORGANISM=reference["common"]
    )
    sorted = "{PIPELINE}/{SAMPLE}.{ORGANISM}.sorted.bam".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"], ORGANISM=reference["common"]
    )

    script.write(
        """
#
# sort and mark duplicates
#
if [[ ! -f {SORTED} || ! -f {SORTED}.bai ]]; then
    logthis "${{yellow}}Sorting and marking duplicates${{reset}}"

    {BIN}/bamsormadup \\
        SO=coordinate \\
        threads={THREADS} \\
        level=6 \\
        tmpfile={TEMP}/{SAMPLE} \\
        inputformat=bam \\
        indexfilename={SORTED}.bai \\
        M={STATS}/{SAMPLE}.{ORGANISM}.duplication_metrics <{ALIGNED} >{SORTED}

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
            ALIGNED=aligned,
            SORTED=sorted,
            ORGANISM=reference["common"],
        )
    )


def extractUmappedReads(script: TextIOWrapper, options: OptionsDict, reference):
    pipeline = options["pipeline"]
    sample = options["sample"]

    aligned = "{PIPELINE}/{SAMPLE}.{ORGANISM}.aligned.bam".format(
        PIPELINE=options["pipeline"], SAMPLE=options["sample"], ORGANISM=reference["common"]
    )

    script.write(
        """
#
# extract unmapped reads
#
if [[ ! -f {PIPELINE}/{SAMPLE}.{ORGANISM}_unmapped.fastq ]]; then
    logthis "${{yellow}}Extracting unmapped reads${{reset}}"

    samtools fastq -N -f 4 \\
        -0 {PIPELINE}/{SAMPLE}.{ORGANISM}_unmapped_other.fastq \\
        -s {PIPELINE}/{SAMPLE}.{ORGANISM}_unmapped_singleton.fastq \\
        -1 {PIPELINE}/{SAMPLE}.{ORGANISM}_unmapped_R1.fastq \\
        -2 {PIPELINE}/{SAMPLE}.{ORGANISM}_unmapped_R2.fastq \\
        {ALIGNED}

    logthis "${{yellow}}Unmapped read extraction completed${{reset}}"
else
    logthis "Unmapped reads extracted to initial FASTQ, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            PIPELINE=pipeline,
            ALIGNED=aligned,
            ORGANISM=reference["common"],
        )
    )


def sortAndExtractUnmapped(script: TextIOWrapper, options: OptionsDict, reference):
    processUnmapped = options["processUnmapped"]
    alignOnly = options["alignOnly"]

    sortAlignedAndMappedData(script, options, reference)
    generateDepth(script, options, reference)

    if processUnmapped == True:
        extractUmappedReads(script, options, reference)

    if alignOnly == True:
        script.write(
            """

# align-only flag set
logthis "align-only set, exiting pipeline early"
exit
"""
        )


def callVariants(script: TextIOWrapper, options: OptionsDict, reference):
    sample = options["sample"]
    pipeline = options["pipeline"]
    referenceAssembly = reference["assembly"]
    bin = options["bin"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using freebayes${{reset}}"

    # calling variants using freebays
    freebayes \\
        --fasta-reference {REFERENCE_ASSEMBLY} \\
        --ploidy 1 \\
        --max-complex-gap 18 \\
        --use-duplicate-reads \\
        {PIPELINE}/{SAMPLE}.{ORGANISM}.sorted.bam >{PIPELINE}/{SAMPLE}.{ORGANISM}.tmp.vcf

    logthis "Normalizing called variants"
    bcftools norm \\
        --atomize \\
        --fasta-ref {REFERENCE_ASSEMBLY} \\
        -o {PIPELINE}/{SAMPLE}.{ORGANISM}.norm.bcf \\
        -Ou \\
        {PIPELINE}/{SAMPLE}.{ORGANISM}.tmp.vcf

    logthis "Filtering called variants"
    bcftools filter -i "QUAL > 100" \\
        -o {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf \\
        -Ou \\
        {PIPELINE}/{SAMPLE}.{ORGANISM}.norm.bcf

    logthis "Indexing normalized BCF and calling consensus"
    bcftools index --force {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf

    bcftools consensus \\
        --sample {SAMPLE} \\
        --fasta-ref {REFERENCE_ASSEMBLY} \\
        {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf \\
       | sed -E 's/^>([A-Z0-9\\.]+)/>{SAMPLE} | \\1/' >{PIPELINE}/{SAMPLE}.{ORGANISM}.consensus.fa

    logthis "Converting freebayes temporary BCF to VCF"
    bcftools convert \\
        -Ov \\
        {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf \\
        -o {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf

    bgzip --force --keep {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf
    tabix --force -p vcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz

    logthis "Cleaning up temporary freebayes output"
    rm -f {PIPELINE}/{SAMPLE}.{ORGANISM}.tmp.vcf
    rm -f {PIPELINE}/{SAMPLE}.{ORGANISM}.norm.bcf
    rm -f {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf
    rm -f {PIPELINE}/{SAMPLE}.{ORGANISM}.filtered.bcf.csi

    logthis "${{green}}freebayes variant calling completed${{reset}}"
else
    logthis "Variants already called via bcftools for {PIPELINE}/{SAMPLE}.{ORGANISM}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            BIN=bin,
            REFERENCE_ASSEMBLY=referenceAssembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            ORGANISM=reference["common"],
        )
    )


def runVariantPipeline(script: TextIOWrapper, options: OptionsDict, reference):
    script.write("\n")

    callVariants(script, options, reference)
    script.write("\n")


def doVariantQC(script: TextIOWrapper, options: OptionsDict, reference):
    referenceAssembly = reference["assembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    checks = []

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.frq ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --freq2 --out {STATS}/{SAMPLE}.{ORGANISM} --max-alleles 2 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.idepth ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --depth --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.ldepth.mean ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --site-mean-depth --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.lqual ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --site-quality --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.imiss ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --missing-indv --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.lmiss ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --missing-site --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.het ]]; then vcftools --gzvcf {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz --het --out {STATS}/{SAMPLE}.{ORGANISM} 2>/dev/null; fi' \\\n""".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -d {STATS}/{SAMPLE}.{ORGANISM}_bcfstats ]]; then bcftools stats --fasta-ref {REFERENCE_ASSEMBLY} {PIPELINE}/{SAMPLE}.{ORGANISM}.unannotated.vcf.gz > {STATS}/{SAMPLE}.{ORGANISM}.chk; fi' \\\n""".format(
            REFERENCE_ASSEMBLY=referenceAssembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    return checks


def doAlignmentQC(script: TextIOWrapper, options: OptionsDict, reference):
    referenceAssembly = reference["assembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]
    doPicardQc = options["doPicardQc"]
    bin = options["bin"]

    sorted = "{PIPELINE}/{SAMPLE}.{ORGANISM}.sorted.bam".format(
        PIPELINE=options["pipeline"],
        SAMPLE=options["sample"],
        ORGANISM=reference["common"],
    )

    checks = []

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.flagstat.txt ]]; then samtools flagstat -@ 8 {SORTED} >{STATS}/{SAMPLE}.{ORGANISM}.flagstat.txt; fi' \\\n""".format(
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
            ORGANISM=reference["common"],
        )
    )

    if doPicardQc == True:
        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.alignment_metrics.txt ]]; then java -jar {BIN}/picard.jar CollectAlignmentSummaryMetrics --VERBOSITY ERROR -R {REFERENCE_ASSEMBLY} -I {SORTED} -O {STATS}/{SAMPLE}.{ORGANISM}.alignment_metrics.txt; fi' \\\n""".format(
                REFERENCE_ASSEMBLY=referenceAssembly,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
                SORTED=sorted,
                BIN=bin,
                ORGANISM=reference["common"],
            )
        )

        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then java -jar {BIN}/picard.jar CollectGcBiasMetrics --VERBOSITY ERROR -R {REFERENCE_ASSEMBLY} -I {SORTED} -O {STATS}/{SAMPLE}.{ORGANISM}.gc_bias_metrics.txt -CHART {STATS}/{SAMPLE}.{ORGANISM}.gc_bias_metrics.pdf -S {STATS}/{SAMPLE}.{ORGANISM}.gc_bias_summary.txt; fi' \\\n""".format(
                REFERENCE_ASSEMBLY=referenceAssembly,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
                SORTED=sorted,
                BIN=bin,
                ORGANISM=reference["common"],
            )
        )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.samstats ]]; then samtools stats -@ 8 -r {REFERENCE_ASSEMBLY} {SORTED} >{STATS}/{SAMPLE}.{ORGANISM}.samstats; fi' \\\n""".format(
            REFERENCE_ASSEMBLY=referenceAssembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.samidx ]]; then samtools idxstats -@ 8 {SORTED} >{STATS}/{SAMPLE}.{ORGANISM}.samidx; fi' \\\n""".format(
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
            ORGANISM=reference["common"],
        )
    )

    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.samtools.coverage ]]; then samtools coverage -d 0 --reference {REFERENCE_ASSEMBLY} {SORTED} >{STATS}/{SAMPLE}.{ORGANISM}.samtools.coverage; fi' \\\n""".format(
            REFERENCE_ASSEMBLY=referenceAssembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
            ORGANISM=reference["common"],
        )
    )

    # todo: fix this, the reference names are broken
    checks.append(
        """'if [[ ! -f {STATS}/{SAMPLE}.{ORGANISM}.bedtools.coverage ]]; then samtools view -bq 30 -F 1284 {SORTED} | bedtools genomecov -d -ibam stdin | awk "\\$2 % 100 == 0 {{print \\$1,\\$2,\\$3}}" >{STATS}/{SAMPLE}.{ORGANISM}.bedtools.coverage; fi' \\\n""".format(
            REFERENCE_ASSEMBLY=referenceAssembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
            ORGANISM=reference["common"],
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


def doQualityControl(script: TextIOWrapper, options: OptionsDict, reference):
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    plot_vcfstats = (
        """plot-vcfstats --prefix {STATS}/{SAMPLE}.{ORGANISM}_bcfstats {STATS}/{SAMPLE}.{ORGANISM}.chk""".format(
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    # this doesn't have a test, it's fast enough that we can afford to run it
    plot_bamstats = (
        """plot-bamstats --prefix {STATS}/{SAMPLE}.{ORGANISM}_samstats/ {STATS}/{SAMPLE}.{ORGANISM}.samstats""".format(
            SAMPLE=sample,
            STATS=stats,
            ORGANISM=reference["common"],
        )
    )

    alignment_checks = doAlignmentQC(script, options, reference)
    variant_checks = doVariantQC(script, options, reference)

    cmd = ""
    if options["doAlignmentQc"] == True or options["doVariantQc"] == True:
        cmd = "parallel --joblog {PIPELINE}/{SAMPLE}.{ORGANISM}.qc.log ::: \\\n".format(
            PIPELINE=pipeline,
            SAMPLE=sample,
            ORGANISM=reference["common"],
        )

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


def doFastQc(script: TextIOWrapper, options: OptionsDict, filenames):
    if len(filenames) == 2:
        script.write(
            """
if [[ ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.html || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.html ]]; then 
    fastqc --svg --threads 2 --outdir {STATS} --noextract {R1} {R2}
fi
""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                R1=filenames[0],
                R2=filenames[1],
            )
        )
    else:
        script.write(
            """
if [[ ! -f {STATS}/{SAMPLE}.fastqc.zip ]]; then
    fastqc --svg --threads 1 --outdir {STATS} --noextract {FILENAME}
fi
""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                FILENAME=filenames[0],
            )
        )


def generateDepth(script: TextIOWrapper, options: OptionsDict, reference):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# calculate depth by position
#
if [[ ! -f {PIPELINE}/{SAMPLE}.{ORGANISM}.depth.gz ]]; then
    logthis "${{yellow}}Calculating depth by position${{reset}}"

    samtools depth -@ 8 -aa -a -J {PIPELINE}/{SAMPLE}.{ORGANISM}.sorted.bam | gzip >{PIPELINE}/{SAMPLE}.{ORGANISM}.depth.gz

    logthis "${{yellow}}Depth calculation complete${{reset}}"
else
    logthis "Depth calculation already complete, ${{green}}skipping${{reset}}"
fi
    """.format(
            SAMPLE=sample,
            PIPELINE=pipeline,
            ORGANISM=reference["common"],
        )
    )


def commonPipeline(
    script: TextIOWrapper,
    options: OptionsDict,
    instrument: str,
    filenames,
    references,
    preprocess,
    align,
):
    preprocess(script, options)
    options["_instrument"] = instrument

    for requestReference in options["references"]:
        reference = references[requestReference]

        align(script, options, reference)
        sortAndExtractUnmapped(script, options, reference)
        runVariantPipeline(script, options, reference)

        if options["runQc"] == True:
            doQualityControl(script, options, reference)

    if options["runQc"] == True and options["skipFastQc"] == False:
        doFastQc(script, options, filenames)

    if options["runQc"] == True and options["doMultiQc"] == True:
        runMultiQC(script, options)
