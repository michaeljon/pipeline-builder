from io import TextIOWrapper

from bio_types import *


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

    {BIN}/bamsormadup \\
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


def extractUmappedReads(script: TextIOWrapper, options: OptionsDict):
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# extract unmapped reads
#
if [[ ! -f {PIPELINE}/{SAMPLE}_unmapped.fastq ]]; then
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


def sortAndExtractUnmapped(script: TextIOWrapper, options: OptionsDict):
    processUnmapped = options["processUnmapped"]
    alignOnly = options["alignOnly"]

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


def callVariants(script: TextIOWrapper, options: OptionsDict):
    sample = options["sample"]
    pipeline = options["pipeline"]
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    bin = options["bin"]

    script.write(
        """
# call variants
if [[ ! -f {PIPELINE}/{SAMPLE}.unannotated.vcf.gz ]]; then
    logthis "${{yellow}}Calling variants using freebayes${{reset}}"

    # calling variants using freebays
    freebayes \\
        --fasta-reference {REFERENCE}/{ASSEMBLY}.fna \\
        --ploidy 1 \\
        {PIPELINE}/{SAMPLE}.sorted.bam >{PIPELINE}/{SAMPLE}.tmp.vcf

    logthis "Normalizing called variants"
    bcftools norm \\
        --atomize \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        -o {PIPELINE}/{SAMPLE}.norm.bcf \\
        -Ou \\
        {PIPELINE}/{SAMPLE}.tmp.vcf

    logthis "Filtering called variants"
    bcftools filter -i "QUAL > 10" \\
        -o {PIPELINE}/{SAMPLE}.filtered.bcf \\
        -Ou \\
        {PIPELINE}/{SAMPLE}.norm.bcf

    logthis "Indexing normalized BCF and calling consensus"
    bcftools index --force {PIPELINE}/{SAMPLE}.filtered.bcf

    bcftools consensus \\
        --sample {SAMPLE} \\
        --fasta-ref {REFERENCE}/{ASSEMBLY}.fna \\
        -o {PIPELINE}/{SAMPLE}.consensus.fa \\
        {PIPELINE}/{SAMPLE}.filtered.bcf

    logthis "Converting freebayes temporary BCF to VCF"
    bcftools convert \\
        -Ov \\
        {PIPELINE}/{SAMPLE}.filtered.bcf \\
        -o {PIPELINE}/{SAMPLE}.unannotated.vcf

    bgzip --force --keep {PIPELINE}/{SAMPLE}.unannotated.vcf
    tabix --force -p vcf {PIPELINE}/{SAMPLE}.unannotated.vcf.gz

    logthis "Cleaning up temporary freebayes output"
    rm -f {PIPELINE}/{SAMPLE}.tmp.vcf
    rm -f {PIPELINE}/{SAMPLE}.norm.bcf
    rm -f {PIPELINE}/{SAMPLE}.filtered.bcf
    rm -f {PIPELINE}/{SAMPLE}.filtered.bcf.csi

    logthis "${{green}}freebayes variant calling completed${{reset}}"
else
    logthis "Variants already called via bcftools for {PIPELINE}/{SAMPLE}.sorted.bam, ${{green}}skipping${{reset}}"
fi
""".format(
            BIN=bin,
            REFERENCE=reference,
            ASSEMBLY=assembly,
            PIPELINE=pipeline,
            SAMPLE=sample,
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

    callVariants(script, options)
    if options["__canAnnotateVariants"] == True:
        annotate(script, options)
    if options["__canAssignClades"] == True:
        assignClade(script, options)

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


def doAlignmentQC(script: TextIOWrapper, filenames, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]
    threads = options["cores"]
    skipFastQc = options["skipFastQc"]
    doPicardQc = options["doPicardQc"]
    bin = options["bin"]

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
            """'if [[ ! -f {STATS}/{SAMPLE}.alignment_metrics.txt ]]; then java -jar {BIN}/picard.jar CollectAlignmentSummaryMetrics --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.alignment_metrics.txt; fi' \\\n""".format(
                REFERENCE=reference,
                ASSEMBLY=assembly,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
                SORTED=sorted,
                BIN=bin,
            )
        )

        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}.gc_bias_metrics.txt || ! -f {STATS}/{SAMPLE}.gc_bias_metrics.pdf || ! -f {STATS}/{SAMPLE}.gc_bias_summary.txt ]]; then java -jar {BIN}/picard.jar CollectGcBiasMetrics --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.gc_bias_metrics.txt -CHART {STATS}/{SAMPLE}.gc_bias_metrics.pdf -S {STATS}/{SAMPLE}.gc_bias_summary.txt; fi' \\\n""".format(
                REFERENCE=reference,
                ASSEMBLY=assembly,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
                SORTED=sorted,
                BIN=bin,
            )
        )

        checks.append(
            """'if [[ ! -f {STATS}/{SAMPLE}.wgs_metrics.txt ]]; then java -jar {BIN}/picard.jar CollectWgsMetrics --VERBOSITY ERROR -R {REFERENCE}/{ASSEMBLY}.fna -I {SORTED} -O {STATS}/{SAMPLE}.wgs_metrics.txt --MINIMUM_BASE_QUALITY 20 --MINIMUM_MAPPING_QUALITY 20 --COVERAGE_CAP 10000 --READ_LENGTH 151 --INTERVALS {REFERENCE}/{ASSEMBLY}_autosomal.interval_list --USE_FAST_ALGORITHM --INCLUDE_BQ_HISTOGRAM; fi' \\\n""".format(
                REFERENCE=reference,
                ASSEMBLY=assembly,
                PIPELINE=pipeline,
                SAMPLE=sample,
                STATS=stats,
                THREADS=threads,
                SORTED=sorted,
                BIN=bin,
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
        """'if [[ ! -f {STATS}/{SAMPLE}.bedtools.coverage ]]; then samtools view -bq 30 -F 1284 {SORTED} | bedtools genomecov -d -ibam stdin | awk "\\$2 % 100 == 0 {{print \\$1,\\$2,\\$3}}" | sed "s/AF304460.1/hcov_229e/;s/AY597011.2/hcov_hku1/;s/AY567487.2/hcov_nl63/;s/AY585228.1/hcov_oc43/;s/MN908947.3/sars_cov_2/;s/NC_038311.1/hrv_a/;s/NC_038312.1/hrv_b/;s/NC_038878.1/hrv_c/" >{STATS}/{SAMPLE}.bedtools.coverage; fi' \\\n""".format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            STATS=stats,
            SORTED=sorted,
        )
    )

    if skipFastQc == False:
        if len(filenames) == 2:
            checks.append(
                """'if [[ ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R1.trimmed_fastqc.html || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.zip || ! -f {STATS}/{SAMPLE}_R2.trimmed_fastqc.html ]]; then fastqc --svg --threads 2 --outdir {STATS} --noextract {R1} {R2}; fi' \\\n""".format(
                    SAMPLE=sample,
                    STATS=stats,
                    R1=filenames[0],
                    R2=filenames[1],
                )
            )
        else:
            checks.append(
                """'if [[ ! -f {STATS}/{SAMPLE}.fastqc.zip ]]; then fastqc --svg --threads 1 --outdir {STATS} --noextract {FILENAME}; fi' \\\n""".format(
                    SAMPLE=sample,
                    STATS=stats,
                    FILENAME=filenames[0],
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


def doQualityControl(script: TextIOWrapper, options: OptionsDict, filenames: List[str]):
    pipeline = options["pipeline"]
    sample = options["sample"]
    stats = options["stats"]

    plot_vcfstats = """plot-vcfstats --prefix {STATS}/{SAMPLE}_bcfstats {STATS}/{SAMPLE}.chk""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    # this doesn't have a test, it's fast enough that we can afford to run it
    plot_bamstats = """plot-bamstats --prefix {STATS}/{SAMPLE}_samstats/ {STATS}/{SAMPLE}.samstats""".format(
        SAMPLE=sample,
        STATS=stats,
    )

    alignment_checks = doAlignmentQC(script, filenames, options)
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


def generateDepth(script: TextIOWrapper, options: OptionsDict):
    reference = options["reference"]
    assembly = options["referenceAssembly"]
    pipeline = options["pipeline"]
    sample = options["sample"]

    script.write(
        """
#
# calculate depth by position
#
if [[ ! -f {PIPELINE}/{SAMPLE}.depth.gz ]]; then
    logthis "${{yellow}}Calculating depth by position${{reset}}"

    samtools depth -@ 8 -aa -a -J --reference {REFERENCE}/{ASSEMBLY}.fna {PIPELINE}/{SAMPLE}.sorted.bam | gzip >{PIPELINE}/{SAMPLE}.depth.gz

    logthis "${{yellow}}Depth calculation complete${{reset}}"
else
    logthis "Depth calculation already complete, ${{green}}skipping${{reset}}"
fi
    """.format(
            REFERENCE=reference,
            ASSEMBLY=assembly,
            SAMPLE=sample,
            PIPELINE=pipeline,
        )
    )


def commonPipeline(
    script: TextIOWrapper,
    options: OptionsDict,
    filenames,
    preprocessAndAlign,
):
    preprocessAndAlign(script, options)
    sortAndExtractUnmapped(script, options)
    runVariantPipeline(script, options)

    # we'll wait here to make sure all the background stuff is done before we
    # run multiqc and cleanup
    script.write('logthis "${green}Done with front-end processing${reset}"\n')

    if options["runQc"] == True:
        doQualityControl(script, options, filenames)
        script.write('logthis "${green}Done with back-end processing${reset}"\n')
