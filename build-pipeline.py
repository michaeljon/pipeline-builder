import json
import optparse
from datetime import datetime
from math import ceil
from os.path import exists, expandvars

segmentSize = 50_000_000
lastBlockMax = segmentSize * 0.25


def loadIntervals():
    with open("chromosomeSizes.json", "r") as file:
        chromosomeSizes = json.load(file)
        return chromosomeSizes


def computeIntervals():
    intervals = []

    with open("chromosomeSizes.json", "r") as file:
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


def getFileNames(sample, pipeline, trimmed):
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


def genTrimmer(script, r1, r2, pipeline, stats, runQC):
    script.write("#\n")
    script.write("# Trim reads\n")
    script.write("#\n")

    if runQC == True:
        script.write(
            """trim_galore \\
    --illumina \\
    --cores 72 \\
    --output_dir {PIPELINE} \\
    --fastqc_args "--outdir {STATS} --noextract" \\
    {R1} \\
    {R2}
""".format(
                R1=r1, R2=r2, PIPELINE=pipeline, STATS=stats
            )
        )
    else:
        script.write(
            """trim_galore \\
    --illumina \\
    --cores 72 \\
    --output_dir {PIPELINE} \\
    {R1} \\
    {R2}
""".format(
                R1=r1, R2=r2, PIPELINE=pipeline
            )
        )


def genBWA(script, r1, r2, reference, sample, working, stats, threads, output):
    script.write("#\n")
    script.write("# Align, sort, and mark duplicates\n")
    script.write("#\n")

    script.write(
        """bwa-mem2 mem -t {THREADS} \\
    {REFERENCE}/Homo_sapiens_assembly38.fasta \\
    {R1} \\
    {R2} \\
    -Y \\
    -R "@RG\\tID:{SAMPLE}\\tPL:ILLUMINA\\tPU:MJS.SEQUENCER.7\\tLB:{SAMPLE}\\tSM:{SAMPLE}" |
bamsormadup \\
    SO=coordinate \\
    threads={THREADS} \\
    level=6 \\
    inputformat=sam \\
    indexfilename={SORTED}.bai \\
    M={STATS}/{SAMPLE}.duplication_metrics >{SORTED}

""".format(
            R1=r1,
            R2=r2,
            REFERENCE=reference,
            SAMPLE=sample,
            WORKING=working,
            STATS=stats,
            THREADS=threads,
            SORTED=output,
        )
    )


def splinter(script, bam, sorted, interval):
    script.write("        # Splinter our interval off the main file\n")
    script.write(
        """        samtools view -@ 4 -bh {SORTED} {INTERVAL} >{BAM}\n""".format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )
    script.write(
        """        samtools index -@ 4 {BAM}\n""".format(
            SORTED=sorted, INTERVAL=interval, BAM=bam
        )
    )


def genBQSR(script, reference, interval, bam, bqsr):
    script.write(
        """
        # run base quality score recalibration - build the bqsr table
        gatk BaseRecalibrator \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BAM} \\
            -O {BQSR}.table \\
            --preserve-qscores-less-than 6 \\
            --known-sites {REFERENCE}/Homo_sapiens_assembly38.dbsnp138.vcf \\
            --known-sites {REFERENCE}/Homo_sapiens_assembly38.known_indels.vcf \\
            --known-sites {REFERENCE}/Mills_and_1000G_gold_standard.indels.hg38.vcf \\
            -L {INTERVAL}

        # run base quality score recalibration - apply the bqsr table
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
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr
        )
    )


def callVariants(script, reference, interval, bam, bqsr, vcf):
    script.write(
        """
        # call variants
        gatk HaplotypeCaller \\
            -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
            -I {BQSR} \\
            -O {VCF} \\
            --pairHMM AVX_LOGLESS_CACHING_OMP \\
            --native-pair-hmm-threads 8 \\
            -L {INTERVAL}
""".format(
            REFERENCE=reference, INTERVAL=interval, BAM=bam, BQSR=bqsr, VCF=vcf
        )
    )


def filterSNPs(script, reference, vcf, snps, filtered):
    script.write(
        """
        # pull snps out of out called variants and annotate them
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
""".format(
            REFERENCE=reference, VCF=vcf, SNPS=snps, FILTERED=filtered
        )
    )


def filterINDELs(script, reference, vcf, indels, filtered):
    script.write(
        """
        # pull indels out of out called variants and annotate them
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
""".format(
            REFERENCE=reference, VCF=vcf, INDELS=indels, FILTERED=filtered
        )
    )


def filterVariants(script, reference, vcf):
    snps = vcf.replace(".vcf", ".snps.vcf")
    filtered = vcf.replace(".vcf", ".snps.filtered.vcf")
    filterSNPs(script, reference, vcf, snps, filtered)

    indels = vcf.replace(".vcf", ".indels.vcf")
    filtered = vcf.replace(".vcf", ".indels.filtered.vcf")
    filterINDELs(script, reference, vcf, indels, filtered)


def annotate(script, reference, vep, pipeline, vcf, filtered, annotated, summary):
    input = filtered.split("/")[-1]
    output = annotated.split("/")[-1]

    script.write(
        """
        # remove these if they're leftover from a previous run
        rm -f {ANNOTATED}
        rm -f {SUMMARY}

        # run vep in a docker container to annotate the vcf
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
""".format(
            VEP=vep,
            PIPELINE=pipeline,
            REFERENCE=reference,
            VCF=vcf,
            ANNOTATED=annotated,
            FILTERED=filtered,
            SUMMARY=summary,
            INPUT=input,
            OUTPUT=output,
        )
    )


def annotateVariants(script, reference, working, pipeline, vcf):
    vep = "{WORKING}/vep_data".format(WORKING=working)

    filtered = vcf.replace(".vcf", ".snps.filtered.vcf")
    annotated = vcf.replace(".vcf", ".snps.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".snps.filtered.vcf_summary.html")
    annotate(script, reference, vep, pipeline, vcf, filtered, annotated, summary)

    filtered = vcf.replace(".vcf", ".indels.filtered.vcf")
    annotated = vcf.replace(".vcf", ".indels.filtered.annotated.vcf")
    summary = vcf.replace(".vcf", ".indels.filtered.vcf_summary.html")
    annotate(script, reference, vep, pipeline, vcf, filtered, annotated, summary)


def runIntervals(script, working, reference, pipeline, prefix, sorted):
    intervals = computeIntervals()

    for interval in intervals:
        bam = """{PREFIX}.{INTERVAL}.bam""".format(PREFIX=prefix, INTERVAL=interval[1])
        bqsr = """{PREFIX}.{INTERVAL}_bqsr.bam""".format(
            PREFIX=prefix, INTERVAL=interval[1]
        )
        vcf = """{PREFIX}.{INTERVAL}.vcf""".format(PREFIX=prefix, INTERVAL=interval[1])

        script.write("    #\n")
        script.write("    # Run interval {INTERVAL}\n".format(INTERVAL=interval[0]))
        script.write("    #\n")
        script.write("    (\n")

        splinter(script, bam, sorted, interval[0])
        genBQSR(script, reference, interval[0], bam, bqsr)
        callVariants(script, reference, interval[0], bam, bqsr, vcf)
        filterVariants(script, reference, vcf)
        annotateVariants(script, reference, working, pipeline, vcf)

        script.write("    )&\n")
        script.write("\n")

    script.write("wait\n")


def merge(script, reference, pipeline, sample):
    script.write(
        """
#
# Create final VCF file(s)
# 
/bin/ls -1 {PIPELINE}/*.indels.filtered.annotated.vcf >{PIPELINE}/merge.indels.list
/bin/ls -1 {PIPELINE}/*.snps.filtered.annotated.vcf >{PIPELINE}/merge.snps.list

gatk MergeVcfs \\
    -I {PIPELINE}/merge.indels.list \\
    -O {PIPELINE}/{SAMPLE}.indels.final.vcf &

gatk MergeVcfs \\
    -I {PIPELINE}/merge.snps.list \\
    -O {PIPELINE}/{SAMPLE}.snps.final.vcf &
wait

gatk MergeVcfs \\
    -I {PIPELINE}/{SAMPLE}.snps.final.vcf \\
    -I {PIPELINE}/{SAMPLE}.indels.final.vcf \\
    -O {PIPELINE}/{SAMPLE}.final.unfiltered.vcf

gatk SelectVariants \\
    -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
    -V {PIPELINE}/{SAMPLE}.final.unfiltered.vcf \\
    -O {PIPELINE}/{SAMPLE}.final.filtered.vcf \\
    --exclude-filtered
""".format(
            REFERENCE=reference, PIPELINE=pipeline, SAMPLE=sample
        )
    )


def runQC(script, reference, pipeline, sample, stats, threads):
    script.write(
        """
#
# RUN QC process
# 
samtools merge -@ 72 -o {PIPELINE}/{SAMPLE}.merged.bqsr.bam {PIPELINE}/*bqsr*.bam
samtools index -@ 72 {PIPELINE}/{SAMPLE}.merged.bqsr.bam

samtools flagstat \\
    {PIPELINE}/{SAMPLE}.merged.bqsr.bam >{STATS}/{SAMPLE}.bqsr.flagstat.txt &

gatk CollectInsertSizeMetrics \\
    -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    -O {STATS}/{SAMPLE}.bqsr.insert_metrics.txt \\
    -H {STATS}/{SAMPLE}.bqsr.insert_metrics.pdf \\
    -M 0.5 &

gatk CollectAlignmentSummaryMetrics \\
    -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
    -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    -O {STATS}/{SAMPLE}.bqsr.alignment_metrics.txt &

gatk CollectGcBiasMetrics \\
    -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
    -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    -O {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.txt \\
    -CHART {STATS}/{SAMPLE}.bqsr.gc_bias_metrics.pdf \\
    -S {STATS}/{SAMPLE}.bqsr.gc_bias_summary.txt &

gatk CollectWgsMetrics \\
    -R {REFERENCE}/Homo_sapiens_assembly38.fasta \\
    -I {PIPELINE}/{SAMPLE}.merged.bqsr.bam \\
    -O {STATS}/{SAMPLE}.metrics.txt \\
    --READ_LENGTH 151 \\
    -INTERVALS {REFERENCE}/ref_genome_autosomal.interval_list \\
    --USE_FAST_ALGORITHM \\
    --INCLUDE_BQ_HISTOGRAM &

fastqc \\
    --threads={THREADS} \\
    --outdir {STATS} \\
    {PIPELINE}/{SAMPLE}.merged.bqsr.bam &
wait
""".format(
            REFERENCE=reference,
            PIPELINE=pipeline,
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
        )
    )


def runInputQC(script, pipeline, sample, stats, threads):
    filenames = getFileNames(sample, pipeline, False)

    script.write(
        """
#
# RUN QC process on input files
# 
fastqc \\
    --threads={THREADS} \\
    --outdir {STATS} \\
    --noextract \\
    {R1} &

fastqc \\
    --threads={THREADS} \\
    --outdir {STATS} \\
    --noextract \\
    {R2} &
wait
""".format(
            R1=filenames[0],
            R2=filenames[1],
            SAMPLE=sample,
            STATS=stats,
            THREADS=threads,
        )
    )


def runMultiQC(script, stats):
    script.write(
        """
#
# Run MultiQC across everything
# 
(
    mkdir -p {STATS}/qc
    cd {STATS}/qc
    multiqc -f {STATS}
)    

""".format(
            STATS=stats
        )
    )


def cleanup(script, prefix, pipeline, sample, byInterval):
    intervals = computeIntervals()

    script.write("#\n")
    script.write("# Clean up all intermediate interval files\n")
    script.write("#\n")

    if byInterval == True:
        for interval in intervals:
            bam = """{PREFIX}.{INTERVAL}.bam""".format(
                PREFIX=prefix, INTERVAL=interval[1]
            )
            bqsr = """{PREFIX}.{INTERVAL}_bqsr.bam""".format(
                PREFIX=prefix, INTERVAL=interval[1]
            )
            vcf = """{PREFIX}.{INTERVAL}.vcf""".format(
                PREFIX=prefix, INTERVAL=interval[1]
            )

            script.write("rm -f {BAM}\n".format(BAM=bam))
            script.write("rm -f {BAM}.bai\n".format(BAM=bam))

            script.write("\n")
            script.write("rm -f {BQSR}\n".format(BQSR=bqsr))
            script.write("rm -f {BQSR}.bai\n".format(BQSR=bqsr))
            script.write("rm -f {BQSR}.table\n".format(BQSR=bqsr))

            script.write("\n")
            script.write("rm -f {VCF}\n".format(VCF=vcf))

            for type in ["snps", "indels"]:
                script.write("\n")
                script.write(
                    "rm -f {VCF}\n".format(
                        VCF=vcf.replace(".vcf", ".{TYPE}.vcf").format(TYPE=type)
                    )
                )
                script.write(
                    "rm -f {VCF}\n".format(
                        VCF=vcf.replace(".vcf", ".{TYPE}.filtered.vcf").format(
                            TYPE=type
                        )
                    )
                )
                script.write(
                    "rm -f {VCF}\n".format(
                        VCF=vcf.replace(
                            ".vcf", ".{TYPE}.filtered.annotated.vcf"
                        ).format(TYPE=type)
                    )
                )
                script.write(
                    "rm -f {VCF}\n".format(
                        VCF=vcf.replace(
                            ".vcf", ".{TYPE}.filtered.vcf_summary.html"
                        ).format(TYPE=type)
                    )
                )

            script.write("\n")
    else:
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
    script.write("# Assumed chromosome sizes\n")
    intervals = loadIntervals()
    for c in intervals.keys():
        script.write("#   {CHROME} = {SIZE}\n".format(CHROME=c, SIZE=intervals[c]))

    script.write("#\n")
    script.write("# Split parameters\n")
    script.write("#   segmentSize = {P}\n".format(P=segmentSize))
    script.write("#   lastBlockMax = {P}\n".format(P=lastBlockMax))


def writeVersions(script):
    script.write("#\n")
    script.write("# This will write version numbers of tools here...")
    script.write("#\n")


def main():
    parser = optparse.OptionParser()
    parser.set_defaults(doQC=False, trim=True, cleanIntermediateFiles=True)
    parser.add_option(
        "-s",
        "--sample",
        action="store",
        metavar="SAMPLE",
        dest="sample",
        help="short name of sample, e.g. DPZw_k file must be in <WORKING>/pipeline/<sample>_R[12].fastq.gz",
    )
    parser.add_option(
        "-q",
        "--skip-qc",
        action="store_false",
        dest="doQC",
        default=True,
        help="Skip QC process on input and output files",
    )
    parser.add_option(
        "-t",
        "--trim",
        action="store_true",
        dest="trim",
        default=False,
        help="Run trim-galore on the input FASTQ",
    )
    parser.add_option(
        "-c",
        "--clean",
        action="store_true",
        dest="cleanIntermediateFiles",
        default=False,
        help="Clean up the mess we make",
    )
    parser.add_option(
        "-w",
        "--work-dir",
        action="store",
        metavar="WORKING_DIR",
        dest="working",
        help="Working directory, e.g. base for $WORKING/pipeline, $WORKING/stats",
    )
    parser.add_option(
        "-r",
        "--reference-dir",
        action="store",
        metavar="REFERENCE_DIR",
        dest="reference",
        help="Location of the reference genome files",
    )
    parser.add_option(
        "-p",
        "--pipeline-dir",
        action="store",
        metavar="PIPELINE_DIR",
        dest="pipeline",
        help="Location of R1 and R2 files",
    )
    parser.add_option(
        "-o",
        "--stats-dir",
        action="store",
        metavar="STATS_DIR",
        dest="stats",
        help="Destination for statistics and QC files",
    )
    parser.add_option(
        "-b",
        "--bin-dir",
        action="store",
        metavar="BIN_DIR",
        dest="bin",
        default="$HOME/bin",
        help="Install location of all tooling",
    )
    parser.add_option(
        "-j",
        "--script",
        action="store",
        metavar="SHELL_SCRIPT",
        dest="script",
        default="pipeline-runner",
        help="Filename of bash shell to create",
    )

    (opts, args) = parser.parse_args()

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

    for opt in ["working", "reference", "pipeline", "stats", "bin", "script"]:
        options[opt] = expandvars(options[opt])

    filenames = getFileNames(options["sample"], options["pipeline"], False)

    if (
        exists(expandvars(filenames[0])) == False
        or exists(expandvars(filenames[1])) == False
    ):
        print(
            "Unable to locate the R1 or R2 files at {FILENAMES}".format(
                FILENAMES=filenames
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

        # assume we're not trimming, this gets the original R1/R2
        filenames = getFileNames(options["sample"], options["pipeline"], False)

        if options["trim"] == True:
            genTrimmer(
                script,
                filenames[0],
                filenames[1],
                options["pipeline"],
                "${HOME}/stats",
                options["doQC"],
            )

            # trimming writes new files under different names, this will grab those
            filenames = getFileNames(
                options["sample"], options["pipeline"], options["trim"]
            )

        genBWA(
            script,
            filenames[0],
            filenames[1],
            options["reference"],
            options["sample"],
            options["working"],
            options["stats"],
            72,
            sorted,
        )
        runIntervals(
            script,
            options["working"],
            options["reference"],
            options["pipeline"],
            prefix,
            sorted,
        )
        merge(script, options["reference"], options["pipeline"], options["sample"])

        if options["doQC"]:
            runQC(
                script,
                options["reference"],
                options["pipeline"],
                options["sample"],
                options["stats"],
                24,
            )

            runInputQC(
                script, options["pipeline"], options["sample"], options["stats"], 24
            )

            runMultiQC(script, options["stats"])

        if options["cleanIntermediateFiles"] == True:
            cleanup(script, prefix, options["pipeline"], options["sample"], False)

        script.write(
            """echo "Done processing {SAMPLE}\\n\\tstats in {STATS}\\n\\tVCFs in {PIPELINE}"
""".format(
                SAMPLE=options["sample"],
                STATS=options["stats"],
                PIPELINE=options["pipeline"],
            )
        )


if __name__ == "__main__":
    main()
