#!/usr/bin/env bash

set -e

SAMPLE=
RUNQC=false
REMOVE_INTERMEDIATES=false
SKIP_PROCESSING=false
TRIM=true

THREADS=$(grep 'processor' /proc/cpuinfo | wc -l)
WORKING=$HOME
cd ${WORKING}/pipeline

usage() {
    echo "usage: $0 -s <sample> [-n] [-q] [-c] [-l]"
    echo " -s short name of sample, e.g. DPZw_k"
    echo "       file must be in ${WORKING}/pipeline/<sample>_R[12].fastq.gz"
    echo " -l skip trimming (via TrimGalore), default false"
    echo " -n skip processing data, default false"
    echo " -c clean up intermediate files, default false"
    echo " -q run QC process without output in ${WORKING}/stats, default false"
}

if [ "${1}" == "--help" ]; then
    usage
    exit 0
fi

while getopts "s:nqct" opt; do
    case $opt in
    s)
        SAMPLE=$OPTARG
        ;;
    n)
        SKIP_PROCESSING=true
        ;;
    q)
        RUNQC=true
        ;;
    t)
        TRIM=false
        ;;
    c)
        REMOVE_INTERMEDIATES=true
        ;;
    :)
        echo "Missing argument for ${OPTARG}"
        usage
        exit 1
        ;;
    esac
done

echo "SAMPLE=$SAMPLE"
echo "REMOVE_INTERMEDIATES=$REMOVE_INTERMEDIATES"
echo "RUNQC=$RUNQC"
echo "SKIP_PROCESSING=$SKIP_PROCESSING"

if [ -z "${SAMPLE}" ]; then
    echo "Sample name is a required parameter and must be present in ${WORKING}/pipeline"
    usage
    exit 2
fi

if [ ! -f ${WORKING}/pipeline/${SAMPLE}_R1.fastq.gz ]; then
    echo "Missing sample R1 in ${WORKING}/pipeline/${SAMPLE}_R1.fastq.gz"
    exit 3
fi

if [ ! -f ${WORKING}/pipeline/${SAMPLE}_R2.fastq.gz ]; then
    echo "Missing sample R1 in ${WORKING}/pipeline/${SAMPLE}_R2.fastq.gz"
    exit 3
fi

SORTED=${WORKING}/pipeline/${SAMPLE}.sorted.bam
CHR=${SORTED%.bam}
UNSORTED=${WORKING}/pipeline/${SAMPLE}
VCF=$(basename $UNSORTED)

# comment this out to disable
GATK_OPTIONS="--gatk-config-file ${WORKING}/gatk.properties"
export PATH=${WORKING}/bin/FastQC:${WORKING}/bin/TrimGalore-0.6.7:${WORKING}/bin/gatk-4.2.3.0:${WORKING}/bin:$PATH

if [ "${SKIP_PROCESSING}" == "false" ]; then
    R1=${WORKING}/pipeline/${SAMPLE}_R1.fastq.gz
    R2=${WORKING}/pipeline/${SAMPLE}_R2.fastq.gz

    if [ "${TRIM}" == "true" ]; then
        if [ "${RUNQC}" == "true" ]; then
            trim_galore \
                --illumina \
                --cores 72 \
                --output_dir ${WORKING}/pipeline \
                --fastqc_args "--outdir ${WORKING}/stats --noextract" \
                ${R1} ${R2}
        else
            trim_galore --illumina --cores 72 --output_dir ${WORKING}/pipeline ${R1} ${R2}
        fi

        # trimming generates new output files, we'll use those
        R1=${WORKING}/pipeline/${SAMPLE}_R1_trimmed.fq.gz
        R2=${WORKING}/pipeline/${SAMPLE}_R2_trimmed.fq.gz
    fi

    ${WORKING}/bin/bwa-mem2 mem -t ${THREADS} \
        ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
        ${R1} ${R2} \
        -Y \
        -R "@RG\tID:${SAMPLE}\tPL:ILLUMINA\tPU:MJS.SEQUENCER.7\tLB:${SAMPLE}\tSM:${SAMPLE}" |
        ${WORKING}/bin/bamsormadup \
            SO=coordinate \
            threads=${THREADS} \
            level=6 \
            inputformat=sam \
            indexfilename=${SORTED}.bai \
            M=${WORKING}/stats/${SAMPLE}.duplication_metrics >${SORTED}

    for i in {1..22} X Y M; do
        (
            # split the region from the aligned file
            ${WORKING}/bin/samtools view -@ 4 -bh ${SORTED} chr${i} >${CHR}_chr${i}.bam
            ${WORKING}/bin/samtools index -@ 4 ${CHR}_chr${i}.bam

            # run base quality score recalibration - build the bqsr table
            ${WORKING}/bin/gatk-4.2.3.0/gatk BaseRecalibrator ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -I ${CHR}_chr${i}.bam \
                -O ${CHR}_bqsr.chr${i}.table \
                --preserve-qscores-less-than 6 \
                --known-sites ${WORKING}/reference/Homo_sapiens_assembly38.dbsnp138.vcf \
                --known-sites ${WORKING}/reference/Homo_sapiens_assembly38.known_indels.vcf \
                --known-sites ${WORKING}/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf \
                -L chr${i}

            # run base quality score recalibration - apply the bqsr table
            ${WORKING}/bin/gatk-4.2.3.0/gatk ApplyBQSR ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -I ${CHR}_chr${i}.bam \
                -O ${CHR}_bqsr.chr${i}.bam \
                --preserve-qscores-less-than 6 \
                --static-quantized-quals 10 \
                --static-quantized-quals 20 \
                --static-quantized-quals 30 \
                --bqsr-recal-file ${CHR}_bqsr.chr${i}.table \
                -L chr${i}

            # index that file, this is our target for getting vcf
            ${WORKING}/bin/samtools index -@ 4 ${CHR}_bqsr.chr${i}.bam

            ${WORKING}/bin/gatk-4.2.3.0/gatk HaplotypeCaller ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -I ${CHR}_bqsr.chr${i}.bam \
                -O ${CHR}_variants_chr${i}.vcf \
                --bam-output ${CHR}_variants_chr${i}.bam \
                --pairHMM AVX_LOGLESS_CACHING_OMP \
                --native-pair-hmm-threads 8 \
                -L chr${i}

            # index that file, this is our target for getting vcf
            ${WORKING}/bin/samtools index -@ 4 ${CHR}_variants_chr${i}.bam

            # pull snps out of out called variants and annotate them
            ${WORKING}/bin/gatk-4.2.3.0/gatk SelectVariants ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -V ${CHR}_variants_chr${i}.vcf \
                -select-type SNP \
                -O ${CHR}_variants_chr${i}.snps.vcf

            ${WORKING}/bin/gatk-4.2.3.0/gatk VariantFiltration ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -V ${CHR}_variants_chr${i}.snps.vcf \
                --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
                --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
                --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
                --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
                --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
                -O ${CHR}_variants_chr${i}.snps.filtered.vcf

            rm -f ${WORKING}/pipeline/${VCF}.annotated_variants_chr${i}.snps.filtered.vcf
            rm -f ${WORKING}/pipeline/${VCF}.annotated_variants_chr${i}.snps.filtered.vcf_summary.html

            sudo docker run \
                -v ${WORKING}/vep_data:/opt/vep/.vep:Z \
                -v ${WORKING}/pipeline:/opt/vep/.vep/input:Z \
                -v ${WORKING}/pipeline:/opt/vep/.vep/output:Z \
                -v ${WORKING}/reference:/opt/vep/.vep/reference:Z \
                ensemblorg/ensembl-vep \
                ./vep --cache --format vcf --merged --offline --use_given_ref --vcf --verbose \
                --fasta /opt/vep/.vep/reference/Homo_sapiens_assembly38.fasta \
                --input_file /opt/vep/.vep/input/${VCF}.sorted_variants_chr${i}.snps.filtered.vcf \
                --output_file /opt/vep/.vep/output/${VCF}.annotated_variants_chr${i}.snps.filtered.vcf

            # pull indels out of out called variants and annotate them
            ${WORKING}/bin/gatk-4.2.3.0/gatk SelectVariants ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -V ${CHR}_variants_chr${i}.vcf \
                -select-type INDEL \
                -O ${CHR}_variants_chr${i}.indels.vcf

            ${WORKING}/bin/gatk-4.2.3.0/gatk VariantFiltration ${GATK_OPTIONS} \
                -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
                -V ${CHR}_variants_chr${i}.indels.vcf \
                --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
                --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \
                --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS_lt_n20" \
                -O ${CHR}_variants_chr${i}.indels.filtered.vcf

            rm -f ${WORKING}/pipeline/${VCF}.annotated_variants_chr${i}.indels.filtered.vcf
            rm -f ${WORKING}/pipeline/${VCF}.annotated_variants_chr${i}.indels.filtered.vcf_summary.html

            sudo docker run \
                -v ${WORKING}/vep_data:/opt/vep/.vep:Z \
                -v ${WORKING}/pipeline:/opt/vep/.vep/input:Z \
                -v ${WORKING}/pipeline:/opt/vep/.vep/output:Z \
                -v ${WORKING}/reference:/opt/vep/.vep/reference:Z \
                ensemblorg/ensembl-vep \
                ./vep --cache --format vcf --merged --offline --use_given_ref --vcf --verbose \
                --fasta /opt/vep/.vep/reference/Homo_sapiens_assembly38.fasta \
                --input_file /opt/vep/.vep/input/${VCF}.sorted_variants_chr${i}.indels.filtered.vcf \
                --output_file /opt/vep/.vep/output/${VCF}.annotated_variants_chr${i}.indels.filtered.vcf
        ) &
    done
    wait

    /bin/ls -1 ${WORKING}/pipeline/${VCF}.annotated_variants_chr[0-9XYM]*.indels.filtered.vcf >${WORKING}/pipeline/merge.indels.list
    /bin/ls -1 ${WORKING}/pipeline/${VCF}.annotated_variants_chr[0-9XYM]*.snps.filtered.vcf >${WORKING}/pipeline/merge.snps.list

    ${WORKING}/bin/gatk-4.2.3.0/gatk MergeVcfs \
        -I ${WORKING}/pipeline/merge.indels.list \
        -O ${WORKING}/pipeline/${SAMPLE}.indels.final.vcf &

    ${WORKING}/bin/gatk-4.2.3.0/gatk MergeVcfs \
        -I ${WORKING}/pipeline/merge.snps.list \
        -O ${WORKING}/pipeline/${SAMPLE}.snps.final.vcf &
    wait

    ${WORKING}/bin/gatk-4.2.3.0/gatk MergeVcfs \
        -I ${WORKING}/pipeline/${SAMPLE}.snps.final.vcf \
        -I ${WORKING}/pipeline/${SAMPLE}.indels.final.vcf \
        -O ${WORKING}/pipeline/${SAMPLE}.final.unfiltered.vcf

    ${WORKING}/bin/gatk-4.2.3.0/gatk SelectVariants ${GATK_OPTIONS} \
        -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
        -V ${WORKING}/pipeline/${SAMPLE}.final.unfiltered.vcf \
        -O ${WORKING}/pipeline/${SAMPLE}.final.filtered.vcf \
        --exclude-filtered

    ${WORKING}/bin/samtools merge -@ 72 -o ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam ${WORKING}/pipeline/*bqsr*.bam
    ${WORKING}/bin/samtools index -@ 72 ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam
fi

if [ "${RUNQC}" == "true" ]; then
    ${WORKING}/bin/samtools flagstat \
        ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam >${WORKING}/stats/${SAMPLE}.bqsr.flagstat.txt &

    ${WORKING}/bin/gatk-4.2.3.0/gatk CollectInsertSizeMetrics \
        -I ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam \
        -O ${WORKING}/stats/${SAMPLE}.bqsr.insert_metrics.txt \
        -H ${WORKING}/stats/${SAMPLE}.bqsr.insert_metrics.pdf \
        -M 0.5 &

    ${WORKING}/bin/gatk-4.2.3.0/gatk CollectAlignmentSummaryMetrics \
        -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
        -I ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam \
        -O ${WORKING}/stats/${SAMPLE}.bqsr.alignment_metrics.txt &

    ${WORKING}/bin/gatk-4.2.3.0/gatk CollectGcBiasMetrics \
        -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
        -I ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam \
        -O ${WORKING}/stats/${SAMPLE}.bqsr.gc_bias_metrics.txt \
        -CHART ${WORKING}/stats/${SAMPLE}.bqsr.gc_bias_metrics.pdf \
        -S ${WORKING}/stats/${SAMPLE}.bqsr.gc_bias_summary.txt &

    ${WORKING}/bin/gatk-4.2.3.0/gatk CollectWgsMetrics \
        -R ${WORKING}/reference/Homo_sapiens_assembly38.fasta \
        -I ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam \
        -O ${WORKING}/stats/${SAMPLE}.metrics.txt \
        --READ_LENGTH 151 \
        -INTERVALS ${WORKING}/reference/ref_genome_autosomal.interval_list \
        --USE_FAST_ALGORITHM \
        --INCLUDE_BQ_HISTOGRAM &

    ${WORKING}/bin/FastQC/fastqc \
        --threads=24 \
        --outdir ${WORKING}/stats \
        ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam &

    if [ "${TRIM}" == "false" ]; then
        ${WORKING}/bin/FastQC/fastqc \
            --threads=24 \
            --outdir ${WORKING}/stats \
            --noextract \
            ${WORKING}/pipeline/${SAMPLE}_R1.fastq.gz &

        ${WORKING}/bin/FastQC/fastqc \
            --threads=24 \
            --outdir ${WORKING}/stats \
            --noextract \
            ${WORKING}/pipeline/${SAMPLE}_R2.fastq.gz &
    fi
    wait

    mkdir -p ${WORKING}/stats/qc
    cd ${WORKING}/stats/qc
    multiqc -f ${WORKING}/stats
fi

if [ "${REMOVE_INTERMEDIATES}" == "true" ]; then
    # clean up intermediate bam and index files
    rm -fv ${WORKING}/pipeline/${SAMPLE}.sorted.bam
    rm -fv ${WORKING}/pipeline/${SAMPLE}.sorted.bam.bai

    rm -fv ${WORKING}/pipeline/${SAMPLE}.bam
    rm -fv ${WORKING}/pipeline/${SAMPLE}.bam.bai

    rm -fv ${WORKING}/pipeline/*table
    rm -fv ${WORKING}/pipeline/${SAMPLE}.annotated_variants_chr*
    rm -fv ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam
    rm -fv ${WORKING}/pipeline/${SAMPLE}.merged.bqsr.bam.bai
    rm -fv ${WORKING}/pipeline/${SAMPLE}.sorted_bqsr*
    rm -fv ${WORKING}/pipeline/${SAMPLE}.sorted_chr*
    rm -fv ${WORKING}/pipeline/${SAMPLE}.sorted_variants_chr*

    rm -fv ${WORKING}/pipeline/${SAMPLE}.snps.*
    rm -fv ${WORKING}/pipeline/${SAMPLE}.indels.*
fi
