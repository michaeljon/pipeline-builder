# trim full FASTQ
Read1 before filtering:
total reads: 362785659
total bases: 54780634509
Q20 bases: 53100481029(96.9329%)
Q30 bases: 50523200788(92.2282%)

Read2 before filtering:
total reads: 362785659
total bases: 54780634509
Q20 bases: 52790550441(96.3672%)
Q30 bases: 49880186036(91.0544%)

Read1 after filtering:
total reads: 358596782
total bases: 53363494387
Q20 bases: 51934220305(97.3216%)
Q30 bases: 49503908066(92.7674%)

Read2 after filtering:
total reads: 358596782
total bases: 53361353861
Q20 bases: 51681222660(96.8514%)
Q30 bases: 48909638642(91.6574%)

Filtering result:
reads passed filter: 717193564
reads failed due to low quality: 8091020
reads failed due to too many N: 70080
reads failed due to too short: 216654
reads with adapter trimmed: 34922311
bases trimmed due to adapters: 1443700393

Duplication rate: 11.019%

Insert size peak (evaluated by paired-end reads): 271

JSON report: /home/ubuntu/pipeline/DPZw_k/DPZw_k-trim-fastp.json
HTML report: /home/ubuntu/pipeline/DPZw_k/DPZw_k-trim-fastp.html

fastp --in1 /home/ubuntu/stats/FASTQ/DPZw_k_R1.fastq.gz --in2 /home/ubuntu/stats/FASTQ/DPZw_k_R2.fastq.gz --out1=/home/ubuntu/pipeline/parts/DPZw_k_trimmed_R1.fastq.gz --out2=/home/ubuntu/pipeline/parts/DPZw_k_trimmed_R2.fastq.gz --verbose --thread 16 --detect_adapter_for_pe -j /home/ubuntu/pipeline/DPZw_k/DPZw_k-trim-fastp.json -h /home/ubuntu/pipeline/DPZw_k/DPZw_k-trim-fastp.html
fastp v0.23.2, time used: 1253 seconds

real	20m53.116s
user	202m53.450s
sys	    3m47.720s


# split FASTQ into 64 segments
Read1 before filtering:
total reads: 358596782
total bases: 53363494387
Q20 bases: 51934220305(97.3216%)
Q30 bases: 49503908066(92.7674%)

Read2 before filtering:
total reads: 358596782
total bases: 53361353861
Q20 bases: 51681222660(96.8514%)
Q30 bases: 48909638642(91.6574%)

Read1 after filtering:
total reads: 358596782
total bases: 53363494387
Q20 bases: 51934220305(97.3216%)
Q30 bases: 49503908066(92.7674%)

Read2 after filtering:
total reads: 358596782
total bases: 53361353861
Q20 bases: 51681222660(96.8514%)
Q30 bases: 48909638642(91.6574%)

Filtering result:
reads passed filter: 717193564
reads failed due to low quality: 0
reads failed due to too many N: 0

Duplication rate: 11.2502%

Insert size peak (evaluated by paired-end reads): 271

JSON report: /home/ubuntu/pipeline/DPZw_k/DPZw_k-split-fastp.json
HTML report: /home/ubuntu/pipeline/DPZw_k/DPZw_k-split-fastp.html

fastp --in1=/home/ubuntu/pipeline/parts/DPZw_k_trimmed_R1.fastq.gz --in2=/home/ubuntu/pipeline/parts/DPZw_k_trimmed_R2.fastq.gz --out1=/home/ubuntu/pipeline/parts/DPZw_k_R1.fastq.gz --out2=/home/ubuntu/pipeline/parts/DPZw_k_R2.fastq.gz --verbose --thread 16 --split 64 --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering -j /home/ubuntu/pipeline/DPZw_k/DPZw_k-split-fastp.json -h /home/ubuntu/pipeline/DPZw_k/DPZw_k-split-fastp.html
fastp v0.23.2, time used: 533 seconds

real	8m52.947s
user	116m35.489s
sys	    1m43.029s

# align segnments (representative per segment)
No. of OMP threads: 72
Processor is running @3000.478564 MHz
Runtime profile:

	Time taken for main_mem function: 88.10 sec

	IO times (sec) :
	Reading IO time (reads) avg: 46.27, (46.27, 46.27)
	Writing IO time (SAM) avg: 17.21, (17.21, 17.21)
	Reading IO time (Reference Genome) avg: 4.40, (4.40, 4.40)
	Index read time avg: 6.17, (6.17, 6.17)

	Overall time (sec) (Excluding Index reading time):
	PROCESS() (Total compute time + (read + SAM) IO time) : 76.68
	MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 56.96, (56.96, 56.96)

	 SAM Processing time (sec):
	--WORKER_SAM avg: 13.14, (13.14, 13.14)

	Kernels' compute time (sec):
	Total kernel (smem+sal+bsw) time avg: 42.50, (42.50, 42.50)
		SMEM compute avg: 14.23, (14.58, 13.88)
		SAL compute avg: 8.73, (9.11, 8.48)
				MEM_SA avg: 3.79, (3.96, 3.64)

		BSW time, avg: 15.26, (15.41, 15.03)

Important parameter settings:
	BATCH_SIZE: 512
	MAX_SEQ_LEN_REF: 256
	MAX_SEQ_LEN_QER: 128
	MAX_SEQ_LEN8: 128
	SEEDS_PER_READ: 500
	SIMD_WIDTH8 X: 64
	SIMD_WIDTH16 X: 32
	AVG_SEEDS_PER_READ: 64

real	1m29.761s
user	55m5.413s
sys	    2m1.258s

# sorting and duplicate marking (representative per segment)
real	2m24.407s
user	131m34.471s
sys	    0m43.528s

# compress sam (gz) approx same size as aligned bam (per segment)
size (bytes)    component (segment 0020 example)
994983924	    DPZw_k.0020.aligned.sam.gz
1373	        DPZw_k.0020.duplication_metrics
1071385832	    DPZw_k.0020.sorted.bam
7501696	        DPZw_k.0020.sorted.bam.bai

# all sam.gz files can be deleted after sorting and dup marking
# and can happen on a per file basis

# combined size of aligned sams (gzipped)
63621796264

# combined size of sorted and marked bams (plus index files)
68557961160 (bams)
479536680 (index files)

# samtools merge and convert to bam
real	19m34.043s
user	163m40.812s
sys	    2m59.965s

# resulting size of aligned and sorted bam
41007349002 (bam)
5084787 (index)

# the next operation here is based on creating intervals from the bam
