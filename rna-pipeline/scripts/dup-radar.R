# Command line argument processing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: dupRadar.R <input.bam> <annotation.gtf> <strandDirection:0=unstranded/1=forward/2=reverse> <paired/single> <nbThreads>", call. = FALSE)
}
input_bam <- args[1]
annotation_gtf <- args[2]
stranded <- as.numeric(args[3])
paired_end <- if (args[4] == "paired") TRUE else FALSE
threads <- as.numeric(args[5])

bamRegex <- "(.+)\\.bam$"

if (!(grepl(bamRegex, input_bam) && file.exists(input_bam) && (!file.info(input_bam)$isdir))) stop("First argument '<input.bam>' must be an existing file (not a directory) with '.bam' extension...")
if (!(file.exists(annotation_gtf) && (!file.info(annotation_gtf)$isdir))) stop("Second argument '<annotation.gtf>' must be an existing file (and not a directory)...")
if (is.na(stranded) || (!(stranded %in% (0:2)))) stop("Third argument <strandDirection> must be a numeric value in 0(unstranded)/1(forward)/2(reverse)...")
if (is.na(threads) || (threads <= 0)) stop("Fifth argument <nbThreads> must be a strictly positive numeric value...")

# Remove bam file extension to generate basename
input_bam_basename <- gsub(bamRegex, "\\1", input_bam)

# Debug messages (stderr)
message("Input bam      (Arg 1): ", input_bam)
message("Input gtf      (Arg 2): ", annotation_gtf)
message("Strandness     (Arg 3): ", c("unstranded", "forward", "reverse")[stranded + 1])
message("paired/single  (Arg 4): ", ifelse(paired_end, "paired", "single"))
message("Nb threads     (Arg 5): ", threads)
message("Output basename       : ", input_bam_basename)


# Load / install packages
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!require("dupRadar")) {
    BiocManager::install("dupRadar")
}

library("dupRadar")

if (!require("parallel")) {
    install.packages("parallel", dependencies = TRUE, repos = "http://cloud.r-project.org/")
    library("parallel")
}

# Duplicate stats
dm <- analyzeDuprates(input_bam, annotation_gtf, stranded, paired_end, threads)
write.table(dm, file = paste(input_bam_basename, "_dupMatrix.txt", sep = ""), quote = F, row.name = F, sep = "\t")

# 2D density scatter plot
pdf(paste0(input_bam_basename, "_duprateExpDens.pdf"))
duprateExpDensPlot(DupMat = dm)
title("Density scatter plot")
mtext(input_bam_basename, side = 3)
dev.off()

fit <- duprateExpFit(DupMat = dm)
cat(
    paste("- dupRadar Int (duprate at low read counts):", fit$intercept),
    paste("- dupRadar Sl (progression of the duplication rate):", fit$slope),
    fill = TRUE, labels = input_bam_basename,
    file = paste0(input_bam_basename, "_intercept_slope.txt"), append = FALSE
)

# Create a multiqc file dupInt
sample_name <- gsub("Aligned.sortedByCoord.out.markDups", "", input_bam_basename)
line <- "#id: DupInt
#plot_type: 'generalstats'
#pconfig:
#    dupRadar_intercept:
#        title: 'dupInt'
#        namespace: 'DupRadar'
#        description: 'Intercept value from DupRadar'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'
#        format: '{:.2f}%'
Sample dupRadar_intercept"

write(line, file = paste0(input_bam_basename, "_dup_intercept_mqc.txt"), append = TRUE)
write(paste(sample_name, fit$intercept), file = paste0(input_bam_basename, "_dup_intercept_mqc.txt"), append = TRUE)

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y <- 100 * predict(fit$glm, data.frame(x = curve_x), type = "response")

# Remove all of the infinite values
infs <- which(curve_x %in% c(-Inf, Inf))
curve_x <- curve_x[-infs]
curve_y <- curve_y[-infs]

# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]

# Convert x values back to real counts
curve_x <- 10^curve_x

# Write to file
line <- "#id: DupRadar
#section_name: 'DupRadar'
#section_href: 'https://bioconductor.org/packages/release/bioc/html/dupRadar.html'
#description: \"provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.
#    This plot shows the general linear models - a summary of the gene duplication distributions. \"
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xLog: True
#    xlab: 'expression (reads/kbp)'
#    ylab: '% duplicate reads'
#    ymax: 100
#    ymin: 0
#    tt_label: '<b>{point.x:.1f} reads/kbp</b>: {point.y:,.2f}% duplicates'
#    xPlotLines:
#        - color: 'green'
#          dashStyle: 'LongDash'
#          label:
#                style: {color: 'green'}
#                text: '0.5 RPKM'
#                verticalAlign: 'bottom'
#                y: -65
#          value: 0.5
#          width: 1
#        - color: 'red'
#          dashStyle: 'LongDash'
#          label:
#                style: {color: 'red'}
#                text: '1 read/bp'
#                verticalAlign: 'bottom'
#                y: -65
#          value: 1000
#          width: 1"

write(line, file = paste0(input_bam_basename, "_duprateExpDensCurve_mqc.txt"), append = TRUE)
write.table(
    cbind(curve_x, curve_y),
    file = paste0(input_bam_basename, "_duprateExpDensCurve_mqc.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,
)

# Distribution of expression box plot
pdf(paste0(input_bam_basename, "_duprateExpBoxplot.pdf"))
duprateExpBoxplot(DupMat = dm)
title("Percent Duplication by Expression")
mtext(input_bam_basename, side = 3)
dev.off()

# Distribution of RPK values per gene
pdf(paste0(input_bam_basename, "_expressionHist.pdf"))
expressionHist(DupMat = dm)
title("Distribution of RPK values per gene")
mtext(input_bam_basename, side = 3)
dev.off()

# Print sessioninfo to standard out
print(input_bam_basename)
citation("dupRadar")
sessionInfo()
