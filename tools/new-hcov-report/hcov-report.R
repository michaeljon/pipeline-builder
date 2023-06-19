suppressPackageStartupMessages(
  library(
    ggplot2,
    quietly = TRUE, verbose = FALSE
  )
)
suppressPackageStartupMessages(
  library(
    ggthemes,
    quietly = TRUE, verbose = FALSE
  )
)
suppressPackageStartupMessages(
  library(
    tidyverse,
    quietly = TRUE, verbose = FALSE
  )
)
suppressPackageStartupMessages(
  library(
    dplyr,
    quietly = TRUE, verbose = FALSE
  )
)

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop(
    "Usage: Rscript --vanilla gene_coverage.R <sample> <organism> <input-file> <output-folder>",
    call. = FALSE
  )
}

sample <- args[1]
organism <- args[2]
input <- args[3]
outpath <- args[4]

output <- paste0(outpath, "/", sample, "-", organism, ".pdf")

min_median_depth <- 0
min_coverage <- 95

regions <- read_csv("./hcov-regions.csv", show_col_types = FALSE)

gene_list <- regions[regions$organism == organism, c("gene", "start", "stop")]
organism_data <- read_tsv(input, col_names = TRUE, show_col_types = FALSE)

df <- organism_data[, c("position", "depth")]
non_zero <- nrow(filter(df, df$depth != 0))
all <- nrow(df)
coverage <- non_zero / all * 100.0

if (sum(df$depth) == 0) {
  print(paste(
    "Not processing ", sample, " for ", organism, ". No reads found.",
    sep = ""
  ))
  quit(status = 100)
}

if (coverage <= min_coverage) {
  print(paste(
    "Not processing ", sample, " for ", organism, ". Breadth ", as.integer(coverage), "% too low.",
    sep = ""
  ))
  quit(status = 200)
}

pdf(
  output,
  onefile = TRUE,
  title = paste(
    "Feature breadth / depth for", sample, ":", organism
  ),
  paper = "us"
)

median_depth <- df %>%
  pull(depth) %>%
  median(na.rm = FALSE)

mean_depth <- df %>%
  pull(depth) %>%
  mean(na.rm = FALSE)

min_depth <- df %>%
  pull(depth) %>%
  sort() %>%
  min(na.rm = FALSE)

max_depth <- df %>%
  pull(depth) %>%
  sort() %>%
  max(na.rm = FALSE)

if (!is.na(median_depth) && median_depth >= min_median_depth) {
  partitions <- c()
  labels <- c()
  positions <- c()

  for (g in 1:nrow(gene_list)) {
    gene <- gene_list[g, ]

    partitions <- append(partitions, as.integer(gene$start))
    position <- if (g %% 2 == 0) {
      as.integer(max_depth * 0.95)
    } else {
      as.integer(0)
    }

    positions <- append(
      positions,
      position
    )

    label <- gene$gene
    if (label == "3-prime-utr") {
      label <- "3'"
    } else if (label == "5-prime-utr") {
      label <- "5'"
    }

    labels <- append(labels, label)

    if (g == nrow(gene_list)) {
      partitions <- append(partitions, as.integer(gene$stop))
      labels <- append(labels, "*")
      positions <- append(positions, as.integer(max_depth * 0.95))
    }
  }

  suppressWarnings(
    print(
      ggplot(data = df, mapping = aes(x = position, y = depth)) +
        ylim(0, NA) +
        xlab("Position in organism") +
        ylab("Depth") +
        geom_line(color = "darkgray") +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
        geom_vline(
          xintercept = partitions,
          color = "orange",
          linetype = "solid",
          linewidth = 0.25
        ) +
        geom_hline(
          yintercept = median_depth,
          color = "darkblue",
          linetype = "solid",
          linewidth = 1.0
        ) +
        geom_hline(
          yintercept = mean_depth,
          color = "darkgreen",
          linetype = "solid",
          linewidth = 1.0
        ) +
        theme_hc(base_size = 16) +
        labs(
          title = sprintf("Sample: %s", sample),
          subtitle = sprintf("Organism: %s", organism),
          caption = paste(
            "Min. depth: ", min_depth, "\n",
            "Median: ", as.integer(median_depth), "\n",
            "Mean: ", as.integer(mean_depth), "\n",
            sep = ""
          )
        ) +
        annotate(
          "text",
          x = partitions,
          y = positions,
          label = labels,
          angle = 45,
          size = 3
        )
    )
  )
}


for (g in 1:nrow(gene_list)) {
  gene <- gene_list[g, ]

  df <- organism_data[
    (organism_data$position >= gene$start & organism_data$position <= gene$stop),
    c("position", "depth")
  ]

  if (sum(df$depth) > 0) {
    median_depth <- df %>%
      pull(depth) %>%
      median(na.rm = FALSE)

    mean_depth <- df %>%
      pull(depth) %>%
      mean(na.rm = FALSE)

    min_depth <- df %>%
      pull(depth) %>%
      min(na.rm = FALSE)

    if (!is.na(median_depth) && median_depth >= min_median_depth) {
      suppressWarnings(
        print(
          ggplot(data = df, mapping = aes(x = position, y = depth)) +
            ylim(0, NA) +
            xlab("Base relative to full genome") +
            ylab("Depth") +
            geom_line(color = "darkgray") +
            geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
            geom_hline(
              yintercept = median_depth,
              color = "darkblue",
              linetype = "solid",
              linewidth = 1.0
            ) +
            geom_hline(
              yintercept = mean_depth,
              color = "darkgreen",
              linetype = "solid",
              linewidth = 1.0
            ) +
            theme_hc(base_size = 16) +
            labs(
              title = sprintf("Feature: %s", gene$gene),
              subtitle = sprintf("Organism: %s", organism),
              caption = paste(
                "Min. depth: ", min_depth, "\n",
                "Median: ", as.integer(median_depth), "\n",
                "Mean: ", as.integer(mean_depth), "\n",
                "Width: ", (gene$stop - gene$start + 1),
                sep = ""
              )
            )
        )
      )
    }
  }
}

graphics.off()
