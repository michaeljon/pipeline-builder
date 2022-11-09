library(ggplot2)
library(ggthemes)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

for (arg in args) {
  input <- sprintf("%s.tsv", arg)
  output <- sprintf("%s.png", arg)

  df <- read_tsv(input, col_names = TRUE, show_col_types = FALSE)
  sample <- df[1, 1]
  organism <- df[1, 2]
  gene <- df[1, 3]

  median_depth <- df %>%
    pull(depth) %>%
    sort() %>%
    median(na.rm = FALSE)

  if (median_depth >= 50) {
    png(output)

    plt <- ggplot(data = df, mapping = aes(x = pos_in_gene, y = depth)) +
      ylim(0, NA) +
      xlab(sprintf("Position in gene %s", gene)) +
      ylab("Depth") +
      geom_hline(
        yintercept = median_depth,
        color = "purple",
        linetype = "dotted",
        size = 1.5
      ) +
      geom_line(color = "darkgray") +
      geom_smooth() +
      theme_hc() +
      labs(
        title = sprintf("Depth by position in gene %s for %s", gene, sample),
        subtitle = sprintf("Organism: %s", organism)
      )

    print(plt)
    dev.off()
  }
}

q()
