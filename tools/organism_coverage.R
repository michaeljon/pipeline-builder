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

  median_depth <- df %>%
    pull(depth) %>%
    sort() %>%
    median(na.rm = FALSE)

  if (median_depth > 5) {
    png(output)

    plt <- ggplot(data = df, mapping = aes(x = position, y = depth)) +
      ylim(0, NA) +
      xlab(sprintf("Position in organism %s", organism)) +
      ylab("Depth") +
      geom_hline(
        yintercept = median_depth,
        color = "purple",
        linetype = "dotted",
        linewidth = 1.5
      ) +
      geom_line(color = "darkgray") +
      geom_smooth() +
      theme_hc(base_size = 16) +
      labs(
        title = sample,
        subtitle = sprintf("Organism: %s", organism)
      )

    print(plt)
    dev.off()
  }
}

q()
