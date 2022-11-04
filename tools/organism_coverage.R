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

  if (median_depth > 50) {
    png(output)

    plt <- ggplot(data = df, mapping = aes(x = position, y = depth)) +
      ylim(0, NA) +
      xlab(sprintf("Position in organism %s", organism)) +
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
        title = sprintf("Depth by position in organism for %s", sample),
        subtitle = sprintf("Organism: %s", organism)
      )

    print(plt)
    dev.off()
  }
}

<<<<<<< HEAD
q()
=======
q()
>>>>>>> 8cdaa679289a4db107c19664625145afc6d96a98
