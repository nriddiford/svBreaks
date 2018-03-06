#' A function to plot the size distribution of SV types

#' sizeDist
#' Plot the size distribution for different types of structual variants
#' @keywords size
#' @import tidyverse
#' @import RColorBrewer
#' @export

sizeDist <- function() {
  bp_data <- getData()
  bp_data <- bp_data %>%
    dplyr::filter(type != "TRA", type != "BND", bp_no != "bp2") %>%
    mutate(length = ifelse(length == 0, 0.01, length)) %>%
    mutate(type = ifelse(type == "TANDUP", "DUP", as.character(type))) %>%
    droplevels()

  fillCols <- setCols(bp_data, "genotype", fill = "Y")
  cols <- setCols(bp_data, "genotype", fill = "N")

  bp_data <- transform(bp_data, genotype = reorder(genotype, length))

  p <- ggplot(bp_data, aes(genotype, length))
  p <- p + geom_violin(aes(fill = genotype))
  p <- p + geom_jitter(width = .1, height = .1, size = 0.2, alpha = 0.05)

  p <- p + scale_y_log10("Length (Kb)")
  p <- p + scale_x_discrete(labels = c("S_n", "G_r", "G_p", "S_t"))
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + facet_wrap(~type, nrow = 3)
  p <- p + fillCols
  p <- p + cols


  sizedistOut <- paste("sizeDist.png")
  cat("Writing file", sizedistOut, "\n")
  ggsave(paste("plots/", sizedistOut, sep = ""), width = 15, height = 10)

  p
}
