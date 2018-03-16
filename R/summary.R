# A collection of summary functions
# for plotting breakpints per sample

#' bpStats
#' Get some basic stats for breakpoints
#' @keywords stats
#' @import tidyverse
#' @export
#' @param colSample Set to <samplename> to colour that sample red in the plot

bpStats <- function(colSample=NA) {
  bp_data <- getData()

  sampleSvs <- bp_data %>%
    dplyr::filter(bp_no == "bp1") %>%
    group_by(sample, genotype) %>%
    tally() %>%
    mutate(som_count = ifelse(genotype=='somatic_tumour', n, 1))

  exclude_samples <- c("A373R1", "A373R7", "A512R17", "A373R11")

  if (!is.na(colSample)) {
    sampleSvs$colour <- ifelse(sampleSvs$sample %in% exclude_samples, "#C72424FE", "grey37")
  }
  else {
    sampleSvs$colour <- "grey37"
  }

  p <- ggplot(sampleSvs)
  p <- p + geom_bar(aes(fct_reorder(sample, som_count), n, fill = colour), stat = "identity")
  # p <- p + scale_y_continuous("Number of variants", limits = c(0, 20), expand = c(0.01, 0.01), breaks=seq(0,20,by=2))
  p <- p + scale_y_continuous("Number of variants", expand = c(0.01, 0.01))

  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
      axis.title.x = element_blank()
    )
  p <- p + scale_fill_identity()
  p <- p + facet_wrap(~genotype, ncol = 1, scales = "free_y")


  sampleSVs <- paste("SVs_sample.png")
  cat("Writing file", sampleSVs, "\n")
  ggsave(paste("plots/", sampleSVs, sep = ""), width = 20, height = 10)

  cat("sample", "SVs", sep = "\t", "\n")
  bp_data <- dplyr::filter(bp_data, genotype == 'somatic_tumour')
  rank <- sort(table(bp_data$sample), decreasing = TRUE)
  rank <- as.array(rank)

  total <- 0

  scores <- list()
  for (i in 1:nrow(rank)) {
    cat(names(rank[i]), rank[i], sep = "\t", "\n")
    total <- total + rank[i]
    scores[i] <- rank[i]
  }
  cat("--------------", "\n")
  scores <- unlist(scores)

  mean <- as.integer(mean(scores))
  med <- as.integer(median(scores))

  cat("total", total, sep = "\t", "\n")
  cat("samples", nrow(rank), sep = "\t", "\n")

  cat("--------------", "\n")
  cat("mean", mean, sep = "\t", "\n")
  cat("median", med, sep = "\t", "\n")

  cat("\n")
  dels <- nrow(dplyr::filter(bp_data, type == "DEL"))
  invs <- nrow(dplyr::filter(bp_data, type == "INV"))
  dups <- nrow(dplyr::filter(bp_data, type == "DUP"))
  tandups <- nrow(dplyr::filter(bp_data, type == "TANDUP"))
  tra <- nrow(dplyr::filter(bp_data, type == "TRA"))

  cat("Dels ", dels, sep = "", "\n")
  cat("Invs ", invs, sep = "", "\n")
  cat("Dups ", dups, sep = "", "\n")
  cat("Tandups ", tandups, sep = "", "\n")
  cat("Tra ", tra, sep = "", "\n")

  p
}


#' bpFeatures
#' Get some basic stats for breakpoints
#' @keywords stats
#' @import tidyverse
#' @export
#' @param notch - Set to 1 to exclude Notch events

bpFeatures <- function(notch=0) {
  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }
  else {
    bp_data <- getData()
    ext <- ".pdf"
  }

  # To condense exon counts into "exon"
  bp_data$feature <- as.factor(gsub("_.*", "", bp_data$feature))

  bp_data <- dplyr::filter(bp_data, feature != "pseudogene" & feature != "snoRNA" & feature != "stop" & feature != "tRNA")
  bp_data <- droplevels(bp_data)
  # Reoders descending
  bp_data$feature <- factor(bp_data$feature, levels = names(sort(table(bp_data$feature), decreasing = TRUE)))

  cols <- setCols(bp_data, "feature")

  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(feature, fill = genotype, group = genotype), stat = "count", position = position_dodge())
  # p<-p + cols
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous(expand = c(0.01, 0.01))
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + facet_wrap(~genotype)

  features_outfile <- paste("Breakpoints_features_count", ext, sep = "")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep = ""), width = 20, height = 10)

  p
}

#' svsbySample
#' Plot the number of structural breakpoints per sample
#' @keywords samples
#' @import tidyverse
#' @export

svsbySample <- function() {
  bp_data <- getData()
  bp_data <- dplyr::filter(bp_data, bp_no != "bp2")
  bp_data$sample <- factor(bp_data$sample)

  sampleSvs <- bp_data %>%
    group_by(sample, genotype) %>%
    tally() %>%
    ungroup() %>%
    arrange(-n)


  p <- ggplot(sampleSvs)
  p <- p + geom_bar(aes(sample, n, fill = genotype), stat = "identity")
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous("Number of SVs", expand = c(0.01, 0.01))
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  p
}
