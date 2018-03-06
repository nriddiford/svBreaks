# Functions related to calculating enrichment
# of svs in genomic features

# bpFeatureEnrichment
#'
#' Calculate the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import tidyverse
#' @export

bpFeatureEnrichment <- function(features=system.file("extdata", "genomic_features.txt", package="svBreaks"), genome_length=118274340, print=NA) {
  genome_features <- read.delim(features, header = T)
  bp_data <- getData()

  mutCount <- nrow(bp_data)

  # To condense exon counts into "exon"
  bp_data$feature <- as.factor(gsub("exon_.*", "exon", bp_data$feature))

  classCount <- table(bp_data$feature)
  classLengths <- setNames(as.list(genome_features$length), genome_features$feature)

  fun <- function(f) {
    # Calculate the fraction of geneome occupied by each feature
    featureFraction <- classLengths[[f]] / genome_length

    # How many times should we expect to see this feature hit in our bp_data (given number of obs. and fraction)?
    featureExpect <- (mutCount * featureFraction)

    # observed/expected
    fc <- classCount[[f]] / featureExpect
    fc <- round(fc, digits = 1)
    featureExpect <- round(featureExpect, digits = 3)

    # Binomial test
    if (!is.null(classLengths[[f]])) {
      if (classCount[f] >= featureExpect) {
        stat <- binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "greater")
        test <- "enrichment"
      }
      else {
        stat <- binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "less")
        test <- "depletion"
      }
      sig_val <- "F"
      if (stat$p.value <= 0.05) {
        sig_val <- "T"
      }
      p_val <- format.pval(stat$p.value, digits = 3, eps = 0.0001)
      Log2FC <- log2(fc)
      # Log2FC<-round(Log2FC, 1)
      list(feature = f, observed = classCount[f], expected = featureExpect, Log2FC = Log2FC, test = test, sig = sig_val, p_val = p_val)
    }

  }

  enriched <- lapply(levels(bp_data$feature), fun)
  enriched <- do.call(rbind, enriched)
  featuresFC <- as.data.frame(enriched)
  # Sort by FC value
  featuresFC <- dplyr::arrange(featuresFC, desc(abs(as.numeric(Log2FC))))
  featuresFC$Log2FC <- round(as.numeric(featuresFC$Log2FC), 1)
  featuresFC$expected <- round(as.numeric(featuresFC$expected), 1)

  if (!is.na(print)) {
    cat("printing")
    first.step <- lapply(featuresFC, unlist)
    second.step <- as.data.frame(first.step, stringsAsFactors = F)

    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))

    feat_enrichment_table <- paste("feature_enrichment_table.tiff")
    # ggexport(filename = paste("plots/", "ex_", feat_enrichment_table, sep=""),
    #          width = 480, height = 480, pointsize = 12, res = 250,
    #          verbose = TRUE)
    ggsave(paste("plots/", feat_enrichment_table, sep = ""), width = 5.2, height = (nrow(featuresFC) / 3), dpi = 300)
  }

  else {
    return(featuresFC)
  }
}

# bpFeatureEnrichmentPlot
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import tidyverse
#' @export
#'
bpFeatureEnrichmentPlot <- function() {
  feature_enrichment <- bpFeatureEnrichment()
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -Log2FC))

  feature_enrichment <- dplyr::filter(feature_enrichment, observed >= 3)
  feature_enrichment <- droplevels(feature_enrichment)

  p <- ggplot(feature_enrichment)
  p <- p + geom_bar(aes(feature, Log2FC, fill = as.character(test)), stat = "identity")
  p <- p + guides(fill = FALSE)
  p <- p + ylim(-3, 3)
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  # p <- p + facet_wrap(~genotype)

  feat_plot <- paste("feat_plot.pdf")
  cat("Writing file", feat_plot, "\n")
  ggsave(paste("plots/", feat_plot, sep = ""), width = 5, height = 10)
  p
}
