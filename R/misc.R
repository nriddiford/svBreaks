
#' Misc functions

#' svTypes
#'
#' Plot counts of different svTypes
#' @import tidyverse
#' @export

svTypes <- function(notch=0, object=NA) {
  if (is.na(object)) {
    object <- "type"
  }

  if (notch) {
    bp_data <- getData()
    bp_data <- dplyr::filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)

    ext <- "_Notch.pdf"
  }
  else {
    bp_data <- getData()
    ext <- ".pdf"
  }

  bp_data$type <- ifelse(bp_data$type == "BND", "INV", as.character(bp_data$type))
  bp_data <- droplevels(bp_data)

  cols <- setCols(bp_data, "type")

  # Reorder by count
  bp_data$type <- factor(bp_data$type, levels = names(sort(table(bp_data$type), decreasing = TRUE)))

  if (object == "sample") {
    # Reorder by count
    bp_data$sample <- factor(bp_data$sample, levels = names(sort(table(bp_data$sample), decreasing = TRUE)))
  }

  # Only take bp1 for each event
  bp_data <- dplyr::filter(bp_data, bp_no != "bp2")
  bp_data <- droplevels(bp_data)

  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(get(object), fill = type))
  p <- p + cols
  p <- p + cleanTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 40), axis.title = element_text(size = 90)
    )
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous("Number of calls", expand = c(0.01, 0.01))
  p <- p + facet_wrap(~genotype)
  # p<-p + coord_flip()
  # p<-p + scale_y_reverse()

  types_outfile <- paste("sv_types_by_", object, ext, sep = "")
  cat("Writing file", types_outfile, "\n")
  ggsave(paste("plots/", types_outfile, sep = ""), width = 22, height = 22)

  p
}

#' featureDensity
#'
#' Plot density of several feature tyopes across the genome
#' @import dplyr
#' @import ggplot2
#' @export

featureDensity <- function(feature_file1 = system.file("extdata", "g4_positions.txt", package="svBreaks"), feature1 = 'G4',
                           feature_file2 = NA, feature2 = NA,
                           feature_file3 = NA, feature3 = NA) {

  n = 1

  file_list = list()

  f1 <- read.delim(feature_file1, header = F)
  f1$type <- feature1
  colnames(f1) <- c("chrom", "pos", "type")
  file_list[[n]] <- f1

  if (!is.na(feature_file2)){
    n = n + 1
    f2 <- read.delim(feature_file2, header = F)
    f2$type <- feature2
    colnames(f2) <- c("chrom", "pos", "type")
    file_list[[n]]  <- f2
  }
  if (!is.na(feature_file3)){
    n = n + 1
    f3 <- read.delim(feature_file3, header = F)
    f3$type <- feature3
    colnames(f3) <- c("chrom", "pos", "type")
    file_list[[n]]  <- f3
  }

  files <- do.call(rbind, file_list)
  rownames(files) <- NULL

  chroms = c("2L", "2R", "3L", "3R", "X")
  locations <- files %>%
    mutate(type = as.factor(type)) %>%
    mutate(pos = as.numeric(pos/1000000)) %>%
    filter(chrom %in% chroms) %>%
    arrange(chrom, pos) %>%
    droplevels()


  p <- ggplot(locations)
  p <- p + geom_density(aes(pos, fill = type), alpha = 0.4)
  p <- p + geom_rug(aes(pos, colour = type), sides = "tb", alpha = 0.05)
  p <- p + facet_wrap(~chrom, scale = "free_x", nrow=length(levels(locations$chrom)))
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, max(locations$pos), by = 1))
  p <- p + slideTheme() +
    theme(axis.text.y = element_blank())

  featurePlot <- paste("feature_dist.png")
  cat("Writing file", featurePlot, "\n")
  ggsave(paste("plots/", featurePlot, sep = ""), width = 20, height = 10)

  p
}

#' typeLen
#'
#' Plot the length of different sv types
#' @import tidyverse
#' @export

typeLen <- function(size_threshold = 1, notch=0) {
  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }
  else {
    bp_data <- getData()
    ext <- ".pdf"
  }

  cols <- setCols(bp_data, "type")

  # Only take bp1 for each event
  bp_data <- dplyr::filter(bp_data, type != "TRA", type != "BND", bp_no != "bp2")

  bp_data$length <- (bp_data$length / 1000)

  if (is.na(size_threshold)) {
    size_threshold <- max(bp_data$length)
  }

  if (size_threshold <= 1) {
    breaks <- 0.1
  }
  else {
    breaks <- 1
  }

  p <- ggplot(bp_data, aes(length))
  p <- p + geom_density(aes(fill = type), alpha = 0.4)
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"))
  p <- p + scale_x_continuous("Size in Mb", expand = c(0, 0), breaks = seq(0, size_threshold, by = breaks), limits = c(0, (size_threshold + 0.1)))
  p <- p + scale_y_continuous(expand = c(0, 0))
  p <- p + cols

  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep = ""), width = 20, height = 10)

  p
}

#' typeLenCount
#'
#' Plot the length of different sv types as counts
#' @import tidyverse
#' @export

typeLenCount <- function(size_threshold = 1, notch=0) {
  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }
  else {
    bp_data <- getData()
    ext <- ".pdf"
  }

  cols <- setCols(bp_data, "type")

  # Only take bp1 for each event
  bp_data <- dplyr::filter(bp_data, type != "TRA", type != "BND", bp_no != "bp2")

  bp_data$length <- (bp_data$length / 1000)

  if (is.na(size_threshold)) {
    size_threshold <- max(bp_data$length)
  }

  if (size_threshold <= 1) {
    breaks <- 0.1
  }
  else {
    breaks <- 1
  }

  p <- ggplot(bp_data, aes(length))
  p <- p + geom_histogram(aes(length, ..count.., fill = type), colour = "black", binwidth = 0.05, position = "dodge")
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"))
  p <- p + scale_x_continuous("Size in Mb", expand = c(0, 0), breaks = seq(0, size_threshold, by = breaks), limits = c(0, (size_threshold + 0.1)))
  p <- p + scale_y_continuous(expand = c(0, 0))
  p <- p + geom_density(aes(fill = type), alpha = 0.4, colour = NA)
  p <- p + cols

  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep = ""), width = 20, height = 10)

  p
}


#' genomeHits
#'
#' Plot distribution of brekpoints across the genome
#' @import tidyverse
#' @export

genomeHits <- function(notch=0) {
  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }
  else {
    bp_data <- getData()
    ext <- ".pdf"
  }

  p <- ggplot(bp_data)
  p <- p + geom_point(aes(bp / 1000000, sample, colour = sample, shape = type, size = 0.5), alpha = 0.7)
  p <- p + guides(color = FALSE, size = FALSE)
  p <- p + cleanTheme() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 20),
      strip.text.x = element_text(size = 15),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, 33, by = 1), limits = c(0, 33), expand = c(0.01, 0.01))

  sv_gen_dist <- paste("bp_gen.dist", ext, sep = "")
  cat("Writing file", sv_gen_dist, "\n")
  ggsave(paste("plots/", sv_gen_dist, sep = ""), width = 20, height = 10)

  p
}



#' bpGenAll
#'
#' Plot distribution of brekpoints across the genome
#' @import tidyverse
#' @export

bpGenAll <- function(object=NA, notch=0) {
  bp_data <- getData()
  ext <- ".pdf"
  if (is.na(object)) {
    object <- "type"
    cols <- setCols(bp_data, "type")
  }

  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }

  cat("Plotting SVs by", object, "\n")

  p <- ggplot(bp_data)
  p <- p + geom_histogram(aes(bp / 1000000, fill = get(object)), binwidth = 0.1, alpha = 0.8)
  p <- p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, 33, by = 1), limits = c(0, 33), expand = c(0.01, 0.01))
  p <- p + scale_y_continuous("Number of Breakpoints", expand = c(0.01, 0.01))
  p <- p + cleanTheme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 20),
      strip.text.x = element_text(size = 15)
    )

  if (object == "type") {
    p <- p + cols
  }

  chrom_outfile <- paste("Breakpoints_chroms_by_", object, ext, sep = "")
  cat("Writing file", chrom_outfile, "\n")
  ggsave(paste("plots/", chrom_outfile, sep = ""), width = 20, height = 10)

  p
}

#' bpChromDist
#'
#' Plot distribution of brekpoints across the genome by chromosome
#' @import tidyverse
#' @export

bpChromDist <- function(object=NA, notch=0) {
  bp_data <- getData()
  ext <- ".pdf"

  if (is.na(object)) {
    object <- "type"
    cols <- setCols(bp_data, "type")
  }

  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }

  chromosomes <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
  lengths <- c(23513712, 25286936, 28110227, 32079331, 23542271, 3667352, 1348131)

  karyotype <- setNames(as.list(lengths), chromosomes)

  for (c in chromosomes) {
    len <- karyotype[[c]]
    len <- len / 1000000

    cat("Chrom", c, "length:", len, sep = " ", "\n")

    per_chrom <- dplyr::filter(bp_data, chrom == c)

    p <- ggplot(per_chrom)
    p <- p + geom_histogram(aes(bp / 1000000, fill = get(object)), binwidth = 0.1, alpha = 0.8)
    p <- p + scale_x_continuous("Mbs", breaks = seq(0, len, by = 1), limits = c(0, len + 0.1), expand = c(0.01, 0.01))
    p <- p + scale_y_continuous("Number of Breakpoints", limits = c(0, 35), expand = c(0.01, 0.01))
    p <- p + cleanTheme() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 20)
      )
    p <- p + ggtitle(paste("Chromosome: ", c))

    if (object == "type") {
      p <- p + cols
    }

    per_chrom <- paste("Breakpoints_on_", c, "_by_", object, ext, sep = "")
    cat("Writing file", per_chrom, "\n")
    ggsave(paste("plots/", per_chrom, sep = ""), width = 20, height = 10)
  }
}
