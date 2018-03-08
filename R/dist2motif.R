# Functions to calculate the distance
# from each breakpoint to user-provided loci (e.g. TSS)

#' generateData
#' Prepare data for dist2motif
#' @keywords simulate
#' @import tidyverse
#' @export

generateData <- function(sim=NA){
  real_data <- getData()
  # real_data <- notchFilt(keep=0)
  real_data <- real_data %>%
    dplyr::filter(chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" ) %>%
    dplyr::rename(pos = bp) %>%
    droplevels()

  if (!is.na(sim)) {
    # for (i in (1:iterations)){
    # cat("Running iteration", i, "\n")
    simByChrom <- list()

    for (c in levels(real_data$chrom)){
      hitCount <- nrow(real_data[real_data$chrom== c,])
      cat(paste("Simulating", hitCount, "breakpoints on chromosome", c), "\n")
      bp_data <- dtableShuffle(nSites = hitCount, byChrom = c)
      simByChrom[[c]] <- bp_data
      # simByChrom[[]]$iteration <- i
    }
    bp_data <- as.data.frame(do.call(rbind, simByChrom))
    rownames(bp_data) <- NULL
    return(bp_data)
  } else{
    return(real_data)
  }
}


#' dist2Motif
#' Calculate the distance from each breakpoint to closest motif
#' @keywords motif
#' @import tidyverse
#' @export

dist2Motif <- function(feature_file = system.file("extdata", "tss_locations.txt", package="svBreaks"), sim=NA, print=0, send=0, feature="tss") {
  bp_data <- generateData(sim=sim)

  feature <- paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep = "")

  if (feature == "Promoter") {
    feature_locations <- getPromoter()
    cat("Getting gene promoter locations...\n")
  } else {
    feature_locations <- read.delim(feature_file, header = F)
    cat("Reading in file:", feature_file, sep = " ", "\n")
  }

  cat("Calculating distances to", feature, sep = " ", "\n")

  # feature_locations <- read.delim(feature_file, header = F)
  colnames(feature_locations) <- c("chrom", "pos")

  feature_locations$pos <- as.integer(feature_locations$pos)

  # Will throw error if SVs don't exist on a chrom...

  # Removes chroms with fewer than 10 observations
  svCount <- table(bp_data$chrom)
  bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 10]))
  bp_data <- droplevels(bp_data)

  feature_locations <- subset(feature_locations, chrom %in% levels(bp_data$chrom))
  feature_locations <- droplevels(feature_locations)


  fun2 <- function(p) {
    index <- which.min(abs(tss_df$pos - p))
    closestTss <- tss_df$pos[index]
    chrom <- as.character(tss_df$chrom[index])
    gene <- as.character(tss_df$gene[index])
    dist <- (p - closestTss)
    list(p, closestTss, dist, chrom, gene)
  }

  l <- list()

  for (c in levels(bp_data$chrom)) {
    df <- dplyr::filter(bp_data, chrom == c)
    tss_df <- dplyr::filter(feature_locations, chrom == c)
    dist2tss <- lapply(df$pos, fun2)
    dist2tss <- do.call(rbind, dist2tss)
    dist2tss <- as.data.frame(dist2tss)

    colnames(dist2tss) <- c("bp", "closest_tss", "min_dist", "chrom", "closest_gene")
    dist2tss$min_dist <- as.numeric(dist2tss$min_dist)
    l[[c]] <- dist2tss
  }

  dist2tss <- do.call(rbind, l)
  dist2tss <- as.data.frame(dist2tss)
  dist2tss$chrom <- as.character(dist2tss$chrom)

  dist2tss <- arrange(dist2tss, (abs(min_dist)))

  if (send == 1) {
    return(dist2tss)
  } else {
    p <- ggplot(dist2tss)
    p <- p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
    p <- p + scale_x_continuous(
      paste("Distance to", feature, "(Kb)", sep = " "),
      limits = c(-10000, 10000),
      breaks = c(-10000, -1000, 1000, 10000),
      expand = c(.0005, .0005),
      labels = c("-10", "-1", "1", "10")
    )

    p <- p + scale_y_continuous("Density")
    p <- p + geom_vline(xintercept = 0, colour = "black", linetype = "dotted")
    p <- p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
    p <- p + geom_rug(aes(min_dist, colour = chrom))
    p <- p + slideTheme() +
      theme(
        strip.text = element_text(size = 20),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(4, "lines")
      )

    # p <- p + facet_wrap(~chrom, ncol = 7, scales = "free_y")

    if (is.na(sim)) {
      distout <- paste("bp", feature, "dist.pdf", sep = "")
    } else {
      distout <- paste("bp", feature, "dist_sim.pdf", sep = "")
    }

    cat("Writing file", distout, "\n")
    ggsave(paste("plots/", distout, sep = ""), width = 20, height = 10)

    p
  }
}

# distOverlay
#'
#' Calculate the distance from each breakpoint to closest motif
#' Overlay the same number of random simulated breakpoints
#' @keywords motif
#' @import tidyverse
#' @export

distOverlay <- function(feature_file=system.file("extdata", "tss_locations.txt", package="svBreaks"), feature="tss", lim=10, all=NA) {
  feature <- paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep = "")

  if (feature == "promoter") {
    real_data <- dist2Motif(send = 1, feature = feature)
    sim_data <- dist2Motif(feature = feature, sim = 1, send = 1)
  }
  else {
    real_data <- dist2Motif(feature_file = feature_file, send = 1, feature = feature)
    sim_data <- dist2Motif(feature_file = feature_file, feature = feature, sim = 1, send = 1)
  }

  real_data$Source <- "Real"
  sim_data$Source <- "Sim"

  sim_data <- dplyr::filter(sim_data, chrom != "Y", chrom != 4)
  sim_data <- droplevels(sim_data)
  real_data <- dplyr::filter(real_data, chrom != "Y", chrom != 4)
  real_data <- droplevels(real_data)

  colours <- c("#E7B800", "#00AFBB")

  scale <- "(Kb)"
  if (lim == 0.1) {
    cat("Setting limits to -+100bp\n")
    lims <- c(-100, 100)
    brks <- c(-100, -10, 10, 100)
    expnd <- c(.0005, .0005)
    labs <- c("-100", "-10", "10", "100")
    scale <- "(bp)"
  }

  else if (lim == 0.5) {
    cat("Setting limits to -+0.5kb\n")
    lims <- c(-500, 500)
    brks <- c(-500, -100, 100, 500)
    expnd <- c(.0005, .0005)
    labs <- c("-500", "-100", "100", "500")
    scale <- "(bp)"
  }

  else if (lim == 1) {
    cat("Setting limits to -+1kb\n")
    lims <- c(-1000, 1000)
    brks <- c(-1000, 1000)
    expnd <- c(.0005, .0005)
    labs <- c("-1", "1")
  }
  else if (lim == 5) {
    cat("Setting limits to -+1kb\n")
    lims <- c(-5000, 5000)
    brks <- c(-5000, 5000)
    expnd <- c(.0005, .0005)
    labs <- c("-5", "5")
  }
  else {
    cat("Setting limits to -+10kb\n")
    lims <- c(-10000, 10000)
    brks <- c(-10000, -1000, 1000, 10000)
    expnd <- c(.0005, .0005)
    labs <- c("-10", "-1", "1", "10")
  }

  p <- ggplot()
  p <- p + geom_density(data = real_data, aes(min_dist, fill = Source), alpha = 0.4)
  p <- p + geom_density(data = sim_data, aes(min_dist, fill = Source), alpha = 0.4)
  if (is.na(all)) {
    p <- p + facet_wrap(~chrom, ncol = 2)
  }

  p <- p + scale_x_continuous(
    paste("Distance to", feature, scale, sep = " "),
    limits = lims,
    breaks = brks,
    expand = expnd,
    labels = labs
  )
  p <- p + scale_y_continuous("Density")
  p <- p + geom_vline(xintercept = 0, colour = "black", linetype = "dotted")

  p <- p + geom_rug(data = real_data, aes(min_dist, colour = Source), sides = "b")
  p <- p + geom_rug(data = sim_data, aes(min_dist, colour = Source), sides = "t")

  p <- p + scale_fill_manual(values = colours)
  p <- p + scale_colour_manual(values = colours)

  p <- p + slideTheme() +
    theme(
      strip.text = element_text(size = 20),
      axis.text.y = element_blank(),
      legend.position = "top"
    )

  overlay <- paste("bp", feature, "dist_overlay.pdf", sep = "")
  cat("Writing file", overlay, "\n")
  ggsave(paste("plots/", overlay, sep = ""), width = 20, height = 10)

  p
}


#' bpSim
#'
#' Generate simulated SV breakpoints acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random breakpoints to generate [Default nrow(bp_data)]
#' @import GenomicRanges
#' @import rtracklayer
#' @keywords sim
#' @export

bpSim <- function(intervals=system.file("extdata", "intervals.bed", package="svBreaks"), N=1000, write=F) {
  intFile <- import.bed(intervals)
  space <- sum(width(intFile))
  positions <- sample(c(1:space), N)
  cat("Simulating", N, "breakpoints", sep = " ", "\n")
  new_b <- GRanges(
    seqnames = as.character(rep(seqnames(intFile), width(intFile))),
    ranges = IRanges(start = unlist(mapply(seq, from = start(intFile), to = end(intFile))), width = 1)
  )
  bedOut <- new_b[positions]
  if (write) {
    export.bed(new_b[positions], "data/simulatedBPs.bed")
  }
  remove(new_b)
  return(data.frame(bedOut))
}

#' dtableShuffle
#'
#' Generate simulated SV breakpoints acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random breakpoints to generate [Default nrow(bp_data)]
#' @import data.table
#' @keywords sim
#' @export

dtableShuffle <- function(nSites = 1e3, byChrom = NA, iterations = 2 ){
  intervals=system.file("extdata", "intervals.bed", package="svBreaks")
  bed <- fread(intervals)
  colnames(bed) <- c("chrom", "start", "end", "NA")

  bed <- bed %>%
    select(chrom:end) %>%
    mutate(size = (end-start)+1)  %>%
    as.data.table()

  if(is.na(byChrom)){
    # Randomly sample bed file rows, proportional to the length of each range
    simulatedBps <- bed[sample(.N, size=nSites, replace=TRUE, prob=bed$size)]

    # Randomly sample uniformly within each chosen range
    simulatedBps[, position := sample(start:end, size=1), by=1:dim(simulatedBps)[1]]
  } else {
    bed <- bed %>%
      filter(chrom == byChrom) %>%
      as.data.table()

    # Randomly sample bed file rows, proportional to the length of each range
    simulatedBps <- bed[sample(.N, size=nSites, replace=TRUE, prob=bed$size)]

    # Randomly sample uniformly within each chosen range
    simulatedBps[, position := sample(start:end, size=1), by=1:dim(simulatedBps)[1]]
  }


  simulatedBps <- simulatedBps %>%
    mutate(pos = position) %>%
    select(chrom, pos)  %>%
    mutate(chrom = as.factor(chrom))

  return(simulatedBps)
}

