# Functions to calculate the distance
# from each breakpoint to user-provided loci (e.g. TSS)

#' generateData
#' Prepare data for dist2motif
#' @keywords simulate
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @export

generateData <- function(..., breakpoints=NA, sim=NA){
  
  if(is.na(breakpoints)){
    real_data <- getData(..., genotype=='somatic_tumour', !sample %in% c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"))
    # real_data <- notchFilt(keep=0)
    real_data <- real_data %>%
      dplyr::filter(chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" ) %>%
      dplyr::mutate(pos = bp) %>%
      dplyr::select(chrom, pos) %>%
      droplevels()
  } else{
    real_data <- read.table(breakpoints, header = F)
    if(is.null(real_data$V3)){
      real_data$V3 <- real_data$V2 + 2
    }
    colnames(real_data) <- c("chrom", "start", "end")
    real_data <- real_data %>% 
      dplyr::filter(chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" ) %>%
      dplyr::mutate(pos = (end+start)/2) %>%
      dplyr::select(chrom, pos) %>% 
      droplevels()
    
  }
  
  if (!is.na(sim)) {
    byIteration <- list()
    #run each iteration
    for (i in 1:sim){
      cat("Running simulation", i, "\n")
      simByChrom <- list()

      for (c in levels(real_data$chrom)){
        hitCount <- nrow(real_data[real_data$chrom== c,])
        hitCount <- (hitCount*10)
        if (i == 1){
          cat(paste("Simulating", hitCount, "breakpoints on chromosome", c), "\n")
        }
        bp_data <- bpSim(nSites = hitCount, byChrom = c)
        bp_data$iteration <- i
        simByChrom[[c]] <- bp_data
      }
      result <- as.data.frame(do.call(rbind, simByChrom))
      rownames(result) <- NULL
      byIteration[[i]] <- result
    }

    #combine each iteration into one data frame
    # final <- dplyr::bind_rows(byIteration)
    final <- as.data.frame(do.call(rbind, byIteration))
    final$iteration <- as.factor(final$iteration)

    return(final)
  } else{
    cat("Using real data", "\n")
    real_data$iteration <- as.factor(1)
    return(real_data)
  }
}


#' dist2Motif
#' Calculate the distance from each breakpoint to closest motif
#' @keywords motif
#' @import ggplot2
#' @import dplyr
#' @import RColorBrewer
#' @export

dist2Motif <- function(..., breakpoints=NA,feature_file = system.file("extdata", "tss_locations.txt", package="svBreaks"), sim=NA, print=0, send=0, feature="tss") {
  # df : chrom pos iteration
  bp_data <- generateData(..., breakpoints=breakpoints, sim=sim)

  feature <- paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep = "")

  if (feature == "Promoter") {
    feature_locations <- getPromoter()
    cat("Getting gene promoter locations...\n")
  } else {
    feature_locations <- read.table(feature_file, header = F)
    cat("Reading in file:", feature_file, sep = " ", "\n")
  }

  cat("Calculating distances to", feature, sep = " ", "\n")
  
  # Not a good idea if the file supplies single nucleotide data as will shift by 1 ->
  if(is.null(feature_locations$V3)){
    feature_locations$V3 <- feature_locations$V2 + 2
  }
  colnames(feature_locations) <- c("chrom", "start", "end")
  
  feature_locations <- feature_locations %>% 
    dplyr::mutate(pos = as.integer((end+start)/2)) %>%
    dplyr::select(chrom, pos)
# 
#   feature_locations <- feature_locations[,c(1,2)]
# 
#   colnames(feature_locations) <- c("chrom", "pos")

  # feature_locations$pos <- as.integer(feature_locations$pos)

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
    dist <- (p - closestTss)
    list(p, closestTss, dist, chrom)
  }

  byIteration <- list()

  for (i in levels(bp_data$iteration)){
    byChrom <- list()
    df1 <- dplyr::filter(bp_data, iteration == i)

    for (c in levels(bp_data$chrom)) {
      df <- dplyr::filter(df1, chrom == c)
      tss_df <- dplyr::filter(feature_locations, chrom == c)

      dist2tss <- lapply(df$pos, fun2)
      dist2tss <- do.call(rbind, dist2tss)

      new <- data.frame(matrix(unlist(dist2tss), nrow=nrow(df)))
      new$iteration <- i
      colnames(new) <- c("bp", "closest_tss", "min_dist", "chrom", "iteration")

      byChrom[[c]] <- new
    }
    perIter <- do.call(rbind, byChrom)
    byIteration[[i]] <- perIter
  }
  final <- do.call(rbind, byIteration)
  rownames(final) <- NULL
  final$iteration <- as.factor(final$iteration)
  final$chrom <- as.character(final$chrom)
  final$min_dist <- as.numeric(as.character(final$min_dist))
  dist2tss <- final

  if (send == 1) {
    return(dist2tss)
  } else {
    # here we just want to make a plot for the simulated data so just show for iteration 1
    dist2tss <- dplyr::filter(dist2tss, iteration==1)
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
    # p <- p + facet_grid(chrom ~ iteration, )
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

    if (is.na(sim)) {
      distout <- paste("bp", feature, "dist.pdf", sep = "")
    } else {
      distout <- paste("bp", feature, "dist_sim.pdf", sep = "")
    }

    cat("Writing file", distout, "\n")
    ggsave(paste("plots/", distout, sep = ""), width = 10, height = 10)

    p
  }
}

# distOverlay
#'
#' Calculate the distance from each breakpoint to closest motif
#' Overlay the same number of random simulated breakpoints
#' @keywords motif
#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import RColorBrewer
#' @export

distOverlay <- function(..., breakpoints = NA, feature_file=system.file("extdata", "tss_locations.txt", package="svBreaks"), feature="tss", lim=10, byChrom=NA, n=10) {
  feature <- paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep = "")

  scaleFactor <- lim*1000

  if (feature == "promoter") {
    real_data <- dist2Motif(send = 1, feature = feature)
    sim_data <- dist2Motif(feature = feature, sim = n, send = 1)
  }
  else {
    real_data <- dist2Motif(..., breakpoints=breakpoints, feature_file = feature_file, send = 1, feature = feature)
    sim_data <- dist2Motif(..., feature_file = feature_file, feature = feature, sim = n, send = 1)
  }

  real_data$Source <- "Real"
  sim_data$Source <- "Sim"

  dummy_iterations <- list()
  for (i in levels(sim_data$iteration)){
    real_data$iteration <- as.factor(i)
    dummy_iterations[[i]] <- real_data
  }
  real_data <- do.call(rbind, dummy_iterations)
  rownames(real_data) <- NULL

  real_data$iteration <- factor(real_data$iteration, levels = 1:n)
  sim_data$iteration <- factor(sim_data$iteration, levels = 1:n)

  # Perform significance testing
  combined <- simSig(r = real_data, s = sim_data, max_dist = scaleFactor)

  colours <- c("#E7B800", "#00AFBB")

  cat("Setting limits to -+", lim, "Kb\n", sep='')

  scale <- "(Kb)"

  lims <- c(as.numeric(paste("-", scaleFactor, sep = '')), scaleFactor)
  brks <- c(as.numeric(paste("-", scaleFactor, sep = '')),
            as.numeric(paste("-", scaleFactor/10, sep = '')),
            scaleFactor/10,
            scaleFactor)
  labs <- as.character(brks/1000)
  expnd <- c(.0005, .0005)


  # pp <- ggboxplot(combined, x = "Source", y = "min_dist",
  #                 color = "Source", palette = "jco",
  #                 add = "jitter",
  #                 facet.by = c("iteration", "chrom"))
  #  # Use only p.format as label. Remove method name.
  # pp + stat_compare_means(label = "p.format",
  #                         method = "t.test")
  # pp

  # pp <- ggdensity(combined, x = "min_dist",
  #                 color = "Source",
  #                 fill = "Source",
  #                 add = "median",
  #                 facet.by = "iteration",
  #                 rug = TRUE,
  #                 palette = c("#00AFBB", "#E7B800"))
  # pp
  #

  p <- ggplot(combined)
  p <- p + geom_density(aes(min_dist, fill = Source, group = Source), alpha = 0.4, adjust=0.5)
  
  # p <- p + geom_freqpoly(aes(min_dist, y = ..count../sum(..count..), group = Source, colour = Source))

  # p <- ggplot(combined, aes(sample = min_dist, colour = Source, group = Source))
  # p <- p + stat_qq()
  # 
  
  p <- p + scale_x_continuous(
    paste("Distance to", feature, scale, sep = " "),
    limits = lims,
    breaks = brks,
    expand = expnd,
    labels = labs
  )
  p <- p + scale_y_continuous("Density")
  p <- p + geom_vline(xintercept = 0, colour = "black", linetype = "dotted")

  p <- p + geom_rug(data = combined[combined$Source=='Real',], aes(min_dist, colour = Source), sides = "b")
  p <- p + geom_rug(data = combined[combined$Source=='Sim',], aes(min_dist, colour = Source), sides = "t")

  p <- p + scale_fill_manual(values = colours)
  p <- p + scale_colour_manual(values = colours)

  if (!is.na(byChrom)) {
    p <- p + facet_wrap(iteration~chrom, ncol = length(levels(as.factor(combined$chrom))))
  } else {
    p <- p + facet_wrap(~iteration)
  }

  p <- p + slideTheme() +
    theme(
      strip.text = element_text(size = 20),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12),
      legend.position = "top"
    )

  if(n>=20){
    p <- p + theme(
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.8, "lines")
      # axis.text = element_text(size = 12)
    )
    # p <- p + facet_wrap(~iteration, ncol = 5)
  }

  overlay <- paste("bp", feature, "_dist_overlay_", lim, "kb.png", sep = "")
  cat("Writing file", overlay, "\n")
  ggsave(paste("plots/", overlay, sep = ""), width = 20, height = 10)

  p
  # pp
}



#' simSig
#'
#' Calculate the median distance between breakpoints (per-chrom) and feature
#' and perform test
#' @param real_data Dataframe containing real data (as produced by generateData())
#' @import dplyr
#' @import ggplot2
#' @import RColorBrewer
#' @import broom
#' @import car
#' @import PerformanceAnalytics
#' @keywords sim
#' @export

simSig <- function(r, s, test=NA, max_dist=5000){

  # if(!(is.na(test))){
  #   ####
  #   # dummy data
  #   # set.seed(42)
  #   dummyReal <- data.frame(Source=rep('Real'), min_dist=sample(-1000:1000, 1000, replace=TRUE), iteration=rep(1:10))
  #   dummySim <- data.frame(Source=rep('Sim'), min_dist=sample(-2000:2000, 1000, replace=TRUE), iteration=rep(1:10))
  #   dummyCombined <- suppressWarnings(dplyr::full_join(dummyReal, dummySim))
  # 
  #   dummyCombined$Source <- as.factor(dummyCombined$Source)
  #   dummyCombined$iteration <- as.factor(dummyCombined$iteration)
  # 
  #   dummyCoc <- dummyCombined %>%
  #     group_by(iteration) %>%
  #     do(tidy(wilcox.test(min_dist ~ Source, exact=FALSE, data = . ))) %>%
  # 
  #     arrange(p.value) %>%
  #     ungroup()
  # 
  #   dummyCombined %>%
  #     group_by(iteration) %>%
  #     do(tidy(leveneTest(dummyCombined$min_dist, dummyCombined$Source, center = 'mean', data = .))) %>%
  # 
  #     arrange(p.value) %>%
  #     ungroup()
  # 
  # 
  #   p <- ggplot(dummyCombined)
  #   p <- p + geom_density(aes(min_dist, fill = Source, group = Source), alpha = 0.4)
  #   p
  # 
  #   return(dummyCombined)
  # }
  #####

  simulated <- s %>%
    group_by(chrom, iteration) %>%
    dplyr::mutate( count = n(),
            median = median(min_dist),
            mean = mean(min_dist),
            sd = sd(min_dist),
            Source = factor(Source)) %>%
    filter(abs(min_dist) <= max_dist ) %>%
    ungroup()

  real <- r %>%
    dplyr::group_by(chrom, iteration) %>%
    dplyr::mutate( count = n(),
                   median = median(min_dist),
                   mean = mean(min_dist),
                   sd = sd(min_dist),
                   Source = factor(Source)) %>%
    filter(abs(min_dist) <= max_dist ) %>%
    ungroup()

  combined <- suppressWarnings(dplyr::full_join(real, simulated))
  combined$Source <- as.factor(combined$Source)

  combined$iteration <- as.factor(combined$iteration)
  pVals = list()
  for(i in levels(combined$iteration)){
    df <- filter(combined, iteration==i)
    rl <- filter(df, Source == "Real")
    sm <- filter(df, Source == "Sim")
    result1 <- suppressWarnings(ks.test(rl$min_dist, sm$min_dist))
    ksPval <- round(result1$p.value, 4)
    

    result2 <- leveneTest(df$min_dist, df$Source, center='median')
    result3 <- bartlett.test(df$min_dist, df$Source)
    bPval <- round(result3$p.value, 4)
    lPval <- round(result2$`Pr(>F)`[1], 4)
    rmed <- round(median(rl$min_dist)/1000, 2)
    smed <- round(median(sm$min_dist)/1000, 2)
    rsd <- round(sd(rl$min_dist)/1000, 2)
    ssd <- round(sd(sm$min_dist)/1000, 2)
    rKurtosis <- round(kurtosis(rl$min_dist), 2)
    sKurtosis <- round(kurtosis(sm$min_dist), 2)
    rSkew <- round(skewness(rl$min_dist), 2)
    sSkew <- round(skewness(sm$min_dist), 2)
    fStat <- var.test(min_dist ~ Source , df, alternative = "two.sided")
    fRatio <- round(fStat$statistic, 2)
    fStat <- round(fStat$p.value, 4)
    
    
    sig <- ifelse(lPval <= 0.001, "***",
                      ifelse(lPval <= 0.01, "**",
                             ifelse(lPval <= 0.05, "*", "")))

    vals <- data.frame(iteration = i,
                       KS = ksPval,
                       Levenes = lPval,
                       # Bartlett = bPval,
                       Fstat_ratio = fRatio,
                       Fstat = fStat,
                       real_median = rmed,
                       sim_median = smed,
                       real_sd = rsd,
                       sim_sd = ssd,
                       real_kurtosis = rKurtosis,
                       sim_kurtosis = sKurtosis,
                       real_skew = rSkew,
                       sim_skew = sSkew,
                       sig = sig)
    pVals[[i]] <- vals

  }

  pVals_df <- do.call(rbind, pVals)
  rownames(pVals_df) <- NULL
  pVals_df <- pVals_df %>%
    arrange(Levenes, KS)
  
  print(pVals_df, row.names = FALSE)

  ## Boxplot per chrom

  colours <- c("#E7B800", "#00AFBB")

  # qqnorm(combined$min_dist);qqline(combined$min_dist, col = 2)

  p <- ggplot(combined)
  p <- p + geom_boxplot(aes(chrom, min_dist, fill = Source), alpha = 0.6)
  p <- p + scale_y_continuous("Distance", limits=c(-5000, 5000))
  p <- p + facet_wrap(~iteration, ncol = 2)
  p <- p + scale_fill_manual(values = colours)

  #p
  return(combined)
}

#' bpSim
#'
#' Generate simulated SV breakpoints acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random breakpoints to generate [Default nrow(bp_data)]
#' @import dplyr
#' @importFrom data.table fread as.data.table
#' @keywords sim
#' @export

bpSim <- function(nSites = 1e3, byChrom = NA, iterations = 10 ){
  intervals=system.file("extdata", "intervals.bed", package="svBreaks")
  bed <- data.table::fread(intervals)
  colnames(bed) <- c("chrom", "start", "end", "NA")

  bed <- bed %>%
    dplyr::select(chrom:end) %>%
    dplyr::mutate(size = (end-start)+1)  %>%
    as.data.table() %>%
    droplevels()

  if(is.na(byChrom)){
    # Randomly sample bed file rows, proportional to the length of each range
    simulatedBps <- bed[sample(.N, size=nSites, replace=TRUE, prob=bed$size)]

    # Randomly sample uniformly within each chosen range
    simulatedBps[, position := sample(start:end, size=1), by=1:dim(simulatedBps)[1]]
  } else {
    bed <- bed %>%
      dplyr::filter(chrom == byChrom) %>%
      as.data.table()

    # Randomly sample bed file rows, proportional to the length of each range
    simulatedBps <- bed[sample(.N, size=nSites, replace=TRUE, prob=bed$size)]

    # Randomly sample uniformly within each chosen range
    simulatedBps[, position := sample(start:end, size=1), by=1:dim(simulatedBps)[1]]
  }

  simulatedBps <- simulatedBps %>%
    dplyr::mutate(pos = position) %>%
    dplyr::select(chrom, pos)  %>%
    dplyr::mutate(chrom = as.factor(chrom))

  return(simulatedBps)
}

