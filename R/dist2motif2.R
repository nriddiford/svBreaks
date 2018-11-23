#' dist2Motif2
#' Calculate the distance from each breakpoint to closest motif in a directory of files
#' @keywords motif
#' @import ggplot2 dplyr tidyr RColorBrewer
#' @export
dist2motif2 <- function(..., breakpoints=NA, feature_file = NA, featureDir = 'rawdata/features/', sim=NA, keep=NULL, position = 'centre') {
  
  bp_data <- generateData(..., confidence=='precise', breakpoints=breakpoints, sim=sim, keep=keep)
  
  cat("Calculating distances to", position, 'of regions', sep = " ", "\n")
  
  svCount <- table(bp_data$chrom)
  bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 5]))
  # bp_data <- droplevels(bp_data)

  minDist <- function(p) {
    index <- which.min(abs(tss_df$pos - p))
    closestTss <- tss_df$pos[index]
    chrom <- as.character(tss_df$chrom[index])
    dist <- (p - closestTss)
    list(p, closestTss, dist, chrom)
  }
  
  scores <- list()

  fileNames <- dir(featureDir, pattern = ".bed")
  # cat("Analysing all files in directory:", bedFiles, "\n")
  for (i in 1:length(fileNames)){
    filename <- basename(tools::file_path_sans_ext(fileNames[i]))
    parts <- unlist(strsplit(filename, split = '\\.'))
    feature <- parts[1]
    
    cat("Analysing file:", fileNames[i], 'with feature:', feature, "\n")
    
    feature_locations <- read.table(paste(featureDir, fileNames[i], sep='/'), header = F)
    feature_locations <- feature_locations[,c(1,2,3)]
    colnames(feature_locations) <- c("chrom", "start", "end")
    
    # fCount <- table(feature_locations$chrom)
    # 
    # bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 5]))
    # 
    
    feature_locations <- feature_locations %>% 
      dplyr::filter(chrom %in% levels(bp_data$chrom))
    
  
    if(position == 'centre'){
      feature_locations <- feature_locations %>% 
        dplyr::mutate(end = as.integer(((end+start)/2)+1)) %>%
        dplyr::mutate(pos = as.integer(end-1)) %>%
        dplyr::select(chrom, pos)
    } else if(position == 'edge'){
      feature_locations <- feature_locations %>% 
        tidyr::gather(c, pos, start:end, factor_key=TRUE) %>% 
        dplyr::select(chrom, pos)
    }
    byIteration <- list()
    for (j in levels(bp_data$iteration)){
      byChrom <- list()
      df1 <- dplyr::filter(bp_data, iteration == j)
      for (c in levels(bp_data$chrom)) {
        df <- dplyr::filter(df1, chrom == c)
        tss_df <- dplyr::filter(feature_locations, chrom == c)
        dist2tss <- lapply(df$pos, minDist)
        dist2tss <- do.call(rbind, dist2tss)
        new <- data.frame(matrix(unlist(dist2tss), nrow=nrow(df)))
        new$iteration <- j
        new$feature <- as.factor(feature)
        colnames(new) <- c("bp", "closest_tss", "min_dist", "chrom", "iteration", "feature")
        byChrom[[c]] <- new
      }
      perIter <- do.call(rbind, byChrom)
      byIteration[[j]] <- perIter
    }
    dist2feat <- do.call(rbind, byIteration)
    scores[[i]] <- dist2feat
  }
  
  final <- do.call(rbind, scores)
  rownames(final) <- NULL
  final$iteration <- as.factor(final$iteration)
  final$chrom <- as.character(final$chrom)
  final$min_dist <- as.numeric(as.character(final$min_dist))
  

  return(final)
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
distOverlay2 <- function(..., breakpoints = NA, featureDir = 'rawdata/features/', from='bps', lim=2.5, n=2, plot = TRUE, keep=NULL, position = 'centre') {
  scaleFactor <- lim*1000
  real_data <- dist2motif2(..., breakpoints = breakpoints, featureDir = featureDir, keep=keep, position = position)
  sim_data <- dist2motif2(..., featureDir = featureDir, sim = n, position = position)
  
  real_data$Source <- as.factor("Real")
  sim_data$Source <- as.factor("Sim")
  
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
  pVals_and_df <- simSig2(r = real_data, s = sim_data, max_dist = scaleFactor)
  
  combined <- pVals_and_df[[1]]
  pVals <- pVals_and_df[[2]]
  
  if(plot==T){
    print(plotdistanceOverlay2(..., d=combined, from=from, facetPlot=FALSE, byChrom=byChrom, lim=lim, n=n, position=position ))
    print(pVals)
  }else{
    print(pVals)
    return(list(combined, pVals))
  }
}


#' plotdistanceOverlay
#'
#' Plot the distance overlay 
#' @param d Dataframe containing combined real + sim data (d <- distOverlay())
#' @import dplyr ggplot2 RColorBrewer scales colorspace cowplot
#' @keywords distance
#' @export
plotdistanceOverlay2 <- function(..., d, from='bps', lim=2.5, n=2, position='centre', histo=FALSE, binWidth = 500){
  grDevices::pdf(NULL)
  
  scaleFactor <- lim*1000
  scale <- "(Kb)"
  
  lims <- c(as.numeric(paste("-", scaleFactor, sep = '')), scaleFactor)
  brks <- c(as.numeric(paste("-", scaleFactor, sep = '')),
            as.numeric(paste("-", scaleFactor/10, sep = '')),
            scaleFactor/10,
            scaleFactor)
  labs <- as.character(brks/1000)
  expnd <- c(0, 0)
  
  new <- d %>% 
    mutate(iteration = as.factor(ifelse(Source=='Real', 0, iteration)))
    
  real_fill <- '#3D9DEB'
  iterFill <- colorspace::rainbow_hcl(n)
  
  colours <- c(real_fill, iterFill)

  plts <- list()
  for (i in 1:(length(levels(new$feature)))){
    d <- new %>% 
      filter(feature == levels(new$feature)[i])
    
    p <- ggplot(d)
    if(histo) {
      p <- p + geom_histogram(data=d[d$Source=="Sim",], aes(min_dist, fill = Source, group = iteration), alpha = 0.1, binwidth = binWidth,  position="identity")
      p <- p + geom_histogram(data=d[d$Source=="Real",], aes(min_dist, fill = Source, group = iteration), alpha = 0.5, binwidth = binWidth, position="identity")
      p <- p + scale_fill_manual(values=colours)
      p <- p + scale_y_continuous(paste("Count per", binWidth, "bp bins"))
    } else {
      p <- p + geom_line(data=d[d$Source=="Real",], aes(min_dist, colour = iteration), size=2, stat='density')
      p <- p + geom_line(aes(min_dist, group = interaction(iteration, Source), colour = iteration), alpha = 0.7, size=1, stat='density')
      p <- p + scale_color_manual(values=colours)
    }
    
    p <- p + scale_x_continuous(
      limits = lims,
      breaks = brks,
      expand = expnd,
      labels = labs
    )
    p <- p + 
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.text.x = element_text(size = 16),
        axis.line.y = element_line(color = "black", size = 0.5),
        plot.title = element_text(size=22, hjust = 0.5)
      )
    p <- p + labs(title = paste(d$feature, "\n", position))
    plts[[i]] <- p
  }
  cat("Plotting", length(levels(new$feature)), "plots", "\n")
  grDevices::dev.off()
  cowplot::plot_grid(plotlist=plts)
}


simSig2 <- function(r, s, test=NA, max_dist=5000){
  cat("Calculating descriptive statistics\n")
  arrange_data <- function(x){
    x <- x %>% 
      group_by(iteration, feature) %>%
      dplyr::mutate( count = n(),
                     median = median(min_dist),
                     mean = mean(min_dist),
                     sd = sd(min_dist),
                     Source = Source) %>%
      dplyr::filter(abs(min_dist) <= max_dist ) %>%
      ungroup()
    return(x)
  }
  simulated <- arrange_data(s)
  real <- arrange_data(r)
  
  combined <- suppressWarnings(dplyr::full_join(real, simulated))
  combined$Source <- as.factor(combined$Source)
  
  simbyFeat = list()
  for (f in levels(combined$feature)){
    pVals = list()
    c <- dplyr::filter(combined, feature==f)
    for(i in levels(c$iteration)){
      df <- dplyr::filter(c, iteration==i)
      rl <- dplyr::filter(df, Source == "Real")
      sm <- dplyr::filter(df, Source == "Sim")
      result1 <- tryCatch(suppressWarnings(ks.test(rl$min_dist, sm$min_dist)), error=function(err) NA)
      result1 <- suppressWarnings(ks.test(rl$min_dist, sm$min_dist))
      ksPval <- round(result1$p.value, 4)
      
      result2 <- car::leveneTest(df$min_dist, df$Source, center='median')
      result3 <- stats::bartlett.test(df$min_dist, df$Source)
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
      # fStat <- var.test(min_dist ~ Source , df, alternative = "two.sided")
      # fRatio <- round(fStat$statistic, 2)
      # fStat <- round(fStat$p.value, 4)
      
      sig <- ifelse(lPval <= 0.001, "***",
                    ifelse(lPval <= 0.01, "**",
                           ifelse(lPval <= 0.05, "*", "")))
      
      vals <- data.frame(iteration = i,
                         feature = f,
                         KS = ksPval,
                         Levenes = lPval,
                         # Bartlett = bPval,
                         # Fstat_ratio = fRatio,
                         # Fstat = fStat,
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
  simbyFeat[[f]] <- pVals_df
  }
  
  combined_sig_vals <- do.call(rbind, simbyFeat)
  
  rownames(combined_sig_vals) <- NULL
  combined_sig_vals <- combined_sig_vals %>%
    arrange(Levenes, KS)
  
  # print(pVals_df, row.names = FALSE)
  
  ## Boxplot per chrom
  
  # colours <- c("#E7B800", "#00AFBB")
  # cat("Plotting qq plot of min distances\n")
  # qqnorm(combined$min_dist)
  # qqline(combined$min_dist, col = 2)
  
  # p <- ggplot(combined)
  # p <- p + geom_boxplot(aes(chrom, min_dist, fill = Source), alpha = 0.6)
  # p <- p + scale_y_continuous("Distance", limits=c(-5000, 5000))
  # p <- p + facet_wrap(~iteration, ncol = 2)
  # p <- p + scale_fill_manual(values = colours)
  
  # p
  return(list(combined, combined_sig_vals))
}

