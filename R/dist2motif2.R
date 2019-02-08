# Functions to calculate the distance
# from each breakpoint to user-provided loci (e.g. TSS)

#' generateData2
#' Prepare data for dist2motif
#' @keywords simulate
#' @import ggplot2
#' @import dplyr
#' @import colorspace
#' @import RColorBrewer
#' @export
generateData2 <- function(..., df, breakpoints, sim=FALSE, chroms=c('2L', '2R', '3L', '3R', '4', 'X', 'Y')){
  if(missing(df) && missing(breakpoints)) stop("\n[!] Must provide either a df or bed file! Exiting.")

  if(!missing(df)) {
    cat("Reading data from df\n")
    real_data <- df %>% 
      dplyr::filter(...) %>% 
      dplyr::mutate(pos = bp) %>%
      dplyr::select(chrom, pos)
    
  } else if(!missing(breakpoints)) {
      cat("Reading data from bed file\n")
      real_data <- read.table(breakpoints, header = F)
      if(ncol(real_data)>=3){
        # real_data <- real_data[c(1,2)]
        # real_data$end <- real_data[2] + 2
        real_data <- real_data[c(1,2,3)]
        colnames(real_data) <- c("chrom", "start", "end")
        real_data <- real_data %>% 
          tidyr::gather(stend, val, -chrom) %>% 
          dplyr::mutate(pos=val) %>% 
          dplyr::select(chrom, pos)
      } else if(ncol(real_data)==2) {
        colnames(real_data) <- c("chrom", "pos")
      } else stop("Badly formatted bed file. Exiting")
  }
  real_data <- real_data %>% 
    dplyr::filter(chrom %in% chroms) %>% 
    droplevels()
  
  # if(missing(breakpoints)){
  #   cat("Using data from svBreaks::getData()\n")
  #   real_data <- getData(...) %>%
  #     dplyr::filter(chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" ) %>%
  #     dplyr::mutate(pos = bp) %>%
  #     dplyr::select(chrom, pos) %>%
  #     droplevels()
  # } else{
  #   cat("Using data from ", breakpoints, "\n")
  #   real_data <- read.table(breakpoints, header = F)
  #   if(is.null(real_data$V3)){
  #     real_data$V3 <- real_data$V2 + 2
  #   }
  #   colnames(real_data) <- c("chrom", "start", "end")
  #   real_data <- real_data %>% 
  #     dplyr::filter(chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" ) %>%
  #     dplyr::mutate(pos = (end+start)/2) %>%
  #     dplyr::select(chrom, pos) %>% 
  #     droplevels()
  # }
  
  if(sim) {
    byIteration <- list()
    #run each iteration
    for (i in 1:sim){
      cat("Running simulation", i, "of", sim, "\n")
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
    real_data$iteration <- as.factor(1)
    return(real_data)
  }
}


#' dist2Motif2
#' Calculate the distance from each breakpoint to closest motif in a directory of files
#' @keywords motif
#' @import ggplot2 dplyr tidyr RColorBrewer
#' @export
dist2motif2 <- function(..., df, breakpoints, feature_file, featureDir = system.file("extdata", "features", package="svBreaks"), chroms=c('2L', '2R', '3L', '3R', '4', 'X', 'Y'), sim=FALSE, position = 'centre') {
  if(missing(df) && missing(breakpoints)) stop("\n[!] Must provide either a df or bed file! Exiting.")
  if(!dir.exists(featureDir)){
    cat("Not a dir")
  }
  if(length(grep(list.files(featureDir), pattern = '.bed', value=T)) == 0) stop("\n[!] Must provide a directory containing bed files! Exiting.") else print(grep(list.files(featureDir), pattern = '.bed', value=T))
  
  bp_data <- svBreaks::generateData2(..., df=df, breakpoints=breakpoints, sim=sim, chroms=chroms)

  cat("Calculating distances to", position, 'of regions', sep = " ", "\n")
 
  # svCount <- table(bp_data$chrom)
  # bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 5]))
  # bp_data <- droplevels(bp_data)

  minDist <- function(p) {
    index <- which.min(abs(tss_df$pos - p))
    closestTss <- tss_df$pos[index]
    chrom <- as.character(tss_df$chrom[index])
    dist <- (p - closestTss)
    list(p, closestTss, dist, chrom)
  }

  scores <- list()

  add_col <- function(x){
    if(is.null(x$V3)){
      x$V3 <- x$V2 + 2
    }
    return(x)
  }
  
  if(!missing(feature_file)){
    fileNames <- list(feature_file)
  } else{
    fileNames <- list.files(featureDir, pattern = ".bed")
  }
  # cat("Analysing all files in directory:", bedFiles, "\n")
  for (i in 1:length(fileNames)){
    filename <- basename(tools::file_path_sans_ext(fileNames[i]))
    parts <- unlist(strsplit(filename, split = '\\.'))
    feature <- parts[1]

    cat("Analysing file:", filename, 'with feature:', feature, "\n")

    feature_locations <- read.table(paste(featureDir, fileNames[i], sep='/'), header = F)
    feature_locations <- add_col(feature_locations)
    feature_locations <- feature_locations[,c(1,2,3)]
    colnames(feature_locations) <- c("chrom", "start", "end")

    # fCount <- table(feature_locations$chrom)
    #
    # bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 5])) %>% droplevels()
    #

    feature_locations <- feature_locations %>%
      dplyr::filter(chrom %in% levels(bp_data$chrom)) %>% 
      droplevels()

    if(length(levels(feature_locations$chrom)) < length(levels(bp_data$chrom))){
      cat("Feature:", feature, "not found on all chroms. Trimming real data to match chroms in ", filename, paste0("[",levels(feature_locations$chrom), "]"), "\n")
      bp_data <- bp_data %>% 
        dplyr::filter(chrom %in% levels(feature_locations$chrom)) %>% 
        droplevels()
    }

    if(position == 'centre'){
      feature_locations <- feature_locations %>%
        dplyr::mutate(middle = as.integer(((end+start)/2)+1)) %>%
        dplyr::mutate(pos = as.integer(middle-1)) %>%
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
#' @import dplyr ggplot2 ggpubr RColorBrewer
#' @export
distOverlay2 <- function(..., df, breakpoints, feature_file, featureDir = system.file("extdata", "features", package="svBreaks"),
                         from='bps', chroms=c('2L', '2R', '3L', '3R', '4', 'X', 'Y'),
                         lim=2.5, n=2, plot = TRUE, position = 'centre', in_window=5, write=FALSE, out_dir) {
  
  scaleFactor <- lim*1e3

  real_data <- svBreaks::dist2motif2(..., df=df, breakpoints=breakpoints, feature_file=feature_file, featureDir=featureDir, position=position, chroms=chroms)
  sim_data <- svBreaks::dist2motif2(..., df=df, breakpoints=breakpoints, feature_file=feature_file, featureDir=featureDir, sim=n, position=position, chroms=chroms)
  
  in_range <- plyr::round_any(in_window/100*nrow(real_data), 10)
  
  bps_within_range <- real_data %>% 
    dplyr::group_by(feature) %>% 
    dplyr::mutate(min_count = sum(abs(min_dist)<=scaleFactor)) %>% 
    dplyr::filter(min_count >= in_range) %>% 
    dplyr::ungroup() %>% 
    droplevels()
  
  removed_features <- real_data %>% dplyr::filter(!feature %in% levels(bps_within_range$feature)) %>% droplevels()
  removed_features <- levels(removed_features$feature)
  if(length(removed_features) > 0){
    cat("Dropping features:", paste0("'", removed_features, "'"), "with fewer than 10 hits wihin +/-", lim, "Kb\n")
    if(length(levels(bps_within_range$feature))==0){
      stop("There are no feaures in dir that have >= ", in_window, "breakpoints within specified range ",paste0("(lim=", lim,"Kb) "), "Exiting.")
    }
    sim_data <- sim_data %>% 
      dplyr::filter(feature %in% levels(bps_within_range$feature)) %>% 
      droplevels()
  }
  bps_in_window <- real_data[abs(real_data$min_dist)<=scaleFactor,]
  perc_in_window <- plyr::round_any((nrow(bps_in_window)/nrow(real_data))*100, 1)
  cat("There are ", paste0(nrow(bps_in_window), "/", nrow(real_data), " [", perc_in_window, "%", "]"), "breakpoints within specified range", paste0("(lim=", lim,"Kb) "),"\n")
  
  
  real_data <- bps_within_range
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
  pVals_and_df <- svBreaks::simSig2(r = real_data, s = sim_data, max_dist = scaleFactor)

  combined <- pVals_and_df[[1]] %>% ungroup()
  pVals <- pVals_and_df[[2]] %>% ungroup()
  
  if(write){
    if(missing(out_dir) || !dir.exists(out_dir)){
      stop("Must specify an existing directory to write data to. Exiting")
    }

    bps_in_window <- bps_in_window %>% 
      dplyr::mutate(start = as.numeric(as.character(bp))-1,
                    end = as.numeric(as.character(bp))+1,
                    min_dist = paste0(feature, ":", min_dist)) %>%
      dplyr::filter(iteration == 1) %>% 
      dplyr::select(chrom, bp, end, min_dist, feature) %>% 
      droplevels()
    svBreaks::writeBed(df = bps_in_window, name = paste0("Bps_", lim, "kb_", bps_in_window$feature, ".bed")[1], outDir = out_dir)
    
    for(i in 1:(length(levels(combined$iteration)))){
      df_by_iter <- combined %>% dplyr::filter(iteration == i) %>% droplevels()
      for(s in levels(combined$Source)){
        df_by_source <- df_by_iter %>% 
          dplyr::filter(Source == s) %>% 
          dplyr::mutate(start = as.numeric(as.character(bp))-1,
                        end = as.numeric(as.character(bp))+1,
                        min_dist = paste0(feature, ":", min_dist)) %>% 
          dplyr::select(chrom, start, end, min_dist, Source, iteration, feature) %>% 
          droplevels()
        svBreaks::writeBed(df = df_by_source, name = paste0(df_by_source$Source, "_", df_by_source$feature, "_", df_by_source$iteration, ".bed")[1], outDir = out_dir)
      }
    }
    
    # combined %>%
    #   dplyr::group_by(Source, iteration) %>%
    #   dplyr::mutate(end = as.integer(bp)+2) %>%
    #   dplyr::select(chrom, bp, end, Source, iteration) %>%
    #   svBreaks::writeBed(df= ., name = paste0(.$Source, "_", .$iteration, ".bed")[1], outDir = '~/Desktop')
  }

  if(plot){
    print(plotdistanceOverlay2(..., distances=combined, from=from, facetPlot=FALSE, byChrom=byChrom, lim=lim, n=n, position=position ))
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
plotdistanceOverlay2 <- function(..., distances, from='bps', lim=2.5, n=2, position='centre', histo=FALSE, binWidth = 500){
  grDevices::pdf(NULL)
  
  scaleFactor <- lim*1e3
  scale <- "(Kb)"
  
  lims <- c(as.numeric(paste0("-", scaleFactor)), scaleFactor)
  brks <- c(as.numeric(paste0("-", scaleFactor)),
            0,
            scaleFactor)
  labs <- as.character(brks/1e3)
  expnd <- c(0, 0)
  
  new <- distances %>% 
    dplyr::mutate(iteration = as.factor(ifelse(Source=='Real', 0, iteration))) %>% 
    dplyr::filter(abs(min_dist)<=scaleFactor)
    
  real_fill <- '#3D9DEB'
  iterFill <- colorspace::rainbow_hcl(n)
  
  colours <- c(real_fill, iterFill)
  
  plts <- list()
  for (i in 1:(length(levels(new$feature)))){
    distances <- new %>% 
      dplyr::filter(feature == levels(new$feature)[i])
    p <- ggplot(distances)
    if(histo) {
      
      # distances <- distances %>% 
      #   dplyr::group_by_(feature, iteration) %>% 
      #   dplyr::mutate(newCount = )
      #   
      
      p <- p + geom_histogram(data=distances[distances$Source=="Sim",], aes(min_dist, y=(..count..), fill = Source, group = iteration), alpha = 0.1, binwidth = binWidth,  position="identity")
      p <- p + geom_histogram(data=distances[distances$Source=="Real",], aes(min_dist, y=..count.., fill = Source, group = iteration), alpha = 0.5, binwidth = binWidth, position="identity")
      
      # ggplot(data=distances[distances$Source=="Real",]) + geom_bar(mapping = aes(min_dist, y=..count.., fill = Source, group = iteration), alpha = 0.1, binwidth = binWidth, stat = "bin")
      # ggplot(data=distances[distances$Source=="Sim",]) + geom_bar(mapping = aes(min_dist, y=..count../1e3, fill = Source, group = iteration), alpha = 0.5, binwidth = binWidth, stat = "bin")
      # 
      
      p <- p + scale_fill_manual(values=colours)
      p <- p + scale_y_continuous(paste("Count per", binWidth, "bp bins"))
    } else {
      p <- p + geom_line(data=distances[distances$Source=="Real",], aes(min_dist, colour = iteration), size=2, stat='density')
      # p <- p + geom_density(data=d[d$Source=="Real",], aes(min_dist, y=..scaled..), fill = real_fill, alpha=0.2, adjust=3)
      p <- p + geom_line(aes(min_dist, group = interaction(iteration, Source), colour = iteration), alpha = 0.7, size=1, stat='density')
      
      p <- p + geom_rug(data=distances[distances$Source=="Real",], aes(min_dist, colour = iteration), sides = "b")
      p <- p + geom_rug(data=distances[distances$Source=="Sim",], aes(min_dist, colour = iteration), sides = "t")
      
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
        axis.text.x = element_text(size = 12),
        axis.line.y = element_line(color = "black", size = 0.5),
        plot.title = element_text(size=22, hjust = 0.5)
      )
    p <- p + labs(title = paste(distances$feature[1], "\n", position))
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
      dplyr::group_by(feature, iteration) %>%
      dplyr::mutate(count = n(),
                    median = median(min_dist),
                    mean = mean(min_dist),
                    trimmed_mean = mean(min_dist, trim=0.35),
                    trimmed_median = median(min_dist, trim=0.35),
                    sd = sd(min_dist),
                    Source = Source)
      # This was reducing the data and breaking in the stats loop...
      # dplyr::filter(abs(min_dist) <= max_dist ) %>%
      
    
    # p95 <- quantile(x$min_dist, 0.85)
    # p05 <- quantile(x$min_dist, 0.25)
    # # d <- x %>% dplyr::arrange(-min_dist)
    # trimmed_data <- x[x$min_dist<= p95 & x$min_dist >= p05,] %>% 
    #   dplyr::arrange(-min_dist)
    # x$trimmed_mean <- mean(abs(trimmed_data$min_dist))
    # mean(x$min_dist[x$min_dist<= p95 & x$min_dist >= p05])
    # 
    # mean(x[which(x$min_dist <= p95 & x$min_dist >= p05)])
    # meanTrunc <- mean(numVec[which(numVec <= p95 & numVec >= p05)])
    
    return(x)
  }
  # r <- r %>% dplyr::filter(iteration==1)
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
      # result1 <- suppressWarnings(ks.test(rl$min_dist, sm$min_dist))
      # if(!is.na(result1)){
      ksPval <- round(result1$p.value, 4)
      # }else{
      #   ksPval <- 1
      # }
      
      result2 <- car::leveneTest(df$min_dist, df$Source, center='median')
      result3 <- stats::bartlett.test(df$min_dist, df$Source)
      bPval <- round(result3$p.value, 4)
      lPval <- round(result2$`Pr(>F)`[1], 4)
      rmed <- round(median(rl$min_dist)/1e3, 2)
      smed <- round(median(sm$min_dist)/1e3, 2)
      # rtrim <- round(median(rl$min_dist, trim=0.25)/1000, 2)
      rtrim <- rl$trimmed_median[1]/1e3
      strim <- sm$trimmed_median[1]/1e3
    
      # strim <- round(median(sm$min_dist, trim=0.25)/1000, 2)
      rsd <- round(sd(rl$min_dist)/1e3, 2)
      ssd <- round(sd(sm$min_dist)/1e3, 2)
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
                         real_trim_med = rtrim,
                         sim_trim_med = strim,
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


add_col <- function(x){
  if(is.null(x$V3)){
    x$V3 <- x[,2] + 2
  }
  return(x)
}

writeBed <- function(df, outDir=NULL, name='regions.bed', svBreaks=FALSE){

  if(missing(outDir)){
    outDir <- getwd()
    cat("Writing to", outDir, "\n")
  }
  df <- add_col(df)
  # colnames(df[,c(1,2,3)]) <- c("chrom", "start", "end")
  # names(df)[1:3] <- c("chrom", "start", "end")
  

  df <- df %>%
    filter(as.numeric(df[2]) < as.numeric(df[3])) %>%
    droplevels()

  cat(paste(outDir,name, sep='/'), "\n")

  write.table(df, file = paste(outDir,name, sep='/'), row.names=F, col.names=F, sep="\t", quote = FALSE)
}

