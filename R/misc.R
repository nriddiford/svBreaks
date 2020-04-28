
#' Misc functions

#' transform_types
#'
#' Convert short types into plottable names
#' @export
transform_types <- function(x){
  x <- x %>% 
    dplyr::mutate(type2 = ifelse(type2 == "TRA", "Translocation",
                                 ifelse(type2 == "DEL", "Deletion",
                                        ifelse(type2 == "DUP", "Duplication",
                                               ifelse(type2 == "COMPLEX", "Complex",
                                                      ifelse(type2 == "TANDUP", "Tandem\nDuplication",
                                                             ifelse(type2 == "BND", "Inversion", type2)))))))
  return(x)
}


#' SV colours
#' 
#' @export
sv_colours <- function(){
  return(
    # c("DEL" = '#0073C299',
    #        "COMPLEX" = '#EFC00099',
    #        "DUP" = '#232323',
    #        "BND" = '#868686',
    #        "TRA" = '#CD534C99',
    #        "NA" = "grey")
         # c("DEL" = "#269FBF",
         #   "COMPLEX" = '#85C27C',
         #   "DUP" = '#E4E57E',
         #   "BND" = '#E68B4C',
         #   "TRA" = '#BC5D8F',
         #   "TANDUP" = '#A68BD1',
         #   "NA" = "grey"
         # )
         # 
         c("DEL" = "#269FBF",
           "COMPLEX" = '#85C27C',
           "DUP" = '#A68BD1',
           "BND" = '#E68B4C',
           "TRA" = '#BC5D8F',
           "TANDUP" = '#E4E57E',
           "NA" = "grey"
         )
         
  )
}


#' get sample names
#' 
#' @export
getMissingSamples <- function(..., df=NULL, all_samples_info='~/Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt'){
  all_samples_info <- read.delim(all_samples_info, header = F, stringsAsFactors = T)
  colnames(all_samples_info) <- c("sample", "marius", "paper", "sex", "assay")
  
  names <- all_samples_info %>% 
    dplyr::filter(...) %>% 
    droplevels()
  
  names <- as.list(levels(names$sample))
  missing_samples = list()
  
  r <- 0
  for(i in 1:length(names)){
    if(!(names[[i]] %in% as.list(levels(df$sample)))){
      r <- r + 1
      missing_samples[r] <- names[[i]]
    }
  }
  # print("samples missing from dataframe levels: %s", unlist(missing_samples))
  return(missing_samples) 
}


#' svTypes
#'
#' Plot counts of different svTypes
#' @import tidyverse ggsci
#' @export
svTypes <- function(..., bp_data=NULL, title=expression("Structural variants per sample (" ~italic("Notch")~ "excluded )"), object='type', keep=FALSE, plot=T) {
  if(missing(bp_data)) bp_data <- getData(..., genotype=='somatic_tumour')
  ext <- ".png"
  
  bp_data <- bp_data %>% 
    dplyr::filter(...,
                  bp_no != "bp2") %>% 
    dplyr::group_by(sample, type2) %>% 
    dplyr::summarise(type_count = n()) %>%
    dplyr::mutate(sv_count = sum(type_count)) %>% 
    dplyr::arrange(-sv_count) %>% 
    droplevels() %>% 
    as.data.frame()
  
  # cols <- setCols(df=dat, col=object, set='Pastel1')
  
  missing_samples <- getMissingSamples(..., df = non_notch)
  
  
  missing_samples <- plyr::compact(missing_samples)
  dat <- data.frame(sample = unlist(missing_samples), sv_count = 0, type_count =0, type2 = "NA")
  bp_data <- plyr::join(bp_data, dat, type='full')
  
  order <- levels(fct_reorder(bp_data$sample, bp_data$sv_count))
  colours <- sv_colours()
  
  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(fct_reorder(sample, -sv_count), type_count, colour = type2, fill = type2), alpha=0.6, stat = 'identity')
  p <- p + scale_fill_manual("SV type\n", values = colours)
  p <- p + scale_colour_manual(values = colours)
  p <- p + guides(color = FALSE)
  p <- p + slideTheme() +
    theme(
      axis.title.y = element_blank(),
      # panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      panel.grid.major.x = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
      legend.position = "top",
      axis.text.y = element_text(size=15),
      axis.title.x=element_blank()
    )

  if(plot){
    p
  } else {
    return(list(bp_data, p))
  }
}





#' featureDensity
#'
#' Plot density of several feature tyopes across the genome
#' @import dplyr
#' @import ggplot2
#' @export

featureDensity <- function(feature_file1 = system.file("extdata", "g4_positions.txt", package="svBreaks"), feature1 = 'G4',
                           feature_file2 = NA, feature2 = NA,
                           feature_file3 = NA, feature3 = NA,
                           write=TRUE, chroms = c("2L", "2R", "3L", "3R", "X")) {
  n = 1

  file_list = list()

  f1 <- read.delim(feature_file1, header = F)
  f1 <- f1[,c(1,2,3)]
  if(is.null(f1$V3)){
    f1$V3 <- f1$V2 + 2
  }
  
  f1$type <- feature1
  colnames(f1) <- c("chrom", "start", "end", "type")
  
  f1 <- f1 %>% 
    dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
    dplyr::select(chrom, mid, type)

  file_list[[n]] <- f1

  if (!is.na(feature_file2)){
    n = n + 1
    f2 <- read.delim(feature_file2, header = F)
    f2 <- f2[,c(1,2,3)]
    if(is.null(f2$V3)){
      f2$V3 <- f2$V2 + 2
    }
    
    f2$type <- feature2
    colnames(f2) <- c("chrom", "start", "end", "type")
    
    f2 <- f2 %>% 
      dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
      dplyr::select(chrom, mid, type)
    
    file_list[[n]]  <- f2
  }
  if (!is.na(feature_file3)){
    n = n + 1
    f3 <- read.delim(feature_file3, header = F)
    if(is.null(f3$V3)){
      f3$V3 <- f3$V2 + 2
    }
    
    f3$type <- feature3
    colnames(f3) <- c("chrom", "start", "end", "type")
    
    f3 <- f3 %>% 
      dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
      dplyr::select(chrom, mid, type)
    
    file_list[[n]]  <- f3
  }

  files <- do.call(rbind, file_list)
  rownames(files) <- NULL

  # chromosomes <- data.frame(chroms = c("2L", "2R", "3L", "3R", "X"), lengths = c(23513712, 25286936, 28110227, 32079331, 23542271))
                       
  locations <- files %>%
    dplyr::mutate(type = as.factor(type)) %>%
    dplyr::mutate(pos = as.numeric(mid/1000000)) %>%
    dplyr::filter(chrom %in% chroms) %>%
    dplyr::arrange(chrom, mid) %>%
    droplevels()

  p <- ggplot(locations)
  p <- p + geom_density(aes(pos, fill = type), alpha = 0.4, adjust=0.07)
  # p <- p + geom_rug(aes(pos, colour = type), sides = "tb", alpha = 0.05)
  p <- p + facet_wrap(~chrom~type, scales = "free_y", ncol=1)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, max(locations$pos), by = 5))
  p <- p + geom_rug(data=locations, aes(pos, colour = type), sides = "b")
  p <- p + cleanTheme() +
    theme(axis.text.y = element_blank())

  if(write){
    featurePlot <- paste("feature_dist.png")
    cat("Writing file", featurePlot, "\n")
    ggsave(paste("plots/", featurePlot, sep = ""), width = 20, height = 10)
  }
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
  } else {
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
  } else {
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

typeLenCount <- function(..., bp_data=NULL, size_threshold = 1, notch=0) {
  
  if(missing(bp_data)){
    if (notch) {
      bp_data <- notchFilt()
      ext <- "_excl.N.pdf"
    }
    else {
      bp_data <- getData(...)
      ext <- ".pdf"
    }
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
  p <- p + facet_wrap(~chrom, scales = "free_x", ncol = 2)
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
  p <- p + facet_wrap(~chrom, scales = "free_x", ncol = 2)
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



#' bootStrap
#'
#' Plot reav vs sim data bootstrapped n times
#' @import ggplot2
#' @import dplyr
#' @export

bootStrap <- function(..., feature_file=system.file("extdata", "tss_locations.txt", package="svBreaks"),
                      feature="tss", n=10) {
  
  
  nsim <- function(n, m = 0, s = 1) {
    z <- rnorm(n)
    m + s * ((z - mean(z)) / sd(z))
  }
  
  
  nboot <- function(x, R) {
    n <- length(x)
    m <- mean(x)
    s <- sd(x)
    do.call(rbind,
            lapply(1 : R,
                   function(i) {
                     xx <- sort(nsim(n, m, s))
                     p <- seq_along(x) / n - 0.5 / n
                     data.frame(x = xx, p = p, sim = i)
                   }))
  }
  
  
  
  real_data <- dist2Motif(..., feature_file = feature_file, send = 1, feature = feature)
  sim_data <- dist2Motif(..., feature_file = feature_file, feature = feature, sim = 1, send = 1)
  
  
  real_data <- real_data %>% 
    filter(abs(min_dist) <= 30000 )
    
  
  real <- nboot(real_data$min_dist, n)
  sim <- nboot(sim_data$min_dist, n)
  
  real$source <- 'real'
  sim$source <- 'sim'
  
  p <- ggplot()
  p <- p + geom_density(data=real, aes(x, fill = source), alpha = 0.4, adjust=1)
  p <- p + geom_density(data=sim, aes(x, fill = source), alpha = 0.4, adjust=1)
  
  # p <- p + geom_freqpoly(data = real,
  #               mapping = aes(x, y = ..count../sum(..count..), colour = source))
  # p <- p + geom_freqpoly(data = sim,
  #                        mapping = aes(x, y = ..count../sum(..count..), colour = source))
  
  
  p <- p + scale_x_continuous(
    paste("Distance to", feature, "(Kb)", sep = " "),
    limits = c(-20000, 20000),
    breaks = c(-20000, -2000, 2000, 20000),
    expand = c(.0005, .0005),
    labels = c("-20", "-2", "2", "20")
  )
  
  p <- p + geom_rug(data = real, aes(x, colour = source), sides = "b")
  p <- p + geom_rug(data = sim, aes(x, colour = source), sides = "t")
  
  p <- p + facet_wrap(~sim)
  p
 
   # ggplot() +
   #   geom_line(aes(x = qnorm(p), y = x, group = sim),
   #             color = "gray", data = real)
   # 
  
}

