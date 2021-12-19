#' Run regioneR for to regions in bed format
#' @param regionA RegionA [Required]
#' @param regionB RegionB [Required]
#' @param mappable_regions Bed file with mappable geneome [Required]
#' @param chrom_lengths Lengths for for chromosomes
#' @param slop Add this much slop to each region
#' @param limit2chrom Restrict permutaions to one chromsome [Default: TRUE] 
#' @param chrom Chromosome to run permutations on if limit2chrom=TRUE [Default: "X"] 
#' @param from Restrict to region between from and to
#' @param to Restrict to region between from and to
#' @param n Number of permutations to run
#' @param plot Return plot or not? [Default: TRUE] 
#' @param region Limit to region (chr:start-stop) [Default: NULL] 
#' @import regioneR
#' @export
bpRegioneR <- function(..., 
                       regionA = system.file("extdata", "notch_breakpoint_regions_500_mappable.bed", package="svBreaks"),
                       regionB = system.file("extdata", "motif_1.mappable.bed", package="svBreaks"),
                       mappable_regions = system.file("extdata", "dmel6_mappable.bed", package="svBreaks"),
                       chrom_lengths = system.file("extdata", "chrom.sizes.txt", package="svBreaks"),
                       slop, limit2chrom=TRUE, chrom='X', from, to, n = 100, plot = TRUE,  bins='auto', region=NULL){

  if(missing(region)){
    region_name <- tools::file_path_sans_ext(basename(regionA))
    feature_name <- tools::file_path_sans_ext(basename(regionB))
  }
  genome = read.table(chrom_lengths, header=F, stringsAsFactors = TRUE)
  mappable = read.table(mappable_regions, header=F, stringsAsFactors = TRUE)
  
  mappable <- mappable %>%
    dplyr::filter(V1 %in% levels(genome$V1)) %>% 
    dplyr::select(-V4) %>%
    droplevels()
  
  if(limit2chrom && ( !missing(from) || !missing(to) ) ){
    mappable <- mappable %>%
      dplyr::filter(V1 == chrom,
                    V2 >= (from - 5000),
                    V3 <= (to + 5000)) %>%
      droplevels()
    # genome$V2 = from
    # genome$V3 = to
    cat(paste0("o Genome: ", chrom, ":", from, ":", to), "\n")
  }
  
  test <- read.table(regionA, header=F, stringsAsFactors = TRUE)
  test <- test[,c(1,2,3)]
  test$V1 <- stringr::str_remove(test$V1, 'chr')
  
  test <- test %>% 
    dplyr::mutate(length = V3 - V2) %>% 
    dplyr::filter(...,
                  V1 %in% levels(genome$V1)) %>%
    dplyr::select(-length)
  
  cat(paste0("o Looking for enrichment of ", feature_name, " in ", nrow(test), " regions of ", region_name), "\n")
  cat("o Shuffling feature", feature_name, " ", n, "times\n")
  
  feature <- read.delim(regionB, header = F)
  feature <- feature[,c(1,2,3)]
  
  feature$V1 <- stringr::str_remove(feature$V1, 'chr')
  
  if(!missing(slop)){
    cat("o Expanding ", feature_name, " regions by ", slop, "\n")
    feature <- feature %>%
      dplyr::mutate(V2 = V2 - slop,
                    V3 = V3 + slop)
  }
  
  # mappable_genome <- getGenomeAndMask("dm6", mask=exclude)$genome
  # mappable_genome <- getGenomeAndMask(genome=genome, mask=exclude)
  # mappable_genome <- regioneR::overlapRegions(A = genome, B = mappable, type="BinA") %>% 
  #   dplyr::select(chr, startB, endB)
  
  # pt <- regioneR::overlapPermTest(A=test, B=feature, ntimes=n, genome=mappable_genome,  per.chromosome=TRUE)
  # plot(pt)
  
  # genome <- getGenomeAndMask(genome = 'dm6', mask = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed')
  res <- regioneR::permTest(A=test, B=feature, ntimes=n, randomize.function=regioneR::randomizeRegions,
                           genome=mappable, evaluate.function=regioneR::numOverlaps, per.chromosome=limit2chrom)
  
  d <- data.frame(feature = c('real', rep('shuffled', length(res$`regioneR::numOverlaps`$permuted))),
                  overlaps = c(res$`regioneR::numOverlaps`$observed, res$`regioneR::numOverlaps`$permuted)
  )
  
  title_string <- paste0("Region: ", region_name, "\n",
                         "Feature: ", feature_name, "\n",
                         "shuffles: ", res$`regioneR::numOverlaps`$ntimes, "\n",
                         "observed: ", res$`regioneR::numOverlaps`$observed, "\n",
                         "p-val: ", format.pval(res$`regioneR::numOverlaps`$pval) ,"\n")
  
  if(plot) {
    if(bins == 'auto') bins <- max(d$overlaps)
    print(plot_bpRegioner(d, title = title_string, bins = bins))
  } else{
    print(title_string)
    print(format.pval(res$`regioneR::numOverlaps`$pval))
    return(d)
  }
  
}

#' Plot overlaps from  bpRegioner
#' @export
plot_bpRegioner <- function(df, title=NULL, bins=30, rug=FALSE){
  
  if(missing(title)) title <- ''
  
  grey <- '#333333'
  green <- '#259FBF'
  cols <- c(grey, green)
  
  df$feature = factor(df$feature, levels=c("shuffled","real"), labels=c("shuffled", "real")) 
  
  compare_mean <- df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(med = median(overlaps))
  
  
  p <- ggplot(df, aes(overlaps, fill = feature, colour = feature)) +
    geom_histogram(alpha=0.6, bins = bins, position = "identity") +
    geom_vline(data = compare_mean, aes(xintercept = med, colour = feature),
               linetype = "dashed", size = .7) +
    cleanTheme() +
    theme(panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")) +
    scale_x_continuous("Overlaps", expand = c(0, 0.1)) +
    scale_y_continuous("Count", expand = c(0, 0)) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    ggtitle(title)
  
  print(p)
  
}