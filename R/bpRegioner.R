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
#' @param plot 
#' @param region
#' @param feature
#' 
#' @export
bpRegioneR <- function(..., regionA = '~/Desktop/misc_bed/breakpoints/Notch_CFS/notch_10kb_merged_breakpoints.bed',
                       regionB = '~/GitHub/BardinLab/meme/out/notch_CFS/fimo_out/motif2_merged.bed',
                       mappable_regions = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_mappable.bed',
                       chrom_lengths = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/chrom.sizes.txt',
                       slop, limit2chrom=TRUE, chrom='X', from, to, n = 100, plot = TRUE, region=NULL, feature=NULL){

  if(missing(feature) && missing(region)){
    region_name <- tools::file_path_sans_ext(basename(regionA))
    feature_name <- tools::file_path_sans_ext(basename(regionB))
  }
  genome = read.delim(chrom_lengths, header=F)
  mappable = read.delim(mappable_regions, header=F)
  
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
    cat("Genome starts at", from, "and ends ", to, "\n")
  }
  
  test <- read.delim(regionA, header =F)
  test <- test[,c(1,2,3)]
  test$V1 <- stringr::str_remove(test$V1, 'chr')
  
  test <- test %>% 
    dplyr::mutate(length = V3 - V2) %>% 
    dplyr::filter(...,
                  V1 %in% levels(genome$V1)) %>%
    dplyr::select(-length)
  
  cat("looking for enrichment of ", feature_name, "in ", nrow(test), "regions of ", region_name, "\n")
 
  
  feature <- read.delim(regionB, header = F)
  feature <- feature[,c(1,2,3)]
  
  feature$V1 <- stringr::str_remove(feature$V1, 'chr')
  
  if(!missing(slop)){
    cat("Expanding ", feature_name, " regions by ", slop, "\n")
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
  
  if(plot){
    
    d <- data.frame(feature = c('real', rep('shuffled', length(res$`regioneR::numOverlaps`$permuted))),
                    overlaps = c(res$`regioneR::numOverlaps`$observed, res$`regioneR::numOverlaps`$permuted)
    )
    
    title_string <- paste0("Region: ", region_name, "\n",
                           "Feature: ", feature_name, "\n",
                           "shuffles: ", res$`regioneR::numOverlaps`$ntimes, "\n",
                           "observed: ", res$`regioneR::numOverlaps`$observed, "\n",
                           "p-val: ", format.pval(res$`regioneR::numOverlaps`$pval) ,"\n")
    # library(ggpubr)
    
    gghistogram(d, x = "overlaps",
                add = "mean", rug = TRUE,
                color = "feature", fill = "feature",
                palette = c("#00AFBB", "#E7B800"),
                title = title_string, 
                ggtheme = theme_minimal())
    
  } else{
    return(res)
  }
  
  
  # lz <- localZScore(pt=pt, A=test, B=feature)
  # plot(lz)
}

# feature <- feature %>% 
#   dplyr::rename(chrom = V1,
#                 start = V2,
#                 end   = V3)
# 
# test <- test %>% 
#   dplyr::rename(chrom = V1,
#                 start = V2,
#                 end   = V3)
# 
# r <- data.table(feature)
# b <- data.table(test)
# setkey(r)
# 
# annotated <- as.data.frame(data.table::foverlaps(b, r, by.x = names(b), type = "any", mult = "first", nomatch = NA))
# 
# bpRegions <- annotated %>% 
#   dplyr::mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
#   dplyr::select(chrom, start, end, feature, i.start, i.end)

# # Limit to mappable genome
# pt <- permTest(A=test, ntimes=50, randomize.function=resampleRegions, genome=genome, universe = mappable,
#                evaluate.function=numOverlaps, B=feature, verbose=FALSE)

# 
# bpMeanDistance <- function(test_file = '~/Desktop/misc_bed/breakpoints/all_prec_dels',
#                            ref_file = '~/Desktop/misc_bed/damID/RPA70_damid.bed',
#                            n = 50){
#   
#   test <- read.delim(test_file, header =F)
#   ref <- read.delim(ref_file, header = F)
#   genome = read.delim("~/Documents/Curie/Data/Genomes/Dmel_v6.12/chrom.sizes.txt", header=F)
#   
#   shuff <- randomizeRegions(test, genome=genome)
#   
#   
#   real_dist <- meanDistance(A = test, B=ref)
#   sim_dist <- meanDistance(A = shuff, B=ref)
#   
#   cat("real mean distance: ", real_dist, "\n")
#   cat("sim mean distance: ", sim_dist, "\n")
#   
#   
# }
# 
# 
# 
# pt <- permTest(A=test, B=feature, randomize.function=randomizeRegions,
#                evaluate.function=numOverlaps)
