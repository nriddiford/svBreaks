bpRegioneR <- function(..., regionA = '~/Desktop/misc_bed/breakpoints/Notch_CFS/notch_10kb_merged_breakpoints.bed',
                       regionB = '~/GitHub/BardinLab/meme/out/notch_CFS/fimo_out/motif2_merged.bed',
                       mappable_regions = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_mappable.bed',
                       chrom_lengths = '~/Documents/Curie/Data/Genomes/Dmel_v6.12/chrom.sizes.txt',
                       slop, limit2chrom=TRUE, from, to, n = 100){

  
  genome = read.delim(chrom_lengths, header=F)
  mappable = read.delim(mappable_regions, header=F)
  
  mappable <- mappable %>%
    dplyr::select(-V4) %>%
    dplyr::filter(V1 %in% levels(genome$V1)) %>% 
    droplevels()
  
  if(limit2chrom && (!missing(from) || !missing(to)) ){
    mappable <- mappable %>%
      dplyr::filter(V2 >= (from - 5000),
                    V3 <= (to + 5000)
      ) %>%
      droplevels()
    # genome$V2 = from
    # genome$V3 = to
    cat("Genome starts at", from, "and ends ", to, "\n")
  }
  
  test <- read.delim(regionA, header =F)
  test <- test[,c(1,2,3)]
  
  test <- test %>% 
    dplyr::mutate(length = V3 - V2) %>% 
    dplyr::filter(...,
                  V1 %in% levels(genome$V1)) %>%
    dplyr::select(-length)
  
  cat("looking for enrichment in", nrow(test), "regions\n")
 
  if(!missing(slop)){
    cat("Expanding test regions by ", slop, "\n")
    test <- test %>%
      dplyr::mutate(V2 = V2 - slop,
                    V3 = V3 + slop)
  }
  
  feature <- read.delim(regionB, header = F)
  feature <- feature[,c(1,2,3)]
  
  # mappable_genome <- getGenomeAndMask("dm6", mask=exclude)$genome
  # mappable_genome <- getGenomeAndMask(genome=genome, mask=exclude)
  # mappable_genome <- regioneR::overlapRegions(A = genome, B = mappable, type="BinA") %>% 
  #   dplyr::select(chr, startB, endB)
  
  # pt <- regioneR::overlapPermTest(A=test, B=feature, ntimes=n, genome=mappable_genome,  per.chromosome=TRUE)
  # plot(pt)
  pt <- regioneR::permTest(A=test, B=feature, ntimes=n, randomize.function=regioneR::randomizeRegions,
                           genome=mappable, evaluate.function=regioneR::numOverlaps, per.chromosome=limit2chrom)
  plot(pt)
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
