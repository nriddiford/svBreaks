# bpRegionEnrichment
#'
#' Calculate the enrichment of SVs in genomic regions
#' @keywords enrichment
#' @import tidyverse
#' @import data.table
#' @export

bpRegionEnrichment <- function(bedDir='inst/extdata/bed/mappable', breakpoints=NA, slop=0, plot=TRUE, genome_length=118274340, intersect=FALSE, write=FALSE ){
  if(is.na(breakpoints)){
    breakpoints <- getData(genotype=='somatic_tumour', !sample %in% c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"))
    bps <- breakpoints %>% 
      dplyr::rename(start = bp) %>% 
      dplyr::mutate(end = start+1) %>%
      select(chrom, start, end)
    
  } else{
    bps <- read.delim(breakpoints, header = F)
    colnames(bps) <- c("chrom", "start", "end")
  }
  
  if(is.null(bps$end)){
    bps$end <- bps$start + 1
  }
  
  cat("Expanding regions by", slop, "\n\n")
  
  scores <- list()
  regionFC <- list()
  
  fileNames <- dir(bedDir, pattern = ".bed")
  # cat("Analysing all files in directory:", bedFiles, "\n")
  for (f in fileNames){
    # cat("Analysing file:", f, "\n")
    
    regions <- read.delim(paste(bedDir, f, sep='/'), header = F)
    regions <- regions[,c(1,2,3)]
    colnames(regions) <- c("chrom", "start", "end")
    
    myChroms <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
    
    # This slop wont handle cases where the adjusted regions overlap (they will be counted twice...)
    regions <- regions %>% 
      filter(chrom %in% myChroms) %>% 
      mutate(start = start - slop) %>% 
      mutate(end = end + slop) %>% 
      arrange(chrom) %>% 
      droplevels()
    
    if(intersect){
      cat("Analysing file:", f, "\n")
      merged_regions <- mergeOverlaps(regions, dataframe = T)
      cat("Merged file:", f, "\n")
      
      intersected_regions <- subtractUnmappable(merged_regions, dataframe=T)
      cat("Intersected file:", f, "\n")
      intersected_regions$chr <- as.factor(intersected_regions$chr)
      regions <- intersected_regions %>% 
        dplyr::rename(chrom=chr)
      
      if(write){
        cat("Writing bed file:", f, "\n")
        basename <- tools::file_path_sans_ext(f)
        basename <- paste(basename, '_merged_intersected_', slop, '.bed', sep='')
        writeBed(df=regions, name = basename)
      }

    }
   
    # setkey = d1 for breakpoints, d2 for regions
    r <- data.table(regions)
    b <- data.table(bps)
    setkey(r)
    
    # search for bp overlaps with regions (d2, d1) for breakpoints, d1, d2 for regions
    annotated <- as.data.frame(foverlaps(b, r, by.x = names(b), type = "any", mult = "first", nomatch = NA))
    
    bpRegions <- annotated %>% 
      mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
      select(chrom, start, end, feature, i.start, i.end)
    
    regionSpace <- regions %>% 
      group_by(chrom) %>% 
      summarise(chromSpace = sum(end-start)) 
    
    regionWidth <- sum(regionSpace$chromSpace)
    
    mutCount <- nrow(bps)
    
    regionFraction <- regionWidth / genome_length
    
    if(regionFraction>=1){
      cat("Adding", slop, "to", f, "causes the regions to occupy a larger space than the genome. Setting to 1", "\n")
      regionFraction = 0.999
    }
    
    expectedHits <- (mutCount * regionFraction)
    
    inRegion <- sum(ifelse(bpRegions$feature == 'inside', 1, 0))
    outRegion <- sum(ifelse(bpRegions$feature == 'outside', 1, 0))
    
    expectedHits <- round(expectedHits, digits = 2)
    expectedHits <- ifelse(expectedHits==0, 0.001, expectedHits)
    
    fc <- inRegion / expectedHits
    Log2FC <- log2(fc)
    
    fc <- round(fc, digits = 1)

    cat(regionFraction, f, "\n")
    
    biTest <- function(f){
      # Binomial test
      if (inRegion >= expectedHits) {
        stat <- binom.test(x = inRegion, n = mutCount, p = regionFraction, alternative = "greater")
        test <- "enrichment"
      } else {
        stat <- binom.test(x = inRegion, n = mutCount, p = regionFraction, alternative = "less")
        test <- "depletion"
      }
      
      sig_val <- ifelse(stat$p.value <= 0.001, "***",
                        ifelse(stat$p.value <= 0.01, "**",
                               ifelse(stat$p.value <= 0.05, "*", "")))
      
      sig_val <- ifelse(stat$p.value > 0.05, "-", sig_val)
      
      p_val <- format.pval(stat$p.value, digits = 3, eps = 0.0001)
      filename <- tools::file_path_sans_ext(f)
      list(feature = filename, observed = inRegion, expected = expectedHits, Log2FC = Log2FC, test = test, sig = sig_val, p_val = stat$p.value)
    }
    
    biNomialTest <- lapply(f, biTest)
    regionsTested <- do.call(rbind, biNomialTest)
    regionsTested <- as.data.frame(regionsTested)
    
    # regionsTested$Log2FC <- round(as.numeric(regionsTested$Log2FC), 1)
    regionsTested$Log2FC <- as.numeric(regionsTested$Log2FC)
    
    regionsTested$feature <- as.character(regionsTested$feature)
    regionsTested$observed <- as.numeric(regionsTested$observed)
    regionsTested$test <- as.character(regionsTested$test)
    regionsTested$sig <- as.character(regionsTested$sig)
    regionsTested$p_val <- as.numeric(regionsTested$p_val)
    regionsTested$expected <- round(as.numeric(regionsTested$expected), 2)
    scores[[f]] <- regionsTested
    
  }
  
  #combine each iteration into one data frame
  final <- as.data.frame(do.call(rbind, scores))
  
  final$count <- as.numeric(final$observed) + as.numeric(final$expected)
  final <- final %>%
    filter(count >= 5) %>%
    select(-count) %>% 
    arrange(p_val, -abs(Log2FC)) %>% 
    droplevels()
  
  if(plot){
    cat("Plotting volcano plot", "\n")
    print(Volcano(final))
  }
  # print(final)
  return(final)

}


# Volcano
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import dplyr
#' @import plotly
#' @export

Volcano <- function(df){
  feature_enrichment <- df

  # feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -Log2FC))
  ax <- list(
    size = 25
  )
  
  ti <- list(
    size = 25
  )
  
  plot_ly(data = feature_enrichment,
          x = ~Log2FC,
          y = ~-log10(p_val),
          type = 'scatter',
          mode = 'markers',
          # height = 1200,
          # width = 1000,
          # frame = ~p_val,
          text = ~paste("Feature: ", feature, "\n",
                        "Observed: ", observed, "\n",
                        "Expected: ", expected, "\n",
                        "P-val: ", p_val, "\n"),
          color = ~log10(p_val), size = ~-log10(p_val) ) %>% 
    layout(title = "Depletion/Enrichment of breakpoints in genomic features", titlefont = ti,
           xaxis = list(title="Log2(FC)", titlefont = ax),
           yaxis = list(title="-Log10(p)", titlefont = ax),
           showlegend = FALSE)
}


# mergeOverlaps
#'
#' Implements bedtools merge to merge overlapping regions in bedfile
#' @keywords bed merge
#' @import bedr
#' @export
#' 
mergeOverlaps <- function(f, dataframe=FALSE){
  
  if(!dataframe){
    bedFile <- read.delim(f, header = F)
  } else{
    bedFile <- f
  }
  bedFile <- bedFile[,c(1,2,3)]
  colnames(bedFile) <- c("chr", "start", "end")
  bedFile$chr <- as.character(bedFile$chr)
  bedFile.sort <- bedr.sort.region(bedFile, check.chr = FALSE)
  bedFile.merged <- bedr(engine = "bedtools", input = list(i = bedFile.sort), method = "merge", params = "",  check.chr = FALSE)
  
  return(data.frame(bedFile.merged, row.names = NULL))
}


# subtractUnmappable
#'
#' Implements bedtools subtract to remove unmappable regions in bedfiles
#' @keywords bed subtract
#' @import bedr
#' @export

subtractUnmappable <- function(f, u='~/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_75.bed', dataframe=FALSE){
  if(!dataframe){
    bedFile <- read.delim(f, header = F)
  } else{
    bedFile <- f
  }
  bedFile <- bedFile[,c(1,2,3)]
  colnames(bedFile) <- c("chr", "start", "end")
  bedFile$chr <- as.character(bedFile$chr)
  bedFile.sort <- bedr.sort.region(bedFile, check.chr = FALSE, method = 'lexicographical')
  
  excludeFile <- read.delim(u, header = F)
  excludeFile <- excludeFile[,c(1,2,3)]
  colnames(excludeFile) <- c("chr", "start", "end")
  excludeFile$chr <- as.character(excludeFile$chr)
  excludeFile.sort <- bedr.sort.region(excludeFile, check.chr = FALSE,  method = 'lexicographical')
  
  bedFileregionsExcluded <- bedr.subtract.region(bedFile.sort, excludeFile.sort, remove.whole.feature = FALSE, check.chr = FALSE)
  return(data.frame(bedFileregionsExcluded, row.names = NULL))
}


# writeBed
#'
#' Simple function to write out bedfiels - useful for checking after merging/subtracting...
#' @keywords bed
#' @export

writeBed <- function(df, outDir=getwd(), name='regions.bed'){
  # svs <- svs %>% 
  #   mutate(info = paste(sample, type, feature, feature2, sep = "_")) %>% 
  #   select(chrom, bp, bp2, info)
  cat(paste(outDir,name, sep='/'))
  
  write.table(df, file = paste(outDir,name, sep='/'), row.names=F, col.names=F, sep=" ", quote = FALSE)
}























# bpRegionEnrichmentPlot
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import tidyverse
#' @export

bpRegionEnrichmentPlot <- function(...) {
  feature_enrichment <- bpRegionEnrichment()
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -Log2FC))
  
  feature_enrichment <- droplevels(feature_enrichment)
  
  
  
  p <- ggplot(feature_enrichment)
  p <- p + geom_bar(aes(feature, Log2FC, fill = as.character(test)), stat = "identity")
  p <- p + guides(fill = FALSE)
  p <- p + ylim(-2, 2)
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  
  feat_plot <- paste("feat_plot.png")
  cat("Writing file", feat_plot, "\n")
  ggsave(paste("plots/", feat_plot, sep = ""), width = 5, height = 10)
  p
}



delsRegionEnrichment <- function(bedFiles='bed', region=TRUE, breakpoints=NA,  genome_length=118274340 ){
  require(data.table) ## 1.9.3
  if(is.na(breakpoints)){
    breakpoints <- getData(genotype=='somatic_tumour', type=='DEL')
    bps <- breakpoints %>% 
      filter(bp_no == 'bp1') %>% 
      dplyr::rename(start = bp) %>% 
      dplyr::rename(end = bp2) %>%
      # dplyr::rename(chr = bp2) %>%
      mutate(width = (end - start)) %>%
      # mutate(end = ifelse(width < 0, start+1, end)) %>% 
      # mutate(width2 = (end - start)) %>%
      select(chrom, start, end)
  } else{
    bps <- read.delim(breakpoints, header = F)
    colnames(bps) <- c("chrom", "start", "end")
  }
  
  scores <- list()
  regionFC <- list()
  
  fileNames <- dir(bedFiles, pattern = ".bed")
  for (f in fileNames){
    cat(f, "\n")
    regions <- read.delim(paste(bedFiles,f, sep='/'), header = F)
    regions <- regions[,c(1,2,3)]
    colnames(regions) <- c("chrom", "start", "end")
    
    myChroms <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
    
    regions <- regions %>% 
      filter(chrom %in% myChroms) %>% 
      droplevels()
    
    # setkey = r for breakpoints, b for regions
    r <- data.table(regions)
    b <- data.table(bps)
    setkey(r, chrom, start, end)
    
    # search for bp overlaps with regions (b, r) for breakpoints, r, b for regions
    annotated <- as.data.frame(foverlaps(b, r, by.x = names(b), type = "any", mult = "all", nomatch = NA))
    
    bpRegions <- annotated %>% 
      mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
      select(chrom, start, end, feature, i.start, i.end)
    
    # Fraction of genome that is region
    regionSpace <- regions %>% 
      group_by(chrom) %>% 
      summarise(chromSpace = sum(end-start))
    
    regionWidth <- sum(regionSpace$chromSpace)
    regionFraction <- regionWidth / genome_length
    
    # Fraction of genome that is mutated
    mutSpace <- bps %>% 
      group_by(chrom) %>% 
      summarise(chromSpace = sum(end-start))
    
    mutWidth <- sum(mutSpace$chromSpace)
    mutFraction <- mutWidth / genome_length
    
    mutCount <- nrow(bps)
    
    expectedHits <- (mutCount * regionFraction)
    
    inRegion <- sum(ifelse(bpRegions$feature == 'inside', 1, 0))
    outRegion <- sum(ifelse(bpRegions$feature == 'outside', 1, 0))
    
    fc <- inRegion / expectedHits
    fc <- round(fc, digits = 1)
    expectedHits <- round(expectedHits, digits = 1)
    
    
    # f = 'uce'
    biTest <- function(f){
      # Binomial test
      if (inRegion >= expectedHits) {
        stat <- binom.test(x = inRegion, n = mutCount, p = regionFraction, alternative = "greater")
        test <- "enrichment"
      } else {
        stat <- binom.test(x = inRegion, n = mutCount, p = regionFraction, alternative = "less")
        test <- "depletion"
      }
      
      sig_val <- ifelse(stat$p.value <= 0.001, "***",
                        ifelse(stat$p.value <= 0.01, "**",
                               ifelse(stat$p.value <= 0.05, "*", "")))
      
      sig_val <- ifelse(stat$p.value > 0.05, "-", sig_val)
      
      p_val <- format.pval(stat$p.value, digits = 3, eps = 0.0001)
      Log2FC <- log2(fc)
      list(feature = f, observed = inRegion, expected = expectedHits, Log2FC = Log2FC, test = test, sig = sig_val, p_val = p_val)
    }
    
    biNomialTest <- lapply(f, biTest)
    regionsTested <- do.call(rbind, biNomialTest)
    regionsTested <- as.data.frame(regionsTested)
    
    regionsTested$Log2FC <- round(as.numeric(regionsTested$Log2FC), 1)
    regionsTested$expected <- round(as.numeric(regionsTested$expected), 1)
    scores[[f]] <- regionsTested
    
  }
  
  #combine each iteration into one data frame
  # final <- dplyr::bind_rows(byIteration)
  final <- as.data.frame(do.call(rbind, scores))
  final <- final %>% 
    filter(expected >= 5) %>% 
    droplevels()
  
  # print(final)
  
  # return(regionsTested)
  
  
}


#   install_github("favorov/GenometriCorr")

# genoCorr <- function(){
#   library('GenometriCorr')
#   library("rtracklayer")
#   
#   dmel.chrom.length <- c(23513712, 25286936, 28110227, 32079331, 23542271, 3667352, 1348131)
#                          
#   names(dmel.chrom.length) <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
#  
#   my_chr <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
#   
#   pn.area <- 100
#   pn.dist <- 100
#   pn.jacc <- 100
#   a <- as(import(bedFile),"RangedData")
# 
#   bps <- breakpoints %>% 
#     filter(bp_no == 'bp1') %>% 
#     filter(chrom == chrom2) %>% 
#     mutate(width = (bp2 - bp)) %>% 
#     select(chrom, bp, bp2, width)
#   
#   b <- GRanges(seqnames = Rle(bps$chrom), IRanges(bps$bp,
#                                                   bps$bp2))
#   a_vs_b <- GenometriCorrelation(a, b, chromosomes.length = dmel.chrom.length,
#                                  chromosomes.to.proceed = my_chr,
#                                  ecdf.area.permut.number = pn.area,
#                                  mean.distance.permut.number = pn.dist,
#                                  jaccard.measure.permut.number = pn.jacc,
#                                  keep.distributions = TRUE, showProgressBar = TRUE)
#   
#   
#   
#   
#   
#   
#   
# }
#   
