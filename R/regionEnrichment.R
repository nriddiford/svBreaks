# bpRegionEnrichment
#' Calculate the enrichment of SVs in genomic regions
#' @keywords enrichment
#' @import dplyr
#' @import data.table
#' @export
bpRegionEnrichment <- function(..., bedDir='/Users/Nick_curie/Desktop/misc_bed/features', breakpoints=NA, keep=NULL,
                               slop=0, plot=TRUE, genome_length=118274340, intersect=FALSE, outDir=NA, parseName=FALSE, minHits=10){
  if(is.na(breakpoints)){
    breakpoints <- 'svs'
    bps <- getData(...,
                   genotype=='somatic_tumour',
                   confidence == 'precise')
    bps <- bps %>% 
      dplyr::rename(start = bp) %>% 
      dplyr::mutate(end = start+1) %>%
      select(chrom, start, end)
    
  } else if(breakpoints=='notch'){
      Nstart <- 3134870 - 500000
      Nstop <- 3172221 + 500000
      N_window <- Nstop - Nstart
      if(!missing(keep)){
        genome_length <- N_window
      } else {
        genome_length <- genome_length - N_window
      }
      
      cat("X:", Nstart, "-", Nstop, sep='', "\n")
      
      N_filtered <- notchFilt(..., genotype=='somatic_tumour', keep=keep)
      
      bps <- N_filtered %>% 
        dplyr::rename(start = bp) %>% 
        dplyr::mutate(end = start+1) %>%
        select(chrom, start, end)
      
  } else {
    bps <- read.table(breakpoints, header = F)
    
    if(ncol(bps)<3){
      bps$V3 <- bps$V2 + 2
    }
    bps <- bps[,c(1,2,3)]
   
    colnames(bps) <- c("chrom", "start", "end")
  
    bps <- bps %>% 
      dplyr::mutate(end = as.integer(((end+start)/2)+1)) %>%
      dplyr::mutate(start = as.integer(end-1)) %>%
      dplyr::select(chrom, start, end)
  }
  
  cat("Specified genome size:", genome_length, "\n")
  cat("Expanding regions by", slop, "\n\n")
  
  scores <- list()
  regionFC <- list()
  
  fileNames <- dir(bedDir, pattern = ".bed")
  # cat("Analysing all files in directory:", bedFiles, "\n")
  for (f in fileNames){
    # cat("Analysing file:", f, "\n")
    
    if(parseName){
      filename <- basename(tools::file_path_sans_ext(f))
      parts <- unlist(strsplit(filename, split = '_'))
      factor <- parts[1]
      genotype <- parts[2]
      tissue <- parts[3]
      element <- parts[4]
      replicate <- parts[5]
      id <- unlist(strsplit(parts[7], split= "[.]"))[1]
    }else{
      filename <- tools::file_path_sans_ext(f)
      factor <- unlist(strsplit(basename(tools::file_path_sans_ext(f)), split = "\\."))[1]
      genotype <- 'wt'
      tissue <- 'pooled'
      element <- 'misc'
      replicate <- '0'
      id <- 'NA'
    }
    
    regions <- read.table(paste(bedDir, f, sep='/'), header = F)
    regions <- regions[,c(1,2,3)]
    colnames(regions) <- c("chrom", "start", "end")
    
    myChroms <- c("2L", "2R", "3L", "3R", "X", "Y", "4")
    
    # This slop wont handle cases where the adjusted regions overlap (they will be counted twice...)
    regions <- regions %>% 
      filter(chrom %in% myChroms) %>% 
      filter(start < end) %>% 
      mutate(start = start - slop) %>% 
      mutate(end = end + slop) %>% 
      arrange(chrom) %>% 
      droplevels()
    
    if(breakpoints=='notch'){
      if(!missing(keep)){
        regions <- regions %>% 
          filter(chrom == "X", start >= Nstart, end <= Nstop) %>% 
          droplevels()
      } else {
        regions <- regions %>% 
          filter(!(chrom == "X" & start >= Nstart & end <= Nstop) ) %>% 
          droplevels()
      }
    }
    
    if (!nrow(regions)) next

    if(intersect){
      cat("Analysing file:", f, "\n")
      merged_regions <- mergeOverlaps(regions, dataframe = T)
      cat("Merged file:", f, "\n")
      
      intersected_regions <- subtractUnmappable(merged_regions, dataframe=T)
      cat("Intersected file:", f, "\n")
      intersected_regions$chr <- as.factor(intersected_regions$chr)
      regions <- intersected_regions %>% 
        dplyr::rename(chrom=chr)
    }
    
    # if(write){
    #   cat("Writing bed file:", f, "\n")
    #   basename <- tools::file_path_sans_ext(f)
    #   basename <- paste(basename, '_merged_intersected_', slop, '.bed', sep='')
    #   writeBed(df=regions, name = basename)
    # }
   
    # setkey = d1 for breakpoints, d2 for regions
    r <- data.table(regions)
    b <- data.table(bps)
    setkey(r)
    
    # search for bp overlaps with regions (d2, d1) for breakpoints, d1, d2 for regions
    annotated <- as.data.frame(foverlaps(b, r, by.x = names(b), type = "any", mult = "first", nomatch = NA))
    
    bpRegions <- annotated %>% 
      dplyr::mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
      dplyr::select(chrom, start, end, feature, i.start, i.end)
    
    # Replace old write option to now show regions in file that contain breakpoints
    if(!is.na(outDir)){
      cat("Writing bed file of overlaps for :", factor, "\n")
      basename <- factor
      basename <- paste(basename, '_containing_', slop, '.bed', sep='')
      writeBed(df=bpRegions, outDir = outDir, name = basename)
    }
    
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
      
      list(feature = factor, observed = inRegion, expected = expectedHits, Log2FC = Log2FC, test = test, sig = sig_val, p_val = stat$p.value, genotype=genotype, tissue=tissue, element=element, replicate=replicate, id=id, filename=filename)
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
  

  minPval <- min(final$p_val[final$p_val>0])
  
  # final$p_val <- ifelse(final$p_val==0, minPval/abs(final$Log2FC), final$p_val)
  # 
  # final$p_val <- ifelse(final$p_val==0, minPval, final$p_val)

  final <- final %>%
    dplyr::mutate(p_val = ifelse(p_val==0, minPval/abs(Log2FC), p_val)) %>% 
    dplyr::mutate(p_val = ifelse(p_val==0, minPval, p_val)) %>%
    dplyr::mutate(padj = p.adjust(p_val, method = 'hochberg')) %>%
    dplyr::mutate(eScore = round(abs(Log2FC) * -log10(padj),2)) %>% 
    dplyr::filter(count >= minHits) %>%
    dplyr::select(feature:p_val, padj, eScore, genotype, tissue, element, replicate, id, filename, -count) %>% 
    dplyr::arrange(-eScore, padj, -abs(Log2FC)) %>% 
    droplevels()
  
  if(nrow(final)==0){
    cat("Fewer than", minHits, "found in all files. Try adjusting 'minHits' to a lower number", sep=" ", "\n")
  } else { 
    if(plot){
      cat("Plotting volcano plot", "\n")
      print(Volcano(final))
      # print(ggVolcano(df=final))
    }
    return(final)
  }
}


# Volcano
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import dplyr
#' @import plotly
#' @export
Volcano <- function(d, byq=FALSE){
  feature_enrichment <- d
  
  var <- feature_enrichment$p_val
  printVar <- 'P-val'
  if(byq){
    var <- var <- feature_enrichment$padj
    printVar <- 'Adj-p'
  }
  
  # minPval <- min(feature_enrichment$p_val[feature_enrichment$p_val>0])
  # 
  # feature_enrichment$p_val <- ifelse(feature_enrichment$p_val==0, minPval/abs(feature_enrichment$Log2FC), feature_enrichment$p_val)
  # 
  # feature_enrichment$p_val <- ifelse(feature_enrichment$p_val==0, minPval, feature_enrichment$p_val)
  
  
  maxLog2 <- max(abs(feature_enrichment$Log2FC[is.finite(feature_enrichment$Log2FC)]))
  maxLog2 <- as.numeric(round_any(maxLog2, 1, ceiling))
  
  ax <- list(
    size = 25
  )
  
  ti <- list(
    size = 25
  )
  
  p <- plot_ly(data = feature_enrichment,
          x = ~Log2FC,
          y = ~-log10(var),
          type = 'scatter',
          # showlegend = FALSE,
          mode = 'markers',
          # height = 1200,
          # width = 1000,
          # frame = ~p_val,
          text = ~paste("Feature: ", feature, "\n",
                        "Observed: ", observed, "\n",
                        "Expected: ", expected, "\n",
                        "P-val: ", p_val, "\n",
                        "Adj. P-val: ", padj, "\n",
                        "eScore: ", eScore, "\n"),
          color = ~log10(var),
          colors = "Spectral",
          size = ~-log10(var)
          ) %>% 
    layout(
           xaxis = list(title="Log2(FC)", titlefont = ax, range = c(-maxLog2, maxLog2)),
           yaxis = list(title=paste("-Log10(", printVar, ")", sep=''), titlefont = ax)
           )
  p
}

# ggVolcano
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import plyr dplyr ggplot2
#' @import RColorBrewer
#' @export

ggVolcano <- function(..., df){
  feature_enrichment <- df
  
  arguments <- list(...)
  
  if(!is.null(arguments$slop)){
    slop <- arguments$slop
    title <- paste("Enrichment/depletion of genomic features\nfor breakpoints +/-", slop, 'bps')
    outFile <- paste("regionEnrichment_", slop,'.png', sep='')
  } else {
    title <- "Enrichment/depletion of genomic features\nfor breakpoints"
    outFile <- "regionEnrichment.png"
  }
  
  maxLog2 <- max(abs(feature_enrichment$Log2FC))
  maxLog2 <- round_any(maxLog2, 1, ceiling)

  p <- ggplot(feature_enrichment)
  p <- p + geom_point(aes(Log2FC, -log10(p_val), size = -log10(p_val), colour = sig))
  # p <- p + scale_color_hue(l=60, c=50)
  p <- p + scale_color_brewer(palette="Spectral")
  # p <- p + scale_color_gradientn(colours = rainbow(100))
  p <- p + guides(size = FALSE, colour = FALSE) 
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      legend.text=element_text(size=rel(1.5))
      # legend.key.size = unit(3,"line"),
      # legend.title=element_text(size=rel(2))
    )
  p <- p + scale_x_continuous(limits=c(-maxLog2, maxLog2))
  p <- p + ggtitle(title)
  p
  
}


# plotMonteCarlo
#'
#' Plot the result of a monte carlo simulation (n shuffles of feature y)
#' Expects dataframe returned from bpRegionEnrichment()
#' @keywords expected frequency
#' @import plyr dplyr ggplot2
#' @import RColorBrewer
#' @export

plotMonteCarlo <- function(x, bindWith = 10){
  
  x <- x %>% 
    mutate(fillCol = ifelse(grepl("shuff", feature), "grey37", "#C72424FE"))

  ggplot(x, aes(observed, fill = fillCol)) +
  geom_histogram(binwidth = bindWith) +
  cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 30)
      ) + 
  scale_y_continuous("Frequency", expand = c(0,0)) +
  scale_x_continuous("Number of intersections", expand = c(0,0)) +
  
  scale_fill_identity()

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
  bedFile$start <- ifelse(bedFile$start < 0, 0, bedFile$start)
  
  bedFile.sort <- bedr.sort.region(bedFile, check.chr = FALSE, method = 'lexicographical', check.valid = FALSE)
  
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
#' @import dplyr
#' @export

writeBed <- function(df, outDir=getwd(), name='regions.bed', svBreaks=FALSE){
  if(svBreaks){
    df <- df %>%
      dplyr::filter(bp_no=='bp1') %>% 
      dplyr::mutate(info = paste(sample, type, feature, feature2, sep = "_")) %>%
      dplyr::rename(start = bp) %>% 
      dplyr::rename(end = bp2) %>% 
      dplyr::select(chrom, start, end, info)
  } else{
    colnames(df[,c(1,2,3)]) <- c("chrom", "start", "end")
    names(df)[1:3] <- c("chrom", "start", "end")
  }
  
  df <- df %>% 
    dplyr::filter(as.numeric(start) < as.numeric(end)) %>% 
    droplevels()
  
  cat(paste(outDir,name, sep='/'), "\n")
  
  write.table(df, file = paste(outDir,name, sep='/'), row.names=F, col.names=F, sep="\t", quote = FALSE)
}


# bpRegionEnrichmentPlot
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import plyr
#' @import dplyr ggplot2
#' @export

bpRegionEnrichmentPlot <- function(...) {
  feature_enrichment <- bpRegionEnrichment(..., plot=F)
  arguments <- list(...)
  
  if(!is.null(arguments$slop)){
    slop <- arguments$slop
    title <- paste("Enrichment/depletion of genomic features\nfor breakpoints +/-", slop, 'bps')
    outFile <- paste("regionEnrichment_", slop,'.png', sep='')
  } else {
    title <- "Enrichment/depletion of genomic features\nfor breakpoints"
    outFile <- "regionEnrichment.png"
  }
  
  
  feature_enrichment <- feature_enrichment %>% 
    dplyr::mutate(feature = as.character(feature)) %>% 
    dplyr::mutate(sig = ifelse(sig=='-', '',sig)) %>% 
    transform(feature = reorder(feature, -Log2FC))
  
  
  maxLog2 <- max(abs(feature_enrichment$Log2FC))
  maxLog2 <- round_any(maxLog2, 1, ceiling)
  
  adjVal <- (maxLog2*1.5)
  
  print(feature_enrichment)

  p <- ggplot(feature_enrichment)
  p <- p + geom_bar(aes(feature, Log2FC, fill = as.character(test)), stat = "identity")
  p <- p + guides(fill = FALSE)
  p <- p + ylim(-maxLog2, maxLog2)
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  p <- p + geom_text(aes(feature, maxLog2, label=sig, vjust=-adjVal),size=5)
  p <- p + ggtitle(title)
 
  
  cat("Writing file", outFile, "\n")
  ggsave(paste("plots/", outFile, sep = ""), width = nrow(feature_enrichment)/2, height = 10)
  p
}

# delsRegionEnrichment
#'
#' asdas
#' @keywords regions overlap
#' @import data.table
#' @export
delsRegionEnrichment <- function(bedFiles='bed', region=TRUE, breakpoints=NA,  genome_length=118274340 ){
  if(is.na(breakpoints)){
    breakpoints <- getData(genotype=='somatic_tumour', type=='DEL')
    bps <- breakpoints %>% 
      dplyr::filter(bp_no == 'bp1') %>% 
      dplyr::rename(start = bp) %>% 
      dplyr::rename(end = bp2) %>%
      dplyr::mutate(width = (end - start)) %>%
      dplyr::mutate(end = ifelse(width <= 0, start+1, end)) %>% 
      # mutate(width2 = (end - start)) %>%
      dplyr::select(chrom, start, end)
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
      dplyr::filter(chrom %in% myChroms) %>% 
      droplevels()
    
    # setkey = r for breakpoints, b for regions
    r <- data.table(regions)
    b <- data.table(bps)
    setkey(r, chrom, start, end)
    
    # search for bp overlaps with regions (b, r) for breakpoints, r, b for regions
    annotated <- as.data.frame(foverlaps(b, r, by.x = names(b), type = "any", mult = "all", nomatch = NA))
    
    bpRegions <- annotated %>% 
      dplyr::mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
      dplyr::select(chrom, start, end, feature, i.start, i.end)
    
    # Fraction of genome that is region
    regionSpace <- regions %>% 
      dplyr::group_by(chrom) %>% 
      dplyr::summarise(chromSpace = sum(end-start))
    
    regionWidth <- sum(regionSpace$chromSpace)
    regionFraction <- regionWidth / genome_length
    
    # Fraction of genome that is mutated
    mutSpace <- bps %>% 
      dplyr::group_by(chrom) %>% 
      dplyr::summarise(chromSpace = sum(end-start))
    
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
    dplyr::filter(expected >= 5) %>% 
    droplevels()
  
  # print(final)
  
  # return(regionsTested)
}
