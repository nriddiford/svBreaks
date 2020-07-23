# bpRegionEnrichment
#' Calculate the enrichment of SVs in genomic regions
#' @keywords enrichment
#' @import dplyr
#' @importFrom data.table as.data.table
#' @export
bpRegionEnrichment <- function(..., bp_data, dataType='svBreaks', bed_file, bedDir='/Users/Nick_curie/Desktop/misc_bed/features', chroms=c('2L', '2R', '3L', '3R', '4', 'X', 'Y'), restrict=TRUE,
                               slop=0, plot=TRUE, genome_length=118274340, intersect=FALSE, outDir, parseName=FALSE, minHits=10, multiple=TRUE){
  
  if(missing(bp_data) && missing(bed_file)) stop("\n[!] Must provide either a df or bed file! Exiting.")
  
  if(!missing(bp_data) & dataType=='svBreaks'){
    cat("Reading from dataframe\n")
    # breakpoints <- 'svs'
    bps <- bp_data %>% 
      dplyr::filter(..., 
                    confidence == 'precise') %>% 
      dplyr::rename(start = bp) %>% 
      dplyr::mutate(end = start+1) %>% 
                    # end = as.integer(((end+start)/2)+1),
                    # start = as.integer(end-1)) %>%
      dplyr::select(chrom, start, end)
  } else if(dataType=='mutationProfiles'){
    cat("Reading SNVs from dataframe\n")
    bps <- bp_data %>% 
      dplyr::filter(...) %>%
      dplyr::mutate(start = pos,
                    end = pos+1) %>%
      dplyr::select(chrom, start, end)
  } else if(!missing(bed_file)){
    cat("Reading breakpoints from: ", bed_file, "\n")
    bps <- read.table(bed_file, header = F)
    
    if(ncol(bps)<3){
      bps$V3 <- bps$V2 + 2
    }
    bps <- bps[,c(1,2,3)]
   
    colnames(bps) <- c("chrom", "start", "end")
    bps <- bps %>% 
    dplyr::mutate(length = end - start) %>% 
      dplyr::filter(...,
                    chrom %in% chroms) %>%
      dplyr::select(-length) %>% 
      droplevels()
    
    if(median(bps$end-bps$start) > 20){
      cat(paste0("The median distance between breakpoints is large: [", median(bps$end-bps$start), " bps].", "This file is not suitable for this analysis\n"))
    } 
  }
  
  cat("Specified genome size:", genome_length, "\n")
  cat("Expanding regions by", slop, "bps\n\n")
  
  scores <- list()
  regionFC <- list()
  
  fileNames <- dir(bedDir, pattern = ".bed")
  
  if(restrict && length(chroms) == 1){
    genome_start <- min(bps$start[bps$chrom %in% chroms]) - 5000
    genome_start<- ifelse(genome_start < 0, 0, genome_start)
    genome_end <- max(bps$end[bps$chrom %in% chroms]) + 5000
    cat("Restricting regions to chroms:", chroms, "\n")
    cat("Restricting regions to locus: ", genome_start, "-", genome_end, "\n", sep="")
  }
  
  i <- 0
  for (f in fileNames){
    if(parseName){
      filename <- basename(tools::file_path_sans_ext(f))
      parts <- unlist(strsplit(filename, split = '_'))
      factor <- parts[1]
      genotype <- parts[2]
      tissue <- parts[3]
      element <- parts[4]
      replicate <- parts[5]
      id <- unlist(strsplit(parts[7], split= "[.]"))[1]
    } else{
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
    
    # This slop wont handle cases where the adjusted regions overlap (they will be counted twice...)
    regions <- regions %>% 
      dplyr::filter(chrom %in% chroms,
                    start < end) %>% 
      dplyr::mutate(start = start - slop,
                    end = end + slop) %>% 
      dplyr::arrange(chrom) %>% 
      droplevels()
    
    if(restrict && length(chroms) == 1){
      regions <- regions %>% 
        dplyr::filter(start >= genome_start,
                      end <= genome_end) %>% 
        dplyr::arrange(chrom) %>% 
        droplevels()
      
      total_length <- bps %>% 
        dplyr::group_by(chrom) %>% 
        dplyr::summarise(total = max(end) - min(start)) %>% 
        dplyr::summarise(sum(total))
      
      # TODO This doesn't account for unmappable (i.e. the functional genome can be larger than the mappable genome)
      
      genome_length <- as.numeric(genome_end - genome_start)
      cat("Functional genome length: ", genome_length, "\n")
    }
    
    if (!nrow(regions)) next

    if(intersect){
      cat("Analysing file:", f, "\n")
      merged_regions <- svBreaks::mergeOverlaps(regions, dataframe = T)
      cat("Merged file:", f, "\n")
      
      intersected_regions <- svBreaks::subtractUnmappable(merged_regions, dataframe=T)
      cat("Intersected file:", f, "\n")
      intersected_regions$chr <- as.factor(intersected_regions$chr)
      regions <- intersected_regions %>% 
        dplyr::rename(chrom=chr)
    }
    
    # setkey = d1 for breakpoints, d2 for regions
    if(missing(bed_file)){
      r <- data.table(regions)
      b <- data.table(bps)
    } else {
      r <- data.table(bps)
      b <- data.table(regions)
    }
  
    setkey(r)
    
    overlap_parm <- ifelse(multiple, "all", "first")
    
    # search for bp overlaps with regions (d2, d1) for breakpoints, d1, d2 for regions
    annotated <- as.data.frame(data.table::foverlaps(b, r, by.x = names(b), type = "any", mult = overlap_parm, nomatch = NA))
    
    bpRegions <- annotated %>% 
      dplyr::mutate(feature = ifelse(is.na(start), 'outside', 'inside')) %>% 
      dplyr::select(chrom, start, end, feature, i.start, i.end)
    
    # Replace old write option to now show regions in file that contain breakpoints
    if(!missing(outDir)){
      cat("Writing bed file of overlaps for :", factor, "\n")
      if(slop == 0) slop = '' else slop = paste0("_", slop)
      basename <- factor
      basename <- paste0(basename, '_containing', slop, '.bed')
      svBreaks::writeBed(df=bpRegions, outDir = outDir, name = basename)
    }
    
    regionSpace <- regions %>% 
      dplyr::group_by(chrom) %>% 
      dplyr::summarise(chromSpace = sum(end-start)) 
    
    regionWidth <- sum(regionSpace$chromSpace)
    mutCount <- nrow(bps)
    if(!mutCount) next
    regionFraction <- regionWidth / genome_length
    
    if(regionFraction>=1){
      cat("Adding", slop, "to", f, "causes the regions to occupy a larger space than the genome. Setting to 1", "\n")
      regionFraction = 0.999
    }
    
    # Shouldn't this be multiplying breakpoint count by 2? ** 12.2.19 ** -> NO
    expectedHits <- (mutCount * regionFraction)
    
    inRegion <- sum(ifelse(bpRegions$feature == 'inside', 1, 0))
    outRegion <- sum(ifelse(bpRegions$feature == 'outside', 1, 0))
    
    expectedHits <- round(expectedHits, digits = 2)
    expectedHits <- ifelse(expectedHits==0, 0.001, expectedHits)

    fc <- inRegion / expectedHits
    Log2FC <- log2(fc)
    
    fc <- round(fc, digits = 1)

    # cat(regionFraction, f, "\n")
    
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
    
    df <- data.frame(matrix(unlist(regionsTested), nrow=nrow(regionsTested), byrow=T), stringsAsFactors = F)
    colnames(df) <- colnames(regionsTested)
    
    
    regionsTested <- df %>% 
      dplyr::mutate(Log2FC = round(as.numeric(Log2FC), 2),
                    observed = as.numeric(observed),
                    p_val = as.numeric(p_val),
                    expected = round(as.numeric(expected), 2))
  
    # regionsTested$Log2FC <- round(as.numeric(regionsTested$Log2FC), 1)
    # regionsTested$Log2FC <- as.numeric(regionsTested$Log2FC)
    
    # regionsTested$feature <- as.character(regionsTested$feature)
    # regionsTested$observed <- as.numeric(regionsTested$observed)
    # regionsTested$test <- as.character(regionsTested$test)
    # regionsTested$sig <- as.character(regionsTested$sig)
    # regionsTested$p_val <- as.numeric(regionsTested$p_val)
    # regionsTested$expected <- round(as.numeric(regionsTested$expected), 2)
    scores[[f]] <- regionsTested
    i <- i + 1
    if(i %% 10 == 0) cat("Processed ", i, "files\n")
  }
  
  #combine each iteration into one data frame
  final <- as.data.frame(do.call(rbind, scores))
  
  # final$count <- as.numeric(final$observed) + as.numeric(final$expected)
  

  minPval <- min(final$p_val[final$p_val>0])
  
  # final$p_val <- ifelse(final$p_val==0, minPval/abs(final$Log2FC), final$p_val)
  # 
  # final$p_val <- ifelse(final$p_val==0, minPval, final$p_val)
  
  final <- final %>%
    dplyr::mutate(count = observed + expected) %>% 
    dplyr::mutate(p_val = ifelse(p_val==0, minPval/abs(Log2FC), p_val)) %>% 
    dplyr::mutate(p_val = ifelse(p_val==0, minPval, p_val)) %>%
    # dplyr::mutate(padj = p.adjust(p_val, method = 'hochberg')) %>%
    dplyr::mutate(padj = p.adjust(p_val, method = 'hochberg')) %>%
    dplyr::mutate(sig = ifelse(padj <= 0.001, "***",
                                   ifelse(padj <= 0.01, "**",
                                          ifelse(padj <= 0.05, "*", "-")))) %>% 

    dplyr::mutate(eScore = round(abs(Log2FC) * -log10(padj),2)) %>% 
    # dplyr::mutate(eScore = round((abs(Log2FC) * 2) * -log10(padj),2)) %>% 
    dplyr::filter(count >= minHits) %>%
    dplyr::select(feature:p_val, padj, eScore, genotype, tissue, element, replicate, id, filename, -count) %>% 
    dplyr::arrange(-eScore, padj, -abs(Log2FC)) %>% 
    droplevels()
  
  
  if(nrow(final)==0){
    cat("Fewer than", minHits, "found in all files. Try adjusting 'minHits' to a lower number", sep=" ", "\n")
    return(final)
  }
  if(plot){
    cat("Plotting volcano plot", "\n")
    # print(Volcano(final))
    print(ggVolcano(df=final))
    }
  return(final)
}


# Volcano
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import dplyr
#' @importFrom plotly plot_ly
#' @export
Volcano <- function(d, byq=FALSE){
  feature_enrichment <- d
  
  var <- feature_enrichment$p_val
  printVar <- 'P-val'
  if(byq){
    var <- var <- feature_enrichment$padj
    printVar <- 'Adj-p'
  }
  
  maxLog2 <- max(abs(feature_enrichment$Log2FC[is.finite(feature_enrichment$Log2FC)]))
  maxLog2 <- as.numeric(round_any(maxLog2, 1, ceiling))
  
  ax <- list(size = 25)
  ti <- list(size = 25)
  
  p <- plotly::plot_ly(data = feature_enrichment,
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
          size = ~-log10(var)) %>% 
    layout(xaxis = list(title="Log2(FC)", titlefont = ax, range = c(-maxLog2, maxLog2)),
           yaxis = list(title=paste("-Log10(", printVar, ")", sep=''), titlefont = ax))
  p
}

# ggVolcano
#'
#' Plot the enrichment of SVs in genomic features
#' @keywords enrichment
#' @import dplyr ggplot2 ggrepel
#' @import RColorBrewer
#' @export
ggVolcano <- function(..., df, escore_threshold = 3){
  if(!nrow(df)) return(ggplot())
  red <- "#FC4E07"
  blue <- "#00AFBB"
  yellow <- "#E7B800"
  grey <- "grey"
  
  # feature_enrichment <- df
  
  # arguments <- list(...)
  # 
  # if(!is.null(arguments$slop)){
  #   slop <- arguments$slop
  #   title <- paste("Enrichment/depletion of genomic features\nfor breakpoints +/-", slop, 'bps')
  #   outFile <- paste("regionEnrichment_", slop,'.png', sep='')
  # } else {
  #   title <- "Enrichment/depletion of genomic features\nfor breakpoints"
  #   outFile <- "regionEnrichment.png"
  # }
  
  maxLog2 <- max(abs(df$Log2FC))
  maxLog2 <- round_any(maxLog2, 1, ceiling)
  
  df$colour <- ifelse(df$eScore >= escore_threshold, 'yes', 'no')
  df$label <- ifelse(df$eScore >= escore_threshold, df$feature, '')
  
  df$colour = factor(df$colour, levels=c("yes","no"), labels=c("***","ns")) 
  
  cols = c(red, grey)
  if(!nrow(df[df$eScore > escore_threshold,])) cols = c(grey)
  
  p <- ggplot(df, aes(Log2FC, -log10(padj)))
  p <- p + geom_point(aes(colour = colour), size=3, alpha=0.7)
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_colour_manual(values = cols)
  p <- p + geom_text_repel(aes(label = label), inherit.aes = TRUE)
  p <- p + guides(size = FALSE, colour = FALSE) 
  p <- p + cleanTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.title =element_text(size=rel(1.2))
    )
  p <- p + scale_y_continuous("-Log10(padj)")
  p <- p + scale_x_continuous("Log2(FC)", limits=c(-maxLog2, maxLog2), breaks=seq(-maxLog2, maxLog2, by=1))
  p
}

# theme(
#   panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
#   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12),
#   axis.text.y = element_text(size=12),
#   axis.title.x = element_blank(),
#   axis.title.y =element_text(size=15)

# plotMonteCarlo
#'
#' Plot the result of a monte carlo simulation (n shuffles of feature y)
#' Expects dataframe returned from bpRegionEnrichment()
#' @keywords expected frequency
#' @import dplyr ggplot2
#' @import RColorBrewer
#' @export

plotMonteCarlo <- function(x, bindWith = 10){
  
  x <- x %>% 
    dplyr::mutate(fillCol = ifelse(grepl("shuff", feature), "grey37", "#C72424FE"))

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
      dplyr::mutate(end = ifelse(type %in% c("TRA", "COMPLEX_TRA"), start+1, end)) %>% 
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
#' @importFrom plyr round_any
#' @import dplyr ggplot2
#' @export

bpRegionEnrichmentPlot <- function(..., feature_enrichment=NULL, bedDir=NULL, bp_data=NULL) {
  if(!missing(bedDir) && !missing(bp_data)) feature_enrichment <- bpRegionEnrichment(..., bp_data=bp_data, bedDir=bedDir, plot=F)
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
  maxLog2 <- plyr::round_any(maxLog2, 1, ceiling)
  
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
