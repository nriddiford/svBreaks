# A collection of summary functions
# for plotting breakpints per sample

#' bpStats
#' Get some basic stats for breakpoints
#' @keywords stats
#' @import dplyr
#' @import ggplot2
#' @importFrom forcats fct_reorder
#' @export
#' @param colSample Set to <samplename> to colour that sample red in the plot

bpStats <- function(..., bp_data=NULL, colSample=NA, write=FALSE) {
  if(missing(bp_data)){
    bp_data <- getData(...)
  }
  
  blueBar <- '#3B8FC7'
  
  sampleSvs <- bp_data %>%
    dplyr::filter(...,
                  bp_no == "bp1") %>%
    dplyr::group_by(sample, genotype) %>%
    dplyr::tally() %>%
    dplyr::mutate(som_count = ifelse(genotype=='somatic_tumour', n/2, 0))
  
  if (!is.na(colSample)) {
    sampleSvs$colour <- ifelse(sampleSvs$sample %in% exclude_samples, "#C72424FE", "grey37")
  } else {
    sampleSvs$colour <- blueBar
  }

  p <- ggplot(sampleSvs)
  p <- p + geom_bar(aes(forcats::fct_reorder(sample, -som_count), n, fill = colour), stat = "identity")
  # p <- p + scale_y_continuous("Number of variants", limits = c(0, 20), expand = c(0.01, 0.01), breaks=seq(0,20,by=2))
  p <- p + scale_y_continuous("Number of variants", expand = c(0.01, 0.01))

  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 30),
      axis.title.x = element_blank()
    )
  p <- p + scale_fill_identity()
  # p <- p + facet_wrap(~genotype, ncol = 1, scales = "free_y")

  if(write){
    sampleSVs <- paste("SVs_sample.png")
    cat("Writing file", sampleSVs, "\n")
    ggsave(paste("plots/", sampleSVs, sep = ""), width = 20, height = 10)
  }
  
  cat("sample", "SVs", sep = "\t", "\n")
  bp_data <- bp_data %>% 
    dplyr::filter(...,
                  genotype == 'somatic_tumour',
                  bp_no == "bp1") %>% 
    droplevels()
  
  rank <- sort(table(bp_data$sample), decreasing = TRUE)
  rank <- as.array(rank)

  total <- 0

  scores <- list()
  for (i in 1:nrow(rank)) {
    cat(names(rank[i]), rank[i], sep = "\t", "\n")
    total <- total + rank[i]
    scores[i] <- rank[i]
  }
  cat("--------------", "\n")
  scores <- unlist(scores)

  mean <- as.integer(mean(scores))
  med <- as.integer(median(scores))

  cat("total", total, sep = "\t", "\n")
  cat("samples", nrow(rank), sep = "\t", "\n")

  cat("--------------", "\n")
  cat("mean", mean, sep = "\t", "\n")
  cat("median", med, sep = "\t", "\n")

  cat("\n")
  dels <- nrow(dplyr::filter(bp_data, type == "DEL"))
  invs <- nrow(dplyr::filter(bp_data, type == "INV"))
  dups <- nrow(dplyr::filter(bp_data, type == "DUP"))
  tandups <- nrow(dplyr::filter(bp_data, type == "TANDUP"))
  tra <- nrow(dplyr::filter(bp_data, type == "TRA"))

  cat("Dels ", dels, sep = "", "\n")
  cat("Invs ", invs, sep = "", "\n")
  cat("Dups ", dups, sep = "", "\n")
  cat("Tandups ", tandups, sep = "", "\n")
  cat("Tra ", tra, sep = "", "\n")

  p
}


#' bpFeatures
#' Get some basic stats for breakpoints
#' @param notch - Set to 1 to exclude Notch events
#' @keywords stats
#' @import dplyr ggplot2
#' @export
bpFeatures <- function(..., notch=0) {
  if (notch) {
    bp_data <- notchFilt()
    ext <- "_excl.N.pdf"
  }
  else {
    bp_data <- getData(..., !sample %in% c("A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"))
    ext <- ".pdf"
  }

  # To condense exon counts into "exon"
  bp_data$feature <- as.factor(gsub("_.*", "", bp_data$feature))

  bp_data <- dplyr::filter(bp_data, feature != "pseudogene" & feature != "snoRNA" & feature != "stop" & feature != "tRNA")
  bp_data <- droplevels(bp_data)
  # Reoders descending
  bp_data$feature <- factor(bp_data$feature, levels = names(sort(table(bp_data$feature), decreasing = TRUE)))

  # cols <- setCols(bp_data, "feature")

  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(feature, fill = feature), stat = "count", position = position_dodge())
  # p<-p + cols
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous(expand = c(0.01, 0.01))

  features_outfile <- paste("Breakpoints_features_count", ext, sep = "")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep = ""), width = 20, height = 10)

  p
}


#' svsbySample
#' Plot the number of structural breakpoints per sample
#' @keywords samples
#' @import tidyverse
#' @importFrom forcats fct_reorder
#' @export
svsbySample <- function(..., bp_data=NULL, colSample=NA) {
  # excludedSamples <- c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
  if(missing(bp_data)){
    bp_data <- getData(...)
  }
  
  sample_names <- levels(bp_data$sample)
  
  sampleSvs <- bp_data %>%
    dplyr::filter(bp_no != "bp2") %>% 
    dplyr::filter(genotype == 'somatic_tumour') %>% 
    droplevels() %>% 
    dplyr::group_by(sample) %>%
    dplyr::tally()
  
  missing_samples = list()
  
  for(i in 1:length(sample_names)){
    if(!(sample_names[i] %in% levels(sampleSvs$sample))){
      missing_samples[i] <- sample_names[i]
    }
  }
  
  missing_samples <- compact(missing_samples)
  
  dat <- data.frame(sample = c(levels(sampleSvs$sample), paste(missing_samples)), count = (c(sampleSvs$n, rep(0, length(missing_samples)) ))) 
                    
                    
  if(!is.na(colSample)) {
    dat$colour <- ifelse(dat$sample == colSample, "#C72424FE", "#636161FE")
  } else {
    dat$colour <- "grey37"
  }
  
  
  p <- ggplot(dat)
  p <- p + geom_bar(aes(forcats::fct_reorder(sample, -count), count, fill = colour), stat = "identity")
  p <- p + scale_x_discrete(expand = c(0.01, 0.01))
  p <- p + scale_y_continuous("Number of SVs", expand = c(0.01, 0.01), limits = c(0, max(dat$count)), breaks = seq(0,max(dat$count),by=5))
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=18)
    )
  p <- p + scale_fill_identity()
  
  # svsSample <- 'svsBySample.png'
  # cat("Writing file", svsSample, "\n")
  # ggsave(paste("plots/", svsSample, sep = ""), width = 10, height = 10)
  
  print(p)
  dat <- dat %>%
    dplyr::select(-colour) %>% 
    dplyr::arrange(-count)
  
  return(dat)
}

#' genesbySample
#' Plot the number of structural breakpoints per sample
#' @keywords genes
#' @import dplyr 
#' @importFrom forcats fct_reorder
#' @export
genesbySample <- function(..., affected_genes = "inst/extdata/all_genes_filtered.txt", genome_length=118274340) {
  
  # bp_data <- getData(..., genotype=='somatic_tumour', !sample %in% excluded_samples)
  
  allGenes <- read.delim(affected_genes, header = T)
  colnames(allGenes) <- c("event", "sample", "genotype", "type", "af", "chrom", "gene")
  
  allGenes <- allGenes %>% 
    dplyr::filter(...,
                  genotype=='somatic_tumour') %>%
    dplyr::mutate(cell_fraction = ifelse(chrom %in% c('X', 'Y'), af,
                                         ifelse(af*2>1, 1, af*2))) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(affected_genes = n()) %>% 
    dplyr::arrange(-affected_genes) %>% 
    droplevels()
  
  p <- ggplot(allGenes, aes(forcats::fct_reorder(sample, -affected_genes), affected_genes))
  p <- p + geom_bar(stat='identity')
  p <- p + scale_y_continuous("Number of genes", expand = c(0.01, 0.01))
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=40),
      
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=18)
    )

  print(p)
  
  return(allGenes)
}

complexEvents <- function(..., all_samples = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_samples.txt'){
  all_data <- read.delim(all_samples, header = T)

  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  complex_count <- all_data %>% 
    dplyr::filter(
                  is.na(status)) %>% 
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>% 
    dplyr::group_by(sample, event) %>% 
    dplyr::mutate(newCol = as.character(ifelse(any(geneIn('N', affected_genes)),'Notch', 'Other'))) %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome1 == "X", allele_frequency,
                                   ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>% 
    dplyr::group_by(sample, event) %>% 
    dplyr::mutate(complex_count = ifelse(any(type %in% c("COMPLEX_DEL", "COMPLEX_BND", "COMPLEX_DUP")), count(type), 0)) %>% 
    dplyr::filter(complex_count >= 2) %>% 
    dplyr::tally() %>% 
    dplyr::filter(n>=2) %>% 
    dplyr::ungroup() %>% 
    # dplyr::select(sample, newCol, allele_frequency, cell_fraction, affected_genes, complex_count, n) %>% 
    droplevels()
    
    
   p <- ggplot(complex_count)
   p <- p + geom_bar(aes(forcats::fct_reorder(sample, -n), complex_count), alpha = 0.7, stat = "identity")
   p
    
}

length_by_sample <- function(..., bp_data=NULL){
  if(missing(bp_data)){
    bp_data <- getData(...)
  }
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  df <- bp_data %>% 
    dplyr::filter(!str_detect(type, 'TRA')) %>% 
    dplyr::filter(type %in% c('DUP', 'DEL'),
                  sample != 'D106R25') %>% 
    dplyr::distinct(sample, event, type, .keep_all=TRUE) %>% 
    dplyr::group_by(sex) %>% 
    # dplyr::filter(n()>=2) %>% 
    dplyr::mutate(av_length = sum(length)/n())

  p <- ggplot(df, aes(forcats::fct_reorder(sex, -av_length), length+1))
  p <- p + geom_violin(aes(fill=sex), alpha=.8)
  p <- p + geom_jitter(position=position_jitter(0.2), size = .8)
  p <- p + scale_y_log10("SV length (Kb)")
  p <- p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          # axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5),
          axis.text = element_text(size=20), axis.title = element_text(size=20),
          axis.title.x = element_blank()
    )
  p <- p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5)
  
  p <- p + facet_wrap(~type)
  p
  
}

  