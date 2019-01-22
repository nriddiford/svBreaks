#' Functions related to Notch
#' notchFilt
#'
#' Function to filter in/out events affecting a locus (in this case the Drosophila Notch locus)
#' @param infile File to process [Required]
#' @keywords parse
#' @import tidyverse
#' @export
#' @return Dataframe
#'
notchFilt <- function(..., keep=NULL, start=2700000, stop=3400000) {
  excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  
  bp_data <- getData(..., !sample %in% excluded_samples)
  if(!missing(keep)){
    cat("Selecting for bps in Notch\n")
    notchIn <- bp_data %>%
      dplyr::filter(chrom == "X" & chrom2 == "X") %>%
      dplyr::filter(bp >= start && bp < stop) %>% 
      dplyr::filter(bp2 < stop) %>% 
      droplevels()
    return(notchIn)
  } else{
    cat("Excluding bps in Notch\n")
    noNotch <- bp_data %>%
      filter(!(chrom == "X" & bp >= start & bp < stop & bp2 <= stop)) %>%
      filter(!gene %in% c('N', 'dnc', 'kirre'), !gene2 %in% c('N', 'dnc', 'kirre'))  %>%
      droplevels()
    return(noNotch)
  }
}

#' Print barplot showing samples with SVs affecting specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @param plot Show barplot [Default: TRUE]
#' @import stringr
#' @export
geneHit <- function(..., all_samples = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_samples.txt', filter_gene = "N", plot = TRUE) {
  all_data <- read.delim(all_samples, header = T)
  # excluded_samples <- c("A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  # excluded_samples <- c("D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  # excluded_samples <- c()
    geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  gene_hits <- all_data %>% 
    dplyr::filter(...,
                  !type %in% c('COMPLEX_TRA'),
                  is.na(status)) %>% 
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>% 
    dplyr::group_by(sample, event) %>% 
    dplyr::filter(any(geneIn(filter_gene, affected_genes))) %>% 
    # dplyr::distinct(type, .keep_all = TRUE) %>% 
    # group_by(type) %>%  
    dplyr::mutate(start = as.integer(ifelse(type %in% c("COMPLEX_DEL", "COMPLEX_BND", "COMPLEX_TANDUP", "COMPLEX_DUP") && chromosome1 == "X", min(bp1), bp1))) %>% 
    dplyr::mutate(end = as.integer(ifelse(type %in% c("COMPLEX_DEL", "COMPLEX_BND", "COMPLEX_TANDUP", "COMPLEX_DUP") && chromosome2 == "X", max(bp2), bp2))) %>% 
    dplyr::mutate(length = ifelse(type %in% c('TRA', 'COMPLEX_TRA'), 2, (end-start)/1e3),
                  event_length = sum(length)) %>% 
    dplyr::mutate(type2 = as.character(ifelse(type %in% c("COMPLEX_DEL", "COMPLEX_BND", "COMPLEX_TRA", "COMPLEX_DUP", "COMPLEX_TANDUP"), 'COMPLEX', as.character(type)))) %>% 
    dplyr::mutate(event_count = n_distinct(event)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(n = n_distinct(event)) %>% 
    dplyr::distinct(event, .keep_all = TRUE) %>%
    dplyr::ungroup() %>% 
    dplyr::arrange(sample, -allele_frequency) %>% 
    dplyr::mutate(sample_mod = case_when(n >1 ~ make.unique(as.character(sample)), 
                                  TRUE ~ as.character(sample))) %>% 
    # dplyr::mutate(sample_mod = as.character(ifelse(event_count == 1, as.character(sample), paste(sample, event_count, sep = '.')))) %>% 
    dplyr::arrange(desc(length)) %>% 
    dplyr::select(sample, sample_mod, event, chromosome1, chromosome2, type, type2, bp1, start, bp2, end, length, event_length, allele_frequency, cn, affected_genes) %>% 
    droplevels()
  
  sample_names <- all_data %>% 
    dplyr::filter(...) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::distinct(sample) %>% 
    droplevels()
  
  sample_names <- levels(sample_names$sample)
  missing_samples = list()
   
  for(i in 1:length(sample_names)){
    if(!(sample_names[i] %in% levels(gene_hits$sample))){
      missing_samples[i] <- sample_names[i]
    }
  }
   
  missing_samples <- plyr::compact(missing_samples)
  
  dat <- data.frame(sample_mod = unlist(missing_samples), length = 0, type2 = "NA")
  
  all_samples <- plyr::join(gene_hits, dat, type='full')
  
  # dat <- data.frame(sample = c(levels(as.factor(gene_hits$sample_mod)), paste(missing_samples)),
  #                   length = c(gene_hits$length, rep(0, length(missing_samples)))
  # )
  
  if(plot){
    all_samples$star <- ifelse(all_samples$length==0, "X", '')
    
    p <- ggplot(all_samples, aes(fct_reorder(sample_mod, -length), length, fill = fct_reorder(type2, -length), label = star))
    p <- p + geom_bar(alpha = 0.7, stat = "identity")
    p <- p + scale_y_continuous("Length (Kb)", expand = c(0, 0.5))
    p <- p + cleanTheme() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
        axis.title.x = element_blank()
      )
    # p <- p + scale_alpha_continuous(range = c(0.1, 1))
    p <- p + ggtitle(paste("Samples with SVs affecting", filter_gene))
    p <- p + scale_fill_jco()
    p <- p + geom_text(vjust=-5)
    # p <- p + geom_text(data = dat, aes(x=sample_mod, y=50),  label = "X")
    
    # p <- p + coord_flip()
    # p <- p + scale_y_reverse()
    
    p
  } else return(gene_hits)
}



#' Tally type of sv events in specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @export
tally_hits <- function(..., all_samples = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_samples.txt', filter_gene = "N", plot=FALSE, freq=FALSE) {
  nhits <- geneHit(..., plot=F)
  tally <- nhits %>% 
    dplyr::mutate(class = type2) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::arrange(sample, -allele_frequency) %>% 
    dplyr::distinct(sample, .keep_all = TRUE) %>% 
    dplyr::group_by(class) %>%
    dplyr::summarise(count = n())
  
  tally$frequency <- round(tally$count/sum(tally$count)*100,1)
  
  if(plot){
    p <- ggplot(tally)
    if(freq)
      p <- p + geom_bar(aes(fct_reorder(class, -frequency), frequency, fill = fct_reorder(class, -frequency)), alpha = 0.7, stat='identity')
    else
      p <- p + geom_bar(aes(fct_reorder(class, -count), count, fill = fct_reorder(class, -count)), alpha = 0.7, stat='identity')
    p <- p + scale_fill_jco()
    p <- p + cleanTheme() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
        axis.title.x = element_blank()
      )
    p
  } else {
    tally <- tally %>% 
      dplyr::mutate(TOTAL = sum(count))
    return(tally)
  }
}


#' notchHits
#' Plot the breakpint density over Notch
#' @keywords Notch
#' @import tidyverse
#' @export
notchHits <- function(infile = "inst/extdata/Notch_hits.txt") {
  # notch_data <- read.delim(infile, header = T)
  # excluded_samples <- c("A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  
  notch_data <- geneHit(..., plot=F)
  notch_data <- notch_data %>% 
    dplyr::mutate(bp1 = ifelse(chromosome1 == "X", bp1, bp2-2),
                  bp2 = ifelse(chromosome2 == "X", bp2, bp1+2)) %>% 
    dplyr::mutate(sampleax = as.numeric(sample),
                  bp1 = bp1 / 1e6,
                  bp2 = bp2 / 1e6)
   
  # refGene<-read.delim(refgene_file, header=F)
  # colnames(refGene)<-c('Tno',	'id',	'chrom',	'strand',	'Tstart',	'Tstop',	'cdsstart',	'cdsStop',	'exonNo',	'Estart',	'Estop',	'score',	'altname',	'CDSstartStat',	'CDSEndStat',	'exonFrames')
  # # Tid	id	chrom	strand	Tstart	Tstop	cdsstart	cdsStop	exon#	Estart	Estop	score	altname	CDSstartStat	CDSEndStat	exonFrames
  # gene='N'
  # refGene<-dplyr::filter(refGene, altname == gene)
  # refGene<-droplevels(refGene)
  #
  # wStart<-refGene$Tstart - 1000
  # wEnd<-refGene$stop + 1000

  # bp_data<-getData()

  # bp_data<-reshape(bp_data, idvar=c("event", "sample"), timevar="bp_no", direction="wide")
  # #bp_data<-dplyr::filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)
  # bp_data<-dplyr::select(bp_data, sample, chrom.bp1, bp.bp1, gene.bp1, bp.bp2,gene.bp2, type.bp1)
  # colnames(bp_data)<-c("sample", "chrom1", "bp1", "gene1", "bp2", "gene2", "type")
  
  # cols<-setCols(bp_data, "type", set="Set1")

  # if (bp_data$type == "TRA") {
  #   bp_data$bp1 <- bp_data$bp2 - 100
  # }


  p <- ggplot(notch_data)
  p <- p + geom_rect(aes(xmin=bp1, xmax=bp2, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill=type2), color="black", alpha=0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp1, xmax = bp1 + 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp2, xmax = bp2 - 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)

  p <- p + guides(color = FALSE, size = FALSE, sampleax = FALSE, type2 = FALSE)

  p <- p + scale_y_continuous("Sample", expand = c(0.01, 0), breaks = seq(levels(as.factor(notch_data$sampleax))), labels = levels(notch_data$sample))
  p <- p + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(2.7, 3.4, by = 0.05), limits = c(2.70, 3.4))
  p <- p + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)
  

  p <- p + slideTheme() +
    theme(
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
      legend.position = "top",
      axis.text.y = element_blank()
    )

  # p <- p + scale_fill_identity()
  
  dummy = data.frame(x = (seq(2.7,3.4,by=0.2)))
  dummy$y <- 0
  p2 <- ggplot(dummy, aes(x,y))
  p2 <- p2 + scale_y_discrete(expand = c(0, 0))
  p2 <- p2 + annotate("rect", xmin = 2.740000, xmax = 3.134532, ymin = -1.5, ymax = 0, alpha = .4, fill="#CFAEAEFE")
  p2 <- p2 + annotate("rect", xmin = 3.134870, xmax = 3.172221, ymin = -1.5, ymax = 0, alpha = .4, fill="#8FBD80FE")
  p2 <- p2 + annotate("rect", xmin = 3.176440, xmax = 3.334000, ymin = -1.5, ymax = 0, alpha = .4, fill="#A9D0DEFE")
  p2 <- p2 + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)
  p2 <- p2 + slideTheme()+
    theme(axis.text = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank())
          # p<-p + cols

  # Nhits <- paste("Notch_hits.pdf")
  # cat("Writing file", Nhits, "\n")
  # ggsave(paste("plots/", Nhits, sep=""), width = 10, height = 10)
  #

  # Density plot
  # cols2<-setCols(bp_data, "type", fill='Y', set="Set1")

  long_notch <- gather(notch_data, bp, bp1,bp2)

  p3 <- ggplot(notch_data)
  p3 <- p3 + geom_density(aes(bp / 1000000, fill = "royalblue4"), alpha = 0.6)
  p3 <- p3 + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(2.7, 3.4, by = 0.05), limits = c(2.70, 3.4))
  p3 <- p3 + scale_y_continuous("Density", expand = c(0, 0))
  p3 <- p3 + guides(colour = FALSE)
  p3 <- p3 + geom_rug(data = notch_data, aes(bp / 1000000))
  p3 <- p3 + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)

  p3 <- p3 + slideTheme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      axis.title.y = element_blank()
      # strip.text = element_text(size=10)
    )

  # p2 <- p2 + cols2
  p3 <- p3 + scale_fill_identity()

  # p2 <- p2 + facet_wrap(~type, nrow = 3)

  combined_plots <- ggpubr::ggarrange(
    p, p2, p3,
    labels = c("A", "B", "C"),
    ncol = 1, nrow = 3,
    heights = c(5,1,5),
    align = 'v'
  )

  NhitsDen <- paste("Notch_hits_density.png")
  cat("Writing file", NhitsDen, "\n")
  ggsave(paste("plots/", NhitsDen, sep = ""), width = 30, height = 20)

  combined_plots
}


#' notchDels
#' Plot the size of deletions affecting Notch
#' @keywords Notch
#' @import dplyr ggplot2 ggsci forcats
#' @export
notchDels <- function(infile = "inst/extdata/Notch_hits.txt") {
  d <- read.delim(infile, header = T)
  excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  
  notch_data <- d %>% 
    dplyr::rename(length = length.Kb.) %>%
    dplyr::filter(!sample %in% excluded_samples,
                  is.na(status)) %>% 
    dplyr::mutate(length = ifelse(type == 'TRA', 2, length)) %>% 
    dplyr::select(sample, event, source, type, chromosome1, bp1, chromosome2, bp2, length, affected_genes) %>% 
    dplyr::mutate(sampleax = as.numeric(sample)) %>% 
    droplevels()
  
  # notch_data <- notchFilt(keep=1)
  bp_data <- getData()
  sample_names <- levels(bp_data$sample)
  
  # Problem with events not being clustered properly for CNV-Seq
  notch_data <- notch_data %>%
    dplyr::group_by(sample) %>%
    # dplyr::filter(!duplicated(event)) %>%
    dplyr::mutate(allLength = sum(length)) %>% 
    dplyr::ungroup()


  missing_samples = list()
  
  for(i in 1:length(sample_names)){
    if(!(sample_names[i] %in% levels(notch_data$sample))){
      missing_samples[i] <- sample_names[i]
    }
  }
  
  missing_samples <- plyr::compact(missing_samples)
  dat <- data.frame(sample = c(levels(notch_data$sample), paste(missing_samples)),
                    length = c(notch_data$length, rep(0, length(missing_samples)))
                    ) 
  

  p <- ggplot(notch_data)
  p <- p + geom_bar(aes(fct_reorder(sample, -allLength), length, fill = type), alpha = 0.7, stat = "identity")
  p <- p + scale_y_continuous("Length (Kb)")
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank()
    )
  p <- p + scale_fill_jco()
  p
  dels_out <- paste("NotchDels.png")
  cat("Writing file", dels_out, "\n")
  ggsave(paste("plots/", dels_out, sep = ""), width = 30, height = 15)

  p
}
