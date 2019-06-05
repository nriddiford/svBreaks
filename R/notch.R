#' Print barplot showing samples with SVs affecting specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @param plot Show barplot [Default: TRUE]
#' @export
geneHit <- function(..., all_samples, filter_gene = "N", plot = TRUE) {
  if(missing(all_samples)) stop("\n[!] Must a file containing data for all samples (e.g. 'all_samples.txt'! Exiting.")
  
  all_data <- read.delim(all_samples, header = T)
 
    geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  gene_hits <- all_data %>% 
    dplyr::filter(...,
                  type != "-",
                  !status %in% c('F', 'aF'),
                  !type %in% c('COMPLEX_TRA')
                  ) %>% 
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>% 
    dplyr::group_by(sample, event) %>% 
    dplyr::filter(any(geneIn(filter_gene, affected_genes))) %>% 
    # dplyr::distinct(type, .keep_all = TRUE) %>% 
    # group_by(type) %>%  
    dplyr::mutate(start = as.integer(ifelse(stringr::str_detect(type, 'COMPLEX') && chromosome1 == "X", min(bp1), bp1))) %>% 
    dplyr::mutate(end = as.integer(ifelse(stringr::str_detect(type, 'COMPLEX') && chromosome2 == "X", max(bp2), bp2))) %>% 
    dplyr::mutate(length = ifelse(type %in% c('TRA', 'COMPLEX_TRA'), 2, (end-start)/1e3),
                  event_length = sum(length)) %>% 
    dplyr::mutate(type2 = as.character(ifelse(stringr::str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>% 
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
  
  if(plot){
    all_samples$star <- ifelse(all_samples$length==0, "X", '')
  
    # cols <- svBreaks::setCols(all_samples, col = 'type2', set = "Blues")
    # all_samples <- all_samples %>% 
    #   dplyr::mutate(colour = ifelse(type2 == "DEL", primary, 
    #                                 ifelse(type2 == "COMPLEX", secondary,
    #                                        ifelse(type2 == "BND", tertiary, 
    #                                               ifelse(type2 == "TRA", quaternary, 'red')
    #                                               )
    #                                        )
    #                                 )
    #                 )
    
    colours = sv_colours()
    
    
    p <- ggplot(all_samples, aes(fct_reorder(sample_mod, -length), length, fill = fct_reorder(type2, -length), label = star))
    p <- p + geom_bar(alpha = 0.8, stat = "identity")
    p <- p + scale_y_continuous("Length (Kb)", expand = c(0, 0.5))
    p <- p + cleanTheme() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
        axis.title.x = element_blank()
      )
    p <- p + ggtitle(paste("Samples with SVs affecting", filter_gene))
    p <- p + scale_fill_manual("SV type\n", values = colours)
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
tally_hits <- function(..., all_samples, bp_data, filter_gene = "N", plot=FALSE, freq=FALSE) {
  if(missing(all_samples)) stop("\n[!] Must specify a file containing data for all samples (e.g. 'all_samples.txt'! Exiting.")
  
  if(missing(bp_data)){
    nhits <- geneHit(..., all_samples=all_samples, plot=F)
  } else {
    nhits <- bp_data
  }
  tally <- nhits %>% 
    dplyr::mutate(class = type2) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::arrange(sample, -allele_frequency) %>% 
    dplyr::distinct(sample, event, .keep_all = TRUE) %>% # added event 7.3.19
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
notchHits <- function(..., all_samples, filter_gene = "N", show_samples=FALSE, bp_density=FALSE, from=2.7, to=3.4, ticks = 50) {
  if(missing(all_samples)) stop("\n[!] Must a file containing data for all samples (e.g. 'all_samples.txt'! Exiting.")
  
  notch_data <- svBreaks::geneHit(..., plot=F, all_samples=all_samples, filter_gene=filter_gene)
  notch_data <- notch_data %>% 
    dplyr::mutate(bp1 = ifelse(chromosome1 == "X", bp1, bp2-100),
                  bp2 = ifelse(chromosome2 == "X", bp2, bp1+100),
                  type2 = factor(type2)) %>% 
    dplyr::mutate(sampleax = as.numeric(sample),
                  bp1 = bp1 / 1e6,
                  bp2 = bp2 / 1e6) %>% 
    dplyr::mutate(bp1 = ifelse(bp1 <= from, from, bp1),
                  bp2 = ifelse(bp2 >= to, to, bp2)
                  ) %>% 
    dplyr::arrange(-length)
  
  
  colours = sv_colours()
  
  p <- ggplot(notch_data)
  p <- p + geom_rect(aes(xmin=bp1, xmax=bp2, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill=type2), color="black", alpha=0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp1, xmax = bp1 + 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp2, xmax = bp2 - 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)

  p <- p + guides(color = FALSE, size = FALSE, sampleax = FALSE, type2 = FALSE)

  p <- p + scale_y_continuous("Sample", expand = c(0.01, 0.01), breaks = seq(levels(as.factor(notch_data$sampleax))), labels = levels(notch_data$sample))
  p <- p + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(from, to, by = ticks/1e3), limits = c(from, to))
  p <- p + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)
  
  p <- p + slideTheme() +
    theme(
      axis.title.y = element_blank(),
      # panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      panel.grid.major.x = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
      legend.position = "top",
      axis.text.y = element_blank(),
      axis.title.x=element_blank()
    )
  if(show_samples){
    p <- p + theme(axis.text.y = element_text(size=20))
  }
  
  p <- p + scale_fill_manual("SV type\n", values = colours)
  
  
  if(bp_density){

    long_notch_hits <- tidyr::gather(notch_data, bp_n, breakpoint, bp1, bp2, factor_key=TRUE) %>% 
      dplyr::select(sample:chromosome1, type2, sampleax, bp_n, breakpoint)
  
    blueBar <- '#3B8FC7'
    
    p3 <- ggplot(long_notch_hits)
    p3 <- p3 + geom_density(aes(breakpoint, fill = blueBar), alpha = 0.6)
    p3 <- p3 + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(2.7, 3.4, by = 0.05), limits = c(2.70, 3.4))
    p3 <- p3 + scale_y_continuous("Density", expand = c(0, 0))
    p3 <- p3 + guides(colour = FALSE)
    p3 <- p3 + geom_rug(data = long_notch_hits, aes(breakpoint))
    p3 <- p3 + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)

    p3 <- p3 + slideTheme() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
        )

    p3 <- p3 + scale_fill_identity("SV type\n")

    combined_plots <- ggpubr::ggarrange(
    p, p3,
    labels = c("A", "B"),
    ncol = 1, nrow = 2,
    heights = c(7,3.5),
    align = 'v'
    )
    combined_plots
  } else p
  
}






# ######

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
      dplyr::filter(!(chrom == "X" & bp >= start & bp < stop & bp2 <= stop)) %>%
      dplyr::filter(!gene %in% c('N', 'dnc', 'kirre'), !gene2 %in% c('N', 'dnc', 'kirre'))  %>%
      droplevels()
    return(noNotch)
  }
}


#' SV colours
#' 
sv_colours <- function(){
  return(c("DEL" = '#0073C299',
           "COMPLEX" = '#EFC00099',
           "BND" = '#86868699',
           "TRA" = '#CD534C99',
           "NA" = "grey")
         )
}
