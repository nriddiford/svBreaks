#' Print barplot showing samples with SVs affecting specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @param plot Show barplot [Default: TRUE]
#' @importFrom plyr compact join
#' @export
geneHit <- function(..., all_samples, drivers=c("N"), all_samples_df=NULL, plot=TRUE, combined_plot=FALSE, show_sample=TRUE) {
  if(missing(all_samples) && missing(all_samples_df)) stop("\n[!] Must a file containing data for all samples (e.g. 'all_samples.txt'! Exiting.")
  
  if(!missing(all_samples_df)) {
    all_data <- all_samples_df
  } else {
    all_data <- read.delim(all_samples, header = T) 
  }
  all_data <- all_data %>% 
    dplyr::mutate(allele_frequency = as.double(as.character(allele_frequency)))
  
  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  ## This is working 15.1.20
  gene_hits <- all_data %>%
    dplyr::filter(
      ...,
      type != "-",
      !status %in% c('F', 'aF'),
      allele_frequency >= 0.1, # hack to exclude large BND in R33
      !type %in% c('COMPLEX_TRA')
    ) %>%
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>%
    dplyr::group_by(sample, event) %>%
    dplyr::filter(any(geneIn(drivers, affected_genes))) %>% 
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
    dplyr::mutate(sample_mod = sample) %>%
    
    # dplyr::mutate(sample_mod = case_when(n >1 ~ make.unique(as.character(sample)),
    #                                      TRUE ~ as.character(sample))) %>%
    dplyr::arrange(desc(length)) %>%
    dplyr::select(sample, sample_mod, event, chromosome1, chromosome2, type, type2, bp1, start, bp2, end, length, event_length, allele_frequency, cn, affected_genes) %>%
    droplevels()
  
  print(gene_hits)
  
  sample_names <- all_data %>% 
    dplyr::filter(...) %>%
    dplyr::group_by(sample) %>% 
    dplyr::distinct(sample) %>% 
    droplevels()
  
  sample_names <- levels(sample_names$sample)
  
  missing_samples = list()
  if(!combined_plot){ 
    for(i in 1:length(sample_names)){
      if(!(sample_names[i] %in% levels(gene_hits$sample))){
        missing_samples[i] <- sample_names[i]
      }
    }
    
    missing_samples <- plyr::compact(missing_samples)
    dat <- data.frame(sample_mod = unlist(missing_samples), length = 0, type2 = "NA")
    all_samples <- plyr::join(gene_hits, dat, type='full')
  } else{
    all_samples <- gene_hits
  }
  
  if(plot){
    all_samples$star <- ifelse(all_samples$length==0, "X", '')
    colours = sv_colours()
    if(combined_plot){
      p <- ggplot(all_samples, aes(fct_reorder(sample_mod, length), length, fill = type2, colour = type2, label = star))
    } else {
      p <- ggplot(all_samples, aes(fct_reorder(sample_mod, -length), length, fill = type2, colour = type2, label = star))
    }
    p <- p + geom_bar(alpha = 0.7, stat = "identity")
    p <- p + guides(colour = FALSE)
    p <- p + scale_y_continuous("Length (kb)", expand = c(0, 0.5))
    p <- p + cleanTheme() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=18),
        axis.title.x = element_blank()
      )
    # p <- p + ggtitle("Structural Variants affecting Notch")
    p <- p + scale_fill_manual("SV type\n", values = colours)
    p <- p + scale_colour_manual(values = colours)
    
    p <- p + geom_text(vjust=-5)
    if(!show_sample) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())
    
    p
  } else return(gene_hits)
}


#' Tally type of sv events in specified gene
#' @param infile File to process [Required]
#' @param filter_gene The gene of interest [Default: N]
#' @importFrom ggsci scale_fill_jco
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
    p <- p + ggsci::scale_fill_jco()
    p <- p + s
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
#' @import dplyr tidyr
#' @importFrom cowplot plot_grid
#' @export
notchHits <- function(..., all_samples, all_samples_df=NULL, drivers = c("N"), show_samples=FALSE, barchart=FALSE, bp_density=FALSE, adjust=0.3, from=2.7, to=3.5, ticks = 50) {
  if(missing(all_samples) && missing(all_samples_df)) stop("\n[!] Must a file containing data for all samples (e.g. 'all_samples.txt'! Exiting.")
  

  if(!missing(all_samples_df)) {
    notch_data <- svBreaks::geneHit(..., plot=F, all_samples_df=all_samples_df, drivers=drivers)
  } else {
    notch_data <- svBreaks::geneHit(..., plot=F, all_samples=all_samples, drivers=drivers)
  }

  # if(barchart) p1 <- svBreaks::geneHit(..., all_samples_df=all_samples_df, show_sample = T, combined_plot=barchart, drivers=drivers, plot=T)

  notch_data <- notch_data %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(length_adj = as.double(length + runif(1, 0, .1))) %>% 
    dplyr::mutate(bp1 = start / 1e6,
                  bp2 = end / 1e6) %>% 
    ungroup()
  
  notch_data <- notch_data %>% 
    dplyr::mutate(bp1 = ifelse(chromosome1 == "X", bp1, bp2-(1/1e4)),
                  bp2 = ifelse(chromosome2 == "X", bp2, bp1+(1/1e4)),
                  type2 = factor(type2)) %>% 
    dplyr::mutate(sampleax = as.numeric(sample),
                  rank = dense_rank(length_adj)) %>% 
    dplyr::mutate(bp1 = ifelse(bp1 <= from, from, bp1),
                  bp2 = ifelse(bp2 >= to, to, bp2)
                  ) %>% 
    dplyr::arrange(-length)
  
  colours = sv_colours()
  
  # notch_data$sampleax <- fct_reorder(as.factor(notch_data$sampleax), notch_data$length)
  p <- ggplot(notch_data)
  p <- p + geom_rect(aes(xmin=bp1, xmax=bp2, ymin=(rank-0.5),ymax=(rank+0.5),fill=type2, color=type2), alpha=0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp1, xmax = bp1 + 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)
  # p <- p + geom_rect(data = notch_data, aes(xmin = bp2, xmax = bp2 - 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)

  p <- p + guides(color = FALSE, size = FALSE, sampleax = FALSE, type2 = FALSE)

  p <- p + scale_y_continuous( expand = c(0.01, 0.01), breaks = seq(notch_data$rank), labels = levels(fct_reorder(notch_data$sample_mod, notch_data$rank)))
  p <- p + scale_x_continuous("", expand = c(0, 0), breaks = seq(from, to, by = ticks/1e3), limits = c(from, to))
  p <- p + geom_vline(xintercept = 3.135669, linetype = "solid", size = .3)
  p <- p + geom_hline(yintercept = notch_data$rank, linetype = "dotted", size = 0.2)
  
  
  p <- p + cleanTheme() +
    theme(
      axis.title.y = element_blank(),
      # panel.grid.major.x = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # panel.grid.major = element_line(color = "grey80", size = 0.2, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20),
      legend.position = "top",
      axis.text.y = element_blank(),
      axis.title.x=element_blank()
    )
  
  if(show_samples) p <- p + theme(axis.text.y = element_text(size=10))
  
  p <- p + scale_fill_manual("SV type\n", values = colours)
  p <- p + scale_colour_manual(values = colours)
  
  if(bp_density){
    p <- p + theme(axis.text.x = element_blank())

    long_notch_hits <- tidyr::gather(notch_data, bp_n, breakpoint, bp1, bp2, factor_key=TRUE) %>% 
      dplyr::filter(breakpoint > from, breakpoint < to) %>% 
      dplyr::select(sample:chromosome1, type2, sampleax, bp_n, breakpoint)
    
    blueBar <- '#3B8FC7'
    
    p3 <- ggplot(long_notch_hits)
    p3 <- p3 + stat_density(aes(breakpoint), alpha = 0.6, adjust = adjust)
    p3 <- p3 + scale_x_continuous("", expand = c(0, 0), breaks = seq(from, to, by = ticks/1e3), limits = c(from, to))
    p3 <- p3 + scale_y_continuous("Density", expand = c(0, 0))
    p3 <- p3 + guides(colour = FALSE)
    p3 <- p3 + geom_rug(data = long_notch_hits, aes(breakpoint))
    p3 <- p3 + geom_vline(xintercept = 3.135669, linetype = "solid", size = .3)
    p3 <- p3 + cleanTheme() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
        axis.title.x = element_blank(),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
        )

    regions <- p + theme(plot.margin = unit(c(.5, .5, .1, .5), "cm"))
    den <- p3 + theme(plot.margin = unit(c(.1, .5, .5, .5), "cm"))
    cowplot::plot_grid(regions, den, align = 'v', ncol = 1, rel_heights = c(8,2))
    
  } else {
    print(p)
  }
  
}


#' Functions related to Notch
#' notchFilt
#'
#' Function to filter in/out events affecting a locus (in this case the Drosophila Notch locus)
#' @param infile File to process [Required]
#' @keywords parse
#' @import dplyr
#' @export
#' @return Dataframe
notchFilt <- function(..., keep=NULL, start=2700000, stop=3400000) {
  excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01","D050R03", "D050R05", "D050R07-1", "D050R07-2","D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  
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


