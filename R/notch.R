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
notchFilt <- function(..., keep=FALSE, start=2700000, stop=3400000) {
  bp_data <- getData(..., !sample %in% c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"), genotype=='somatic_tumour')
  if(keep) {
    cat("Selecting for bps in Notch\n")
    notchIn <- bp_data %>%
      filter(chrom == "X" & bp >= start & bp2 <= stop) %>%
      droplevels()
    return(notchIn)
  } else{
    cat("Excluding bps in Notch\n")
    noNotch <- bp_data %>%
      filter(!(chrom == "X" & bp >= start & bp2 <= stop)) %>%
      filter(!gene %in% c('N', 'dnc', 'kirre'), !gene2 %in% c('N', 'dnc', 'kirre'))  %>%
      droplevels()
    return(noNotch)
  }
}


#' notchHits
#' Plot the breakpint density over Notch
#' @keywords Notch
#' @import tidyverse
#' @export
notchHits <- function(infile = "data/Notch_hits.txt") {
  bp_data <- read.delim(infile, header = T)
  bp_data <- dplyr::select(bp_data, "sample", "event", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "length.Kb.", "affected_genes")
  bp_data <- dplyr::rename(bp_data, length = length.Kb.)
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
  bp_data <- dplyr::filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17")

  bp_data$type <- ifelse(bp_data$type == "BND", "INV", as.character(bp_data$type))
  bp_data <- droplevels(bp_data)

  bp_data$type <- as.factor(bp_data$type)
  bp_data$sampleax <- as.numeric(bp_data$sample)


  # cols<-setCols(bp_data, "type", set="Set1")

  if (bp_data$type == "TRA") {
    bp_data$bp1 <- bp_data$bp2 - 100
  }

  bp_data$bp1 <- bp_data$bp1 / 1000000
  bp_data$bp2 <- bp_data$bp2 / 1000000

  p <- ggplot(bp_data)
  # p<-p + geom_rect(data=bp_data, aes(xmin=bp1, xmax=bp2, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill='royalblue4'), color="black", alpha=0.6)
  p <- p + geom_rect(data = bp_data, aes(xmin = bp1, xmax = bp1 + 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)
  p <- p + geom_rect(data = bp_data, aes(xmin = bp2, xmax = bp2 - 0.001, ymin = (as.numeric(sampleax - 0.5)), ymax = (as.numeric(sampleax + 0.5)), fill = "royalblue4"), color = "black", alpha = 0.6)

  p <- p + guides(color = FALSE, size = FALSE, sampleax = FALSE, type = FALSE)

  p <- p + scale_y_continuous("Sample", expand = c(0, 0), breaks = seq(levels(as.factor(bp_data$sampleax))), labels = levels(bp_data$sample))
  p <- p + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(2.7, 3.4, by = 0.05), limits = c(2.70, 3.4))

  p <- p + slideTheme() +
    theme(
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      axis.text.y = element_text(size = 20)
    )

  p <- p + scale_fill_identity()


  p <- p + annotate("rect", xmin = 2.740000, xmax = 3.134532, ymin = -1.5, ymax = 0, alpha = .4, fill = "slategray1")
  p <- p + annotate("rect", xmin = 3.134870, xmax = 3.172221, ymin = -1.5, ymax = 0, alpha = .4, fill = "slategray3")
  p <- p + annotate("rect", xmin = 3.176440, xmax = 3.334000, ymin = -1.5, ymax = 0, alpha = .4, fill = "slategrey")
  p <- p + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)
  # p<-p + cols

  # Nhits <- paste("Notch_hits.pdf")
  # cat("Writing file", Nhits, "\n")
  # ggsave(paste("plots/", Nhits, sep=""), width = 10, height = 10)
  #

  # Density plot

  bp_data <- getData()
  bp_data <- dplyr::filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)

  bp_data$type <- ifelse(bp_data$type == "BND", "INV", as.character(bp_data$type))

  # cols2<-setCols(bp_data, "type", fill='Y', set="Set1")

  bp_data <- dplyr::filter(bp_data, type != "TRA")
  bp_data <- droplevels(bp_data)

  p2 <- ggplot(bp_data)
  p2 <- p2 + geom_density(aes(bp / 1000000, fill = "royalblue4"), alpha = 0.6)
  p2 <- p2 + scale_x_continuous("Mbs", expand = c(0, 0), breaks = seq(2.7, 3.4, by = 0.05), limits = c(2.70, 3.4))
  p2 <- p2 + scale_y_continuous("Density", expand = c(0, 0))
  p2 <- p2 + guides(colour = FALSE)
  p2 <- p2 + geom_rug(data = bp_data, aes(bp / 1000000))
  p2 <- p2 + geom_vline(xintercept = 3.135669, linetype = "dotted", size = 1)

  p2 <- p2 + slideTheme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      axis.title.y = element_blank()
      # strip.text = element_text(size=10)
    )

  # p2 <- p2 + cols2
  p2 <- p2 + scale_fill_identity()

  # p2 <- p2 + facet_wrap(~type, nrow = 3)

  combined_plots <- ggarrange(
    p, p2,
    labels = c("A", "B"),
    ncol = 1, nrow = 2
  )

  NhitsDen <- paste("Notch_hits_density.pdf")
  cat("Writing file", NhitsDen, "\n")
  ggsave(paste("plots/", NhitsDen, sep = ""), width = 30, height = 20)

  combined_plots
}


#' notchDels
#' Plot the size of deletions affecting Notch
#' @keywords Notch
#' @import dplyr ggplot2 ggsci forcats
#' @export
notchDels <- function() {
  # excludedSamples <- c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
  
  notch_data <- notchFilt(keep=1)
  bp_data <- getData()
  sample_names <- levels(bp_data$sample)
  
  # Problem with events not being clustered properly for CNV-Seq
  notch_data <- notch_data %>%
    # dplyr::filter(type != 'TRA') %>%
    dplyr::filter(bp_no == 'bp1') %>%
    group_by(sample, event) %>%
    dplyr::filter(!duplicated(event)) %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(allLength = sum(length))
  
  missing_samples = list()
  
  for(i in 1:length(sample_names)){
    if(!(sample_names[i] %in% levels(notch_data$sample))){
      missing_samples[i] <- sample_names[i]
    }
  }
  
  missing_samples <- compact(missing_samples)
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
