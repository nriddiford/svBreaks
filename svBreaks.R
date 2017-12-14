list.of.packages <- c('ggplot2', 'plyr', 'dplyr', 'RColorBrewer', 'ggpubr', 'rtracklayer')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat('Installing missing packages...\n')
  install.packages(new.packages)
}
cat('Silently loading packages...')
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(rtracklayer))


#set.seed(42)

## @knitr getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords parse 
#' @import dplyr
#' @export
#' @return Dataframe
#' 
getData <- function(infile = "data/all_bps_filtered.txt", gene_lengths_file="data/gene_lengths.txt", expression_data='data/isc_genes_rnaSeq.csv'){
  # bp_data<-read.delim(infile, header = F)
  
  ###
  
  germline<-read.delim('../svParser/germline/summary/merged/all_bps_filtered.txt')
  colnames(germline) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  germline$genotype <- 'germline'
  somatic<-read.delim('../svParser/filtered/summary/merged/all_bps_filtered.txt')
  colnames(somatic) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  somatic$genotype <- 'somatic'
  
  bp_data <- rbind(somatic,germline)
  
  ###
  
  # colnames(bp_data) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  
  gene_lengths<-read.delim(gene_lengths_file, header = T)
  
  gene_lengths<-gene_lengths[,c("gene", "id")]  
  
  bp_data <- join(bp_data,gene_lengths,"gene", type = 'left')
  
  # Read in tissue specific expression data
  seq_data<-read.csv(header = F, expression_data)
  colnames(seq_data)<-c('id', 'fpkm')
  
  bp_data <- join(bp_data,seq_data,"id", type = 'left')
  
  bp_data$fpkm <- ifelse(is.na(bp_data$fpkm), 0, round(bp_data$fpkm, 1))

  # Filter for genes expressed in RNA-Seq data
  # bp_data<-filter(bp_data, fpkm > 0.1)
  
  # Filter for genes NOT expressed in RNA-Seq data
  # cat("Filtering out expressed genes\n")
  # bp_data<-filter(bp_data, fpkm == 0)
  
  #filter on chroms
  bp_data<-filter(bp_data, chrom != "211000022280116")
  
  #filter out samples
  bp_data<-filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  
  # filter on genotype
  # bp_data<-filter(bp_data, genotype == 'somatic')
  
  # Filter for old/new data
  # bp_data <- filter(bp_data, !grepl("^A|H", sample))
  
  # Filter for SV type
  # bp_data <- filter(bp_data, type == "DEL")
  
  bp_data<-droplevels(bp_data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(bp_data)
}

## exNotch

exNotch <- function(){
  cat("Excluding bps in Notch\n")
  bp_data<-getData()
  bp_data<-filter(bp_data, !(chrom == "X" & bp >= 2700000 & bp <= 3400000))
  bp_data<-droplevels(bp_data)
  return(bp_data)
}

## cleanTheme

cleanTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=20),
    axis.title = element_text(size=30)
  )
}

slideTheme <- function(base_size = 25){
  theme(
    plot.title = element_text(hjust = 0.5, size = 50),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=30),
    axis.title = element_text(size=50),
    strip.text = element_text(size=25)
  )
}


## setCols

setCols <- function(df, col, fill='Y',set="Pastel2"){
  names<-levels(as.factor(df[[col]]))
  names<-sort(names)
  cat("Setting colour levels:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, set)
  names(mycols) <- names
  fillScale <- scale_fill_manual(name = col,values = mycols)
  colScale <- scale_colour_manual(name = col,values = mycols)
  
  if(fill == 'Y') return(fillScale)
  if(fill == 'N') return(colScale)
}

## bpStats


bpStats <- function(colSample=NA){
  bp_data<-getData()
  bp_data<-filter(bp_data, bp_no =='bp1')
  bp_data<-droplevels(bp_data)

  sampleSvs <- bp_data %>%
    group_by(sample) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    transform(sample = reorder(sample, -n))
  
  # sampleSvs <- transform(count, sample = reorder(sample, -n))
  if(!is.na(colSample)){
    sampleSvs$colour<-ifelse(sampleSvs$sample == colSample, '#FF3333', 'grey37' )
  }
  else{
    sampleSvs$colour<-'grey37'
  }
  
  p <- ggplot(sampleSvs)
  p <- p + geom_histogram(aes(sample, n, fill=colour), stat='identity')
  # p <- p + scale_y_continuous("Number of variants", limits = c(0, 20), expand = c(0.01, 0.01), breaks=seq(0,20,by=2))
  p <- p + scale_y_continuous("Number of variants", expand = c(0.01, 0.01))
  
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=20),
      axis.title.x=element_blank()
    )
  p <- p + scale_fill_identity()
  
  
  sampleSVs<-paste("SVs_sample.png")
  cat("Writing file", sampleSVs, "\n")
  ggsave(paste("plots/", sampleSVs, sep=""), width = 10, height = 10)
  
  cat("sample", "SVs", sep='\t', "\n")
  rank<-sort(table(bp_data$sample), decreasing = TRUE)
  rank<-as.array(rank)
  
  total=0
  
  scores=list()
  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
    total<-total + rank[i]
    scores[i]<-rank[i]
  }
  cat('--------------', '\n')
  scores<-unlist(scores)
  
  mean<-as.integer(mean(scores))
  med<-as.integer(median(scores))
  
  cat('total', total, sep='\t', '\n')
  cat('samples', nrow(rank), sep='\t', '\n')
  
  cat('--------------', '\n')
  cat('mean', mean, sep='\t', '\n')
  cat('median', med, sep='\t', '\n')
  
  cat('\n')
  dels<-nrow(filter(bp_data, type == "DEL" ))
  invs<-nrow(filter(bp_data, type == "INV" ))
  dups<-nrow(filter(bp_data, type == "DUP" ))
  tandups<-nrow(filter(bp_data, type == "TANDUP" ))
  tra<-nrow(filter(bp_data, type == "TRA" ))
  
  cat("Dels ", dels,  sep='', '\n')
  cat("Invs ", invs,  sep='', '\n')
  cat("Dups ", dups,  sep='', '\n')
  cat("Tandups ", tandups,  sep='', '\n')
  cat("Tra ", tra,  sep='', '\n')
  
  p
  
  
}


## bpFeatures

bpFeatures <- function(notch=0){
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  
  # To condense exon counts into "exon"
  bp_data$feature<-as.factor(gsub("_.*", "", bp_data$feature))
  
  # Reoders descending
  bp_data$feature<-factor(bp_data$feature, levels = names(sort(table(bp_data$feature), decreasing = TRUE)))
  
  cols<-setCols(bp_data, "feature")
  
  p<-ggplot(bp_data)
  p<-p + geom_bar(aes(feature, fill = genotype, group=reverse(genotype)),stat="count", position=position_dodge())
  #p<-p + cols
  p<-p + slideTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  p<-p + scale_fill_brewer(palette="Paired")
  # p <- p + facet_wrap(~genotype)
  
  features_outfile<-paste("Breakpoints_features_count", ext, sep = "")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)
  
  p
}



sizeDist <- function(){
  bp_data<-getData()
  
  bp_data<-filter(bp_data, type != "TRA", type != "BND", bp_no != "bp2")
  
  bp_data$length <- ifelse(bp_data$length == 0, 0.01, bp_data$length)
  bp_data$length <- (bp_data$length*1000)
  bp_data<-droplevels(bp_data)
  fillCols<-setCols(bp_data, "type", fill='Y')
  cols<-setCols(bp_data, "type", fill='N')
  
  
  bp_data <- transform(bp_data, type = reorder(type, length))
  
  p <- ggplot(bp_data, aes(type, length))
  p <- p + geom_violin(aes(fill=type))
  p <- p + geom_jitter(width=.1, height=.1)
  
  p <- p + scale_y_log10("Length (Bp)") 
  p <- p + slideTheme() +
    theme(
      axis.title.x=element_blank(),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + facet_wrap(~genotype)
  p <- p + fillCols
  p <- p + cols 
  
  
  sizedistOut<-paste("sizeDist.png")
  cat("Writing file", sizedistOut, "\n")
  ggsave(paste("plots/", sizedistOut, sep=""), width = 15, height = 10)
  
  p
  
}



## notchHits


notchHits <- function(infile = "data/Notch_hits.txt"){
  bp_data<-read.delim(infile, header =T)
  bp_data<-dplyr::select(bp_data, "sample", "event", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2","length.Kb.","affected_genes")
  bp_data<-dplyr::rename(bp_data, length = length.Kb.) 
  # refGene<-read.delim(refgene_file, header=F)
  # colnames(refGene)<-c('Tno',	'id',	'chrom',	'strand',	'Tstart',	'Tstop',	'cdsstart',	'cdsStop',	'exonNo',	'Estart',	'Estop',	'score',	'altname',	'CDSstartStat',	'CDSEndStat',	'exonFrames')
  # # Tid	id	chrom	strand	Tstart	Tstop	cdsstart	cdsStop	exon#	Estart	Estop	score	altname	CDSstartStat	CDSEndStat	exonFrames
  # gene='N'
  # refGene<-filter(refGene, altname == gene)
  # refGene<-droplevels(refGene)
  # 
  # wStart<-refGene$Tstart - 1000
  # wEnd<-refGene$stop + 1000
  
  # bp_data<-getData()

  # bp_data<-reshape(bp_data, idvar=c("event", "sample"), timevar="bp_no", direction="wide")
  # #bp_data<-filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)
  # bp_data<-dplyr::select(bp_data, sample, chrom.bp1, bp.bp1, gene.bp1, bp.bp2,gene.bp2, type.bp1)
  # colnames(bp_data)<-c("sample", "chrom1", "bp1", "gene1", "bp2", "gene2", "type")
  bp_data<-filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )

  bp_data$type <- ifelse(bp_data$type=="BND", "INV", as.character(bp_data$type))
  bp_data<-droplevels(bp_data)
  
  bp_data$type <- as.factor(bp_data$type)
  bp_data$sampleax <- as.numeric(bp_data$sample)
  

  # cols<-setCols(bp_data, "type", set="Set1")  
  
  if(bp_data$type == "TRA"){
    bp_data$bp1 = bp_data$bp2-100
  }
  
  bp_data$bp1<-bp_data$bp1/1000000
  bp_data$bp2<-bp_data$bp2/1000000
  
  p<-ggplot(bp_data)
  # p<-p + geom_rect(data=bp_data, aes(xmin=bp1, xmax=bp2, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill='royalblue4'), color="black", alpha=0.6)
  p<-p + geom_rect(data=bp_data, aes(xmin=bp1, xmax=bp1+0.001, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill='royalblue4'), color="black", alpha=0.6)
  p<-p + geom_rect(data=bp_data, aes(xmin=bp2, xmax=bp2-0.001, ymin=(as.numeric(sampleax-0.5)),ymax=(as.numeric(sampleax+0.5)),fill='royalblue4'), color="black", alpha=0.6)
  
  p<-p + guides(color = FALSE, size = FALSE, sampleax = FALSE, type=FALSE)
  
  p<-p + scale_y_continuous("Sample",expand = c(0,0), breaks = seq(levels(as.factor(bp_data$sampleax))), labels=levels(bp_data$sample))
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(2.7,3.4,by=0.05), limits=c(2.70, 3.4))
  
  p<-p + slideTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.position="top",
          axis.text.y = element_text(size = 20)
    )
  
  p<-p + scale_fill_identity()

  
  p<-p + annotate("rect", xmin=2.740000, xmax=3.134532, ymin=-1.5, ymax=0, alpha=.4, fill="slategray1")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=-1.5, ymax=0, alpha=.4, fill="slategray3")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.334000, ymin=-1.5, ymax=0, alpha=.4, fill="slategrey")
  p<-p + geom_vline(xintercept = 3.135669, linetype="dotted", size = 1)
  # p<-p + cols
  
  # Nhits <- paste("Notch_hits.pdf")
  # cat("Writing file", Nhits, "\n")
  # ggsave(paste("plots/", Nhits, sep=""), width = 10, height = 10)
  # 
  
  # Density plot

  bp_data<-getData()
  bp_data<-filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)
  
  bp_data$type <- ifelse(bp_data$type=="BND", "INV", as.character(bp_data$type))
  
  # cols2<-setCols(bp_data, "type", fill='Y', set="Set1")  
  
  bp_data <- filter(bp_data, type != 'TRA')
  bp_data<-droplevels(bp_data)
  
  p2<-ggplot(bp_data)
  p2<-p2 + geom_density(aes(bp/1000000, fill='royalblue4'), alpha = 0.6)
  p2<-p2 + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(2.7,3.4,by=0.05), limits=c(2.70, 3.4))
  p2<-p2 + scale_y_continuous("Density", expand = c(0,0))
  p2<-p2 + guides(colour = FALSE)
  p2<-p2 + geom_rug(data=bp_data,aes(bp/1000000))
  p2<-p2 + geom_vline(xintercept = 3.135669, linetype="dotted",size = 1)
  
  p2<-p2 + slideTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position="top",
          axis.title.y=element_blank()
          # strip.text = element_text(size=10)
          )
  
  # p2 <- p2 + cols2
  p2<-p2 + scale_fill_identity()
  
  # p2 <- p2 + facet_wrap(~type, nrow = 3)
  
  combined_plots <- ggarrange(p, p2, 
                              labels = c("A", "B"),
                              ncol = 1, nrow = 2)
  
  NhitsDen <- paste("Notch_hits_density.pdf")
  cat("Writing file", NhitsDen, "\n")
  ggsave(paste("plots/", NhitsDen, sep=""), width = 30, height = 20)
  
  combined_plots
}

## notchDels

notchDels <- function(infile = "data/Notch_hits.txt"){
  bp_data<-read.delim(infile, header =T)
  bp_data<-dplyr::select(bp_data, "sample", "event", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2","length.Kb.","affected_genes")
  bp_data<-dplyr::rename(bp_data, length = length.Kb.) 
  
  # bp_data<-filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  
  
  bp_data<-bp_data %>%
    group_by(sample) %>%
    slice(which.max(length))
  
  bp_data<-droplevels(bp_data)
  
  bp_data$colour<-ifelse( bp_data$sample == "A373R1" | bp_data$sample ==  "A373R7" | bp_data$sample == "A512R17", '#FF3333', 'gray37' )
  
  bp_data <- transform(bp_data, sample = reorder(sample, -length))
  
  p<-ggplot(bp_data)
  p<-p + geom_bar(aes(sample,length, fill=colour),stat="identity")
  #p<-p + scale_y_discrete(expand = c(0.01,0.01), breaks=seq(0,500,by=100))
  p<-p + slideTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1)
    )

  p<-p + scale_fill_identity()
  dels_out<-paste("NotchDels.pdf")
  cat("Writing file", dels_out, "\n")
  ggsave(paste("plots/", dels_out, sep=""), width = 30, height = 15)
  
  p
}


## bpFeatureEnrichment

bpFeatureEnrichment <- function(features='data/genomic_features.txt', genome_length=118274340, print=NA){
  genome_features<-read.delim(features, header = T)
  bp_data<-getData()
  
  mutCount<-nrow(bp_data)
  
  # To condense exon counts into "exon"
  bp_data$feature<-as.factor(gsub("exon_.*", "exon", bp_data$feature))
  
  classCount<-table(bp_data$feature)
  classLengths<-setNames(as.list(genome_features$length), genome_features$feature)
  
  fun <- function(f) {
    # Calculate the fraction of geneome occupied by each feature
    featureFraction<-classLengths[[f]]/genome_length
    
    # How many times should we expect to see this feature hit in our bp_data (given number of obs. and fraction)?
    featureExpect<-(mutCount*featureFraction)
    
    # observed/expected 
    fc<-classCount[[f]]/featureExpect
    fc<-round(fc,digits=1)
    featureExpect<-round(featureExpect,digits=3)
    
    # Binomial test
    if(!is.null(classLengths[[f]])){
      if(classCount[f] >= featureExpect){
        stat<-binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "greater")
        test<-"enrichment"
      }
      else{
        stat<-binom.test(x = classCount[f], n = mutCount, p = featureFraction, alternative = "less")
        test<-"depletion"
      }
      sig_val<-'F'
      if(stat$p.value <= 0.05){ sig_val<-'T'}
      p_val<-format.pval(stat$p.value, digits = 3, eps=0.0001)
      Log2FC<-log2(fc)
      # Log2FC<-round(Log2FC, 1)
      list(feature = f, observed = classCount[f], expected = featureExpect, Log2FC = Log2FC, test = test, sig = sig_val, p_val = p_val)
    }
  }
  
  enriched<-lapply(levels(bp_data$feature), fun)
  enriched<-do.call(rbind, enriched)
  featuresFC<-as.data.frame(enriched)
  # Sort by FC value
  featuresFC<-dplyr::arrange(featuresFC,desc(abs(as.numeric(Log2FC))))
  featuresFC$Log2FC<-round(as.numeric(featuresFC$Log2FC), 1)
  featuresFC$expected<-round(as.numeric(featuresFC$expected), 1)
  
  if(!is.na(print)){
    cat("printing")
    first.step <- lapply(featuresFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)

    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))
    
    feat_enrichment_table <- paste("feature_enrichment_table.tiff")
    # ggexport(filename = paste("plots/", "ex_", feat_enrichment_table, sep=""),
    #          width = 480, height = 480, pointsize = 12, res = 250,
    #          verbose = TRUE)
    ggsave(paste("plots/", feat_enrichment_table, sep=""), width = 5.2, height = (nrow(featuresFC)/3), dpi=300)
  }
  
  else{ return(featuresFC) }
}

bpFeatureEnrichmentPlot <- function() {
  feature_enrichment<-bpFeatureEnrichment()
  
  
  #feature_enrichment$Log2FC <- log2(as.numeric(feature_enrichment$fc))
  
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  #feature_enrichment$fc <- as.numeric(feature_enrichment$fc)
  
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -Log2FC))
  
  feature_enrichment <- filter(feature_enrichment, observed >= 3)
  feature_enrichment <- droplevels(feature_enrichment)
  
  p<-ggplot(feature_enrichment)
  p<-p + geom_bar(aes(feature, Log2FC, fill = as.character(test)), stat="identity")
  p<-p + guides(fill=FALSE)
  p<-p + ylim(-2,2)
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1)
    )
  
  feat_plot <- paste("feat_plot.pdf")
  cat("Writing file", feat_plot, "\n")
  ggsave(paste("plots/", feat_plot, sep=""), width = 5, height = 10)
  p
  
}
## bpGeneEnrichment

bpGeneEnrichment <- function(gene_lengths="data/gene_lengths.txt", n=3, genome_length=118274340, print=NA){
  
  cat("Showing genes hit at least", n, "times", "\n")
  gene_lengths<-read.delim(gene_lengths, header = T)
  # bp_data<-read.delim('data/all_samples.txt',header=T)
  bp_data<-getData()
  bp_data<-filter(bp_data, gene != "intergenic")
  
  bp_count<-nrow(bp_data)
  
  hit_genes<-table(bp_data$gene)
  
  # hit_genes<-table(bp_data[!duplicated(bp_data),]$gene)
  
  genes<-setNames(as.list(gene_lengths$length), gene_lengths$gene)
  expressed_genes<-setNames(as.list(bp_data$fpkm), bp_data$gene)
  
  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction<-genes[[g]]/genome_length
    
    # How many times should we expect to see this gene hit in our bp_data (given number of obs. and fraction)?
    gene_expect<-bp_count*(genefraction)
    
    # observed/expected 
    fc<-hit_genes[[g]]/gene_expect
    fc<-round(fc,digits=1)
    log2FC = log2(fc)
    
    gene_expect<-round(gene_expect,digits=3)
    list(gene = g, length = genes[[g]], fpkm = expressed_genes[[g]], observed = hit_genes[g], expected = gene_expect, fc = fc, log2FC = log2FC)
  }
  
  enriched<-lapply(levels(bp_data$gene), fun)
  enriched<-do.call(rbind, enriched)
  genesFC<-as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC<-filter(genesFC, observed >= n)
  # Sort by FC value
  genesFC<-dplyr::arrange(genesFC,desc(as.integer(fc)))
  genesFC$expected<-round(as.numeric(genesFC$expected),digits=2)
  genesFC$log2FC<-round(as.numeric(genesFC$log2FC),digits=1)
  
  if(!is.na(print)){
    cat("printing")
    first.step <- lapply(genesFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)
    
    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))
    
    gene_enrichment_table <- paste("gene_enrichment_table.tiff")
    ggsave(paste("plots/", gene_enrichment_table, sep=""), width = 5, height = (nrow(genesFC)/3),dpi=300)
  }
  
  else{ return(genesFC) }
}

## bpGeneEnrichmentPlot

bpGeneEnrichmentPlot <- function(n=2) {
  gene_enrichment<-bpGeneEnrichment(n=n)
  
  gene_enrichment$Log2FC <- log2(as.numeric(gene_enrichment$fc))
  
  gene_enrichment$gene <- as.character(gene_enrichment$gene)
  gene_enrichment$fc <- as.numeric(gene_enrichment$fc)
  
  gene_enrichment <- transform(gene_enrichment, gene = reorder(gene, -fc))
  
  gene_enrichment$test <- ifelse(gene_enrichment$Log2FC>=0, "enriched", "depleted")
  
  gene_enrichment <- filter(gene_enrichment, observed >= 2)
  gene_enrichment<-droplevels(gene_enrichment)
  
  # highlightedGene <- filter(gene_enrichment, gene == "N")
  # highlightedGene <- droplevels(highlightedGene)
  
  p<-ggplot(gene_enrichment)
  p<-p + geom_bar(aes(gene, Log2FC, fill = as.character(test)), stat="identity")
 # p<-p + geom_bar(data=highlightedGene, aes(gene, Log2FC, fill="red"), colour="black", stat="identity")
  p<-p + guides(fill=FALSE)
  p<-p + slideTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=40), axis.title = element_text(size=90)
          
    )
  
  gene_enrichment_plot <- paste("gene_enrichment.pdf")
  cat("Writing file", gene_enrichment_plot, "\n")
  ggsave(paste("plots/", gene_enrichment_plot, sep=""), width = 30, height = 10)
  p
  
}


getPromoter <- function(gene_lengths_in="data/gene_lengths.txt"){
  gene_lengths<-read.delim(gene_lengths_in, header = T)
  
  gene_lengths$promoter<-ifelse(gene_lengths$start<gene_lengths$end,
                                gene_lengths$start- 1500,
                                gene_lengths$end + 1500)
  
  
  gene_lengths<-gene_lengths[,c("chrom", "promoter")]
  colnames(gene_lengths)<-NULL
  return(gene_lengths)
}


### for all genes affected by SVs ...
bpAllGenes <- function(gene_lengths="data/gene_lengths.txt", n=3, genome_length=118274340){
  cat("Showing genes affected by a SV at least", n, "times", "\n")
  gene_lengths<-read.delim(gene_lengths, header = T)
  allGenes<-read.delim('data/all_genes_filtered.txt',header=F)
  colnames(allGenes) <- c("event", "sample", "type", "chrom", "gene")

  allGenes<-filter(allGenes, sample != "A373R1" & sample != "A373R7" & sample != "A512R17", gene != '-')


  geneCount<-nrow(allGenes)
  hit_genes<-table(allGenes$gene)

  genes  <- setNames(as.list(gene_lengths$length), gene_lengths$gene)
  chroms <- setNames(as.list(gene_lengths$chrom), gene_lengths$gene)
  
  # expressed_genes<-setNames(as.list(bp_data$fpkm), bp_data$gene)
  
  
  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction<-genes[[g]]/genome_length
    
    # How many times should we expect to see this gene hit in our bp_data (given number of obs. and fraction)?
    gene_expect<-geneCount*(genefraction)
    
    # observed/expected 
    fc<-hit_genes[[g]]/gene_expect
    fc<-round(fc,digits=1)
    log2FC = log2(fc)
    
    gene_expect<-round(gene_expect,digits=3)
    list(gene = g, length = genes[[g]], chromosome=as.character(chroms[[g]]), observed = hit_genes[g], expected = gene_expect, fc = fc, log2FC = log2FC)
  }
  
  enriched<-lapply(levels(allGenes$gene), fun)
  enriched<-do.call(rbind, enriched)
  genesFC<-as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC<-filter(genesFC, observed >= n)
  # Sort by FC value
  genesFC<-dplyr::arrange(genesFC,desc(as.integer(observed)))
  genesFC$expected<-round(as.numeric(genesFC$expected),digits=2)
  genesFC$log2FC<-round(as.numeric(genesFC$log2FC),digits=1)
  
  
  return(genesFC)
  
}

dist2Motif <- function(feature_file="data/tss_locations.txt",sim=NA, print=0,send=0, feature='tss'){
  if(is.na(sim)){
    bp_data<-getData()
    bp_data<-dplyr::rename(bp_data, pos = bp) 
  }
  
  else{
    cat("Generating simulated bp_data\n")
    hit_count<-nrow(getData())
    bp_data<-bpSim(N=hit_count, write=print)
    colnames(bp_data)<-c("chrom", "pos", "v3", "v4", "v5")
    bp_data<-filter(bp_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    bp_data<-droplevels(bp_data)
  }
  
  feature<-paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep='')
  
  if(feature=='Promoter'){
    feature_locations<-getPromoter()
    cat("Getting gene promoter locations...\n")
  }
  else{ 
    feature_locations<-read.delim(feature_file, header = F)
    cat("Reading in file:", feature_file, sep =' ', "\n")
  }
  
  cat("Calculating distances to", feature, sep=' ', "\n")

  feature_locations<-read.delim(feature_file, header = F)
  colnames(feature_locations)<-c("chrom", "pos")
  
  feature_locations$pos<-as.integer(feature_locations$pos)
  
  # Will throw error if SVs don't exist on a chrom...
  
  # Removes chroms with fewer than 10 observations
  svCount <- table(bp_data$chrom)
  bp_data <- subset(bp_data, chrom %in% names(svCount[svCount >= 10]))
  bp_data<-droplevels(bp_data)
  
  feature_locations <- subset(feature_locations, chrom %in% levels(bp_data$chrom))
  feature_locations<-droplevels(feature_locations)  
  
  
  fun2 <- function(p) {
    index<-which.min(abs(tss_df$pos - p))
    closestTss<-tss_df$pos[index]
    chrom<-as.character(tss_df$chrom[index])
    gene<-as.character(tss_df$gene[index])
    dist<-(p-closestTss)
    list(p, closestTss, dist, chrom, gene)
  }
  
  l <- list()
  
  for (c in levels(bp_data$chrom)){
    df<-filter(bp_data, chrom == c)
    tss_df<-filter(feature_locations, chrom == c)
    dist2tss<-lapply(df$pos, fun2)
    dist2tss<-do.call(rbind, dist2tss)
    dist2tss<-as.data.frame(dist2tss)
    
    colnames(dist2tss)=c("bp", "closest_tss", "min_dist", "chrom", "closest_gene")
    dist2tss$min_dist<-as.numeric(dist2tss$min_dist)
    l[[c]] <- dist2tss
  }
  
  dist2tss<-do.call(rbind, l)
  dist2tss<-as.data.frame(dist2tss)
  dist2tss$chrom<-as.character(dist2tss$chrom)
  
  dist2tss<-arrange(dist2tss,(abs(min_dist)))
  
  if(send==1){
    return(dist2tss)
  }
  else{
    p<-ggplot(dist2tss)
    p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
    p<-p + scale_x_continuous(paste("Distance to", feature, "(Kb)", sep=' '),
    limits=c(-10000, 10000),
    breaks=c(-10000,-1000, 1000, 10000),
    expand = c(.0005, .0005),
    labels=c("-10", "-1", "1", "10") )
    
    p<-p + scale_y_continuous("Density")
    p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
    #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
    p <- p + geom_rug(aes(min_dist, colour=chrom))
    p<-p + slideTheme() +
      theme(strip.text = element_text(size=20),
          legend.position="top")
    
    p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")

    if(is.na(sim)){
      distout<-paste("bp", feature, 'dist.pdf', sep='')
    }
    else{
      distout<-paste("bp", feature, 'dist_sim.pdf', sep='')
    }
    
    cat("Writing file", distout, "\n")
    ggsave(paste("plots/", distout, sep=""), width = 20, height = 10)

  p
  }
}


distOverlay <- function(feature_file="data/tss_locations.txt", feature='tss', lim=10, all=NA){
  feature<-paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep='')

  if(feature=='promoter'){
    real_data<-dist2Motif(send=1, feature=feature)
    sim_data<-dist2Motif(feature=feature, sim=1, send=1)
  }
  else{ 
    real_data<-dist2Motif(feature_file=feature_file, send=1, feature=feature)
    sim_data<-dist2Motif(feature_file=feature_file, feature=feature, sim=1, send=1)
  }

  real_data$Source<-"Real"
  sim_data$Source<-"Sim"
  
  sim_data<-filter(sim_data, chrom != "Y", chrom != 4)
  sim_data<-droplevels(sim_data)
  real_data<-filter(real_data, chrom != "Y", chrom != 4)
  real_data<-droplevels(real_data)
  
  colours<-c( "#E7B800", "#00AFBB")

  scale<-"(Kb)"
  if(lim==0.1){
    cat("Setting limits to -+100bp\n")
    lims=c(-100, 100)
    brks=c(-100, -10, 10, 100)
    expnd = c(.0005, .0005)
    labs=c("-100", "-10", "10", "100")
    scale<-"(bp)"
  } 
  
  else if(lim==0.5){
    cat("Setting limits to -+0.5kb\n")
    lims=c(-500, 500)
    brks=c(-500, -100,100, 500)
    expnd = c(.0005, .0005)
    labs=c("-500", "-100", "100", "500")
    scale<-"(bp)"
  }  
  
  else if(lim==1){
    cat("Setting limits to -+1kb\n")
    lims=c(-1000, 1000)
    brks=c(-1000, 1000)
    expnd = c(.0005, .0005)
    labs=c("-1", "1")
  }  
  else{
    cat("Setting limits to -+10kb\n")
    lims=c(-10000, 10000)
    brks=c(-10000,-1000, 1000, 10000)
    expnd = c(.0005, .0005)
    labs=c("-10", "-1", "1", "10")
  }
  
  p<-ggplot()
  p<-p + geom_density(data=real_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + geom_density(data=sim_data,aes(min_dist, fill = Source), alpha = 0.4)
  if(is.na(all)){
    p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  }

  p<-p + scale_x_continuous(paste("Distance to", feature, scale, sep=' '),
                            limits=lims,
                            breaks=brks,
                            expand=expnd,
                            labels=labs )
  p<-p + scale_y_continuous("Density")
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")

  p <- p + geom_rug(data=real_data,aes(min_dist, colour=Source),sides="b")
  p <- p + geom_rug(data=sim_data,aes(min_dist, colour=Source),sides="t")
  
  p <- p + scale_fill_manual(values=colours)
  p <- p + scale_colour_manual(values=colours)
  
  p<-p + slideTheme() +
    theme(strip.text = element_text(size=20),
          legend.position="top")
  
  overlay<-paste("bp", feature, 'dist_overlay.pdf', sep='')
  cat("Writing file", overlay, "\n")
  ggsave(paste("plots/", overlay, sep=""), width = 20, height = 10)
  
  p
}


featureDensity <- function(){
  # tss_positions<-read.delim("data/tss_locations.txt", header=F)
  # tss_positions$type<-"TSS"
  # colnames(tss_positions)<-c("chrom", "pos", "type")
  g4_positions<-read.delim("data/g4_positions.txt", header=F)
  g4_positions$type<-"G4"
  colnames(g4_positions)<-c("chrom", "pos", "type")
  invR_positions<-read.delim("data/invRepeats.txt", header=F)
  invR_positions$type<-"SIR"
  colnames(invR_positions)<-c("chrom", "pos", "type")
  cru_positions<-read.delim("data/cruciform_positions.txt", header=F)
  cru_positions$type<-"Cru"
  colnames(cru_positions)<-c("chrom", "pos", "type")
  
  locations<-rbind.data.frame(g4_positions,invR_positions,cru_positions)
  
  locations$type<-as.factor(locations$type)
  locations$pos<-as.numeric(locations$pos/1000000)
  locations<-filter(locations, chrom != 'Y', chrom != 4)
  locations<-droplevels(locations)
  
  p<-ggplot(locations)
  p<-p + geom_density(aes(pos, fill=type),alpha = 0.4)
  p <- p + geom_rug(aes(pos, colour=type),sides="b", alpha = 0.05)
  p<-p + facet_wrap(type~chrom, scale = "free_x", ncol = 5)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,max(locations$pos),by=10))
  
  nonBDNA<-paste("nonBform.pdf")
  cat("Writing file", nonBDNA, "\n")
  ggsave(paste("plots/", nonBDNA, sep=""), width = 20, height = 10)
  
  p
  
}







bpinGene <- function(gene_lengths="data/gene_lengths.txt", gene2plot='dnc'){
  gene_lengths<-read.delim(gene_lengths, header = T)
  region<-filter(gene_lengths, gene == gene2plot)
  
  wStart<-(region$start - 10000)
  wEnd<-(region$end + 10000)
  wChrom<-as.character(region$chrom)
  wTss<-suppressWarnings(as.numeric(levels(region$tss))[region$tss])
  bp_data<-getData()
  bp_data<-filter(bp_data, chrom == wChrom & bp >= wStart & bp <= wEnd)
  
  if(nrow(bp_data) == 0){
    stop(paste("There are no svs in", gene2plot, "- Exiting", "\n"))
  }
  
  p<-ggplot(bp_data)
  p<-p + geom_point(aes(bp/1000000, sample, colour = type, size = 1.5), position=position_jitter(width=0.005, height=0))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(round(wStart/1000000, digits = 2),round(wEnd/1000000, digits = 2),by=0.05), limits=c(wStart/1000000, wEnd/1000000))
  p<-p + annotate("rect", xmin=region$start/1000000, xmax=region$end/1000000, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")
  p<-p + geom_vline(xintercept = wTss/1000000, colour="red", alpha=.7, linetype="solid")
  
  p<-p + geom_segment(aes(x = wTss/1000000, y = 0, xend= wTss/1000000, yend = 0.1), colour="red")
  middle<-((wEnd/1000000+wStart/1000000)/2)
  p <- p + annotate("text", x = middle, y = 0.05, label=gene2plot, size=6)
  p<-p + ggtitle(paste("Chromosome:", wChrom))
  
  p
}


#' bpRainfall
#'
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @keywords rainfall
#' @export

bpRainfall <- function(){
  bp_data<-getData()
  
  #bp_data<-filter(bp_data, sample != "A373R11" & sample != 'A373R13')
  distances<-do.call(rbind, lapply(split(bp_data[order(bp_data$chrom, bp_data$bp),], bp_data$chrom[order(bp_data$chrom, bp_data$bp)]),
                                   function(a) 
                                     data.frame(a,
                                                dist=c(diff(a$bp), NA),
                                                logdist = c(log10(diff(a$bp)), NA))
  )
  )
  
  
  distances$logdist[is.infinite(distances$logdist)] <- 0
  distances<-filter(distances, chrom != 4, chrom != "Y")
  
  p<-ggplot(distances)
  p<-p + geom_point(aes(bp/1000000, logdist, colour = sample))
  p <- p + cleanTheme() +
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          strip.text = element_text(size=20)
    )
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,max(distances$bp),by=10))
  p<-p + scale_y_continuous("Genomic Distance")
  
  rainfall_out<-paste("rainfall.pdf")
  cat("Writing file", rainfall_out, "\n")
  ggsave(paste("plots/", rainfall_out, sep=""), width = 15, height = 5)
  
  p
}




#' bpSim
#'
#' Generate simulated SV breakpoints acroos genomic regions (e.g. mappable regions)
#' @param intervals File containing genomic regions within which to simulate SNVs [Default 'data/intervals.bed]
#' @param N Number of random breakpoints to generate [Default nrow(bp_data)]
#' @import GenomicRanges
#' @keywords sim
#' @export

bpSim <- function(intervals="data/intervals.bed", N=1000, write=F){
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(rtracklayer))
  
  intFile <- import.bed(intervals)
  space <- sum(width(intFile))
  positions <- sample(c(1:space), N)
  cat("Simulating", N, "breakpoints", sep = ' ', "\n")
  new_b <- GRanges(seqnames=as.character(rep(seqnames(intFile), width(intFile))),
                   ranges=IRanges(start=unlist(mapply(seq, from=start(intFile), to=end(intFile))), width=1))
  bedOut<-new_b[positions]
  if(write){
    export.bed(new_b[positions], "data/simulatedBPs.bed")
  }
  remove(new_b)
  return(data.frame(bedOut))
}






## svTypes

svTypes<-function(notch=0,object=NA){
  if(is.na(object)){
    object<-'type'
  }
  
  if(notch){
    bp_data<-getData()
    bp_data<-filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)
    
    ext<-'_Notch.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  
  bp_data$type <- ifelse(bp_data$type=="BND", "INV", as.character(bp_data$type))
  bp_data<-droplevels(bp_data)
  
  cols<-setCols(bp_data, "type")
  
  # Reorder by count
  bp_data$type<-factor(bp_data$type, levels = names(sort(table(bp_data$type), decreasing = TRUE)))
  
  if(object == 'sample'){
    # Reorder by count
    bp_data$sample<-factor(bp_data$sample, levels = names(sort(table(bp_data$sample), decreasing = TRUE)))
  }
  
  # Only take bp1 for each event
  bp_data<-filter(bp_data, bp_no != "bp2")
  bp_data<-droplevels(bp_data)
  
  p<-ggplot(bp_data)
  p<-p + geom_bar(aes(get(object), fill=type))
  p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=40), axis.title = element_text(size=90)
    )
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous("Number of calls", expand = c(0.01, 0.01))
  p<-p + facet_wrap(~genotype)
  # p<-p + coord_flip()
  # p<-p + scale_y_reverse()
  
  types_outfile<-paste("sv_types_by_", object, ext, sep = "")
  cat("Writing file", types_outfile, "\n")
  ggsave(paste("plots/", types_outfile, sep=""), width = 22, height = 22)
  
  p
}



























typeLen <- function(size_threshold = 1, notch=0){
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  
  cols<-setCols(bp_data, "type")
  
  # Only take bp1 for each event
  bp_data<-filter(bp_data, type != "TRA", type != "BND", bp_no != "bp2")
  
  bp_data$length<-(bp_data$length/1000)
  
  if(is.na(size_threshold)){
    size_threshold<-max(bp_data$length)
  }
  
  if(size_threshold <= 1){
    breaks<-0.1
  }
  else{
    breaks<-1
  }
  
  p<-ggplot(bp_data, aes(length))
  p<-p + geom_density(aes(fill = type), alpha = 0.4)
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_continuous("Size in Mb", expand = c(0,0), breaks = seq(0,size_threshold,by=breaks), limits=c(0, (size_threshold+0.1)))
  p<-p + scale_y_continuous(expand = c(0,0))
  p<-p + cols
  
  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep=""), width = 20, height = 10)
  
  p
}


genomeHits <- function(notch=0){
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  
  p<-ggplot(bp_data)
  p<-p + geom_point(aes(bp/1000000, sample, colour = sample, shape = type, size = 0.5), alpha = 0.7)
  p<-p + guides(color = FALSE, size = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  
  sv_gen_dist <- paste("bp_gen.dist", ext, sep = "")
  cat("Writing file", sv_gen_dist, "\n")
  ggsave(paste("plots/", sv_gen_dist, sep=""), width = 20, height = 10)
  
  p
}


typeLenCount <- function(size_threshold = 1, notch=0){
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  
  cols<-setCols(bp_data, "type")
  
  # Only take bp1 for each event
  bp_data<-filter(bp_data, type != "TRA", type != "BND", bp_no != "bp2")
  
  bp_data$length<-(bp_data$length/1000)
  
  if(is.na(size_threshold)){
    size_threshold<-max(bp_data$length)
  }
  
  if(size_threshold <= 1){
    breaks<-0.1
  }
  else{
    breaks<-1
  }
  
  p<-ggplot(bp_data, aes(length))
  p<-p + geom_histogram(aes(length, ..count.., fill = type), colour = "black", binwidth = 0.05, position="dodge")
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_continuous("Size in Mb", expand = c(0,0), breaks = seq(0,size_threshold,by=breaks), limits=c(0, (size_threshold+0.1)))
  p<-p + scale_y_continuous(expand = c(0,0))
  p<-p + geom_density(aes(fill = type),alpha=0.4, colour = NA)
  p<-p + cols
  
  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep=""), width = 20, height = 10)
  
  p
}


bpGenAll <- function(object=NA, notch=0){
  bp_data<-getData()
  ext<-'.pdf'
  if(is.na(object)){
    object<-'type'
    cols<-setCols(bp_data, "type")
  }
  
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  
  cat("Plotting SVs by", object, "\n")
  
  p<-ggplot(bp_data)
  p<-p + geom_histogram(aes(bp/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  p<-p + scale_y_continuous("Number of Breakpoints", expand = c(0.01, 0.01))
  p<-p + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )
  
  if (object == 'type'){
    p<-p + cols
  }
  
  chrom_outfile<-paste("Breakpoints_chroms_by_", object, ext, sep = "")
  cat("Writing file", chrom_outfile, "\n")
  ggsave(paste("plots/", chrom_outfile, sep=""), width = 20, height = 10)
  
  p
}


bpChromDist <- function(object=NA, notch=0){
  bp_data<-getData()
  ext<-'.pdf'
  
  if(is.na(object)){
    object<-'type'
    cols<-setCols(bp_data, "type")
  }
  
  if(notch){
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  
  chromosomes<-c("2L", "2R", "3L", "3R", "X", "Y", "4")
  lengths<-c(23513712, 25286936, 28110227, 32079331, 23542271, 3667352, 1348131)
  
  karyotype<-setNames(as.list(lengths), chromosomes)
  
  for (c in chromosomes) {
    
    len<-karyotype[[c]]
    len<-len/1000000
    
    cat("Chrom", c, "length:", len, sep = " ",  "\n")
    
    per_chrom<-filter(bp_data, chrom == c)
    
    p<-ggplot(per_chrom)
    p<-p + geom_histogram(aes(bp/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)
    p<-p + scale_x_continuous("Mbs", breaks = seq(0,len,by=1), limits = c(0, len+0.1),expand = c(0.01, 0.01))
    p<-p + scale_y_continuous("Number of Breakpoints", limits = c(0, 35), expand = c(0.01, 0.01))
    p<-p + cleanTheme() +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            axis.title = element_text(size=20)
      )
    p<-p + ggtitle(paste("Chromosome: ", c))
    
    if (object == 'type'){
      p<-p + cols
    }
    
    per_chrom<-paste("Breakpoints_on_", c, "_by_", object, ext, sep = "")
    cat("Writing file", per_chrom, "\n")
    ggsave(paste("plots/", per_chrom, sep=""), width = 20, height = 10)
  }
}

