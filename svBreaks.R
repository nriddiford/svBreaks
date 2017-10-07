list.of.packages <- c('ggplot2', 'dplyr', 'plyr', 'RColorBrewer', 'ggpubr')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat('Installing missing packages...\n')
  install.packages(new.packages)
}
cat('Silently loading packages...')
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))

set.seed(42)


getData <- function(infile = "data/all_bps_filtered.txt", gene_lengths_file="data/gene_lengths.txt", expression_data='data/isc_genes_rnaSeq.csv'){
  bp_data<-read.delim(infile, header = F)
  colnames(bp_data) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  
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
  # Filter for non expressed genes in RNA-Seq data
  bp_data<-filter(bp_data, fpkm == 0)
  
  #filter on chroms
  #bp_data<-filter(bp_data, chrom != "Y" & chrom != "4")
  
  #filter out samples
  bp_data<-filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  
  # Filter for old/new data
  # bp_data <- filter(bp_data, !grepl("^A|H", sample))
  
  # Filter for SV type
  # bp_data <- filter(bp_data, type == "DEL")
  
  bp_data<-droplevels(bp_data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(bp_data)
}


exNotch <- function(x){
  cat("Excluding bps in Notch\n")
  bp_data<-getData()
  bp_data<-filter(bp_data, !(chrom == "X" & bp >= 3000000 & bp <= 3300000))
  bp_data<-droplevels(bp_data)
  return(bp_data)
}


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


setCols <- function(df, col, fill='Y'){
  names<-levels(as.factor(df[[col]]))
  names<-sort(names)
  cat("Setting colour levels:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, "Set2")
  names(mycols) <- names
  fillScale <- scale_fill_manual(name = col,values = mycols)
  colScale <- scale_colour_manual(name = col,values = mycols)
  
  if(fill == 'Y') return(fillScale)
  if(fill == 'N') return(colScale)
}


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
  p<-p + geom_bar(aes(feature, fill = feature))
  #p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  
  features_outfile<-paste("Breakpoints_features_count", ext, sep = "")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)
  
  p
}


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
  
  bp_data<-filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)
  
  
  bp_data$type <- ifelse(bp_data$type=="BND", "INV", as.character(bp_data$type))
  bp_data<-droplevels(bp_data)
  
  bp_data$type <- as.character(bp_data$type)
  cols<-setCols(bp_data, "type")
  
  # Reorder by count
  bp_data$type<-factor(bp_data$type, levels = names(sort(table(bp_data$type), decreasing = TRUE)))
  
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
          axis.title = element_text(size=20)
    )
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
 # p<-p + coord_flip()
 # p<-p + scale_y_reverse()

  types_outfile<-paste("sv_types_by_", object, ext, sep = "")
  cat("Writing file", types_outfile, "\n")
  ggsave(paste("plots/", types_outfile, sep=""), width = 5, height = 7)
  
  p
}




notchHits <- function(){
  bp_data<-getData()
  bp_data<-filter(bp_data, chrom == "X" & bp >= 2750000 & bp <= 3500000)

  bp_data$type <- ifelse(bp_data$type=="BND", "INV", as.character(bp_data$type))
  bp_data<-droplevels(bp_data)
  
  bp_data$type <- as.factor(bp_data$type)
  
  cols<-setCols(bp_data, "type", fill='N')  
  
  p<-ggplot(bp_data)
  p<-p + geom_jitter(aes(bp/1000000, sample, colour = type, shape=type),size = 3, alpha = 0.9)
  p<-p + guides(color = FALSE, size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          legend.position="top",
          axis.text.y = element_text(size = 10)
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(2.7,3.4,by=0.05), limits=c(2.7, 3.4))
  
  p<-p + annotate("rect", xmin=2.740000, xmax=3.134532, ymin=-1.5, ymax=0, alpha=.2, fill="green")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=-1.5, ymax=0, alpha=.2, fill="skyblue")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.334000, ymin=-1.5, ymax=0, alpha=.2, fill="red")
  p<-p + geom_vline(xintercept = 3.135669, colour="red", linetype="dotted")
  p<-p + cols
  
  # Nhits <- paste("Notch_hits.pdf")
  # cat("Writing file", Nhits, "\n")
  # ggsave(paste("plots/", Nhits, sep=""), width = 10, height = 10)
  # 
  
  # Density plot

  cols2<-setCols(bp_data, "type", fill='Y')  
  
  bp_data <- filter(bp_data, type != 'TRA')
  bp_data<-droplevels(bp_data)
  
  p2<-ggplot(bp_data)
  p2<-p2 + geom_density(aes(bp/1000000, fill = type), alpha = 0.4)
  p2<-p2 + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(2.7,3.4,by=0.05), limits=c(2.70, 3.4))
  p2<-p2 + scale_y_continuous("Density", expand = c(0,0))
  p2<-p2 + guides(colour = FALSE)
  p2<-p2 + geom_rug(data=bp_data,aes(bp/1000000))
  # p2<-p2 + annotate("rect", xmin=2.740000, xmax=3.134532, ymin=-1, ymax=0, alpha=.2, fill="green")
  # p2<-p2 + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=-1, ymax=0, alpha=.2, fill="skyblue")
  # p2<-p2 + annotate("rect", xmin=3.176440, xmax=3.334000, ymin=-1, ymax=0, alpha=.2, fill="red")
  p2<-p2 + geom_vline(xintercept = 3.135669, colour="red", linetype="dotted")
  
  p2<-p2 + cleanTheme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position="top",
          axis.title.y=element_blank(),
          strip.text = element_text(size=10)
          )
  
  p2 <- p2 + cols2
  p2 <- p2 + facet_wrap(~type, nrow = 3)
  
  combined_plots <- ggarrange(p, p2, 
                              labels = c("A", "B"),
                              ncol = 1, nrow = 2)
  
  NhitsDen <- paste("Notch_hits_density.pdf")
  cat("Writing file", NhitsDen, "\n")
  ggsave(paste("plots/", NhitsDen, sep=""), width = 10, height = 10)
  
  combined_plots
}

notchDels <- function(){
  infile = "data/all_bps_filtered.txt"
  bp_data<-read.delim(infile, header = F)
  colnames(bp_data) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  
  bp_data<-filter(bp_data, chrom == "X", bp >= 2750000, bp <= 3400000)
  
  bp_data <- filter(bp_data, type == "DEL" & bp_no == "bp1")
  
  # Make sample names unique to see samples with multiple DELS
  # bp_data$sample <- make.unique(as.character(bp_data$sample))
  # Hack to remove multiple deletions in same sample
  bp_data <- bp_data[!duplicated(bp_data$sample), ]
  bp_data <- transform(bp_data, sample = reorder(sample, -length))
  
  p<-ggplot(bp_data)
  p<-p + geom_bar(aes(sample,length),stat="identity")
  #p<-p + scale_y_discrete(expand = c(0.01,0.01), breaks=seq(0,500,by=100))
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=15))
  
  dels_out<-paste("NotchDels.pdf")
  cat("Writing file", dels_out, "\n")
  ggsave(paste("plots/", dels_out, sep=""), width = 20, height = 10)
  
  p
}



bpFeatureEnrichment <- function(features='data/genomic_features.txt', genome_length=137547960, print=NA){
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
  featuresFC<-arrange(featuresFC,desc(as.numeric(Log2FC)))
  featuresFC$Log2FC<-round(as.numeric(featuresFC$Log2FC), 1)
  featuresFC$expected<-round(as.numeric(featuresFC$expected), 1)
  
  if(!is.na(print)){
    cat("printing")
    first.step <- lapply(featuresFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)

    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))
    
    feat_enrichment_table <- paste("feature_enrichment_table.pdf")
    ggsave(paste("plots/", feat_enrichment_table, sep=""), width = 5.2, height = (nrow(featuresFC)/3))
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

bpGeneEnrichment <- function(gene_lengths="data/gene_lengths.txt", n=3, genome_length=137547960, print=NA){
  gene_lengths<-read.delim(gene_lengths, header = T)
  bp_data<-getData()
  bp_data<-filter(bp_data, gene != "intergenic")
  
  bp_count<-nrow(bp_data)
  
  hit_genes<-table(bp_data$gene)
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
  genesFC<-arrange(genesFC,desc(as.integer(fc)))
  genesFC$expected<-round(as.numeric(genesFC$expected),digits=2)
  genesFC$log2FC<-round(as.numeric(genesFC$log2FC),digits=1)
  
  if(!is.na(print)){
    cat("printing")
    first.step <- lapply(genesFC, unlist) 
    second.step <- as.data.frame(first.step, stringsAsFactors = F)
    
    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))
    
    gene_enrichment_table <- paste("gene_enrichment_table.pdf")
    ggsave(paste("plots/", gene_enrichment_table, sep=""), width = 5, height = (nrow(genesFC)/3))
  }
  
  else{ return(genesFC) }
}

bpGeneEnrichmentPlot <- function() {
  gene_enrichment<-bpGeneEnrichment(n=1)
  
  gene_enrichment$Log2FC <- log2(as.numeric(gene_enrichment$fc))
  
  gene_enrichment$gene <- as.character(gene_enrichment$gene)
  gene_enrichment$fc <- as.numeric(gene_enrichment$fc)
  
  gene_enrichment <- transform(gene_enrichment, gene = reorder(gene, -fc))
  
  gene_enrichment$test <- ifelse(gene_enrichment$Log2FC>=0, "enriched", "depleted")
  
  gene_enrichment <- filter(gene_enrichment, observed >= 5)
  gene_enrichment<-droplevels(gene_enrichment)
  
  highlightedGene <- filter(gene_enrichment, gene == "N")
  highlightedGene <- droplevels(highlightedGene)
  
  p<-ggplot(gene_enrichment)
  p<-p + geom_bar(aes(gene, Log2FC, fill = as.character(test)), stat="identity")
 # p<-p + geom_bar(data=highlightedGene, aes(gene, Log2FC, fill="red"), colour="black", stat="identity")
  p<-p + guides(fill=FALSE)
  p<-p + cleanTheme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 90, hjust=1),
          axis.text = element_text(size=7)
    )
  
  gene_enrichment_plot <- paste("gene_enrichment.pdf")
  cat("Writing file", gene_enrichment_plot, "\n")
  ggsave(paste("plots/", gene_enrichment_plot, sep=""), width = 5, height = 10)
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

dist2Feat <- function(feature_file="data/tss_locations.txt",sim=NA, print=0,send=0, feature='tss'){
  if(is.na(sim)){
    bp_data<-getData()
    bp_data<-rename(bp_data, pos = bp) 
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
    p<-p + cleanTheme() +
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


distOverlay <- function(feature_file="data/tss_locations.txt", feature='tss'){
  feature<-paste(toupper(substr(feature, 1, 1)), substr(feature, 2, nchar(feature)), sep='')

  if(feature=='promoter'){
    real_data<-dist2Feat(send=1, feature=feature)
    sim_data<-dist2Feat(feature=feature, sim=1, send=1)
  }
  else{ 
    real_data<-dist2Feat(feature_file=feature_file, send=1, feature=feature)
    sim_data<-dist2Feat(feature_file=feature_file, feature=feature, sim=1, send=1)
  }

  real_data$Source<-"Real"
  sim_data$Source<-"Sim"
  
  sim_data<-filter(sim_data, chrom != "Y", chrom != 4)
  sim_data<-droplevels(sim_data)
  real_data<-filter(real_data, chrom != "Y", chrom != 4)
  real_data<-droplevels(real_data)
  
  
  colours<-c( "#E7B800", "#00AFBB")
  
  
  p<-ggplot()
  p<-p + geom_density(data=real_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + geom_density(data=sim_data,aes(min_dist, fill = Source), alpha = 0.4)
  # p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  p<-p + scale_x_continuous(paste("Distance to", feature, "(Kb)", sep=' '),
                            limits=c(-10000, 10000),
                            breaks=c(-10000,-1000, 1000, 10000),
                            expand = c(.0005, .0005),
                            labels=c("-10", "-1", "1", "10") )
  p<-p + scale_y_continuous("Density")
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")

  p <- p + geom_rug(data=real_data,aes(min_dist, colour=Source),sides="b")
  p <- p + geom_rug(data=sim_data,aes(min_dist, colour=Source),sides="t")
  
  
  p <- p + scale_fill_manual(values=colours)
  p <- p + scale_colour_manual(values=colours)
  
  p<-p + cleanTheme() +
    theme(strip.text = element_text(size=20),
          legend.position="top")
  
  overlay<-paste("bp", feature, 'dist_overlay.pdf', sep='')
  cat("Writing file", overlay, "\n")
  ggsave(paste("plots/", overlay, sep=""), width = 25, height = 10)
  
  
  
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
  cru_positions<-read.delim("data/cruciform.txt", header=F)
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
  #p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,max(distances$bp),by=10))
  
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






rnaseqCompare <- function(maheva="data/isc_EB_genes_rnaSeq_Maheva.csv", dutta="data/isc_genes_rnaSeq.csv", upper = 1000, lower=1){
  maheva_data<-read.csv(header = F, maheva)
  dutta_data<-read.csv(header = F, dutta)
  
  colnames(maheva_data)<-c('id', 'mfpkm')
  colnames(dutta_data)<-c('id', 'dfpkm')

  expression_data<-join(maheva_data, dutta_data, 'id', type='left')

  expression_data$mfpkm <- ifelse(is.na(expression_data$mfpkm), 0.001, round(expression_data$mfpkm, 1))
  expression_data$dfpkm <- ifelse(is.na(expression_data$dfpkm), 0.001, round(expression_data$dfpkm, 1))
  expression_data$diff<- ifelse(expression_data$mfpkm > expression_data$dfpkm, "mH", "mL")

  
  # colnames(maheva_data)<-c('id', 'fpkm')
  # colnames(dutta_data)<-c('id', 'fpkm')
  # maheva_data$type<-"Maheva"
  # dutta_data$type<-"Dutta"
  # 
  # maheva_data$fpkm <- as.numeric(maheva_data$fpkm)
  # dutta_data$fpkm <- as.numeric(dutta_data$fpkm)
  # 
  # 
  # 
  # # expression_data$fpkm<-as.numeric(expression_data$fpkm)
  # # 
  # # expression_data<-droplevels(expression_data)
  # 
  expression_data<-filter(expression_data, mfpkm <= upper, dfpkm <= upper, dfpkm >= lower, mfpkm >= lower)
  # expression_data <- transform(expression_data, id = reorder(id, -fpkm))

  
  p <- ggplot(expression_data)
  p <- p + geom_point(aes(mfpkm, dfpkm, colour=diff))
  p <- p + scale_x_continuous("Maheva FPKM")
  p <- p + scale_y_continuous("Dutta FPKM")
  p <- p + ggtitle("Gene expression comparison between RNA-Seq experimeents")
  p
  
}










bpTssDist <- function(tss_pos="data/tss_positions.txt",sim=NA, print=0,return=0){
  tss_locations<-read.delim(tss_pos, header = T)
  tss_locations$tss<-as.integer(tss_locations$tss)
  
  if(is.na(sim)){
    bp_data<-getData()
    bp_data<-rename(bp_data, pos = bp) 
  }
  
  else{
    cat("Generating simulated bp_data\n")
    hit_count<-nrow(getData())
    bp_data<-bpSim(N=hit_count, write=print)
    colnames(bp_data)<-c("chrom", "pos", "v3", "v4", "v5")
    bp_data<-filter(bp_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    bp_data<-droplevels(bp_data)
  }
  
  
  
  # Will throw error if SVs don't exist on a chrom...
  
  # Removes chroms with fewer than 20 observations
  svCount <- table(bp_data$chrom)
  bp_data <- subset(bp_data, chrom %in% names(svCount[svCount > 10]))
  bp_data<-droplevels(bp_data)
  
  tss_locations <- subset(tss_locations, chrom %in% levels(bp_data$chrom))
  tss_locations<-droplevels(tss_locations)  
  
  
  fun2 <- function(p) {
    index<-which.min(abs(tss_df$tss - p))
    closestTss<-tss_df$tss[index]
    chrom<-as.character(tss_df$chrom[index])
    gene<-as.character(tss_df$gene[index])
    dist<-(p-closestTss)
    list(p, closestTss, dist, chrom, gene)
  }
  
  l <- list()
  
  for (c in levels(bp_data$chrom)){
    df<-filter(bp_data, chrom == c)
    tss_df<-filter(tss_locations, chrom == c)
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
  
  # Removes chroms with fewer than 20 observations
  # svCount <- table(dist2tss$chrom)
  # dist2tss <- subset(dist2tss, chrom %in% names(svCount[svCount > 10]))
  
  if(return==1){
    return(dist2tss)
  }
  else{
    p<-ggplot(dist2tss)
    p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
    p<-p + scale_x_continuous("Distance to TSS (Kb)",
                              limits=c(-10000, 10000),
                              breaks=c(-10000,-1000, 1000, 10000),
                              expand = c(.0005, .0005),
                              labels=c("-10", "-1", "1", "10") )
    p<-p + scale_y_continuous("Density")
    p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
    #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
    p <- p + geom_rug(aes(min_dist, colour=chrom))
    p<-p + cleanTheme() +
      theme(strip.text = element_text(size=20),
            legend.position="top")
    p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
    
    if(is.na(sim)){
      tssDistout<-paste("bpTSSdist.pdf")
    }
    else{
      tssDistout<-paste("bpTSSdist_sim.pdf")
    }
    cat("Writing file", tssDistout, "\n")
    ggsave(paste("plots/", tssDistout, sep=""), width = 20, height = 10)
    
    p
  }
}

bpTssDistOverlay <- function(){
  real_data<-bpTssDist(return=1)
  real_data$Source<-"Real"
  sim_data<-bpTssDist(sim=1, return=1)
  sim_data$Source<-"Sim"
  
  
  colours<-c( "#E7B800", "#00AFBB")
  
  
  p<-ggplot()
  p<-p + geom_density(data=real_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + geom_density(data=sim_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  p<-p + scale_x_continuous("Distance to TSS (Kb)",
                            limits=c(-10000, 10000),
                            breaks=c(-10000,-1000, 1000, 10000),
                            expand = c(.0005, .0005),
                            labels=c("-10", "-1", "1", "10") )
  p<-p + scale_y_continuous("Density")
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
  
  p <- p + geom_rug(data=real_data,aes(min_dist, colour=Source),sides="b")
  p <- p + geom_rug(data=sim_data,aes(min_dist, colour=Source),sides="t")
  
  
  p <- p + scale_fill_manual(values=colours)
  p <- p + scale_colour_manual(values=colours)
  
  p<-p + cleanTheme() +
    theme(strip.text = element_text(size=20),
          legend.position="top")
  
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  
  tssDistout<-paste("bpTSSdist_overlay.pdf")
  
  cat("Writing file", tssDistout, "\n")
  ggsave(paste("plots/", tssDistout, sep=""), width = 25, height = 10)
  
  
  
  p
  
  
}





bpG4Dist <- function(g4_pos="data/g4_positions.txt",sim=NA, print=0,return=0){
  g4_locations<-read.delim(g4_pos, header = T)
  g4_locations$g4<-as.integer(g4_locations$g4)
  
  if(is.na(sim)){
    bp_data<-getData()
  }
  
  else{
    cat("Generating simulated bp_data\n")
    hit_count<-nrow(getData())
    bp_data<-bpSim(N=hit_count, write=print)
    colnames(bp_data)<-c("chrom", "pos", "v3", "v4", "v5")
    bp_data<-filter(bp_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    bp_data<-droplevels(bp_data)
  }
  
  
  
  # Will throw error if SVs don't exist on a chrom...
  
  # Removes chroms with fewer than 20 observations
  svCount <- table(bp_data$chrom)
  bp_data <- subset(bp_data, chrom %in% names(svCount[svCount > 30]))
  bp_data<-droplevels(bp_data)
  
  g4_locations <- subset(g4_locations, chrom %in% levels(bp_data$chrom))
  g4_locations<-droplevels(g4_locations)  
  
  
  fun2 <- function(p) {
    index<-which.min(abs(g4_df$g4 - p))
    closestTss<-g4_df$g4[index]
    chrom<-as.character(g4_df$chrom[index])
    gene<-as.character(g4_df$gene[index])
    dist<-(p-closestTss)
    list(p, closestTss, dist, chrom, gene)
  }
  
  l <- list()
  
  for (c in levels(bp_data$chrom)){
    df<-filter(bp_data, chrom == c)
    g4_df<-filter(g4_locations, chrom == c)
    dist2g4<-lapply(df$bp, fun2)
    dist2g4<-do.call(rbind, dist2g4)
    dist2g4<-as.data.frame(dist2g4)
    
    colnames(dist2g4)=c("bp", "closest_g4", "min_dist", "chrom", "closest_gene")
    dist2g4$min_dist<-as.numeric(dist2g4$min_dist)
    l[[c]] <- dist2g4
  }
  
  dist2g4<-do.call(rbind, l)
  dist2g4<-as.data.frame(dist2g4)
  dist2g4$chrom<-as.character(dist2g4$chrom)
  
  dist2g4<-arrange(dist2g4,(abs(min_dist)))
  
  # Removes chroms with fewer than 20 observations
  # svCount <- table(dist2g4$chrom)
  # dist2g4 <- subset(dist2g4, chrom %in% names(svCount[svCount > 10]))
  
  if(return==1){
    return(dist2g4)
  }
  else{
    p<-ggplot(dist2g4)
    p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
    p<-p + scale_x_continuous("Distance to G4 (Kb)",
                              limits=c(-10000, 10000),
                              breaks=c(-10000,-1000, 1000, 10000),
                              expand = c(.0005, .0005),
                              labels=c("-10", "-1", "1", "10") )
    p<-p + scale_y_continuous("Density")
    p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
    #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
    p <- p + geom_rug(aes(min_dist, colour=chrom))
    p<-p + cleanTheme() +
      theme(strip.text = element_text(size=20),
            legend.position="top")
    p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
    
    if(is.na(sim)){
      g4Distout<-paste("bpG4dist.pdf")
    }
    else{
      g4Distout<-paste("bpG4dist_sim.pdf")
    }
    cat("Writing file", g4Distout, "\n")
    ggsave(paste("plots/", g4Distout, sep=""), width = 20, height = 10)
    
    p
  }
}



g4DistOverlay <- function(){
  real_data<-bpG4Dist(return=1)
  real_data$Source<-"Real"
  sim_data<-bpTssDist(sim=1, return=1)
  sim_data$Source<-"Sim"
  
  sim_data<-filter(sim_data, chrom != "Y", chrom != 4)
  sim_data<-droplevels(sim_data)
  real_data<-filter(real_data, chrom != "Y", chrom != 4)
  real_data<-droplevels(real_data)
  
  
  colours<-c( "#E7B800", "#00AFBB")
  
  
  p<-ggplot()
  p<-p + geom_density(data=real_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + geom_density(data=sim_data,aes(min_dist, fill = Source), alpha = 0.4)
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  p<-p + scale_x_continuous("Distance to G4 (Kb)",
                            limits=c(-10000, 10000),
                            breaks=c(-10000,-1000, 1000, 10000),
                            expand = c(.0005, .0005),
                            labels=c("-10", "-1", "1", "10") )
  p<-p + scale_y_continuous("Density")
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  #p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
  
  p <- p + geom_rug(data=real_data,aes(min_dist, colour=Source),sides="b")
  p <- p + geom_rug(data=sim_data,aes(min_dist, colour=Source),sides="t")
  
  
  p <- p + scale_fill_manual(values=colours)
  p <- p + scale_colour_manual(values=colours)
  
  p<-p + cleanTheme() +
    theme(strip.text = element_text(size=20),
          legend.position="top")
  
  p<-p + facet_wrap(~chrom, ncol = 3, scales = "free_y")
  
  
  g4Distout<-paste("bpG4dist_overlay.pdf")
  
  cat("Writing file", g4Distout, "\n")
  ggsave(paste("plots/", g4Distout, sep=""), width = 25, height = 10)
  
  
  
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

