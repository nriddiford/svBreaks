library(scales)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


getData <- function(infile = "data/all_bps.txt"){
  bp_data<-read.delim(infile, header = F)
  colnames(bp_data) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")
  
  #filter on chroms
  #bp_data<-filter(bp_data, chrom != "Y" & chrom != "4")
  #filter out samples
  bp_data<-filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
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


setCols <- function(df, col){
  names<-levels(df[[col]])
  cat("Setting colour levles:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, "Set2")
  names(mycols) <- names
  colScale <- scale_fill_manual(name = col,values = mycols)
  return(colScale)
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
    bp_data<-exNotch()
    ext<-'_excl.N.pdf'
  }
  else{
    bp_data<-getData()
    ext<-'.pdf'
  }
  cols<-setCols(bp_data, "type")
  
  # Reorder by count
  bp_data$type<-factor(bp_data$type, levels = names(sort(table(bp_data$type), decreasing = TRUE)))
  
  # Only take bp1 for each event
  bp_data<-filter(bp_data, bp_no != "bp2")
  
  p<-ggplot(bp_data)
  p<-p + geom_bar(aes(get(object), fill = type))
  p<-p + cols
  p<-p + cleanTheme() +
    theme(axis.title.x=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size=20)
    )
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))
  
  types_outfile<-paste("sv_types_by_", object, ext, sep = "")
  cat("Writing file", types_outfile, "\n")
  ggsave(paste("plots/", types_outfile, sep=""), width = 20, height = 10)
  
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


notch_hits <- function(){
  bp_data<-getData()
  bp_data<-filter(bp_data, chrom == "X", bp >= 2750000, bp <= 3400000)
  
  p<-ggplot(bp_data)
  p<-p + geom_point(aes(bp/1000000, sample, colour = sample, shape = type, size = 2))
  p<-p + guides(color = FALSE, size = FALSE, sample = FALSE)
  p<-p + cleanTheme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(3,3.3,by=0.05), limits=c(3, 3.301))
  
  p<-p + annotate("rect", xmin=3.000000, xmax=3.134532, ymin=0, ymax=0.5, alpha=.2, fill="green")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=0, ymax=0.5, alpha=.2, fill="skyblue")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.300000, ymin=0, ymax=0.5, alpha=.2, fill="red")
  
  p
}

notchDels <- function(){
  infile = "data/all_bps_new.txt"
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
          axis.text = element_text(size=20))
  
  dels_out<-paste("NotchDels.pdf")
  cat("Writing file", dels_out, "\n")
  ggsave(paste("plots/", dels_out, sep=""), width = 20, height = 10)
  
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


bpFeatureEnrichment <- function(features='data/genomic_features.txt', genome_length=137547960){
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
      list(feature = f, observed = classCount[f], expected = featureExpect, fc = fc, test = test, sig = sig_val, p_val = p_val)
    }
  }
  
  enriched<-lapply(levels(bp_data$feature), fun)
  enriched<-do.call(rbind, enriched)
  featuresFC<-as.data.frame(enriched)
  # Sort by FC value
  featuresFC<-arrange(featuresFC,desc(as.integer(fc)))
  return(featuresFC)
}

bpEnrichmentPlot <- function() {
  feature_enrichment<-bpFeatureEnrichment()
  
  
  feature_enrichment$Log2FC <- log2(as.numeric(feature_enrichment$fc))
  
  feature_enrichment$feature <- as.character(feature_enrichment$feature)
  feature_enrichment$fc <- as.numeric(feature_enrichment$fc)
  
  feature_enrichment <- transform(feature_enrichment, feature = reorder(feature, -fc))
  
  
  
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

bpGeneEnrichment <- function(gene_lengths="data/gene_lengths.txt", n=3, genome_length=137547960){
  gene_lengths<-read.delim(gene_lengths, header = T)
  bp_data<-getData()
  bp_data<-filter(bp_data, gene != "intergenic")
  
  bp_count<-nrow(bp_data)
  
  hit_genes<-table(bp_data$gene)
  genes<-setNames(as.list(gene_lengths$length), gene_lengths$gene)
  
  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction<-genes[[g]]/genome_length
    
    # How many times should we expect to see this gene hit in our bp_data (given number of obs. and fraction)?
    gene_expect<-bp_count*(genefraction)
    
    # observed/expected 
    fc<-hit_genes[[g]]/gene_expect
    fc<-round(fc,digits=1)
    gene_expect<-round(gene_expect,digits=3)
    list(gene = g, length = genes[[g]], observed = hit_genes[g], expected = gene_expect, fc = fc)
  }
  
  enriched<-lapply(levels(bp_data$gene), fun)
  enriched<-do.call(rbind, enriched)
  genesFC<-as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC<-filter(genesFC, observed >= n)
  # Sort by FC value
  genesFC<-arrange(genesFC,desc(as.integer(fc)))
  return(genesFC)
}


bpTssDist <- function(tss_pos="data/tss_positions.txt",sim=NA, print=0){
  tss_locations<-read.delim(tss_pos, header = T)
  tss_locations$tss<-as.integer(tss_locations$tss)
  
  if(is.na(sim)){
    bp_data<-getData()
  }
  
  else{
    cat("Generating simulated bp_data\n")
    bp_data<-snvSim(N=1000, write=print)
    colnames(bp_data)<-c("chrom", "pos", "v3", "v4", "v5")
    bp_data<-filter(bp_data, chrom == "2L" | chrom == "2R" | chrom == "3L" | chrom == "3R" | chrom == "X" | chrom == "Y" | chrom == "4")
    bp_data<-droplevels(bp_data)
  }
  
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
    dist2tss<-lapply(df$bp, fun2)
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
  snvCount <- table(dist2tss$chrom)
  dist2tss <- subset(dist2tss, chrom %in% names(snvCount[snvCount > 15]))
  
  p<-ggplot(dist2tss)
  p<-p + geom_density(aes(min_dist, fill = chrom), alpha = 0.3)
  p<-p + scale_x_continuous("Distance to TSS (Kb)",
                            limits=c(-100000, 100000),
                            breaks=c(-100000, -10000, -1000, 0, 1000, 10000, 100000),
                            expand = c(.0005, .0005),
                            labels=c("-100", "-10", "-1", 0, "1", "10", "100") )
  p<-p + scale_y_continuous("Density", expand = c(0, 0))
  p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  p<-p + cleanTheme()
  #p<-p + facet_wrap(~chrom, ncol = 2, scales = "free_y")
  p
  
  # p<-ggplot(dist2tss)
  # p<-p + geom_histogram(aes(min_dist, fill = chrom), alpha = 0.6, bins=500)
  # p<-p + scale_x_continuous("Distance to TSS", limits=c(-1000, 1000))
  # p<-p + geom_vline(xintercept = 0, colour="black", linetype="dotted")
  # p
  
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
  
  bp_data<-filter(bp_data, sample != "A373R11" & sample != 'A373R13')
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
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 5)
  #p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,max(distances$bp),by=10))
  
  rainfall_out<-paste("rainfall.pdf")
  cat("Writing file", rainfall_out, "\n")
  ggsave(paste("plots/", rainfall_out, sep=""), width = 15, height = 5)
  
  p
}

