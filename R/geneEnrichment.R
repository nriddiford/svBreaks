# Functions related to calculating enrichment
# of svs in genes

# bpGeneEnrichment
#' Calculate the enrichment of SVs in genes
#' @keywords enrichment
#' @import tidyverse
#' @export
bpGeneEnrichment <- function(..., gene_lengths = system.file("extdata", "gene_lengths.txt", package="svBreaks"), n=3, genome_length=118274340, print=NA) {
  cat("Showing genes hit at least", n, "times", "\n")
  gene_lengths <- read.delim(gene_lengths, header = T)
  # bp_data<-read.delim('data/all_samples.txt',header=T)
  bp_data <- getData(..., gene != "intergenic", confidence == 'precise')
 
  bp_count <- nrow(bp_data)

  hit_genes <- table(bp_data$gene)

  # hit_genes<-table(bp_data[!duplicated(bp_data),]$gene)

  genes <- setNames(as.list(gene_lengths$length), gene_lengths$gene)
  expressed_genes <- setNames(as.list(bp_data$fpkm), bp_data$gene)

  fun <- function(g) {
    # Calculate the fraction of geneome occupied by each gene
    genefraction <- genes[[g]] / genome_length

    # How many times should we expect to see this gene hit in our bp_data (given number of obs. and fraction)?
    gene_expect <- bp_count * (genefraction)

    # observed/expected
    fc <- hit_genes[[g]] / gene_expect
    fc <- round(fc, digits = 1)
    log2FC <- log2(fc)

    gene_expect <- round(gene_expect, digits = 3)
    list(gene = g, length = genes[[g]], fpkm = expressed_genes[[g]], observed = hit_genes[g], expected = gene_expect, fc = fc, log2FC = log2FC)
  }

  enriched <- lapply(levels(bp_data$gene), fun)
  enriched <- do.call(rbind, enriched)
  genesFC <- as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC <- dplyr::filter(genesFC, observed >= n)
  # Sort by FC value
  genesFC <- dplyr::arrange(genesFC, desc(as.integer(fc)))
  genesFC$expected <- round(as.numeric(genesFC$expected), digits = 2)
  genesFC$log2FC <- round(as.numeric(genesFC$log2FC), digits = 1)

  if (!is.na(print)) {
    cat("printing")
    first.step <- lapply(genesFC, unlist)
    second.step <- as.data.frame(first.step, stringsAsFactors = F)

    ggtexttable(second.step, rows = NULL, theme = ttheme("mBlue"))

    gene_enrichment_table <- paste("gene_enrichment_table.tiff")
    ggsave(paste("plots/", gene_enrichment_table, sep = ""), width = 5, height = (nrow(genesFC) / 3), dpi = 300)
  }

  else {
    return(genesFC)
  }
}

# bpGeneEnrichmentPlot
#' Plot the enrichment of SVs in genes
#' @keywords enrichment
#' @import tidyverse
#' @export
#'
bpGeneEnrichmentPlot <- function(n=2) {
  gene_enrichment <- bpGeneEnrichment(n = n)

  gene_enrichment$Log2FC <- log2(as.numeric(gene_enrichment$fc))

  gene_enrichment$gene <- as.character(gene_enrichment$gene)
  gene_enrichment$fc <- as.numeric(gene_enrichment$fc)

  gene_enrichment <- transform(gene_enrichment, gene = reorder(gene, -fc))

  gene_enrichment$test <- ifelse(gene_enrichment$Log2FC >= 0, "#469BDBFE", "#E3464EC6")

  gene_enrichment <- dplyr::filter(gene_enrichment, observed >= 2)
  gene_enrichment <- droplevels(gene_enrichment)

  p <- ggplot(gene_enrichment)
  p <- p + geom_bar(aes(gene, Log2FC, fill = test), stat = "identity")
  p <- p + guides(fill = FALSE)
  p <- p + scale_fill_identity()

  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text = element_text(size = 40), axis.title = element_text(size = 90)
    )

  gene_enrichment_plot <- paste("gene_enrichment.pdf")
  cat("Writing file", gene_enrichment_plot, "\n")
  ggsave(paste("plots/", gene_enrichment_plot, sep = ""), width = 30, height = 10)
  p
}


getPromoter <- function(gene_lengths_in="data/gene_lengths.txt") {
  gene_lengths <- read.delim(gene_lengths_in, header = T)

  gene_lengths$promoter <- ifelse(gene_lengths$start < gene_lengths$end,
                                  gene_lengths$start - 1500,
                                  gene_lengths$end + 1500
  )


  gene_lengths <- gene_lengths[, c("chrom", "promoter")]
  colnames(gene_lengths) <- NULL
  return(gene_lengths)
}


# bpAllGenes
#' Calculate the enrichment of genes affected by SVs (as opposed to breakpoint hits)
#' @keywords enrichment
#' @import tidyverse
#' @export
bpAllGenes <- function(..., gene_lengths_in = system.file("extdata", "gene_lengths.txt", package="svBreaks"),
                       affected_genes = system.file("extdata", "all_genes_filtered.txt", package="svBreaks"),
                       n=3, genome_length=118274340, outfile=NULL) {
  cat("Showing genes affected by a SV at least", n, "times", "\n")
  gene_lengths <- read.delim(gene_lengths_in, header = T)
  allGenes <- read.delim(affected_genes, header = T)
  colnames(allGenes) <- c("event", "sample", "genotype", "type", "allele_frequency", "chromosome", "gene")
  excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  
  allGenes <- allGenes %>% 
    dplyr::filter(genotype=='somatic_tumour', !sample %in% excluded_samples) %>%
    dplyr::mutate(cell_fraction = ifelse(chromosome %in% c('X', 'Y'), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1, allele_frequency*2))) %>% 
    dplyr::filter(...) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::distinct(gene) %>% 
    droplevels()
  
  geneCount <- nrow(allGenes)
  hit_genes <- table(allGenes$gene)

  genes <- setNames(as.list(gene_lengths$length), gene_lengths$gene)
  chroms <- setNames(as.list(gene_lengths$chrom), gene_lengths$gene)

  # expressed_genes<-setNames(as.list(bp_data$fpkm), bp_data$gene)


  fun <- function(g) {
    # Calculate the fraction of genome occupied by each gene
    genefraction <- genes[[g]] / genome_length

    # How many times should we expect to see this gene hit in our bp_data (given number of obs. and fraction)?
    gene_expect <- geneCount * (genefraction)

    # observed/expected
    fc <- hit_genes[[g]] / gene_expect
    
    fc <- round(fc, digits = 1)
    log2FC <- log2(fc)

    gene_expect <- round(gene_expect, digits = 3)
  
    list(gene = g, length = genes[[g]], chromosome = as.character(chroms[[g]]), observed = hit_genes[g], expected = gene_expect, fc = fc, log2FC = log2FC)
  }

  enriched <- lapply(levels(allGenes$gene), fun)
  enriched <- do.call(rbind, enriched)
  genesFC <- as.data.frame(enriched)
  # Filter for genes with few observations
  genesFC <- genesFC %>% 
    dplyr::filter(observed >= n) %>% 
    dplyr::mutate(expected = round(as.numeric(expected), digits = 2)) %>%
    dplyr::mutate(log2FC = round(as.numeric(log2FC), digits = 1)) %>% 
    dplyr::arrange(-as.integer(observed),
                   -log2FC) %>% 
    droplevels()

  geneEnrich <- unnest(genesFC)
  geneEnrich <- geneEnrich %>% 
    dplyr::select(gene, length, chromosome, observed, expected, fc, log2FC)
  
  if(!missing(outfile)){
    write.table(geneEnrich, outfile, sep="\t", row.names=FALSE)
  }

  return(geneEnrich)
}
