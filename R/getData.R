#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords parse breakpoints data
#' @import plyr
#' @import dplyr
#' @export
#' @return Dataframe

getData <- function(...,
                    infile = 'inst/extdata/all_bps_filtered.txt',
                    gene_lengths_file = system.file("extdata", "gene_lengths.txt", package="svBreaks"),
                    expression_data = system.file("extdata", "isc_genes_rnaSeq.csv", package="svBreaks")
                    ) {
  bp_data <- read.delim(infile, header = F)

  cat("Filters applied:\n")
  input_list <- as.list(substitute(list(...)))
  lapply(X=input_list, function(x) {str(x);summary(x)})

  colnames(bp_data) <- c("event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "chrom2",  "bp2", "gene2", "feature2", "type", "length", "af")

  # bp_data$allele_freq <- suppressWarnings(as.numeric(as.character(bp_data$allele_freq)))

  gene_lengths <- read.delim(gene_lengths_file, header = T)
  gene_lengths <- gene_lengths[, c("gene", "id")]

  bp_data <- plyr::join(bp_data, gene_lengths, "gene", type = "left")

  # Read in tissue specific expression data
  seq_data <- read.table(header = F, expression_data)
  colnames(seq_data) <- c("id", "fpkm")
  bp_data <- plyr::join(bp_data, seq_data, "id", type = "left")
  
  # bp_data <- bp_data %>%
  #   mutate(type = as.factor(case_when(
  #     type == 'BND' & chrom == chrom2 ~ 'INV',
  #     type == 'BND' & chrom != chrom2 ~ 'TRA',
  #     TRUE  ~ as.character(type))))
  # 
  bp_data <- bp_data %>%
    mutate(fpkm = ifelse(is.na(fpkm), 0, round(fpkm, 1))) %>%
    mutate(bp = as.numeric(bp)) %>%
    mutate(genotype = factor(genotype)) %>%
    mutate(type = factor(type)) %>% 
    dplyr::filter(...) %>%
    droplevels()
  
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(bp_data)
}
