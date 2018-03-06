#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords parse breakpoints data
#' @import plyr
#' @import tidyverse
#' @export
#' @return Dataframe

getData <- function(infile = system.file("extdata", "all_bps_filtered.txt", package="svBreaks"),
                    gene_lengths_file = system.file("extdata", "gene_lengths.txt", package="svBreaks"),
                    expression_data = system.file("extdata", "isc_genes_rnaSeq.csv", package="svBreaks")
                    ) {
  bp_data <- read.delim(infile, header = F)

  # colnames(bp_data) <- c("event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "type", "length")
  colnames(bp_data) <- c("event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "chrom2",  "bp2", "gene2", "feature2", "type", "length")

  # bp_data$allele_freq <- suppressWarnings(as.numeric(as.character(bp_data$allele_freq)))

  gene_lengths <- read.delim(gene_lengths_file, header = T)
  gene_lengths <- gene_lengths[, c("gene", "id")]

  bp_data <- plyr::join(bp_data, gene_lengths, "gene", type = "left")

  # Read in tissue specific expression data
  seq_data <- read.csv(header = F, expression_data)
  colnames(seq_data) <- c("id", "fpkm")

  bp_data <- plyr::join(bp_data, seq_data, "id", type = "left")
  bp_data$fpkm <- ifelse(is.na(bp_data$fpkm), 0, round(bp_data$fpkm, 1))
  bp_data$bp <- as.numeric(as.character(bp_data$bp))

  # Filter for genes expressed in RNA-Seq data
  # bp_data<-filter(bp_data, fpkm > 0.1)

  # Filter for genes NOT expressed in RNA-Seq data
  # cat("Filtering out expressed genes\n")
  # bp_data<-filter(bp_data, fpkm == 0)

  # filter on chroms
  bp_data <- dplyr::filter(bp_data, chrom != "211000022280116")

  # filter out samples
  bp_data <- dplyr::filter(bp_data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" & sample != "A373R11")

  # filter on genotype
  # bp_data<-dplyr::filter(bp_data, genotype != 'germline_recurrent')
  # bp_data<-dplyr::filter(bp_data, genotype == 'germline_private')
  bp_data<-dplyr::filter(bp_data, genotype == 'somatic_tumour')
  # bp_data<-dplyr::filter(bp_data, genotype == 'somatic_normal')

  # Filter for old/new data
  # bp_data <- filter(bp_data, !grepl("^A|H", sample))

  # Filter for SV type
  # bp_data <- filter(bp_data, type == "DEL")

  bp_data <- droplevels(bp_data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(bp_data)
}


#' notchFilt
#'
#' Function to filter in/out events affecting a locus (in this case the Drosophila Notch locus)
#' @param infile File to process [Required]
#' @keywords parse
#' @import tidyverse
#' @export
#' @return Dataframe
#'
notchFilt <- function(keep=0) {
  bp_data <- getData()
  bp_data <- filter(bp_data, genotype == 'somatic_tumour')
  if (keep) {
    cat("Selecting for bps in Notch\n")
    notchIn <- bp_data %>%
      filter(chrom == "X" & bp >= 2700000 & bp2 <= 3400000) %>%
      # filter(gene != 'N', gene2 != 'N') %>%
      droplevels()
    return(notchIn)
  }
  else{
    cat("Excluding bps in Notch\n")
    noNotch <- bp_data %>%
      filter(!(chrom == "X" & bp >= 2700000 & bp2 <= 3400000)) %>%
      filter(gene != 'N', gene2 != 'N') %>%
      droplevels()
    return(noNotch)
  }
}
