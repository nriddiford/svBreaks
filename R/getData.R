#' getData
#'
#' Function to clean cnv files
#' @param infile File to process [Required]
#' @keywords parse breakpoints data
#' @import plyr dplyr
#' @export
#' @return Dataframe

getData <- function(...,
                    infile = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_bps_mech.txt',
                    attach_info = '../mutationProfiles/data/samples_names_conversion.txt',
                    gene_lengths_file = system.file("extdata", "gene_lengths.txt", package="svBreaks"),
                    expression_data = system.file("extdata", "isc_genes_rnaSeq.csv", package="svBreaks"),
                    exclude = TRUE
                    ) {
  bp_data <- read.delim(infile, header = F)

  cat("Filters applied:\n")
  input_list <- as.list(substitute(list(...)))
  lapply(X=input_list, function(x) {str(x);summary(x)})
  
  colnames(bp_data) <- c("event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "chrom2",  "bp2", "gene2", "feature2", "type", "length", "af", "confidence")
  if(ncol(bp_data) == 18){
    cat("Breakpoints are annotated with mechansims\n")
    colnames(bp_data)[c(17,18)] <- c("microhomology", "mechanism")
  }
  
  if(file.exists(attach_info)){
    name_conversion <- read.delim(attach_info, header=F)
    cat("Attaching assay information to data")
    colnames(name_conversion) <- c("sample", "sample_short", "sex", "assay")
    bp_data <- plyr::join(bp_data, name_conversion, "sample", type = 'left')
  }
  
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

  excluded_samples <- c()
  if(exclude){
    excluded_samples <- c("A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05",
                          "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  }  
  
  bp_data <- bp_data %>%
    dplyr::filter(!sample %in% excluded_samples) %>%
    dplyr::mutate(fpkm = ifelse(is.na(fpkm), 0, round(fpkm, 1))) %>%
    dplyr::mutate(bp = as.numeric(bp)) %>%
    dplyr::mutate(genotype = factor(genotype)) %>%
    dplyr::mutate(type = factor(type)) %>% 
    dplyr::mutate(cell_fraction = ifelse(chrom %in% c('X', 'Y'), af,
                                         ifelse(af*2>1, 1, af*2))) %>%
    dplyr::filter(...) %>%
    droplevels()
  
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(bp_data)
}


#' showSamples
#'
#' A helper function to print the sample names
#' @param infile File to process [Required]
#' @keywords samples
#' @export
showSamples <- function(..., bp_data = NULL){
  if(missing(bp_data)){
    bp_data <- getData(..., exclude = F)
  }
  # print(levels(bp_data$sample))
  dput(as.character(levels(bp_data$sample)))
}
