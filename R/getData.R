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
                    attach_info = system.file("extdata", "samples_names_conversion.txt", package="svBreaks"),
                    gene_lengths_file = system.file("extdata", "gene_lengths.txt", package="svBreaks"),
                    expression_data = system.file("extdata", "isc_genes_rnaSeq.txt", package="svBreaks"),
                    expression_source = 'flygut') {
  
  bp_data <- read.delim(infile) 

  cat("Filters applied:\n")
  input_list <- as.list(substitute(list(...)))
  lapply(X=input_list, function(x) {str(x);summary(x)})
  
  colnames(bp_data) <- c("event", "bp_no", "sample", "genotype", "chrom", "bp", "gene", "feature", "chrom2",  "bp2", "gene2", "feature2", "type", "length", "af", "confidence")
  if(ncol(bp_data) >= 18){
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
  
  if(expression_source == 'flygut'){
    cat("Reading expression data from source: 'flygut [Buchon]'")
    expression_data = system.file("extdata", "Buchon_summary_ISCs.txt", package="svBreaks")
    expression_data <- read.delim(expression_data, header=F)
    colnames(expression_data) <- c('id', 'symbol', 'name', 'isc', 'eb', 'ec', 'ee', 'vm')
    seq_data <- expression_data %>% 
      dplyr::mutate(fpkm = isc) %>% 
      dplyr::select(id, fpkm)
  } else{
    cat("Reading expression data from source: 'Dutta'")
    expression_data = system.file("extdata", "isc_genes_rnaSeq.txt", package="svBreaks")
    seq_data <- read.table(header = F, expression_data)
    colnames(seq_data) <- c("id", "fpkm")
  }
  
  bp_data <- plyr::join(bp_data, seq_data, "id", type = "left")
  
  bp_data <- bp_data %>%
    dplyr::mutate(af = as.double(as.character(af))) %>% 
    dplyr::mutate(fpkm = ifelse(is.na(fpkm), 0, round(as.numeric(fpkm), 1))) %>%
    dplyr::mutate(bp = as.numeric(bp)) %>%
    dplyr::mutate(genotype = factor(genotype)) %>%
    dplyr::mutate(type = factor(type)) %>% 
    dplyr::mutate(type2 = as.character(ifelse(stringr::str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>% 
    dplyr::mutate(type2 = factor(type2)) %>% 
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
