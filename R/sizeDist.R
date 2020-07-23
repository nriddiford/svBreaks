#' sizeDist
#' Plot the size distribution for different types of structual variants
#' @keywords size
#' @import dplyr ggplot2
#' @importFrom forcats fct_reorder fct_recode
#' @import RColorBrewer
#' @import ggsci
#' @export

sizeDist <- function(..., bp_data) {
  if(missing(bp_data)) bp_data <- getData(..., genotype=='somatic_tumour')
  
  bp_data <- bp_data %>%
    dplyr::filter(type != "TRA", type != "BND", bp_no != "bp2") %>%
    dplyr::mutate(length = ifelse(length == 0, 0.01, length)) %>%
    dplyr::mutate(type = ifelse(type == "TANDUP", "DUP", as.character(type))) %>%
    dplyr::add_count(genotype) %>%
    dplyr::mutate(genotype_short  = forcats::fct_recode(genotype, G_r = "germline_recurrent", G_p = "germline_private", S_n = "somatic_normal", S_t = "somatic_tumour"))
  
  p <- ggplot(bp_data, aes(forcats::fct_reorder(genotype_short, -n), length))
  p <- p + geom_violin(aes(fill = genotype), alpha=0.6)
  p <- p + geom_jitter(width = .1, height = .1, size = 0.25, alpha = 0.1)

  p <- p + scale_y_log10("Length (Kb)")
  p <- p + scale_x_discrete()
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + facet_wrap(~type, nrow = 3)
  p <- p + scale_fill_jco()
  
  sizedistOut <- paste("sizeDist.png")
  cat("Writing file", sizedistOut, "\n")
  ggsave(paste("plots/", sizedistOut, sep = ""), width = 15, height = 10)

  p
}


#' genotypeTypeCount
#' Plot the contribution of sv tpyes to toal load per genotype
#' @keywords genotype
#' @import dplyr ggplot2 ggsci
#' @importFrom forcats fct_reorder fct_recode
#' @export
genotypeTypeCount <- function(...) {
  bp_data <- getData(!sample %in% c("A373R1", "A373R7", "A512R17", "A373R11", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9"))
  
  bp_data <- bp_data %>%
    dplyr::filter(type != "TRA", type != "BND", bp_no != "bp2") %>%
    dplyr::mutate(length = ifelse(length == 0, 0.01, length)) %>%
    dplyr::mutate(type = as.character(ifelse(type == "TANDUP", "DUP", as.character(type)))) %>%
    dplyr::mutate(genotype_short  = forcats::fct_recode(genotype, G_r = "germline_recurrent", G_p = "germline_private", S_n = "somatic_normal", S_t = "somatic_tumour")) %>% 
    dplyr::count(genotype_short, type) %>% 
    dplyr::group_by(genotype_short) %>% 
    dplyr::mutate(freq = n / sum(n)) %>%
    dplyr::ungroup()
  
  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(forcats::fct_reorder(genotype_short, freq), freq, fill=type), alpha=0.7, stat = 'identity', position=position_dodge()) 
  p <- p + scale_y_continuous(labels=scales::percent) +
    ylab("Frequency")
  p <- p + slideTheme() +
    theme(
      axis.title.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + scale_fill_jco()
  
  genotypeTypeCount <- paste("genotypeTypeCount.png")
  cat("Writing file", genotypeTypeCount, "\n")
  ggsave(paste("plots/", genotypeTypeCount, sep = ""), width = 10, height = 10)
  
  p
}

