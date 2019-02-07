#' Plot mechanims vs length
#' mechansimSize
#'
#' Function to plot the sizes SVs annotated as being underscored by mechansims
#' @param infile File to process [Required]
#' @keywords parse
#' @import dplyr
#' @export
mechansimSize <- function(..., bp_data=NULL, infile='~/Desktop/parserTest/filtered_231018/summary/merged/all_bps_mech.txt', plot=TRUE) {
  
  if(missing(bp_data)){
    bp_data <- svBreaks::getData(infile=infile)
  }
  
  mechanism_count <- bp_data %>%
    dplyr::filter(mechanism != '-',
                  !stringr::str_detect(type, 'TRA'),
                  bp_no != "bp2") %>%
    dplyr::mutate(type2 = as.character(ifelse(stringr::str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>% 
    dplyr::mutate(length = ifelse(length == 0, 0.01, length)) %>%
    dplyr::add_count(mechanism)

  if(plot){
    p <- ggplot(mechanism_count, aes(fct_reorder(mechanism, -n), length))
    p <- p + geom_violin(aes(fill = mechanism), alpha=0.6)
    p <- p + geom_jitter(width = .1, height = .1, size = 0.5, alpha = 0.5)
    
    p <- p + scale_y_log10("Length (Kb)")
    p <- p + scale_x_discrete()
    p <- p + slideTheme() +
      theme(
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20)
      )
    p <- p + facet_wrap(~type2, nrow = 3)
    p
  } else return(mechanism_count)
}


#' Plot micromology vs type
#' micromology
#'
#' Function to plot the sizes SVs annotated as being underscored by mechansims
#' @param infile File to process [Required]
#' @keywords parse
#' @import dplyr
#' @export
micromologyPlot <- function(..., bp_data=NULL, infile='~/Desktop/parserTest/filtered_231018/summary/merged/all_bps_mech.txt'){
  mech_data <- svBreaks::mechansimSize(..., bp_data=bp_data, plot=F)
  
  mh_sizes <- mech_data %>% 
    dplyr::mutate(mh_length = ifelse(microhomology==0, 0, nchar(as.character(microhomology))))
  
  p <- ggplot(mh_sizes)
  p <- p + geom_bar(aes(mh_length, (..count..)),stat='count')
  p <- p + scale_x_continuous("Microhomology")
  p <- p + scale_y_continuous("Count")
  
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  # p <- p + facet_wrap(~type2, nrow = 2)
  p
}