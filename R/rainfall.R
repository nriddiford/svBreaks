
#' bpRainfall
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @import dplyr
#' @keywords rainfall
#' @export

bpRainfall <- function(..., bp_data = NULL, bed_file=NULL, write = FALSE, chroms= c("2L", "2R", "3L", "3R", "X"), from=NULL, to=NULL, tick_by=10) {
  if(missing(bp_data) && missing(bed_file) ) bp_data <- getData(..., genotype=='somatic_tumour')
  if(!missing(bed_file)){
    bp_data <- read.delim(bed_file, header = F) %>% 
      dplyr::mutate(chrom = V1,
                    bp = V2) %>% 
      dplyr::select(chrom, bp)
      
  }
  
  bp_data <- bp_data %>% 
    dplyr::filter(...) %>% 
    droplevels()
  
  # bp_data <- svBreaks::transform_types(bp_data)
  
  distances <- do.call(rbind, lapply(
    split(bp_data[order(bp_data$chrom, bp_data$bp), ], bp_data$chrom[order(bp_data$chrom, bp_data$bp)]),
    function(a)
      data.frame(
        a,
        dist = c(diff(a$bp), NA),
        logdist = c(log10(diff(a$bp)), NA)
      )
    )
    )
  
  distances$logdist[is.infinite(distances$logdist)] <- 0
  distances <- dplyr::filter(distances, chrom %in% chroms)

  if( !missing(from) && !missing(to) ) {
    distances <- distances %>% 
      dplyr::filter(bp >= (from -5000),
                    bp <= (to + 5000),
                    chrom %in% chroms)
  }
  
  colours = svBreaks::sv_colours()
  
  
  p <- ggplot(distances)
  if(missing(bed_file)) {
    p <- p + geom_point(aes(bp / 1e6, logdist, colour = type2), size=1, alpha=0.8)
  } else{
    p <- p + geom_point(aes(bp / 1e6, logdist))
  }
  p <- p + cleanTheme() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey80", size = 0.8, linetype = "dotted"),
      strip.text = element_text(size = 20)
    )

  p <- p + facet_wrap(~chrom, scales = "free_x", ncol = length(chroms))
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, max(distances$bp), by = tick_by))
  p <- p + scale_y_continuous("Genomic Distance")
  # p <- p + scale_fill_manual("SV type\n", values = colours)
  p <- p + scale_colour_manual(values = colours)
  
  if(write){
    rainfall_out <- paste("rainfall.pdf")
    cat("Writing file", rainfall_out, "\n")
    ggsave(paste("plots/", rainfall_out, sep = ""), width = 15, height = 5)
  }
  p
}
