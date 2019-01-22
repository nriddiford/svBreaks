
#' bpRainfall
#' Plot log10 distances between snvs as rainfall plot
#' @import ggplot2
#' @import dplyr
#' @keywords rainfall
#' @export

bpRainfall <- function(..., bp_data = NULL, write = FALSE) {
  if(missing(bp_data)){
    bp_data <- getData(..., genotype=='somatic_tumour')
  }
  
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
  distances <- dplyr::filter(distances, chrom %in% c("2L", "2R", "3L", "3R", "X"))

  p <- ggplot(distances)
  p <- p + geom_point(aes(bp / 1000000, logdist, colour = -cell_fraction))
  p <- p + cleanTheme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      strip.text = element_text(size = 20)
    )

  p <- p + facet_wrap(~chrom, scales = "free_x", ncol = 5)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, max(distances$bp), by = 10))
  p <- p + scale_y_continuous("Genomic Distance")
  
  if(write){
    rainfall_out <- paste("rainfall.pdf")
    cat("Writing file", rainfall_out, "\n")
    ggsave(paste("plots/", rainfall_out, sep = ""), width = 15, height = 5)
  }
  p
}
