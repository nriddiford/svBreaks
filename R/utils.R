# A collection of theme-related functions

#' Set theme for clean plotting with ggplot2
#' @keywords theme
#' @import ggplot2
#' @export
#'
cleanTheme <- function(base_size = 12) {
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(2, "lines"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 30),
    plot.margin = unit(1:4, "line")
  )
}

## @knitr slideTheme
#'
#' Set theme for slide plotting with ggplot2
#' @keywords theme
#' @import ggplot2
#' @export
#'
slideTheme <- function(base_size = 25) {
  theme(
    plot.title = element_text(hjust = 0.5, size = 50),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(2, "lines"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 50),
    strip.text = element_text(size = 25),
    plot.margin = unit(1:4, "line")
  )
}


## @knitr setCols
#'
#' Set colour scheme for manual colouring
#' @keywords colour scheme
#' @import ggplot2
#' @export
#' @return colourscale object
#'
setCols <- function(df, col, fill="Y", set="Pastel2") {
  names <- levels(as.factor(df[[col]]))
  names <- sort(names)
  cat("Setting colour levels:", names, "\n")
  level_number <- length(names)
  mycols <- brewer.pal(level_number, set)
  names(mycols) <- names
  fillScale <- scale_fill_manual(name = col, values = mycols)
  colScale <- scale_colour_manual(name = col, values = mycols)

  if (fill == "Y") return(fillScale)
  if (fill == "N") return(colScale)
}
