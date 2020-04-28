#' Set theme for clean plotting with ggplot2
#' @keywords theme
#' @export
#'
cleanTheme <- function(base_size = 12) {
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(2, "lines"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 12),
    axis.title= element_text(size = 15),
    strip.text = element_text(size = 15),
    plot.margin = unit(1:4, "line")
  )
}


#' Set theme for slide plotting with ggplot2
#' @keywords theme
#' @export
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

# Modified for black background
#' @keywords theme
#' @export
blackTheme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, colour = 'white'),
    panel.background = element_blank(),
    plot.background = element_rect(colour = "black", fill = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="white", size = 0.5),
    axis.line.y = element_line(color="white", size = 0.5),
    axis.text = element_text(size=20, colour = 'white'),
    # axis.title = element_text(size=20, colour = 'white'),
    axis.ticks = element_line(color = "white"),
    axis.title = element_text(size = 20, color = "white"),
    legend.background = element_rect(color = NA, fill = "black"),
    legend.key = element_rect(color = "white",  fill = "black"),
    legend.text = element_text(size = base_size*0.8, color = "white"),
    legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
    strip.text = element_text(color = "white", size = 25),
    strip.background = element_rect(fill= "#383F4C")
    
  )
}

## setCols
#'
#' Set colour scheme for manual colouring
#' @keywords colour scheme
#' @import ggplot2
#' @export
#' @return colourscale object
#'
setCols <- function(df, col, fill="Y", set="Blues", customVals) {
  names <- levels(as.factor(df[[col]]))
  
  # levels(as.factor(df[with(df, order(length)),][[col]]))
  
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

