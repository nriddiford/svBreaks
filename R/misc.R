
#' Misc functions

#' transform_types
#'
#' Convert short types into plottable names
#' @export
transform_types <- function(x){
  x <- x %>% 
    dplyr::mutate(type2 = ifelse(type2 == "TRA", "Translocation",
                                 ifelse(type2 == "DEL", "Deletion",
                                        ifelse(type2 == "DUP", "Duplication",
                                               ifelse(type2 == "COMPLEX", "Complex",
                                                      ifelse(type2 == "TANDUP", "Tandem\nDuplication",
                                                             ifelse(type2 == "BND", "Inversion", type2)))))))
  return(x)
}


#' SV colours
#' 
#' @export
sv_colours <- function(){
  return(
    c("DEL" = "#67C2EA",
      "COMPLEX" = '#85C385',
      "DUP" = '#3A7FEF',
      "BND" = '#F2E859',
      "TRA" = '#EE715A',
      "TANDUP" = '#E078EF',
      "NA" = "#333333")
    )
}


#' SNV colours
#' 
#' @export
snv_colours <- function(){
  return(
    c("delfrshift" = "#EB5767",
      "insfrshift" = '#F08573',
      "Missense" = '#F1AD32',
      "Nonsense" = '#F4C839',
      "Synonymous" = '#B2AF9D',
      "Essential_Splice" = '#F4C123')
    )
}


#' get sample names
#' @import dplyr
#' @export
getMissingSamples <- function(..., df=NULL, attach_info='~/Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt'){
  all_samples_info <- read.delim(attach_info, header = F, stringsAsFactors = T)
  colnames(all_samples_info) <- c("sample", "marius", "paper", "sex", "assay")
  
  names <- all_samples_info %>% 
    dplyr::filter(...) %>% 
    droplevels()
  
  names <- as.list(levels(names$sample))
  missing_samples = list()
  
  r <- 0
  for(i in 1:length(names)){
    if(!(names[[i]] %in% as.list(levels(df$sample)))){
      r <- r + 1
      missing_samples[r] <- names[[i]]
    }
  }
  # print("samples missing from dataframe levels: %s", unlist(missing_samples))
  return(missing_samples) 
}



#' Swap sample names
#' @import dplyr
#' @importFrom plyr join
#' @export
swapSampleNames <- function(df, attach_info='~/Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt'){
  name_conversion <- read.delim(attach_info, header=F)
  colnames(name_conversion) <- c("sample", "short", "sample_paper", "sex", "assay")
  name_conversion <- name_conversion %>% dplyr::select(-c("short", "sex", "assay"))
  new_df <- plyr::join(df, name_conversion, "sample", type = 'left') %>% 
    dplyr::rename(sample_old = sample,
                  sample = sample_paper) %>%
    dplyr::select(sample, everything()) %>% 
    droplevels()
  return(new_df)
}


#' svTypes
#'
#' Plot counts of different svTypes
#' @import dplyr
#' @importFrom plyr compact join
#' @export
svTypes <- function(..., bp_data=NULL, keep=FALSE, plot=TRUE, add_missing=TRUE, attach_info = '~/Desktop/script_test/mutationProfiles/data/samples_names_conversion.txt'){
  if(missing(bp_data)) bp_data <- getData(..., genotype=='somatic_tumour')
  ext <- ".png"
  
  bp_data <- bp_data %>% 
    dplyr::filter(...,
                  bp_no != "bp2") %>% 
    dplyr::group_by(sample, type2) %>% 
    dplyr::summarise(type_count = n()) %>%
    dplyr::mutate(sv_count = sum(type_count)) %>% 
    dplyr::arrange(-sv_count) %>% 
    droplevels() %>% 
    as.data.frame()
  
  if(add_missing){
    missing_samples <- getMissingSamples(..., df = bp_data, attach_info = attach_info)
    
    missing_samples <- plyr::compact(missing_samples)
    dat <- data.frame(sample = unlist(missing_samples), sv_count = 0, type_count =0, type2 = "NA")
    bp_data <- plyr::join(bp_data, dat, type='full')
  }
  
  order <- levels(fct_reorder(bp_data$sample, bp_data$sv_count))
  colours <- sv_colours()
  
  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(fct_reorder(sample, -sv_count), type_count, colour = type2, fill = type2), alpha=0.6, stat = 'identity')
  p <- p + scale_fill_manual("SV type\n", values = colours)
  p <- p + scale_colour_manual(values = colours)
  p <- p + guides(color = FALSE)
  p <- p + slideTheme() +
    theme(
      axis.title.y = element_blank(),
      # panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15),
      legend.position = "top",
      axis.text.y = element_text(size=15),
      axis.title.x=element_blank()
    )

  if(plot){
    p
  } else {
    return(list(bp_data, p))
  }
}





#' featureDensity
#'
#' Plot density of several feature tyopes across the genome
#' @import dplyr
#' @import ggplot2
#' @export

featureDensity <- function(feature_file1 = system.file("extdata", "g4_positions.txt", package="svBreaks"), feature1 = 'G4',
                           feature_file2 = NA, feature2 = NA,
                           feature_file3 = NA, feature3 = NA,
                           write=TRUE, chroms = c("2L", "2R", "3L", "3R", "X")) {
  n = 1

  file_list = list()

  f1 <- read.delim(feature_file1, header = F)
  f1 <- f1[,c(1,2,3)]
  if(is.null(f1$V3)){
    f1$V3 <- f1$V2 + 2
  }
  
  f1$type <- feature1
  colnames(f1) <- c("chrom", "start", "end", "type")
  
  f1 <- f1 %>% 
    dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
    dplyr::select(chrom, mid, type)

  file_list[[n]] <- f1

  if (!is.na(feature_file2)){
    n = n + 1
    f2 <- read.delim(feature_file2, header = F)
    f2 <- f2[,c(1,2,3)]
    if(is.null(f2$V3)){
      f2$V3 <- f2$V2 + 2
    }
    
    f2$type <- feature2
    colnames(f2) <- c("chrom", "start", "end", "type")
    
    f2 <- f2 %>% 
      dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
      dplyr::select(chrom, mid, type)
    
    file_list[[n]]  <- f2
  }
  if (!is.na(feature_file3)){
    n = n + 1
    f3 <- read.delim(feature_file3, header = F)
    if(is.null(f3$V3)){
      f3$V3 <- f3$V2 + 2
    }
    
    f3$type <- feature3
    colnames(f3) <- c("chrom", "start", "end", "type")
    
    f3 <- f3 %>% 
      dplyr::mutate(mid = as.integer(((end+start)/2)+1)) %>%
      dplyr::select(chrom, mid, type)
    
    file_list[[n]]  <- f3
  }

  files <- do.call(rbind, file_list)
  rownames(files) <- NULL

  # chromosomes <- data.frame(chroms = c("2L", "2R", "3L", "3R", "X"), lengths = c(23513712, 25286936, 28110227, 32079331, 23542271))
                       
  locations <- files %>%
    dplyr::mutate(type = as.factor(type)) %>%
    dplyr::mutate(pos = as.numeric(mid/1000000)) %>%
    dplyr::filter(chrom %in% chroms) %>%
    dplyr::arrange(chrom, mid) %>%
    droplevels()

  p <- ggplot(locations)
  p <- p + geom_density(aes(pos, fill = type), alpha = 0.4, adjust=0.07)
  # p <- p + geom_rug(aes(pos, colour = type), sides = "tb", alpha = 0.05)
  p <- p + facet_wrap(~chrom~type, scales = "free_y", ncol=1)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, max(locations$pos), by = 5))
  p <- p + geom_rug(data=locations, aes(pos, colour = type), sides = "b")
  p <- p + cleanTheme() +
    theme(axis.text.y = element_blank())

  if(write){
    featurePlot <- paste("feature_dist.png")
    cat("Writing file", featurePlot, "\n")
    ggsave(paste("plots/", featurePlot, sep = ""), width = 20, height = 10)
  }
  p
}


#' genomeHits
#'
#' Plot distribution of brekpoints across the genome
#' @import ggplot2
#' @export
genomeHits <- function(bp_data) {
  p <- ggplot(bp_data)
  p <- p + geom_point(aes(bp / 1000000, sample, colour = sample, shape = type, size = 0.5), alpha = 0.7)
  p <- p + guides(color = FALSE, size = FALSE)
  p <- p + cleanTheme() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 20),
      strip.text.x = element_text(size = 15),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted")
    )
  p <- p + facet_wrap(~chrom, scales = "free_x", ncol = 2)
  p <- p + scale_x_continuous("Mbs", breaks = seq(0, 33, by = 1), limits = c(0, 33), expand = c(0.01, 0.01))

  sv_gen_dist <- paste("bp_gen.dist", ext, sep = "")
  cat("Writing file", sv_gen_dist, "\n")
  ggsave(paste("plots/", sv_gen_dist, sep = ""), width = 20, height = 10)

  p
}
