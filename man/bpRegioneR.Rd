% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bpRegioner.R
\name{bpRegioneR}
\alias{bpRegioneR}
\title{Run regioneR for to regions in bed format}
\usage{
bpRegioneR(
  ...,
  regionA = system.file("extdata", "notch_breakpoint_regions_500_mappable.bed", package
    = "svBreaks"),
  regionB = system.file("extdata", "motif_1.mappable.bed", package = "svBreaks"),
  mappable_regions = system.file("extdata", "dmel6_mappable.bed", package = "svBreaks"),
  chrom_lengths = system.file("extdata", "chrom.sizes.txt", package = "svBreaks"),
  slop,
  limit2chrom = TRUE,
  chrom = "X",
  from,
  to,
  n = 100,
  plot = TRUE,
  bins = "auto",
  region = NULL
)
}
\arguments{
\item{regionA}{RegionA [Required]}

\item{regionB}{RegionB [Required]}

\item{mappable_regions}{Bed file with mappable geneome [Required]}

\item{chrom_lengths}{Lengths for for chromosomes}

\item{slop}{Add this much slop to each region}

\item{limit2chrom}{Restrict permutaions to one chromsome [Default: TRUE]}

\item{chrom}{Chromosome to run permutations on if limit2chrom=TRUE [Default: "X"]}

\item{from}{Restrict to region between from and to}

\item{to}{Restrict to region between from and to}

\item{n}{Number of permutations to run}

\item{plot}{Return plot or not? [Default: TRUE]}

\item{region}{Limit to region (chr:start-stop) [Default: NULL]}
}
\description{
Run regioneR for to regions in bed format
}
