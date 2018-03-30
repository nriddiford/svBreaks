# Plotting functions for structural variant breakpoint data

This package contains a collection of functions to explore breakpoint data produced by [svParser](https://github.com/nriddiford/svParser#svparser) i.e. the `all_bps.txt` or `all_bps_filtered.bps` files. Input to svBreaks should contain the following information:

```
event bp_number sample  genotype  chrom1  bp1 gene1 feature1  chrom2  bp2 gene2 feature2 type length(kb)
1	bp_1	A373R11	germline_recurrent	2L	73837	galectin	exon_2	2L	73896	galectin	exon_2	DEL	0.1
```


This tool is under constant development. Please feel free to [contact me](mailto:nick.riddiford@curie.fr), or [raise an issue](https://github.com/nriddiford/svBreaks/issues) if you encounter any problems.


## Installation

Install from github:

```
git clone https://github.com/nriddiford/svBreaks.git
```

Start an R session, and install package:

```{R}
library(devtools)
install_github("nriddiford/svBreaks")
library(svBreaks)
```
