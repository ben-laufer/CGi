# CGi
#### Takes output from [DMRichR](https://github.com/ben-laufer/DMRichR) and performs statistical testing and visualization for CpG annotation enrichments (CpG Islands, CpG Shores, CpG Shelves, and Open Sea).

## Installation

You can install GAT through [conda](https://bioconda.github.io/user/install.html#install-conda): `conda install -c bioconda gat`

## Enrichment Testing

The working directory for the [shell script](GAT_CpG_hyper_hypo.sh) requires:
1) A bed file of DMRs (DMRs.bed) from DMRichR (Extra > GAT).
2) A bed file of background regions (background.bed) from DMRichR (Extra > GAT).
3) A bed file of CpG annotations (hg38_CpG.bed) from DMRichR (Extra > GAT).
4) A isochore file that allows you to correct for GC content. Isochore files can be [downloaded](http://bioinfo2.ugr.es:8080/isochores/maps/) (and some will need to be lifted over), then edited for a [special bed format](https://github.com/AndreasHeger/gat/issues/4) (chr, start, stop, iso), and then sorted. 

Note: All bed files need to have unix and not mac or windows line breaks.

## Data Visualization

Then the enrichment results can be visualized using the [R script](GAT_CpG_hyper_hypo_plot.R).
