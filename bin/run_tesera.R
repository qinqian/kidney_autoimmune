suppressPackageStartupMessages({
    library(tessera)

    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(scCustomize)
    library(sf)
    ## Plotting functions 
    ## Not imported by Tessera
    library(ggplot2)
    library(ggthemes)
    library(viridis)
    library(patchwork)
    library(harmony)
    library(Seurat)
    library(dplyr)
    library(cowplot)    
    library(rcna)
})

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");
