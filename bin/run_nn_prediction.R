set.seed(999)
suppressPackageStartupMessages({
    library(glue)
    library(tessera)
    library(scCustomize)
    ## Downstream analysis in Seurat V5
    library(Seurat)
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
    library(presto)
    library(tibble)
    library(scDotPlot)
    library(SingleR)    
    library(Seurat)
    library(presto)
    library(dplyr)
    library(ggplot2)
    library(argparse)
    library(tidyverse)
})


parser <- ArgumentParser(prog="run_annotation.r", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()

method = args$method
output = args$output
input_data = args$data

#"../xen_seg/res2_short_KPMP_integrate_singlet_merged_annotated_sc.rds")
xen.int <- readRDS(input_data)
neigh = xen.int[["humap_fgraph"]]

neigh.xen = neigh[which(xen.int$tech == 'xenium'),]
neigh.xen = neigh.xen[, which(xen.int$tech == 'sc')]

y_ = xen.int$celltype[which(xen.int$tech == 'sc')]
X = Matrix::sparse.model.matrix(~0+y_)
X = neigh.xen %*% X %>% as.matrix()
X.out=apply(X, 1, which.max)
X.out = gsub("^y_", "", colnames(X)[X.out])
xen.int$celltype[is.na(xen.int$celltype)] = X.out

y_ = xen.int$cell_type[which(xen.int$tech == 'sc')]
X = Matrix::sparse.model.matrix(~0+y_)
X = neigh.xen %*% X %>% as.matrix()
X.out=apply(X, 1, which.max)
X.out = gsub("^y_", "", colnames(X)[X.out])
xen.int$cell_type[is.na(xen.int$cell_type)] = X.out

pdf(glue("{args$output}_Rplots.pdf"), width=16, height=7)
DimPlot_scCustom(xen.int, group.by='celltype', split.by='tech')
DimPlot_scCustom(xen.int, group.by='cell_type', split.by='tech')
DimPlot_scCustom(xen.int, group.by='fine_ids', split.by='tech')
dev.off()

saveRDS(xen.int, glue("output/{args$output}_nn_labels.rds"))

