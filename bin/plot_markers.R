library(Seurat)
library(glue)
library(ggplot2)
library(ggthemes)
library(scCustomize)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

parser <- ArgumentParser(prog="plot_markers.R", description="fetch meta data from others annotation")
parser$add_argument("data", help="rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()
output = args$output
print(args)

adata = readRDS(args$data) 
print(head(adata@meta.data))

pdf(glue("{output}_dimplot.pdf"), width=19, height=8)
DimPlot_scCustom(adata, group.by='orig.ident', split.by='tech', reduction='harmony')
DimPlot_scCustom(adata, group.by='fine_ids', split.by='tech', reduction='humap')
DimPlot_scCustom(adata, group.by='seurat_clusters', reduction='humap')
dev.off()

pdf(glue("{output}_markers.pdf"), width=18, height=8)
FeaturePlot_scCustom(seurat_object = adata, features = "CD4",
                     colors_use = viridis_plasma_dark_high, na_color = "lightgray")
dev.off()
