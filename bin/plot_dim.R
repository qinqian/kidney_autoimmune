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

adata@meta.data$cell_type = gsub("kidney |, human", "", adata@meta.data$cell_type)

pdf(glue("fig/{output}_dimplot1.pdf"), width=18, height=10)
DimPlot_scCustom(adata, group.by='cell_type', label=T, split.by='tech', repel=T) & guides(colour = guide_legend(ncol=5, override.aes = list(size=5))) & theme(legend.title=element_blank(), legend.position='bottom', legend.text=element_text(size=13))
dev.off()

pdf(glue("fig/{output}_dimplot2.pdf"), width=18, height=10)
DimPlot_scCustom(adata, group.by='celltype', label=T, split.by='tech', repel=T) & guides(colour = guide_legend(ncol=9, override.aes = list(size=5))) & theme(legend.title=element_blank(), legend.position='bottom', legend.text=element_text(size=13))
dev.off()

pdf(glue("fig/{output}_markers.pdf"), width=18, height=8)
FeaturePlot_scCustom(seurat_object = adata, features = "PRG4",
                     colors_use = viridis_plasma_dark_high, na_color = "lightgray")
dev.off()
