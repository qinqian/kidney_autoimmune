library(Seurat)
library(glue)
library(ggplot2)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

parser <- ArgumentParser(prog="subset_cell_type.r", description="")
parser$add_argument("data", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()
output = args$output
print(args)

xen = LoadXenium(args$data)
saveRDS(xen@meta.data, glue("{output}_orig_seg_meta.tsv"))

xen <- subset(xen, subset=nCount_Xenium > 0)
xen@meta.data$orig.ident <- output
    
print(head(GetTissueCoordinates(xen, which='centroids')))
spatial_embed = GetTissueCoordinates(xen, which='centroids')[, -3]
x = spatial_embed[, 1]
y = spatial_embed[, 2]
xen@meta.data[, 'x'] = x
xen@meta.data[, 'y'] = y
xen@images[['fov']] <- NULL
xen[['BlankCodeword']] <- NULL
xen[['ControlCodeword']] <- NULL
xen[['ControlProbe']] <- NULL
xen[['GenomicControl']] <- NULL
##cannot delete default
##xen[['Xenium']] <- NULL
saveRDS(xen, glue("{output}_orig_seg.rds"))
