library(Seurat)
library(sf)
library(glue)
library(ggplot2)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

parser <- ArgumentParser(prog="fetch_meta_data.R", description="fetch meta data from others annotation")
parser$add_argument("data", help="rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()
output = args$output
print(args)

adata = readRDS(args$data) 
print(adata)
saveRDS(adata$seur@meta.data, file=glue("{output}_metadata.tsv"))

pdf(glue(output, "_qc", ".pdf"), width=15, height=5)
VlnPlot(adata$seur, features = c("nFeature_RNA", "nCount_RNA"), group.by='orig.ident', ncol = 3, pt.size=0)
dev.off()

