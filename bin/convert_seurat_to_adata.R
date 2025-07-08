library(Seurat)
library(sceasy)
library(glue)
library(ggplot2)
library(argparse)
library(stringr)

parser <- ArgumentParser(prog="convert_seurat_to_adata.R", description="convert seurat data to adata h5ad")
parser$add_argument("data", help="rds")
parser$add_argument("--output", help="output prefix")
args = parser$parse_args()

print(args)
adata = readRDS(args$data) 
print(adata)

use_condaenv("scconverter")
adata[["RNA"]] <- as(adata[["RNA"]], Class="Assay")

print(packageVersion("SeuratObject"))
print(packageVersion("Seurat"))
print(adata)

sceasy::convertFormat(adata, from="seurat", to="anndata", outFile=glue('{args$output}_orig_seg_lognorm.h5ad'),  main_layer="counts", transfer_layers = "data")

