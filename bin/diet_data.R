library(Seurat)
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
print(args)

output = args$output
adata = readRDS(args$data) 

print(adata)
# eat up memories...
adata = JoinLayers(adata)
print(adata)
adata = DietSeurat(adata, layers=c("counts", "data"))

saveRDS(adata, glue("norm/{args$output}_alex_diet.rds"))
