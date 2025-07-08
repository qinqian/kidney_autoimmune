library(Seurat)
library(glue)
library(ggplot2)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)
library(harmony)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(sf)
library(SingleR)
library(patchwork)
##library(ggdendro)
library(sceasy)
library(cowplot)
library(purrr)
library(DelayedArray)
library(CHOIR)


parser <- ArgumentParser(prog="subset_cells", description="")
parser$add_argument("data", help="input xenium directory or rds")
parser$add_argument("--cell", help="cell type")
parser$add_argument("--output", help="output prefix")


args = parser$parse_args()

ref = readRDS(args$data)
print(head(ref@meta.data))

print(dim(ref))

print(colnames(ref@meta.data))
print(table(ref@meta.data$cell_type))
print(args)

#sub.ref = ref[, ref@meta.data$cell_type==args$cell]
sub.ref = ref
print(dim(sub.ref))

sub.ref =
    sub.ref %>%
    NormalizeData(normalization.method = 'LogNormalize', verbose = F, scale.factor = 10000) %>% 
    ScaleData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

#sub.ref <- sub.ref %>%
#    RunPCA(features = VariableFeatures(sub.ref), npcs = 30, verbose = FALSE) %>%
#    harmony::RunHarmony("sample") %>%
#    FindNeighbors(reduction = "harmony") %>%
#    RunUMAP(reduction="harmony", dims=1:30, reduction.key='HUMAP_')

sub.ref = DietSeurat(sub.ref, layers=c("counts","data"))
saveRDS(sub.ref, glue("sub.ref.{args$cell}.rds"))

ref.set = as.SingleCellExperiment(sub.ref)
saveRDS(ref.set, glue("{args$output}.ref.rds"))


