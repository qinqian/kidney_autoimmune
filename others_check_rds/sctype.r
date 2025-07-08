invisible(lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = T))


library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyplots)

kidn_seu_obj = readRDS("../data/miles/shapes_seur_kidney_v5.rds")

write_csv(kidn_seu_obj$seur@meta.data, file="miles_kidney_annotation.tsv")

table(kidn_seu_obj$seur$segmentation_method)

dim(kidn_seu_obj$seur$harmony)

dim(kidn_seu_obj$seur$umap)

pdf("kidney_umap.pdf", width=6, height=6)
DimPlot(kidn_seu_obj$seur, label=T, repel=TRUE)
dev.off()

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");
