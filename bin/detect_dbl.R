#library(DoubletFinder)
library(scater)
library(scuttle)
library(scran)
library(Seurat)
library(Matrix)
library(dplyr)
library(glue)
library(ggplot2)
#library(scCustomize)
library(argparse)
library(stringr)
#library(tidyverse)
#library(patchwork)
#library(purrr)
#library(dplyr)
#library(scuttle)
#library(scran)
#library(scater)
suppressPackageStartupMessages(library(scDblFinder))

set.seed(123)


parser <- ArgumentParser(prog="subset_cell_type.r", description="")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--group", help="group prefix")

args = parser$parse_args()
output = args$output

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu = readRDS(args$data)
print(head(seu@meta.data))
print(colnames(seu@meta.data))

if (c("seur") %in% names(seu)) {
    seu = seu$seur
}
filtered_cells <- seu@meta.data %>% mutate(ids=seq(1:n())) %>% group_by(orig.ident) %>% filter(
                                           nCount_Xenium >= 20 & nFeature_Xenium >= 20 # & nFeature_Xenium <= quantile(nFeature_Xenium, 0.95) & nCount_Xenium <= quantile(nCount_Xenium, 0.95)
                                           ) %>% ungroup()
filtered_cells = rownames(seu@meta.data)[filtered_cells$ids]
print(dim(seu))
seu = seu[,unlist(filtered_cells)]
#CCL5|LGALS3|SOX2-OT
print(grep("(CCL5|LGALS3|SOX2-OT)$", rownames(seu), value=T))
seu = seu[!grepl("(CCL5|LGALS3|SOX2-OT)$", rownames(seu)),]

seu = CreateSeuratObject(counts=seu@assays[['Xenium']]$counts, meta.data=seu@meta.data, project='RA', assay='RNA')
seu <- NormalizeData(seu) |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000)
seu <- ScaleData(seu) |> 
  RunPCA(npcs=10)

##  FindNeighbors() |>
##  FindClusters() |>
##  RunUMAP(dims = 1:10)
##
###DoubletFinder is too slow
##print(dim(seu))
##nExp_poi <- round(0.15*nrow(seu@meta.data))
##seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, sct = FALSE)
##seu@meta.data[, 'DF.class'] = seu@meta.data[, ncol(seu@meta.data)]
##print(table(seu@meta.data[, 'DF.class']))

print('------')
seu.sce <- as.SingleCellExperiment(seu)
print('******')
sce <- scDblFinder(seu.sce, dbr=0.1)
seu$scDblFinder.class = sce$scDblFinder.class
print(table(seu$scDblFinder.class))

sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n=1000)
set.seed(1002)
sce <- runPCA(sce, ncomponents=10, subset_row=hvgs)
sce <- runTSNE(sce, dimred="PCA")
scores <- computeDoubletDensity(sce, subset.row=hvgs)
print(scores)
seu$doublet_density = scores

png(glue("qc/{args$output}_orig_seg_adddbl.png"), width=800, height=800)
plotTSNE(sce, colour_by=I(log1p(scores)))
dev.off()

saveRDS(seu, glue("{args$group}_qc/{args$output}_orig_seg_adddbl.rds"))

