library(Seurat)
library(Matrix)
library(dplyr)
library(glue)
library(ggplot2)
library(scCustomize)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

# run uwot umap
Run_uwot_umap <- function(SeuratObj, min_dist = 0.3, spread = 0.8) {
    HU <- uwot::umap(SeuratObj@reductions$harmony@cell.embeddings, min_dist = min_dist, 
                 spread = spread, ret_extra = 'fgraph', fast_sgd = FALSE)
    colnames(HU$embedding) = c('HUMAP1', 'HUMAP2')
    rownames(HU$fgraph) = colnames(HU$fgraph) = Cells(SeuratObj)
    SeuratObj[['humap']] <- Seurat::CreateDimReducObject(
        embeddings = HU$embedding,
        assay = 'RNA',
        key = 'HUMAP_',
        global = TRUE
    )
    HU_graph <- Seurat::as.Graph(HU$fgraph)
    DefaultAssay(HU_graph) <- DefaultAssay(SeuratObj)
    SeuratObj[['humap_fgraph']] <- HU_graph
    return(SeuratObj)
}


sc.normalize <- function(x, method="lognorm") {
    x <- NormalizeData(x, assay='RNA', 
                       normalization.method='LogNormalize', scale.factor=median(x@meta.data$nCount_RNA))
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
    print(head(x@meta.data))
    x <- ScaleData(x, features=VariableFeatures(x))
    set.seed(9)

    x <- RunPCA(x, npcs = 30, features=VariableFeatures(x), verbose = FALSE) %>%
         harmony::RunHarmony(assay.use="RNA", group.by.vars='orig.ident', plot_convergence = TRUE, max_iter = 20, early_stop = F)

    set.seed(99)

    x <- Run_uwot_umap(x, min_dist = 0.3, spread = 0.8)
    x <- FindClusters(x, graph.name = 'humap_fgraph', resolution = 0.8, verbose = TRUE)
    x
}


parser <- ArgumentParser(prog="normalize_xenium.R", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--method", help="normalization method")

args = parser$parse_args()

method = args$method
output = args$output
input_data = args$data

adata = readRDS(input_data)
adata@meta.data$orig.ident = basename(as.character(adata@meta.data$sample))
print(adata)
cells.pre = adata@meta.data %>% group_by(orig.ident) %>% summarise(n = n(), m1=min(nFeature_RNA), m2=max(nFeature_RNA)) %>% ungroup()
filtered_cells <- adata@meta.data %>% mutate(ids=seq(1:n())) %>% group_by(orig.ident) %>% filter(
                                               #nFeature_RNA <= quantile(nFeature_RNA, 0.98) & nFeature_RNA >= 20
                                               nCount_RNA > 50 & nFeature_RNA > 50
                                               ) %>% ungroup()
cells.post = filtered_cells %>% group_by(orig.ident) %>% summarise(n = n(), m1=min(nFeature_RNA), m2=max(nFeature_RNA)) %>% ungroup()
filtered_cells = rownames(adata@meta.data)[filtered_cells$ids]

print(dim(adata))
adata = adata[,unlist(filtered_cells)]
print(dim(adata))
saveRDS(cells.pre, "cellpre.rds")
saveRDS(cells.post, "cellpost.rds")

adata = sc.normalize(adata)
saveRDS(adata, glue("norm/{args$output}_orig_seg_lognorm.rds"))

