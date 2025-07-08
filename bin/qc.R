library(Seurat)
library(sf)
library(Matrix)
library(ggthemes)
library(dplyr)
library(glue)
library(ggplot2)
library(scCustomize)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

# run uwot umap
Run_uwot_umap <- function(SeuratObj, min_dist = 0.3, spread = 0.8){
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


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_clean(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


sc.normalize <- function(x, method="lognorm") {
    x <- NormalizeData(x, assay='RNA', 
                       normalization.method='LogNormalize', scale.factor=median(x@meta.data$nCount_Xenium))
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    print(head(x@meta.data))
    #x <- ScaleData(x, features=VariableFeatures(x))
    set.seed(9)

    #x <- RunPCA(x, npcs = 30, features=VariableFeatures(x), verbose = FALSE) %>%
    x <- RunGLMPCA(x, npcs = 30, features=VariableFeatures(x), verbose = FALSE, reduction.name='pca') %>%
         harmony::RunHarmony(assay.use="RNA", group.by.vars='orig.ident', plot_convergence = TRUE, max_iter = 10, early_stop = F)

    set.seed(99)

    x <- Run_uwot_umap(x, min_dist = 0.3, spread = 0.8)
    x <- FindClusters(x, graph.name = 'humap_fgraph', resolution = 0.8, verbose = TRUE)
    x
}


parser <- ArgumentParser(prog="normalize_xenium.R", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()

output = args$output
input_data = args$data

adata = readRDS(input_data)
print(head(adata@meta.data))
#adata@meta.data$sample = basename(adata@meta.data$sample)

png(glue("qc/{args$output}_orig_seg_qc.png"), width=1100, height=760)
Stacked_VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA"), group.by='orig.ident', pt.size=0, x_lab_rotate = TRUE)
dev.off()


