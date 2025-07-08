set.seed(999)
suppressPackageStartupMessages({
    library(glue)
    library(tessera)
    library(scCustomize)
    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(sf)
    ## Plotting functions 
    ## Not imported by Tessera
    library(ggplot2)
    library(ggthemes)
    library(viridis)
    library(patchwork)
    library(harmony)
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(presto)
    library(tibble)
    library(scDotPlot)
    library(SingleR)    
    library(Seurat)
    library(presto)
    library(dplyr)
    library(ggplot2)
    library(argparse)
    library(tidyverse)
})


parser <- ArgumentParser(prog="run_annotation.r", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--celltype", help="output prefix", default='cell_type')
parser$add_argument("--method", help="output prefix", default='noncluster')
parser$add_argument("--resolution", help="output prefix", default=0, type="double")

args = parser$parse_args()

output = args$output
input_data = args$data
celltype = args$celltype
method = args$method
resolution = args$resolution

if (length(input_data) < 2) {
    stop("input rds not equal to 2")
}

print(input_data)
#sc
myles.skin = readRDS(input_data[1])
#xenium
sc = readRDS(input_data[2])

overlap.genes = intersect(rownames(myles.skin), rownames(sc))
print(length(overlap.genes))

myles.skin.sub = myles.skin[overlap.genes, ]
sc.sub = subset(sc, features=overlap.genes)

rm(myles.skin)
rm(sc)

sc.sub@meta.data$tech = "xenium"
myles.skin.sub@meta.data$tech = "sc"

spatial <- CreateSeuratObject(sc.sub@assays$RNA, meta.data = sc.sub@meta.data, project = "spatial")

query.set <- as.SingleCellExperiment(spatial)
ref.set   <- as.SingleCellExperiment(myles.skin.sub)
col_data <- colData(ref.set)

print(table(spatial@meta.data$seurat_clusters))
print(sc.sub)

if (resolution > 0) {
   spatial = NormalizeData(spatial, scale.factor=median(spatial$nCount_RNA), normalization.method = "LogNormalize") %>% FindVariableFeatures(nfeatures = 1000) %>% ScaleData() %>% RunPCA(npcs=30) %>% 
    FindNeighbors() %>% FindClusters(resolution=resolution, cluster.name='seurat_clusters', random.seed=999)
}

print(table(spatial@meta.data$seurat_clusters))

if (method == 'cluster') {
    pred.grun2 <- SingleR(test = query.set, ref = ref.set, labels = col_data[, celltype], num.threads=16, de.method="t", method='cluster',  clusters=spatial$seurat_clusters)
} else {
    pred.grun2 <- SingleR(test = query.set, ref = ref.set, labels = col_data[, celltype], num.threads=14, de.method="wilcox", de.n=15)
}

print(table(pred.grun2$labels))
saveRDS(pred.grun2, glue("{args$output}_singleR.rds"))


png(glue("{args$output}_diagnostics.png"), width=1600, height=1000)
print(plotScoreHeatmap(pred.grun2))
dev.off()

png(glue("{args$output}_delta.png"), width=1000, height=800)
plotDeltaDistribution(pred.grun2)
dev.off()


if (method == 'cluster') {
     spatial@meta.data$fine_ids = pred.grun2[match(spatial@meta.data$seurat_clusters, rownames(pred.grun2)),'pruned.labels']
     saveRDS(spatial, glue("{args$output}_xenium_renorm_singleR.rds"))

     merged.obj <- merge(spatial, myles.skin.sub, merge.data=T)
     merged.obj <- JoinLayers(merged.obj)
} else {
    merged.obj <- merge(spatial, myles.skin.sub, merge.data=T)
    merged.obj <- JoinLayers(merged.obj)
    merged.obj@meta.data$fine_ids = pred.grun2[match(rownames(merged.obj@meta.data), rownames(pred.grun2)),'pruned.labels']
}

print(table(merged.obj@meta.data$fine_ids))
print(sum(is.na(merged.obj@meta.data$fine_ids)))

# fill the reference single cell annotation using the same column
merged.obj@meta.data = merged.obj@meta.data %>%
    mutate(fine_ids = ifelse(is.na(fine_ids), !!sym(celltype), fine_ids))
saveRDS(merged.obj, glue("{args$output}_mergedsc_short.rds"))

