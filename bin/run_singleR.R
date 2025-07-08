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

args = parser$parse_args()

output = args$output
input_data = args$data

if (length(input_data) < 2) {
    stop("input rds not equal to 2")
}

print(input_data)
myles.skin = readRDS(input_data[1])
#xenium
sc = readRDS(input_data[2])

overlap.genes = intersect(rownames(myles.skin), rownames(sc))
print(overlap.genes)
print(length(overlap.genes))

myles.skin.sub = myles.skin[overlap.genes, ]
myles.skin.sub@meta.data = myles.skin.sub@meta.data %>% mutate(cluster_name=case_when(
                                                                 grepl('M-0:|M-1:', cluster_name) ~ 'M-01:MERTK+SELENOP+',
                                                                 .default=cluster_name))

print(table(myles.skin.sub@meta.data$cluster_name))
sc.sub = subset(sc, features=overlap.genes)

sc.sub@meta.data$tech = "xenium"
myles.skin.sub@meta.data$tech = "sc"

spatial <- CreateSeuratObject(sc.sub@assays$RNA, meta.data = sc.sub@meta.data, project = "spatial")

##sc.data = exp(LayerData(myles.skin.sub, layer='data'))
##xen.data = exp(LayerData(sc.sub, layer='data'))
##
##print('--------------')
##print(fivenum(colMeans(sc.data)))
##sc.data = t(sc.data)
##sc.data = t(sc.data * median(sc.sub$nCount_RNA))
##
### create an assay using only normalized data
##assay.v5 <- CreateAssay5Object(data = sc.data)
##
### create a Seurat object based on this assay
##myles.skin.sub[['RNA']] <- assay.v5
##print(myles.skin.sub)
##
##print(fivenum(colMeans(sc.data)))
##print(fivenum(colMeans(xen.data)))
##print(dim(xen.data))
##print('aaaaaaaaaaaaaa')

merged.obj <- merge(spatial, myles.skin.sub, merge.data=T)
print(merged.obj)
merged.obj <- JoinLayers(merged.obj)
print(merged.obj)

split.merged = SplitObject(merged.obj, split.by="tech")
query.set <- as.SingleCellExperiment(split.merged$xenium)
ref.set <- as.SingleCellExperiment(split.merged$sc)
pred.grun2 <- SingleR(test = query.set, ref = ref.set, labels = ref.set$cluster_name, num.threads=24, de.method="t")
print(table(pred.grun2$labels))
saveRDS(pred.grun2, glue("{args$output}_singleR.rds"))

png(glue("{args$output}_diagnostics.png"), width=1600, height=1000)
print(plotScoreHeatmap(pred.grun2))
dev.off()

png(glue("{args$output}_delta.png"), width=1000, height=800)
plotDeltaDistribution(pred.grun2)
dev.off()

#merged.obj@meta.data = merged.obj@meta.data %>%
#    mutate(tissue=ifelse(is.na(tissue.type), "xenium", tissue.type))

print(sum(is.na(pred.grun2[match(rownames(merged.obj@meta.data), rownames(pred.grun2)),'pruned.labels'])))
merged.obj@meta.data$fine_ids = pred.grun2[match(rownames(merged.obj@meta.data), rownames(pred.grun2)),'pruned.labels']
merged.obj@meta.data = merged.obj@meta.data %>%
    mutate(fine_ids = ifelse(is.na(fine_ids), cluster_name, fine_ids))
saveRDS(merged.obj, glue("{args$output}_mergedsc.rds"))

