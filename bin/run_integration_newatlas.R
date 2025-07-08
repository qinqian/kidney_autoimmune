library(Seurat)
set.seed(999)
    library(scCustomize)
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
parser$add_argument("--method", help="normalization method")

args = parser$parse_args()

method = args$method
output = args$output
input_data = args$data

## input_data = c("../results/phaseA2_orig_seg_seurat/sub.ref.rds", "../results/phaseA3_subtyping_mono_orig_seg/RA_xenium_adata.myeloid.rds")

print(input_data)

#amp2
if (length(input_data) == 3) {
   myles.skin = readRDS(input_data[1])
   sc = readRDS(input_data[2])
   merged.obj = readRDS(input_data[3])
   print('*********')
   print(dim(merged.obj))
} else if (length(input_data) == 2) {
   myles.skin = readRDS(input_data[1])
   sc = readRDS(input_data[2])
   print(head(sc@meta.data))
   print(head(myles.skin@meta.data))

   myles.skin$tech = 'sc'
   sc$tech = 'xenium'

   overlap.genes = intersect(rownames(myles.skin), rownames(sc))
   print(length(overlap.genes))
   
   myles.skin.sub = myles.skin[overlap.genes, ]
   sc.sub = subset(sc, features=overlap.genes)

   merged.obj = merge(myles.skin.sub, sc.sub)
   merged.obj <- JoinLayers(merged.obj)
} else {
    stop("input rds not right...")
}

print(table(merged.obj@meta.data$tech))
#[1]  5000 27927
#[1]  5000 27927
#    sc xenium
#25578  27927

overlap.genes = intersect(rownames(myles.skin), rownames(sc))
print(length(overlap.genes))

myles.skin.sub = myles.skin[overlap.genes, ]
sc.sub = subset(sc, features=overlap.genes)

sc.sub = FindVariableFeatures(sc.sub, nfeatures=2000)
myles.skin.sub = FindVariableFeatures(myles.skin.sub, nfeatures=2000)
median.count <- mean(median(myles.skin.sub$nCount_RNA), median(sc.sub$nCount_RNA))

VariableFeatures(merged.obj) <- intersect(VariableFeatures(sc.sub), VariableFeatures(myles.skin.sub))

merged.obj <- merged.obj %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = median.count, verbose = F) %>% 
    ScaleData(features=VariableFeatures(merged.obj)) %>% 
    RunPCA(verbose = T)

set.seed(1)
print(head(merged.obj@meta.data, 3))
print(colnames(merged.obj@meta.data))
#print(table(merged.obj@orig.ident))
#merged.seurat <- harmony::RunHarmony(merged.obj, c("tech", "orig.ident"), plot_convergence = TRUE, max_iter = 10, early_stop = F)
print(table(merged.obj@meta.data$tech))
print(head(merged.obj@meta.data$tech))
print(tail(merged.obj@meta.data$tech))
print(sum(is.na(merged.obj@meta.data$tech)))
merged.obj@meta.data$tech = factor(merged.obj@meta.data$tech)

print(head(merged.obj@meta.data))
print(colnames(merged.obj@meta.data))
merged.obj@meta.data = merged.obj@meta.data %>% mutate(sample_ids=ifelse(is.na(orig.ident), SampleID, orig.ident))

print(table(merged.obj@meta.data$sample_ids))

merged.seurat <- harmony::RunHarmony(merged.obj, c("tech", "sample_ids"), plot_convergence = TRUE, max_iter = 20, early_stop = F)

pdf(glue("{args$output}_check_pca_har.pdf"))
DimPlot(merged.seurat, reduction='harmony', group.by='tech')
dev.off()

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

print('------------')
print(merged.seurat)
print(dim(merged.seurat@reductions$harmony@cell.embeddings))
print(dim(merged.seurat))
print(merged.obj)
print(table(merged.obj@meta.data$tech))
merged.seurat <- Run_uwot_umap(merged.seurat, min_dist = 0.3, spread = 0.8)

#merged.seurat <- FindNeighbors(merged.seurat, reduction='harmony')
#set.seed(1)
#merged.seurat <- RunUMAP(merged.seurat, dims=1:30, seed.use=42, reduction='harmony')
#set.seed(1)
#merged.seurat <- FindClusters(merged.seurat, random.seed=42)


#pdf(glue("kpmp/{args$output}_annt_um_split.pdf"), width=18, height=12)
#p1=DimPlot_scCustom(merged.seurat, reduction='humap', group.by='fine_ids', split.by='tech')
#p2=FeaturePlot_scCustom(merged.seurat, features="CCL19", reduction='humap', split.by='tech', na_color = "lightgray", na_cutoff = 1e-06)
#p3=FeaturePlot_scCustom(merged.seurat, features="IFNG", reduction='humap', split.by='tech', na_color = "lightgray", na_cutoff = 1e-06)
#p4=FeaturePlot_scCustom(merged.seurat, features="TGFB1", reduction='humap', split.by='tech', na_color = "lightgray", na_cutoff = 1e-06)
#p5=FeaturePlot_scCustom(merged.seurat, features="TGFB3", reduction='humap', split.by='tech', na_color = "lightgray", na_cutoff = 1e-06)
#print(p1)
#print(p2)
#print(p3)
#print(p4)
#print(p5)
#dev.off()
saveRDS(merged.seurat, glue('kpmp/{args$output}_merged_annotated_sc.rds'))
#
#
#pdf(glue("kpmp/{args$output}_annt_um.pdf"), width=11, height=6)
#print(DimPlot_scCustom(merged.seurat, reduction='humap', group.by='fine_ids', split.by='tech'))
### print(Clustered_DotPlot(merged.seurat, features=c("C1S", "ADAM33", "DCN", "NFATC4", "CXCL12", "PLIN4", "PLIN1", "KRT2", "KRT5", "KRT14", "C3", "CLU", "VWF", "PECAM1", "CLU", "COL6A1", "COL5A1"), group.by="fine_ids"))
#dev.off()
