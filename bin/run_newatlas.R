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

parser <- ArgumentParser(prog="run_annotation.r", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--celltypeA", help="sc which celltype to integrate")
parser$add_argument("--celltypeB", help="xenium which celltype to integrate")
parser$add_argument("--label", help="cell type label for output")

args = parser$parse_args()

celltypeA = args$celltypeA
celltypeB = args$celltypeB
output = args$output
input_data = args$data
label = args$label


newatlas <- readRDS(input_data[1])
xen <- readRDS(input_data[2])$seur
print(table(xen$coarse_ids))

xen.T <- subset(xen, subset = grepl(celltypeB, coarse_ids, ignore.case=TRUE))
newatlas.T = subset(newatlas,  subset = grepl(celltypeA, cell_type, ignore.case=TRUE))


getDE <- function(x, cluster='cell_type', top_n=12) {
    de.celltype <- presto::wilcoxauc(x, group_by=cluster) %>% dplyr::mutate(pct_diff = pct_in - pct_out)
    features <- de.celltype %>% 
         dplyr::group_by(group) %>%
         dplyr::filter(padj < 0.05 & auc > 0.5 & pct_diff > 0.05 & logFC > 0) %>%
         dplyr::arrange(group, dplyr::desc(logFC)) |>
         dplyr::slice_max(logFC, n = top_n) |>
         dplyr::select(group, feature) |>
         dplyr::ungroup() |>
         dplyr::arrange(feature) |>
         dplyr::distinct(feature, .keep_all = TRUE) |>
         dplyr::mutate(gene = feature) %>% dplyr::select(group, gene) %>% 
         dplyr::group_by(group) %>% 
         dplyr::select(group, gene)
    features <- features |> 
        tibble::deframe()
    features
}

get_theme <- function(size=12, angle=45) {
   defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal")
   defined_theme
}

sc.de <- getDE(newatlas.T, "celltype", top_n=500)

markers <- intersect(rownames(xen.T), sc.de)

xen.T.sub = xen.T[markers, ]
newatlas.T.sub = newatlas.T[markers, ]

xen.T.sub$tech = 'xenium'
newatlas.T.sub$tech = 'sc'

xen.meta = xen.T.sub@meta.data %>% select(orig.ident, tech) %>% mutate(celltype='unknown')
sc.meta  = newatlas.T.sub@meta.data %>% select(orig.ident, tech, celltype)

merged.T = CreateSeuratObject(cbind(LayerData(xen.T.sub[['RNA']], 'counts'), newatlas.T.sub[['RNA']]@counts),
                                    meta.data=bind_rows(xen.meta, sc.meta))

VariableFeatures(merged.T) = markers

xen.mergenormalize <- function(x, scale.factor) {
    x <- NormalizeData(x, assay='RNA', scale.factor=scale.factor,
                       normalization.method = "LogNormalize")
    x <- ScaleData(x)
    set.seed(99)
    x <- RunPCA(x, npcs = 30, features=VariableFeatures(x), verbose = FALSE) %>%
         harmony::RunHarmony(assay.use="RNA", group.by.vars=c('tech', 'orig.ident'), dims=1:30,
                             nclust=10, early_stop=F, ncores=4, max_iter=20) %>% 
         FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters(resolution=0.8) %>% 
         RunUMAP(reduction = "harmony", dims = 1:30)
    x
}

merged.scale.factor = mean(median(xen.T.sub$nCount_RNA),
                           median(newatlas.T.sub$nCount_RNA))

median.xen = median(colSums(LayerData(xen.T.sub[['RNA']], "counts")))
median.sc = median(colSums(LayerData(newatlas.T.sub[['RNA']], "counts")))

merged.scale.factor = mean(median.xen, median.sc)

merged.T.proc <- xen.mergenormalize(merged.T, scale.factor=merged.scale.factor)

pdf(glue("kpmp/{output}_{label}_Rplots.pdf"), width=18, height=9)
DimPlot(merged.T.proc, split.by='tech', group.by='celltype', reduction='harmony')+get_theme(size=11)
DimPlot(merged.T.proc, split.by='tech', group.by='celltype', reduction='umap')+get_theme(size=11)
scCustomize::DimPlot_scCustom(merged.T.proc, reduction='umap')+get_theme(size=11)
dev.off()

de.integrate = getDE(subset(merged.T.proc, subset=tech == 'sc'), 'seurat_clusters')

pdf(glue("kpmp/{output}_{label}_Rplots_integrate_clusters_scxen.pdf"), width=13)
DotPlot(subset(merged.T.proc, subset=tech == 'sc'), de.integrate)+get_theme(size=11)
DotPlot(subset(merged.T.proc, subset=tech == 'xenium'), de.integrate)+get_theme(size=11)
dev.off()

saveRDS(merged.T.proc, glue("kpmp/{output}_{label}_newatlas_kidney.rds"))

