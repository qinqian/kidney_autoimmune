library(Seurat)
library(ggplot2)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)
invisible(lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = T))
source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");

parser <- ArgumentParser(prog="load_view", description="")
parser$add_argument("data", help="input xenium directory")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--normalize", action="store_true", help="output prefix")

args = parser$parse_args()

##args = list()
##args$data = "data/kidney/20240803__182820__BWH_20240803_skin_Shruti_kidney/output-XETG00150__0018462__BS22_12012A1__20240803__183643"
##args$normalize = T

xen.obj <- LoadXenium(args$data)

xen.obj <- subset(xen.obj, subset=nCount_Xenium > 0)
xen.obj <- subset(xen.obj, subset=nFeature_Xenium > 20 & nFeature_Xenium <= 5000)

#library(scSHC)
#clusters <- scSHC(xen.obj@assays$Xenium$counts, alpha=0.25, cores=1)
#xen.obj[['shc25']] = clusters[[1]]
#
#clusters <- scSHC(xen.obj@assays$Xenium$counts, alpha=0.05, cores=1)
#xen.obj[['shc05']] = clusters[[1]]


xen.normalize <- function(x) {
    x <- SCTransform(x, assay = "Xenium")
    x <- RunPCA(x, npcs = 30, features = rownames(x))
    x <- FindNeighbors(x, reduction = "pca", dims = 1:30) %>% 
	    FindClusters(resolution = 1.2) %>%
	    FindClusters(resolution = 0.8) %>%
	    FindClusters(resolution = 1.0) %>%
            RunUMAP(dims = 1:30)
    x
}

if (args$normalize) {
    xen.obj <- xen.normalize(xen.obj)
}

xen.obj <- run_sctype(xen.obj, known_tissue_type="Kidney", slot="SCT")
saveRDS(xen.obj, paste0('rds/', args$output, '.rds'))
write_tsv(xen.obj@meta.data, paste0('rds/', args$output, '_meta.tsv'))

pdf(str_c(c('figure/', args$output, "spatial.pdf"), collapse=""))
ImageDimPlot(xen.obj,  group.by="sctype_classification", border.color = "white",size = 0.5, border.size = 0.05)
ImageFeaturePlot(xen.obj, features=c("NOTCH3", "POSTN"), size=0.5, cols=c("white", "red"))
dev.off()

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

pdf(str_c(c('figure/', args$output, "umap.pdf"), collapse=""), width=12)
p.umap = DimPlot(xen.obj)+
    guides(x = axis, y = axis)+
    theme(
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))),
        axis.title = element_text(hjust = 0))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")
p.umap2 = DimPlot(xen.obj, group.by="sctype_classification")+
    guides(x = axis, y = axis)+
    theme(
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))),
        axis.title = element_text(hjust = 0))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")
print(p.umap + p.umap2)
dev.off()


pdf(str_c(c('figure/', args$output, "markers.pdf"), collapse=""))
p.feat = FeaturePlot(xen.obj, features=c("NOTCH3", "POSTN", "COL1A1", "THY1", "COL6A1", "MALAT1"))
print(p.feat)
dev.off()

