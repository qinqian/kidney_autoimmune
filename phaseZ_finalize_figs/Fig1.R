library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
suppressPackageStartupMessages({
    library(tessera)
    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(scCustomize)
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
    library(rcna)
})


fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_blank(), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");

## sopa.orig = readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")

## niche.merge <- readRDS("250711_niches.rds")


obj.merge = readRDS("250721_cells_annotated_lennard.rds")
orig.merge = readRDS("../phaseF_newpipeline/sopa_seg/output/all_KPMP_integrate_singlet_umap_umapnn_labels_umap.rds")

table(orig.merge@meta.data$tech)

meta = read.csv("~/shruti_meta_clean (3).csv")
input_meta = meta[,c('slide_id', 'age', 'sex', 'case_ctrl', 'ICPi',  'malignancy', 'eGFR_base')] %>% arrange(case_ctrl)
#tile_obj$condition = str_trim(meta[match(str_extract(as.character(obj.merge$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

obj.merge@meta.data$cell_label = gsub(" Cell", "", obj.merge@meta.data$lennard_label)
obj.merge@meta.data = obj.merge@meta.data %>% mutate(cell_label = ifelse(cell_label=='Immune', 'Immune (LowQ)', cell_label))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)


## obj.merge@meta.data  = obj.merge@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1)) %>% left_join(input_meta, by=c("sample_id" = "slide_id"))

## obj.merge@meta.data$case_ctrl_num = as.numeric(factor(obj.merge@meta.data$case_ctrl))

## obj.merge <- FindNeighbors(obj.merge, reduction = 'harmony')

cells_to_keep <- colnames(orig.merge)[orig.merge$tech=='xenium']

orig.merge.xen <- subset(orig.merge, cells = cells_to_keep)


# Main cells
main_cells <- colnames(orig.merge.xen)

# Assay cells
assay_cells <- colnames(orig.merge.xen@assays$RNA)  # or whichever assay you're using

# PCA cells
pca_cells <- rownames(orig.merge.xen@reductions$pca@cell.embeddings)

# Graph cells
graph_cells <- colnames(orig.merge.xen[['humap_fgraph']])  # adapt if using another graph

# Active identity names
ident_cells <- names(Idents(orig.merge.xen))

# Check mismatches
length(setdiff(assay_cells, main_cells))
length(setdiff(pca_cells, main_cells))
length(setdiff(graph_cells, main_cells))
length(setdiff(ident_cells, main_cells))


orig.merge.xen@meta.data  = orig.merge.xen@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1)) 
orig.merge.xen@meta.data  = orig.merge.xen@meta.data %>% mutate(case_ctrl=input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl'])
orig.merge.xen@meta.data$case_ctrl_num = as.numeric(factor(str_trim(orig.merge.xen@meta.data$case_ctrl), levels=c("Control", "Case")))

print('---------')
main_cells <- colnames(orig.merge.xen)

# Assay cells
assay_cells <- colnames(orig.merge.xen@assays$RNA)  # or whichever assay you're using

# PCA cells
pca_cells <- rownames(orig.merge.xen@reductions$pca@cell.embeddings)

# Graph cells
graph_cells <- colnames(orig.merge.xen[['humap_fgraph']])  # adapt if using another graph

# Active identity names
ident_cells <- names(Idents(orig.merge.xen))

# Check mismatches
length(setdiff(assay_cells, main_cells))
length(setdiff(pca_cells, main_cells))
length(setdiff(graph_cells, main_cells))
length(setdiff(ident_cells, main_cells))


library(glue)
obj.cna <- association.Seurat(
    ## seurat_object = obj.merge, # somehow lennard's object cannot work, perhaps version of seurat while saving
    seurat_object = orig.merge.xen,
    test_var = 'case_ctrl_num', 
    samplem_key = 'sample_id', 
    graph_use = 'humap_fgraph', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    ## covs = c("age", "sex", "ICPi") ## no covariates to include 
)

names(obj.cna@reductions$cna@misc)


## FeaturePlot_scCustom(obj.cna, features = c('cna_ncorrs'))[[1]] + 
##     scale_color_gradient2_tableau() + 
##     labs(
##         title = 'CNA disease association', color = 'Correlation',
##         subtitle = sprintf('global p=%0.3f', obj.cna@reductions$cna@misc$p)
##     ) + 
pdf("Fig1d.pdf", height=9.5, width=6)
p1=DimPlot_scCustom(obj.merge, group.by="cell_label", label=T, repel=T, seed=99, label.box=T)+ggplot2::theme(legend.position = "none",         axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")
p2=FeaturePlot_scCustom(obj.cna, features = c('cna_ncorrs_fdr10'))[[1]] + 
    scale_color_gradient2(high = "#de2d26", mid = "white", low = "#2c7fb8", midpoint = 0)+
    labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation')+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")+ggplot2::theme(legend.position = "right", axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))
print(p1/p2)
dev.off()

