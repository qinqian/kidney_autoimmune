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

test = obj.merge@meta.data %>% filter(tech=='xenium') %>% janitor::clean_names()
test$condition = str_trim(meta[match(str_extract(as.character(subset(obj.merge, subset=tech=='xenium')$orig.ident), "__(BS.*A[1,2])__(2.+)", group=1), meta$slide_id), 'case_ctrl'])

sample_ids = test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
    mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(sample_ids) %>% pull()
sample_mat = as.matrix(test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
    mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(-sample_ids))
rownames(sample_mat) = str_extract(as.character(sample_ids), "__(BS.*A[1,2])__(2.+)", group=1)

my_sample_col <- data.frame(sample = str_trim(meta$case_ctrl))
row.names(my_sample_col) <- meta$slide_id
my_sample_col = my_sample_col[rownames(sample_mat), ,drop=F]
my_sample_col = my_sample_col %>% mutate(patient=seq(1, 8))

write_tsv(my_sample_col, "my_sample_col.tsv")


sample_mat =  t(sample_mat)
scaled_mat = t(scale(t(sample_mat)))
colnames(scaled_mat) = my_sample_col$patient


my_sample_col = my_sample_col %>% rownames_to_column()

# # Define color function
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5))
ha1 = rowAnnotation(samples = my_sample_col$sample, 
                    col=list(samples=c('Case'='white', 'Control'='black')),
#                    annotation_label_location = "top",  # move label to top
                    annotation_name_gp = gpar(fontsize = 5),  
                    annotation_legend_param = list(
                    samples = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5), direction = "horizontal")))


rownames(scaled_mat)[rownames(scaled_mat)=='Immune Cell'] = 'Immune(Low Q)'

# Main heatmap
ht <- Heatmap(
  t(scaled_mat),
  name = "Expression",
  col = col_fun,
  ## annotation_name_gp = gpar(fontsize = 5),    # annotation label font size
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    direction = "horizontal"
  ),cluster_columns = FALSE,
  column_names_side = "top",
  right_annotation = ha1)

pdf("Fig1c.pdf", height=2.4, width=3)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom")
dev.off()



obj.merge@meta.data$cell_label = gsub(" Cell", "", obj.merge@meta.data$lennard_label)
obj.merge@meta.data = obj.merge@meta.data %>% mutate(cell_label = ifelse(cell_label=='Immune', 'Immune (LowQ)', cell_label))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

pdf("Fig1b.pdf", height=5, width=6)
DimPlot_scCustom(obj.merge, group.by="cell_label", label=T, repel=T, seed=99, label.box=T)+ggplot2::theme(legend.position = "none",         axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")
dev.off()


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
orig.merge.xen@meta.data$case_ctrl_num = as.numeric(factor(orig.merge.xen@meta.data$case_ctrl))

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
pdf("Fig1d.pdf", height=5, width=6)
FeaturePlot_scCustom(obj.cna, features = c('cna_ncorrs_fdr10'))[[1]] + 
    scale_color_gradient2_tableau() + 
    labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation')+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("UMAP1")+ylab("UMAP2")+ggplot2::theme(legend.position = "right", axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))
dev.off()

print(my_sample_col)
write_tsv(my_sample_col, "my_sample_col.tsv")
