library(sccomp)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggrastr)
library(grid)
library(ggsci)

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


xen.seg <- readRDS("../phaseF_newpipeline/xen_seg/comb/RA_orig_seg.rds")


## sopa.orig = readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")
library(ComplexHeatmap)
library(circlize)
library(scico)

sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")

xen.seg@meta.data %>% mutate(sample_id=str_extract(orig.ident, "__(BS\\d*[_-].*)__2024", group=1)) %>% count(sample_id)%>% left_join(sc.niche$obj@meta.data %>%
    mutate(sample_id=str_extract(sample_id, "(BS\\d*[_-].*)__2024", group=1)) %>% count(sample_id), by="sample_id") %>% mutate(group='total')%>%group_by(group) %>% summarise(xenium_seg_total=sum(n.x), sopa_seg_total=sum(n.y))

mean(colSums(xen.seg[['RNA']]$counts))
mean(colSums(sc.niche$obj[['RNA']]$counts))

## sum(colSums(xen.seg$obj[['RNA']]$counts))
## sum(colSums(sc.niche$obj[['RNA']]$counts))

lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
imm.niche <- readRDS("250711_niches.rds")

orig.baysor <- readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")

orig.baysor@meta.data <- orig.baysor@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)
lennard.subtype@meta.data <- lennard.subtype@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)
xy <- Embeddings(orig.baysor, 'spatial')[match(lennard.subtype@meta.data$uniq_id, orig.baysor@meta.data$uniq_id),]

sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(lennard_label=lennard.subtype@meta.data$lennard_label)
sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(tile_label = imm.niche@meta.data[match(sc.niche$obj@meta.data$tile_id, rownames(imm.niche@meta.data)), 'niche_label'])
sc.niche.heatmap<-sc.niche$obj@meta.data %>% count(tile_label, lennard_label) %>% pivot_wider(names_from=lennard_label, values_from=n, values_fill=0) %>% filter(!is.na(tile_label))
sc.niche.heatmap = as.data.frame(sc.niche.heatmap)
rownames(sc.niche.heatmap) = sc.niche.heatmap$tile_label
sc.niche.heatmap = sc.niche.heatmap[,-1]
sc.niche.heatmap <- scale(sc.niche.heatmap)
col_fun <- colorRamp2(c(-2, 0, 2), scico(3, palette = "vik"))  # "vik" is diverging


niche.merge <- readRDS("250711_niches.rds")
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
orig.merge.xen@meta.data  = orig.merge.xen@meta.data %>% mutate(case_ctrl=str_trim(input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl']))
orig.merge.xen@meta.data$case_ctrl_num = as.numeric(factor(str_trim(orig.merge.xen@meta.data$case_ctrl), levels=c("Control", "Case")))

obj.merge@meta.data  = obj.merge@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1)) 
obj.merge@meta.data  = obj.merge@meta.data %>% mutate(case_ctrl=str_trim(input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl']))
obj.merge@meta.data$case_ctrl_num = as.numeric(factor(str_trim(obj.merge@meta.data$case_ctrl), levels=c("Control", "Case")))

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

orig.merge.xen@meta.data$condition = factor(str_trim(orig.merge.xen@meta.data$case_ctrl), levels=c("Control", "Case"))
orig.merge.xen@meta.data$lennard_label = obj.merge@meta.data$lennard_label

library(sccomp)

res = orig.merge.xen |>
    sccomp_estimate(
      formula_composition = ~ condition,
      sample = "sample_id", cell_group = "lennard_label",
      cores = 2, verbose=FALSE
    )
res = res |>  sccomp_test()

niche.merge@meta.data$condition_num = as.numeric(factor(str_trim(niche.merge@meta.data$condition), levels=c("Control", "Case")))
cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
names(cols) <- unique(obj.merge@meta.data$cell_label)

library(ggsci)
library(circlize)
library(scico)

sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
imm.niche <- readRDS("250711_niches.rds")

sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(lennard_label=lennard.subtype@meta.data$lennard_label)
sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(tile_label = imm.niche@meta.data[match(sc.niche$obj@meta.data$tile_id, rownames(imm.niche@meta.data)), 'niche_label'])
sc.niche.heatmap<-sc.niche$obj@meta.data %>% count(tile_label, lennard_label) %>% pivot_wider(names_from=lennard_label, values_from=n, values_fill=0) %>% filter(!is.na(tile_label))
sc.niche.heatmap = as.data.frame(sc.niche.heatmap)
rownames(sc.niche.heatmap) = sc.niche.heatmap$tile_label
sc.niche.heatmap = sc.niche.heatmap[,-1]
sc.niche.heatmap <- scale(sc.niche.heatmap)
col_fun <- colorRamp2(c(-2, 0, 2), scico(3, palette = "vik"))  # "vik" is diverging

sc.niche.heatmap <- sc.niche.heatmap[rownames(sc.niche.heatmap)!="Skeletal Muscle",]

niche_cols = pal_npg("nrc")(9)
names(niche_cols) <- rownames(sc.niche.heatmap)


p.bar2<-niche.merge@meta.data %>% count(condition, niche_label) %>% group_by(condition) %>% mutate(ratio=n/sum(n)) %>%
    ggplot() + geom_bar(aes(x=condition, y=ratio, fill=niche_label), stat='identity') + scale_fill_manual(values=niche_cols) + theme_classic()+  theme(
    legend.title = element_text(size = 6),  # Legend title size
    legend.text  = element_text(size = 5),  # Legend labels size
    legend.key.size = unit(0.4, "cm"),
    legend.position='right',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1.05),
  ) +xlab("")

res = orig.merge.xen |>
    sccomp_estimate(
      formula_composition = ~ condition,
      sample = "sample_id", cell_group = "lennard_label",
      cores = 2, verbose=FALSE
    )
res = res |>  sccomp_test()

case1 <- "BS21-N65682A2"
case2 <- "BS23_52206A2"
cont1 <- "BS22_12012A1"
cont2 <- "BS2_61615A1"

library(scCustomize)
cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
names(cols) <- unique(lennard.subtype@meta.data$lennard_label)

pdf("SuppFig1_part2.pdf", height=6, width=6)
p <- (res |>  plot_1D_intervals()) + theme_bw(base_size=8) + theme(legend.position="bottom") + ylab("Cell type")
print(p|p.bar2)
dev.off()
