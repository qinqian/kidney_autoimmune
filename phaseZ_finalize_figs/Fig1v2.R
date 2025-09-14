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

## sopa.orig = readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")
library(ComplexHeatmap)
library(circlize)
library(scico)

sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
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
ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5))

sc.niche.heatmap <- sc.niche.heatmap[rownames(sc.niche.heatmap)!="Skeletal Muscle",]

niche_cols = pal_npg("nrc")(9)
names(niche_cols) <- rownames(sc.niche.heatmap)

ha1 = rowAnnotation(samples = rownames(sc.niche.heatmap),
                    col=list(samples=niche_cols),
                    annotation_name_gp = gpar(fontsize = 5),    show_legend = FALSE,
                    annotation_legend_param = list(
                    samples = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5), direction = "horizontal")))

ht <- Heatmap(
  sc.niche.heatmap,
  name = "Cell Composition Ratio",
  col = col_fun, width=6.5,
  height=2,
  column_names_gp = gpar(fontsize = 5.4),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5.5),
    labels_gp = gpar(fontsize = 5.5),
    direction = "horizontal"
  ), cluster_columns = T,
  column_names_side = "top", left_annotation = ha1)

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

# Check mismatches
length(setdiff(assay_cells, main_cells))
length(setdiff(pca_cells, main_cells))
length(setdiff(graph_cells, main_cells))
length(setdiff(ident_cells, main_cells))


library(glue)
obj.cna <- association.Seurat(
    seurat_object = orig.merge.xen,
    test_var = 'case_ctrl_num', 
    samplem_key = 'sample_id', 
    graph_use = 'humap_fgraph', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    ## covs = c("age", "sex", "ICPi") ## no covariates to include 
)


## obj.merge@meta.data$condition = factor(str_trim(obj.merge@meta.data$condition), levels=c("Control", "Case"))

## library(sccomp)

## res = orig.merge.xen |>
##     sccomp_estimate(
##       formula_composition = ~ condition,
##       sample = "sample_id", cell_group = "cell_",
##       cores = 2, verbose=FALSE
##     )

## res = res |>  sccomp_test()


niche.merge@meta.data$condition_num = as.numeric(factor(str_trim(niche.merge@meta.data$condition), levels=c("Control", "Case")))

niche.cna <- association.Seurat(
    seurat_object = niche.merge,
    test_var = 'condition_num',
    samplem_key = 'sample_id',
    graph_use = 'RNA_snn',
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    #covs = c("age", "sex", "ICPi") ## no covariates to include
)

cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
names(cols) <- unique(obj.merge@meta.data$cell_label)

p.bar1<-obj.merge@meta.data %>% count(case_ctrl, cell_label) %>% group_by(case_ctrl) %>% mutate(ratio=n/sum(n)) %>%
    ggplot() + geom_bar(aes(x=case_ctrl, y=ratio, fill=cell_label), stat='identity') + scale_fill_manual(values=cols) + theme_classic()+  theme(
    legend.title = element_text(size = 6),  # Legend title size
    legend.text  = element_text(size = 5),  # Legend labels size
    legend.key.size = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1.05),
  )+xlab("")

p.bar2<-niche.merge@meta.data %>% count(condition, niche_label) %>% group_by(condition) %>% mutate(ratio=n/sum(n)) %>%
    ggplot() + geom_bar(aes(x=condition, y=ratio, fill=niche_label), stat='identity') + scale_fill_manual(values=niche_cols) + theme_classic()+  theme(
    legend.title = element_text(size = 6),  # Legend title size
    legend.text  = element_text(size = 5),  # Legend labels size
    legend.key.size = unit(0.4, "cm"),
    legend.position='right',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1.05),
  ) +xlab("")


p11=DimPlot_scCustom(subset(obj.merge, subset = case_ctrl=='Case'), raster = TRUE,        # enable rasterization
  raster.dpi = c(350, 350), group.by="cell_label", seed=99)  +ggplot2::theme(legend.position = "none",         axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + NoLegend()+ggtitle("ICI-AIN")+xlab("")+ylab("")+scale_color_manual(values=cols)
p12=DimPlot_scCustom(subset(obj.merge, subset = case_ctrl=='Control'), raster = TRUE,        # enable rasterization
  raster.dpi = c(350, 350),  group.by="cell_label", seed=99) +ggplot2::theme(legend.position = "none") + NoLegend()+ggtitle("ICI-ATN")+scale_x_continuous(breaks = NULL)+ scale_y_continuous(breaks = NULL)+xlab("")+ylab("")+scale_color_manual(values=cols)

p21=DimPlot_scCustom(subset(subset(niche.merge, condition=='Case'), subset = niche_label != "Skeletal Muscle"),
                    raster = TRUE,        # enable rasterization
                    raster.dpi = c(350, 350),
                    group.by="niche_label", seed=99)+
    guides(x = axis, y = axis)+
    ggplot2::theme(legend.position = "none",
                   axis.line = element_line(arrow = arrow(type = "closed", length = unit(10, 'pt'))))+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + NoLegend() +
    xlab("UMAP1")+ylab("UMAP2") +  scale_color_manual(values=niche_cols) + ggtitle("")

p22=DimPlot_scCustom(subset(subset(niche.merge, condition=='Control'), subset = niche_label != "Skeletal Muscle"),
                     raster = TRUE,        # enable rasterization
                    raster.dpi = c(350, 350),
                    group.by="niche_label", seed=99)+ggplot2::theme(legend.position = "none",
    legend.title = element_text(size = 6),  # Legend title size
    legend.text  = element_text(size = 5),  # Legend labels size
    legend.key.size = unit(0.4, "cm"),
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(3, 'pt')))) + scale_color_manual(values=niche_cols) +xlab("")+ylab("")+    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL)+ggtitle("")+  guides(color = guide_legend(nrow = 3))

p2.niche=FeaturePlot_scCustom(subset(niche.cna, subset = niche_label != "Skeletal Muscle"), features = c('cna_ncorrs_fdr10'), raster=T)[[1]] +
    scale_color_gradient2(high = "#de2d26", mid = "white", low = "#2c7fb8", midpoint = 0)+
    xlab("")+ylab("") + ggplot2::theme(legend.position = "right")+scale_x_continuous(breaks = NULL)+ scale_y_continuous(breaks = NULL) + ggtitle("")
p3=FeaturePlot_scCustom(obj.cna, features = c('cna_ncorrs_fdr10'), raster=T)[[1]] + #
    scale_color_gradient2(high = "#de2d26", mid = "white", low = "#2c7fb8", midpoint = 0)+
    labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation')+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    xlab("")+ylab("")+ ggplot2::theme(legend.position = "right") #+ggtitle("")
## combined1 <- plot_spacer()+p1+p3+ht_plot + p2 + p2.niche + plot_spacer()+ plot_spacer()+plot_layout(ncol = 4, widths=c(4, 2, 2, 5))
## combined1 <- p11+p3+ht_plot + p2 + p2.niche + plot_spacer()+ plot_spacer()+plot_layout(ncol = 4, widths=c(4, 2, 2, 5))




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

common_xlim <- c(0, 5000)
common_ylim <- c(0, 5000)

case1 <- "BS21-N65682A2"
case2 <- "BS23_52206A2"
cont1 <- "BS22_12012A1"
cont2 <- "BS2_61615A1"

library(scCustomize)
cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
names(cols) <- unique(lennard.subtype@meta.data$lennard_label)

patient <- read_tsv("my_sample_col.tsv")
imm.niche@meta.data <- imm.niche@meta.data %>% mutate(patient_id=gsub("__2.*", "", sample_id))
lennard.subtype@meta.data = lennard.subtype@meta.data %>% mutate(patient_id=gsub("output-XETG00150__0018462__|output-XETG00392__0045655__", "", gsub("__2.*", "", sample_ids)))

pdf("Rplots.pdf", width=22, height=4.6)
feature = 'PODXL'
pl0 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, imm.niche$niche_label)%>% filter(patient_id==case1),
            aes(geometry = shape, fill = niche_label), color=NA) + theme_bw(base_size=12) +
    guides(fill = guide_legend(nrow = 3)) +
    coord_sf(expand = FALSE) + NULL+
    scale_fill_d3('category20c') + theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 6),      # legend labels
                                         legend.title = element_text(size = 7))+scale_fill_manual(values=niche_cols) #+theme_void()
pl1 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==case1),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(1700, 2100), ylim=c(900, 1200), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none') + scale_x_continuous(limits = c(1700, 2100)) + scale_y_continuous(limits = c(900, 1200))
feature = 'CD38'
pl2 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==case1),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(3000, 3500), ylim=c(1500, 1800), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none') + scale_x_continuous(limits = c(3000, 3500)) + scale_y_continuous(limits = c(1500, 1800))
pl01 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, imm.niche$niche_label)%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = niche_label), color=NA) + theme_bw(base_size=12) +
    guides(fill = guide_legend(nrow = 3)) +
    coord_sf(expand = FALSE) + NULL+
    scale_fill_d3('category20c') + theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 6),      # legend labels
                                         legend.title = element_text(size = 7))+scale_fill_manual(values=niche_cols) #+theme_void()
pl11 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==cont2),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(1600, 2300), ylim=c(1200, 1600), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none') + scale_x_continuous(limits = c(1600, 2300)) + scale_y_continuous(limits = c(1200, 1600))
feature = 'CD38'
pl21 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==cont2),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='bottom',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(6600, 7200), ylim=c(600, 900), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')  + scale_x_continuous(limits = c(6600, 7200)) + scale_y_continuous(limits = c(600, 900))
combined2 <- pl0+pl1+pl2+plot_spacer()+pl01+pl11+pl21+plot_layout(ncol=7, nrow=1, widths=c(1,1,1,0.01,1,1,1)) #&theme(plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()


pdf("Fig1d.pdf", height=6.4, width=11)
ht_plot <- wrap_elements(full=grid.grabExpr(draw(ht,
                              merge_legend = TRUE, heatmap_legend_side = "bottom",
                              annotation_legend_list = NULL)))
## p_grob  <- ggplotGrob(p.bar2)
combined1 <- p11+p12+p3+p.bar1 +(p21)+(p22)+(p2.niche) + ht_plot +plot_layout(ncol = 4, nrow=2, widths=c(3,3,3,2))
print(combined1)
#print(((combined1)/(combined2))&theme(plot.margin = margin(5, 5, 5, 5))+plot_layout(nrow=2, heights=c(3, 1.2)))
dev.off()

pdf("Rplots.pdf", width=24, height=3)
print(combined2)
dev.off()

pdf("Fig1dsupp.pdf", height=6, width=6)
print(p.bar2)
dev.off()




## meta = read.csv("~/shruti_meta_clean (3).csv")
## input_meta = meta[,c('slide_id', 'age', 'sex', 'case_ctrl', 'ICPi',  'malignancy', 'eGFR_base')] %>% arrange(case_ctrl)
## #tile_obj$condition = str_trim(meta[match(str_extract(as.character(obj.merge$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

## test = obj.merge@meta.data %>% filter(tech=='xenium') %>% janitor::clean_names()
## test$condition = str_trim(meta[match(str_extract(as.character(subset(obj.merge, subset=tech=='xenium')$orig.ident), "__(BS.*A[1,2])__(2.+)", group=1), meta$slide_id), 'case_ctrl'])

## sample_ids = test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
##     mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(sample_ids) %>% pull()
## sample_mat = as.matrix(test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
##     mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(-sample_ids))
## rownames(sample_mat) = str_extract(as.character(sample_ids), "__(BS.*A[1,2])__(2.+)", group=1)

## my_sample_col <- data.frame(sample = str_trim(meta$case_ctrl))
## row.names(my_sample_col) <- meta$slide_id

## my_sample_col <- my_sample_col %>% arrange(sample)

## ## my_sample_col = my_sample_col[rownames(sample_mat), ,drop=F]

## my_sample_col = my_sample_col %>%
##     mutate(group=case_when(
##                sample=="Case" ~ "ICI-AIN",
##                sample=="Control" ~ "ICI-ATN")) %>%
##     mutate(group_simple = case_when(
##                sample=="Case" ~ "AIN",
##                sample=="Control" ~ "ATN")) %>% 
##     mutate(patient=rep(seq(1, 4), 2))  %>%
##     unite("patient_id", c(group_simple, patient), remove=F)

## write_tsv(my_sample_col, "my_sample_col.tsv")


## sample_mat = sample_mat[rownames(my_sample_col),]
## rownames(sample_mat) = my_sample_col$patient_id

## sample_mat <- scale(sample_mat)

## my_sample_col = my_sample_col %>% rownames_to_column()
## my_sample_col = my_sample_col %>% column_to_rownames(var='patient_id')

## # # Define color function
## col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
## ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5))
## ha1 = rowAnnotation(samples = my_sample_col$group, 
##                     col=list(samples=c('ICI-AIN'='gray', 'ICI-ATN'='black')),
## #                    annotation_label_location = "top",  # move label to top
##                     annotation_name_gp = gpar(fontsize = 5),  
##                     annotation_legend_param = list(
##                     samples = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5), direction = "horizontal")))


## colnames(sample_mat)[colnames(sample_mat)=='Immune Cell'] = 'Immune (LowQ)'

## ## cell_label = res %>% filter(!grepl("Intercep", parameter)) %>% arrange(desc(c_effect)) %>% select(cell_label) %>% pull()

## colnames(sample_mat) = gsub(" Cell", "", colnames(sample_mat))

## sample_mat = sample_mat[, cell_label]

## # Desired cell size
## cell_width <- unit(2.5, "mm")   # width per column
## cell_height <- unit(3.5, "mm")   # height per row

## # Main heatmap
## ht2 <- Heatmap(
##   sample_mat,
##   name = "Expression",
##   col = col_fun,
##   width = cell_width * ncol(sample_mat),
##   height = cell_height * nrow(sample_mat),
##   column_names_gp = gpar(fontsize = 8),
##   heatmap_legend_param = list(
##     title_gp = gpar(fontsize = 8),
##     labels_gp = gpar(fontsize = 8),
##     direction = "horizontal"
##   ),cluster_columns = FALSE,
##   column_names_side = "top",
##   right_annotation = ha1)
