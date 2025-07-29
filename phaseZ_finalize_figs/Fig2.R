library(stringr)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(scico)

library(tidyverse)
suppressPackageStartupMessages({
    library(tessera)
    #library(scCustomize)
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
    #library(rcna)
})


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
print(str(sc.niche))

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


ht_opt("heatmap_row_names_gp" = gpar(fontsize = 6))
ht <- Heatmap(
  t(sc.niche.heatmap),
  name = "Cell Composition Ratio",
  col = col_fun,
  ## annotation_name_gp = gpar(fontsize = 5),    # annotation label font size
  column_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6),
    direction = "horizontal"
  ),cluster_columns = FALSE,
  column_names_side = "top")


pdf("Fig2c.pdf", height=3, width=2.5)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom")
dev.off()


imm.niche@meta.data$condition = factor(str_trim(imm.niche@meta.data$condition), levels=c("Control", "Case"))
imm.niche@meta.data$condition_num = as.numeric(factor(str_trim(imm.niche@meta.data$condition), levels=c("Control", "Case")))


library(sccomp)
print(head(imm.niche))
res = imm.niche |>
    sccomp_estimate(
      formula_composition = ~ condition,
      sample = "sample_id", cell_group = "niche_label_fine",
      cores = 2, verbose=FALSE
    )

res = res |>  sccomp_test()
pdf("Fig2e.pdf", height=9, width=4.5)
ht_grob <- grid.grabExpr(draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom"))
# 5. Arrange the two grobs
grob_ggplot = ggplotGrob(res |>  plot_1D_intervals() + theme_bw(base_size=8) + theme(legend.position='bottom') + ylab("Cell type"))
grid.arrange(ht_grob, grob_ggplot, ncol = 1)
dev.off()
