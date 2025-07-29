library(stringr)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(scico)

library(tidyverse)
suppressPackageStartupMessages({
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
    library(rcna)
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


patient <- read_tsv("my_sample_col.tsv")

set.seed(999)
library(glue)
imm.niche@meta.data$condition_num = as.numeric(factor(str_trim(imm.niche@meta.data$condition), levels=c("Control", "Case")))


obj.cna <- association.Seurat(
    ## seurat_object = obj.merge, # somehow lennard's object cannot work, perhaps version of seurat while saving
    seurat_object = imm.niche,
    test_var = 'condition_num', 
    samplem_key = 'sample_id', 
    graph_use = 'RNA_snn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    #covs = c("age", "sex", "ICPi") ## no covariates to include 
)

markers = c("TNXB", "COL5A1", "PDGFRA",
            "CXCL9", "SLAMF7", "CD38", # "IL2RG",
            "CASR", "UMOD", "SIM2",  # "CDH3"
            "HNF4A", "CUBN", "LRP2", #
            "CALB1", "HSD11B2", "SCNN1G",
            "NOTCH3", "CSPG4", "MCAM",
            "HAVCR1", "CDH6", "SOX9",
            "PLA2R1", "PODXL", "WT1")

names(markers) <- c("Fib", "Fib", "Fib",
                    "Immune", "Immune", "Immune",
                    "TAL", "TAL", "TAL",
                    "PT", "PT", "PT",
                    "C-Duct", "C-Duct", "C-Duct",
                    "Vessel", "Vessel", "Vessel",
                    "I-PT", "I-PT", "I-PT",
                    "Gl", "Gl", "Gl"
                    )

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm")
)

print('------')
p2 = DotPlot_scCustom(subset(imm.niche, subset=niche_label != 'Skeletal Muscle'), features=markers, group.by='niche_label') + theme(axis.text.x=element_text(face="bold", angle=90, hjust = 1, vjust=1), axis.text.y=element_text(face="bold")) + xlab("") + ylab("") + get_theme(angle=90, size=5) + theme(legend.box="vertical", legend.margin=margin(), plot.margin = unit(c(0,0,0,0), "cm"), ) + scale_size_continuous(range = c(0.1, 2)) +ggtitle("")

p1 = DimPlot_scCustom(subset(imm.niche, subset=niche_label != 'Skeletal Muscle'), group.by='niche_label', label.size=2, label=T, repel=T, seed=99, label.box=T, raster=F)+ggplot2::theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.text=element_text(size=8), legend.position = "none", axis.line = element_line(arrow = arrow(type = "closed", length = unit(3, 'pt'))))+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")

p3 <- FeaturePlot_scCustom(obj.cna, raster=F, features = c('cna_ncorrs_fdr10'))[[1]] + 
    scale_color_gradient2(high = "#de2d26", mid = "white", low = "#2c7fb8", midpoint = 0)+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")+theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.line = element_line(arrow = arrow(type = "closed", length = unit(3, 'pt'))), legend.text=element_text(size=6), legend.position="bottom", legend.direction="horizontal") 

pdf("Fig2a.pdf", height=4, width=10.5)
print(p1 +p3 + p2 + plot_layout(width=c(1.6, 1.6, 1.8), ncol=3)) 
dev.off()

