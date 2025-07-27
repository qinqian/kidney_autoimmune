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


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


imm.niche <- readRDS("250711_niches.rds")

set.seed(999)
library(glue)
obj.cna <- association.Seurat(
    ## seurat_object = obj.merge, # somehow lennard's object cannot work, perhaps version of seurat while saving
    seurat_object = imm.niche,
    test_var = 'condition_num', 
    samplem_key = 'sample_id', 
    graph_use = 'RNA_snn', 
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    ## covs = c("age", "sex", "ICPi") ## no covariates to include 
)

markers = c("TNXB", "COL5A1", "PDGFRA",
            "CXCL9", "SLAMF7", "CD38", # "IL2RG",
            "CASR", "UMOD", "SIM2",  # "CDH3"
            "HNF4A", "CUBN", "LRP2", #
            "CALB1", "HSD11B2", "SCNN1G",
            "NOTCH3", "CSPG4", "MCAM",
            "HAVCR1", "CDH6", "SOX9",
            "PLA2R1", "PODXL", "WT1"
            )


p2 = DotPlot(subset(imm.niche, subset=niche_label != 'Skeletal Muscle'), features=markers, group.by='niche_label') + theme(axis.text.x=element_text(face="bold", angle=90, hjust = 1, vjust=1), axis.text.y=element_text(face="bold")) + xlab("") + ylab("") + get_theme(angle=90, size=5) + theme(legend.box="vertical", legend.margin=margin(), plot.margin = unit(c(0,0,0,0), "cm"), ) + scale_size_continuous(range = c(0.1, 2)) +ggtitle("")

p1 = (DimPlot_scCustom(subset(imm.niche, subset=niche_label != 'Skeletal Muscle'), group.by='niche_label', label.size=2, label=T, repel=T, seed=99, label.box=T, raster=F)+ggplot2::theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.text=element_text(size=8), legend.position = "none", axis.line = element_line(arrow = arrow(type = "closed", length = unit(3, 'pt'))))+
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2"))

pdf("Fig2a.pdf", height=4, width=10.5)
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "cm")
)
(p1 +
FeaturePlot_scCustom(obj.cna, raster=F, features = c('cna_ncorrs_fdr10'))[[1]] + 
    scale_color_gradient2_tableau() + 
    guides(x = axis, y = axis)+
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")+ggplot2::theme(plot.margin = unit(c(0,0,0,0), "cm"), axis.line = element_line(arrow = arrow(type = "closed", length = unit(3, 'pt'))), legend.text=element_text(size=6), legend.position="bottom", legend.direction="horizontal")  + p2) + plot_layout(width=c(1.6, 1.6, 1.2), ncol=3) #& theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())
dev.off()

library(tmap) 
library(tmap)
library(spdep)
tmap_mode('plot')

tile_sf <- imm.niche@meta.data %>% 
  st_as_sf()
