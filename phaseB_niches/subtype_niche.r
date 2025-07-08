set.seed(999)

suppressPackageStartupMessages({
library(scCustomize)
    library(tessera)
    library(Seurat)
library(sccomp)
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
})

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R")

#0 Fibrosis
#1 Immune
#2 TAL (Thick ascending limb of Loop of Henle)
#3 + #4 Proximal Tubule (there are some subsegments of the PT that are hard to identify from this, I would just stuff this together
#5 Collecting Duct 
#6 Epithelial Cells + Smooth Muscle (Vessels)
#7 Damaged Proximal Tubule
#8 Glomeruli
t.obj = readRDS("tessera/vignettes/ICI_tesera.rds")


pdf("test.pdf", width=12, height=6)
p1 = DimPlot(t.obj[, t.obj@meta.data$condition=='Case  '], reduction = 'umap', group.by = 'seurat_clusters') + scale_color_tableau('Classic 10') 
p2 = FeaturePlot(t.obj[, t.obj@meta.data$condition=='Case  '], feature='CD8A')
p3 = FeaturePlot(t.obj[, t.obj@meta.data$condition=='Case  '], feature='CD3E')
p4 = DimPlot(t.obj[, t.obj@meta.data$condition=='Control'], reduction = 'umap', group.by = 'seurat_clusters') + scale_color_tableau('Classic 10') 
p5 = FeaturePlot(t.obj[, t.obj@meta.data$condition=='Control'], feature='CD8A')
p6 = FeaturePlot(t.obj[, t.obj@meta.data$condition=='Control'], feature='CD3E')
(p1 + p2 + p3) / (p4 + p5 + p6)
dev.off()


levels(t.obj@meta.data$seurat_clusters) = c("Fibrosis", "Immune", "TAL (Thick ascending limb of Loop of Henle)",
                                            "Proximal Tubule", "Proximal Tubule",
                                            "Collecting Duct", "Epithelial Cells + Smooth Muscle (Vessels)", "Damaged Proximal Tubule",
                                            "Glomeruli")


pdf("test2.pdf", width=13.5, height=8)
p1 = DimPlot(t.obj[, t.obj@meta.data$condition=='Case  '], reduction = 'umap', group.by = 'seurat_clusters') + scale_color_tableau('Classic 10') 
p4 = DimPlot(t.obj[, t.obj@meta.data$condition=='Control'], reduction = 'umap', group.by = 'seurat_clusters') + scale_color_tableau('Classic 10') 
print(p1 + p4 + plot_layout(widths=c(1, 1), ncol=2, guides='collect') & theme(legend.position = 'bottom'))
dev.off()


markers = c("TNXB", "COL5A1", "PDGFRA",
            "CXCL9", "SLAMF7", "CD38", # "IL2RG",
            "CASR", "UMOD", "SIM2",  # "CDH3"
            "HNF4A", "CUBN", "LRP2", #
            "CALB1", "HSD11B2", "SCNN1G",
            "NOTCH3", "CSPG4", "MCAM",
            "HAVCR1", "CDH6", "SOX9",
            "PLA2R1", "PODXL", "WT1"
            )

pdf("test3.pdf", width=8, height=11.5)
DotPlot(t.obj, features=markers, group.by='seurat_clusters') + coord_flip() + theme(axis.text.x=element_text(face="bold", angle=45, hjust = 1, vjust=1), axis.text.y=element_text(face="bold"))
dev.off()

pdf("test4.pdf", width=9, height=16)
p1 = FeaturePlot_scCustom(t.obj[, t.obj@meta.data$condition=='Case  '], reduction = 'umap', features=markers, na_color = "lightgray", num_columns = 3, colors_use = viridis_light_high)
print(p1)
dev.off()

t.obj.immune = t.obj[, t.obj@meta.data$seurat_clusters=='Immune']

t.obj.immune =
    t.obj.immune %>% 
    NormalizeData(normalization.method = 'LogNormalize', verbose = F)

VariableFeatures(t.obj.immune) <- split(row.names(t.obj.immune@meta.data), t.obj.immune@meta.data$sample_id) %>% lapply(function(cells_use) {
    t.obj.immune[,cells_use] %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
        VariableFeatures()
    }) %>% unlist %>% unique

t.obj.immune <- t.obj.immune %>%
    ScaleData(verbose = FALSE) %>% 
    RunPCA(features = VariableFeatures(t.obj.immune), npcs = 30, verbose = FALSE) %>%
    RunHarmony("sample_ids") %>% 
    FindNeighbors(reduction = "harmony") %>% 
    FindClusters(resolution = 0.4) %>%
    RunUMAP(reduction="harmony", dims=1:30, reduction.key='HUMAP_')


res =
    t.obj |>
    sccomp_estimate(
      formula_composition = ~ condition,
      .sample = sample_ids, .cell_group = seurat_clusters,
      cores = 1, verbose=FALSE
    )
res = res |>  sccomp_test()

pdf("sccomp_ici_tesera.pdf")
print(res |> sccomp_boxplot(factor = "condition"))
dev.off()


pdf("sccomp_ici_tesera_1d.pdf", width=9.6, height=5)
print((res |>  plot_1D_intervals()) + theme_bw(base_size=15) + ylab("Niche group"))
dev.off()


pdf("niche_immune_subtyping.pdf", width=9.6, height=5)
p1=DimPlot(t.obj.immune)
p2=FeaturePlot(t.obj.immune,features=c("CXCL9", "CXCL10"))
print((p1+p2)+plot_layout(widths=c(2.5, 2, 2)))
dev.off()


t.obj.immune.markers = FindAllMarkers(t.obj.immune)


#s.obj = readRDS("../data/miles/shapes_seur_kidney_v5.rds")
#meta = read.csv("../../xenium/shruti_meta_clean (3).csv")
