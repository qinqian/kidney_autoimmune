library(tidyverse)

library(harmony)
library(Seurat)
library(SingleCellExperiment)
library(sf)
library(SingleR)
library(patchwork)
library(ggdendro)
library(cowplot)
library(sceasy)


ref = "/Users/qian/Documents/active_projects/xenium/data/so_t_nk.rds"
t_nk = readRDS(ref)
### https://github.com/cellgeni/sceasy/issues/82
t_nk[["RNA"]] <- as(t_nk[["RNA"]], "Assay")

## conda install pandas==1.5.3 python==3.10 anndata
## pip install anndata -U
use_condaenv("r")

sceasy::convertFormat(t_nk, from="seurat", to="anndata",
                      outFile='so_t_nk.h5ad')

kidn_seu_obj = readRDS("../data/miles/shapes_seur_kidney_v5.rds")
kidn_seu_obj$seur[["RNA"]] <- as(kidn_seu_obj$seur[["RNA"]], "Assay")


## manually curated on slides
t_nk@meta.data$t_subtype = t_nk@meta.data %>% mutate(test=case_when(
                              grepl("NK$", cluster_name) ~ "NK",
                              grepl("CTL$", cluster_name) ~ "CTL",
                              grepl("Effector", cluster_name) ~ "Effector memory",
                              grepl("CD8", cluster_name) & (!grepl("Effector", cluster_name)) ~ "CD8+ T",
                              grepl("Central|Resident|T18", cluster_name) ~ "CD4+ Central/Resident Memory/Th17-like",
                              .default=cluster_name
                              )) %>% select(test) %>% pull()


t_nk = t_nk[, !grepl("T12", t_nk@meta.data$t_subtype)]
DimPlot(t_nk, group.by='t_subtype')


t_nk.markers = FindAllMarkers(t_nk, only.pos = T, group.by="t_subtype")



seq.top10 <- t_nk.markers %>% group_by(cluster) %>% mutate(pct.d=pct.1-pct.2) %>%
    filter(pct.1>=0.1, p_val_adj<=0.05, abs(avg_log2FC)>=1) %>% slice_min(tibble(-pct.d, -avg_log2FC), n=5)

seq.features.list = list()
for (cl in seq.top10 %>% select(cluster) %>% distinct() %>% pull()) {
    print(cl)
    seq.features.list[[cl]] = unique(seq.top10 %>% filter(cluster==cl) %>% select(gene) %>% pull())
}

png("sanity_check_t_nk_version2.png", width=1800, height=1200)
print(DoHeatmap(object = t_nk, features = unique(seq.top10$gene), label = TRUE, group.by='t_subtype',
                slot="data")+NoLegend())
dev.off()

png("sanity_check_t_nk_dot_version2.png", width=1800, height=1200)
DotPlot(object = t_nk, features = unique(seq.top10$gene), group.by='t_subtype')+theme(axis.text.x=element_text(size=12, angle=90))
## DotPlot(object = t_nk, features = seq.features.list, group.by='cluster_name')+theme(axis.text.x=element_text(size=12, angle=90))
dev.off()

# subset T cells for typing
## kidn_seu_obj_T = kidn_seu_obj$seur[, kidn_seu_obj$seur@meta.data$coarse_ids%in%c("T", "T_dblt")]
kidn_seu_obj_T = kidn_seu_obj$seur[, kidn_seu_obj$seur@meta.data$coarse_ids%in%c("T")]

kidn_seu_obj_T =
    kidn_seu_obj_T %>% 
    NormalizeData(normalization.method = 'LogNormalize', verbose = F)

VariableFeatures(kidn_seu_obj_T) <- split(row.names(kidn_seu_obj_T@meta.data), kidn_seu_obj_T@meta.data$sample_id) %>% lapply(function(cells_use) {
    kidn_seu_obj_T[,cells_use] %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
        VariableFeatures()
    }) %>% unlist %>% unique

kidn_seu_obj_T <- kidn_seu_obj_T %>%
    ScaleData(verbose = FALSE) %>% 
    RunPCA(features = VariableFeatures(kidn_seu_obj_T), npcs = 100, verbose = FALSE) %>%
    RunHarmony("sample_id") %>% 
    FindNeighbors(reduction = "harmony") %>% 
    FindClusters(resolution = 0.4) %>%
    RunUMAP(reduction="harmony", dims=1:100, reduction.key='HUMAP_')

ref.set = as.SingleCellExperiment(t_nk)
sce <- as.SingleCellExperiment(kidn_seu_obj_T)

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$t_subtype, de.method = 'classic', de.n=50, num.threads=12)

write_tsv(as.data.frame(pred.cnts), "T_kidney_lupus_pred_version2.cnts")


lbls.keep <- table(pred.cnts$labels)>10
kidn_seu_obj_T$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')


# diagnostics of singleR results
png("kidney_T_heatmap_version2.png", width=1200, height=600)
plotScoreHeatmap(pred.cnts)
dev.off()

## more diagnostics of delta
png("kidney_T_deltadiag_version2.png", width=1200, height=600)
plotDeltaDistribution(pred.cnts)
dev.off()

png("kidney_T_umap_version2.png", width=1600, height=1200)
p1 = DimPlot(kidn_seu_obj_T, group.by="SingleR.labels", label=T) + NoLegend()
p2 = FeaturePlot(kidn_seu_obj_T, features="FOXP3", label=T)
p3 = FeaturePlot(kidn_seu_obj_T, features="GZMK", label=T)
p4 = DimPlot(t_nk, group.by="t_subtype", label=T) + NoLegend()
print(p1 + p2 + p3 + p4)
dev.off()

png("kidney_T_umap_version2_seurat_clusters.png", width=800, height=500)
p5 = DimPlot(kidn_seu_obj_T, group.by="seurat_clusters", label=T) 
print(p5)
dev.off()

pop.markers = FindAllMarkers(kidn_seu_obj_T, group.by='SingleR.labels', test.use='wilcox', only.pos=T)

top10 <- pop.markers %>% group_by(cluster) %>% mutate(pct.d=pct.1-pct.2) %>%
    filter(pct.1>=0.1, p_val_adj<=0.05, abs(avg_log2FC)>=1) %>% slice_min(tibble(-pct.d, -avg_log2FC), n=5)

png("singleR_kidney_T_marker_heatmap_version2.png", width=1600, height=800)
print(DoHeatmap(object = kidn_seu_obj_T, features = unique(top10$gene), label = TRUE, group.by='SingleR.labels'))+NoLegend()
dev.off()

png("singleR_kidney_T_seqmarker_heatmap_version2.png", width=1600, height=800)
print(DoHeatmap(object = kidn_seu_obj_T, features = unique(seq.top10$gene), label = TRUE, group.by='SingleR.labels'))+NoLegend()
dev.off()

png("singleR_kidney_T_marker_dotplot_version2.png", width=1600, height=800)
print(DotPlot(object = kidn_seu_obj_T, features = unique(top10$gene), group.by='SingleR.labels')+theme(axis.text.x=element_text(size=12, angle=90)))
dev.off()

png("singleR_kidney_T_seqmarker_dotplot_version2.png", width=1600, height=800)
print(DotPlot(object = kidn_seu_obj_T, features = unique(seq.top10$gene), group.by='SingleR.labels')+theme(axis.text.x=element_text(size=12, angle=90)))
dev.off()

sceasy::convertFormat(kidn_seu_obj_T, from="seurat", to="anndata", outFile='shapes_seur_kidney_v5_miles_version2_singleR_T.h5ad', main_layer='data')
saveRDS(kidn_seu_obj_T, "kidn_seu_obj_T.rds")


## TODO map kidn_seu_obj_T t_subtype back to the kidn_seu_obj





## ## cluster dot plot
## kidn_seu_obj_T_dotdata.list = list()
## for (g in top10$gene) {
##     kidn_seu_obj_T_dotdata = data.frame(gene=kidn_seu_obj_T[g, ]@assays$RNA@layers$data, cluster=kidn_seu_obj_T$SingleR.labels)
##     kidn_seu_obj_T_dotdata = kidn_seu_obj_T_dotdata %>% group_by(cluster) %>% summarise(count=mean(gene), cell_ct=sum(gene>0)/n()) %>% mutate(gene=g)
##     kidn_seu_obj_T_dotdata.list[[g]] = kidn_seu_obj_T_dotdata
## }
## kidn_seu_obj_T_dotdata_dfs = bind_rows(kidn_seu_obj_T_dotdata.list)

## # make data square to calculate euclidean distance
## mat <- kidn_seu_obj_T_dotdata_dfs %>% 
##   select(-cell_ct) %>%  # drop unused columns to faciliate widening
##   pivot_wider(names_from = cluster, values_from = count) %>% 
##   data.frame() # make df as tibbles -> matrix annoying
## row.names(mat) <- mat$gene  # put gene in `row`
## mat <- mat[,-1] #drop gene column as now in rows
## clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
## ddgram <- as.dendrogram(clust) # create dendrogram
## ## plot(clust,  hang = -1)
## ggtree_plot <- ggtree::ggtree(ddgram, branch.length = "none")

## # make data square to calculate euclidean distance
## mat <- kidn_seu_obj_T_dotdata_dfs %>% 
##   select(-cell_ct) %>%  # drop unused columns to faciliate widening
##   pivot_wider(names_from = cluster, values_from = count) %>% 
##   data.frame(check.names=F) # make df as tibbles -> matrix annoying

## row.names(mat) <- mat$gene  # put gene in `row`
## mat <- mat[,-1] #drop gene column as now in rows
## v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
## ############ NOTICE THE t() above)

## dotplot = kidn_seu_obj_T_dotdata_dfs %>%
##     filter(count>0, cell_ct>0) %>%
##     mutate(gene=factor(gene, levels = clust$labels[clust$order])) %>% 
##     mutate(cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%  
##     ggplot(aes(x=gene, y = cluster, color = count, size = cell_ct)) + 
##     scale_color_viridis_c() + 
##     geom_point()  + theme(axis.ticks = element_blank()) + cowplot::theme_cowplot() +
##     ylab('') +
##     theme(axis.ticks = element_blank()) +    
##     theme(axis.line  = element_blank()) +
##     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
##     scale_y_discrete(position = "right")

## ddgram_col <- as.dendrogram(v_clust)
## ggtree_plot_col <- ggtree::ggtree(ddgram_col, branch.length = "none") # + aplot::ylim2(dotplot)
## ggtree_plot <- ggtree_plot + ggtree::layout_dendrogram()

## png("test.png", width=1800, height=1200)
## ## ggtree_plot_yset <- ggtree_plot + aplot::ylim2(dotplot)
## ## plot_grid(ggtree_plot_yset, NULL, dotplot, nrow = 1, rel_widths = c(0.8, -0.05, 2), align = 'h')
## ## labels <- ggplot(gene_cluster %>% 
## ##                    mutate(`Cell Type` = Group,
## ##                            cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
## ##                  aes(x = cluster, y = 1, fill = `Cell Type`)) + 
## ##   geom_tile() + 
## ##   scale_fill_brewer(palette = 'Set1') + 
## ##   theme_nothing() +
## ##   xlim2(dotplot)
## ## legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))
## plot_spacer() + plot_spacer() + ggtree_plot +
##   ggtree_plot_col + plot_spacer() + dotplot + 
##   plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.1, 1))
## dev.off()
