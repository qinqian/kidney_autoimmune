##library(celldex)

library(harmony)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(sf)
library(SingleR)
library(patchwork)
library(ggdendro)
library(sceasy)
library(cowplot)


myeloid = readRDS("/Users/qian/Documents/active_projects/xenium/data/so_myeloid.rds")
myeloid[["RNA"]] <- as(myeloid[["RNA"]], "Assay")
kidn_seu_obj = readRDS("../data/miles/shapes_seur_kidney_v5.rds")
kidn_seu_obj$seur[["RNA"]] <- as(kidn_seu_obj$seur[["RNA"]], "Assay")

use_condaenv("r")

sceasy::convertFormat(myeloid, from="seurat", to="anndata",
                      outFile='so_myeloid.h5ad')

myeloid = myeloid[, !grepl("DC15", myeloid@meta.data$cluster_name)]

myeloid.markers = FindAllMarkers(myeloid, only.pos = T, group.by="cluster_name")

seq.top10 <- myeloid.markers %>% group_by(cluster) %>% mutate(pct.d=pct.1-pct.2) %>%
    filter(pct.1>=0.1, p_val_adj<=0.05, abs(avg_log2FC)>=1) %>% slice_min(tibble(-pct.d, -avg_log2FC), n=5)

png("sanity_check_myeloid.png", width=1800, height=1200)
print(DoHeatmap(object = myeloid, features = unique(seq.top10$gene), label = TRUE, group.by='cluster_name',
                slot="data")+NoLegend())
dev.off()

png("sanity_check_myeloid_dot.png", width=1800, height=1200)
DotPlot(object = myeloid, features = unique(seq.top10$gene), group.by='cluster_name')+theme(axis.text.x=element_text(size=12, angle=90))
dev.off()

# subset T cells for typing
kidn_seu_obj_M = kidn_seu_obj$seur[, kidn_seu_obj$seur@meta.data$coarse_ids%in%c("Macrophage")]
kidn_seu_obj_M =
    kidn_seu_obj_M %>% 
    NormalizeData(normalization.method = 'LogNormalize', verbose = F)

VariableFeatures(kidn_seu_obj_M) <- split(row.names(kidn_seu_obj_M@meta.data), kidn_seu_obj_M@meta.data$sample_id) %>% lapply(function(cells_use) {
    kidn_seu_obj_M[,cells_use] %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
        VariableFeatures()
    }) %>% unlist %>% unique

kidn_seu_obj_M <- kidn_seu_obj_M %>%
    ScaleData(verbose = FALSE) %>% 
    RunPCA(features = VariableFeatures(kidn_seu_obj_M), npcs = 100, verbose = FALSE) %>%
    RunHarmony("sample_id") %>% 
    FindNeighbors(reduction = "harmony") %>% 
    FindClusters(resolution = 0.4) %>%
    RunUMAP(reduction="harmony", dims=1:100, reduction.key='HUMAP_')

ref.set = as.SingleCellExperiment(myeloid)
sce <- as.SingleCellExperiment(kidn_seu_obj_M)

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$cluster_name, de.method = 'classic', de.n=50, num.threads=12)
write_tsv(as.data.frame(pred.cnts), "M_kidney_lupus_pred.cnts")


lbls.keep <- table(pred.cnts$labels)>10
kidn_seu_obj_M$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')


# diagnostics of singleR results
png("kidney_M_heatmap.png", width=1200, height=600)
plotScoreHeatmap(pred.cnts)
dev.off()

## more diagnostics of delta
png("kidney_M_deltadiag.png", width=1200, height=600)
plotDeltaDistribution(pred.cnts)
dev.off()

png("kidney_M_umap.png", width=1600, height=1200)
p1 = DimPlot(kidn_seu_obj_M, group.by="SingleR.labels", label=T) + NoLegend()
p2 = FeaturePlot(kidn_seu_obj_M, features="ITGAX", label=T)
p3 = FeaturePlot(kidn_seu_obj_M, features="CXCL9", label=T)
p4 = DimPlot(myeloid, group.by="cluster_name", label=T) + NoLegend()
print(p1 + p2 + p3 + p4)
dev.off()

pop.markers = FindAllMarkers(kidn_seu_obj_M, group.by='SingleR.labels', test.use='wilcox', only.pos=T)

top10 <- pop.markers %>% group_by(cluster) %>% mutate(pct.d=pct.1-pct.2) %>%
    filter(pct.1>=0.1, p_val_adj<=0.05, abs(avg_log2FC)>=1) %>% slice_min(tibble(-pct.d, -avg_log2FC), n=5)

png("singleR_kidney_M_marker_heatmap.png", width=1600, height=800)
print(DoHeatmap(object = kidn_seu_obj_M, features = unique(top10$gene), label = TRUE, group.by='SingleR.labels'))+NoLegend()
dev.off()

png("singleR_kidney_M_seqmarker_heatmap.png", width=1600, height=800)
print(DoHeatmap(object = kidn_seu_obj_M, features = unique(seq.top10$gene), label = TRUE, group.by='SingleR.labels'))+NoLegend()
dev.off()

png("singleR_kidney_M_marker_dotplot.png", width=1600, height=800)
print(DotPlot(object = kidn_seu_obj_M, features = unique(top10$gene), group.by='SingleR.labels')+theme(axis.text.x=element_text(size=12, angle=90)))
dev.off()

png("singleR_kidney_M_seqmarker_dotplot.png", width=1600, height=800)
print(DotPlot(object = kidn_seu_obj_M, features = unique(seq.top10$gene), group.by='SingleR.labels')+theme(axis.text.x=element_text(size=12, angle=90)))
dev.off()

saveRDS(kidn_seu_obj_M, "kidn_seu_obj_M.rds")

sceasy::convertFormat(kidn_seu_obj_M, from="seurat", to="anndata", outFile='shapes_seur_kidney_v5_miles_version2_singleR_M.h5ad', main_layer='data')

## ## cluster dot plot
## kidn_seu_obj_M_dotdata.list = list()
## for (g in top10$gene) {
##     kidn_seu_obj_M_dotdata = data.frame(gene=kidn_seu_obj_M[g, ]@assays$RNA@layers$data, cluster=kidn_seu_obj_M$SingleR.labels)
##     kidn_seu_obj_M_dotdata = kidn_seu_obj_M_dotdata %>% group_by(cluster) %>% summarise(count=mean(gene), cell_ct=sum(gene>0)/n()) %>% mutate(gene=g)
##     kidn_seu_obj_M_dotdata.list[[g]] = kidn_seu_obj_M_dotdata
## }
## kidn_seu_obj_M_dotdata_dfs = bind_rows(kidn_seu_obj_M_dotdata.list)

## # make data square to calculate euclidean distance
## mat <- kidn_seu_obj_M_dotdata_dfs %>% 
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
## mat <- kidn_seu_obj_M_dotdata_dfs %>% 
##   select(-cell_ct) %>%  # drop unused columns to faciliate widening
##   pivot_wider(names_from = cluster, values_from = count) %>% 
##   data.frame(check.names=F) # make df as tibbles -> matrix annoying

## row.names(mat) <- mat$gene  # put gene in `row`
## mat <- mat[,-1] #drop gene column as now in rows
## v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
## ############ NOTICE THE t() above)

## dotplot = kidn_seu_obj_M_dotdata_dfs %>%
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
