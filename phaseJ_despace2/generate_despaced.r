library(Seurat)
library(EnhancedVolcano)
library(DESpace)
library(ggplot2)
library(SpatialExperiment)
library(muSpaData)
library(reshape2)
library(tidyverse)
library(patchwork)
library(splines)
library(edgeR)
library(ggsci)

sopa <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
meta = read.csv("~/shruti_meta_clean (3).csv")
sopa$tile_obj$condition = str_trim(meta[match(str_extract(as.character(sopa$tile_obj$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

sc <- sopa$obj
sc$niche_label <- sopa$tile_obj@meta.data[as.character(sc$tile_id), 'niche_label']
sc$condition = str_trim(meta[match(str_extract(as.character(sc$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

tiles <- sopa$tile_obj
tiles <- tiles[, tiles$niche_label!='Muscle Outlier']
tiles <- as.SingleCellExperiment(tiles)

sc <- sc[,sc$niche_label!='Muscle Outlier']
print(levels(sc$niche_label))

sc@meta.data[,'X'] <- sc@reductions$spatial@cell.embeddings[,1]
sc@meta.data[,'Y'] <- sc@reductions$spatial@cell.embeddings[,2]

sc$condition = factor(sc$condition, levels=c("Control", "Case"))
sc <- as.SingleCellExperiment(sc)



coordinates <- c("X", "Y") # coordinates of cells
spatial_cluster <- 'niche_label' # Banksy spatial clusters
condition_col <- 'condition'       # regeneration time phases
sample_col <- 'sample_id'          # tissue section id
colData(tiles) |> head()

results <- dsp_test(spe = sc,
                    cluster_col = spatial_cluster,
                    sample_col = sample_col,
                    condition_col = condition_col,
                    verbose = TRUE)
head(results$gene_results, 5)

write_tsv(results$gene_results, file="global_despace2.tsv")


pdf("initial_data_tiles.pdf", width=22, height=12)
# The spatial cluster assignments are available in the `colData(spe)`
CD <- colData(sc) |> as.data.frame()
CD %>% 
    ggplot(aes(x = X, y = Y, color=niche_label)) +
#    geom_sf(aes(geometry=shape, fill = factor(niche_label)), color=NA) +
    geom_point(size = 0.25)+
    facet_wrap(~condition+sample_id, ncol=4) +
    theme_void() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(color = NULL, title = "Tessera Spatial Clusters")+ggsci::scale_fill_d3("category20c")
dev.off()

sample_ids <- unique(CD$sample_id)
feature <- results$gene_results$gene[1]

feature <- 'VDR'
plots <- lapply(sample_ids, function(sample_id) {
  # Subset spe for each sample
  spe_j <- sc[, colData(sc)$sample_id == sample_id]
  plot <- FeaturePlot(spe_j, feature,
                      coordinates = coordinates,
                      platform = "Xenium", ncol = 1,
                      diverging = TRUE,
                      point_size = 0.1, legend_exprs = TRUE) + 
    theme(legend.position = "right",
          legend.key.size = unit(0.5, 'cm')) +
    labs(color = "") + ggtitle(sample_id) 
  return(plot)
})

pdf("top1_marker.pdf", width=22, height=12)
combined_plot <- wrap_plots(plots, ncol = 4) + 
    # common legend
    plot_layout(guides = 'collect')
print(combined_plot)
dev.off()

cluster_results <- individual_dsp(sc,
                                  cluster_col = spatial_cluster,
                                  sample_col = sample_col,
                                  condition_col = condition_col, min_counts=10)

names(cluster_results)
write_tsv(cluster_results[["Immune infiltrated interstitial"]], file='Imm_niche_de.tsv')

ress <- lapply(1:20, function(gene_id) {
(feature <- rownames(cluster_results[["Immune infiltrated interstitial"]])[gene_id])
cps <- cpm(results$estimated_y, log = TRUE)
cps_name <- colnames(cps)
mdata <- data.frame(
    log_cpm = cps[feature, ] ,
    Banksy_smooth = factor(sub(".*_", "", cps_name)),
    day = str_trim(meta[match(str_extract(as.character(cps_name), "(.*)__", group=1), meta$slide_id), 'case_ctrl']),
    sample_id = str_extract(as.character(cps_name), "(.*)__", group=1)
)
plt <- ggplot(mdata, aes(x = factor(day), y = log_cpm)) +
    geom_jitter(aes(color = Banksy_smooth), size = 2, width = 0.1) + 
    geom_boxplot(aes(fill = ifelse(Banksy_smooth == "Immune infiltrated interstitial", 
                                   "Immune infiltrated interstitial", "non-Immune niche")), 
                 position = position_dodge(width = 0.8), alpha = 0.5) +
    ## scale_x_discrete(breaks = c(2, 10, 20)) +  
    scale_fill_manual(values = c("#4DAF4A", "grey")) + 
    labs(title = feature, x = "Case control", 
         y = "log-2 counts per million (logCPM)", fill = "",
         color = "Tessera tiles cluster") +
     theme(legend.position = "right")
return(plt)
})

pdf("test_cluster.pdf", width=28, height=16)
print(wrap_plots(ress, ncol=5))
dev.off()

write_tsv(cluster_results[["Fibrogenic interstitial"]], file='fib_niche_de.tsv')

ress <- lapply(1:20, function(gene_id) {
(feature <- rownames(cluster_results[["Fibrogenic interstitial"]])[gene_id])
cps <- cpm(results$estimated_y, log = TRUE)
cps_name <- colnames(cps)
mdata <- data.frame(
    log_cpm = cps[feature, ] ,
    Banksy_smooth = factor(sub(".*_", "", cps_name)),
    day = str_trim(meta[match(str_extract(as.character(cps_name), "(.*)__", group=1), meta$slide_id), 'case_ctrl']),
    sample_id = str_extract(as.character(cps_name), "(.*)__", group=1)
)
plt <- ggplot(mdata, aes(x = factor(day), y = log_cpm)) +
    geom_jitter(aes(color = Banksy_smooth), size = 2, width = 0.1) + 
    geom_boxplot(aes(fill = ifelse(Banksy_smooth == "Fibrogenic interstitial", 
                                   "Fibrogenic interstitial", "Other niche")), 
                 position = position_dodge(width = 0.8), alpha = 0.5) +
    ## scale_x_discrete(breaks = c(2, 10, 20)) +  
    scale_fill_manual(values = c("#4DAF4A", "grey")) + 
    labs(title = feature, x = "Case control", 
         y = "log-2 counts per million (logCPM)", fill = "",
         color = "Tessera tiles cluster") +
     theme(legend.position = "right")
return(plt)
})

pdf("test_cluster2.pdf", width=28, height=16)
print(wrap_plots(ress, ncol=5))
dev.off()


#ress <- lapply(1:10, function(gene_id) {
#(feature <- rownames(cluster_results[["Muscle Outlier"]])[gene_id])
#cps <- cpm(results$estimated_y, log = TRUE)
#cps_name <- colnames(cps)
#mdata <- data.frame(
#    log_cpm = cps[feature, ] ,
#    Banksy_smooth = factor(sub(".*_", "", cps_name)),
#    day = str_trim(meta[match(str_extract(as.character(cps_name), "(.*)__", group=1), meta$slide_id), 'case_ctrl']),
#    sample_id = str_extract(as.character(cps_name), "(.*)__", group=1)
#)
#plt <- ggplot(mdata, aes(x = factor(day), y = log_cpm)) +
#    geom_jitter(aes(color = Banksy_smooth), size = 2, width = 0.1) + 
#    geom_boxplot(aes(fill = ifelse(Banksy_smooth == "Muscle Outlier", 
#                                   "Muscle Outlier", "Other niche")), 
#                 position = position_dodge(width = 0.8), alpha = 0.5) +
#    ## scale_x_discrete(breaks = c(2, 10, 20)) +  
#    scale_fill_manual(values = c("#4DAF4A", "grey")) + 
#    labs(title = feature, x = "Case control", 
#         y = "log-2 counts per million (logCPM)", fill = "",
#         color = "Tessera tiles cluster") +
#     theme(legend.position = "right")
#return(plt)
#})
#pdf("test_cluster3.pdf", width=24, height=12)
#print(wrap_plots(ress, ncol=5))
#dev.off()
