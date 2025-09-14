invisible(lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx"), library, character.only = T))
set.seed(999)

get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

getDE2 <- function(x, cluster='cell_type', top_n=5) {
    de.celltype <- presto::wilcoxauc(x, group_by=cluster) %>% dplyr::mutate(pct_diff = pct_in - pct_out)

    features <- de.celltype %>%
         dplyr::group_by(group) %>%
         dplyr::filter(padj <= 0.05 & auc >= 0.6 & pct_diff >= 0.05 & logFC >= 0.25) %>%
         dplyr::arrange(group, dplyr::desc(logFC)) |>
         dplyr::slice_max(logFC, n = top_n) |>
         dplyr::select(group, feature) |>
         dplyr::ungroup() |>
         dplyr::arrange(feature) |>
         dplyr::distinct(feature, .keep_all = TRUE) |>
         dplyr::mutate(gene = feature) %>% dplyr::select(group, gene) %>%
         dplyr::group_by(group) %>%
         dplyr::select(group, gene)
    features <- features |>
        tibble::deframe()
    features
}

getDE <- function(x, cluster='cell_type', top_n=5) {
    de.celltype <- presto::wilcoxauc(x, group_by=cluster) %>% dplyr::mutate(pct_diff = pct_in - pct_out)

    features <- de.celltype %>%
         dplyr::group_by(group) %>%
         dplyr::filter(padj <= 0.05 & auc >= 0.6 & pct_diff > 0 & logFC > 0.) %>%
         dplyr::arrange(group, dplyr::desc(logFC)) |>
         dplyr::slice_max(logFC, n = top_n) |>
         dplyr::select(group, feature) |>
         dplyr::ungroup() |>
         dplyr::arrange(feature) |>
         #dplyr::distinct(feature, .keep_all = TRUE) |>
         dplyr::mutate(gene = feature) %>% dplyr::select(group, gene) %>%
         dplyr::group_by(group) %>%
         dplyr::select(group, gene)
    #features <- features |>
    #    tibble::deframe()
    features |> mutate(gene_label=paste0(gene, "_", group),
                       gene = gene, 
                       cluster=group)
}

suppressPackageStartupMessages({
    library(glue)
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
    library(presto)
    library(tibble)
    library(scDotPlot)
    library(SingleR)    
    library(Seurat)
    library(presto)
    library(dplyr)
    library(ggplot2)
    library(argparse)
    library(tidyverse)
})

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R");

input_data <- c("../../data/miles/shapes_seur_kidney_v5.rds", "../../phaseF_newpipeline/sopa_seg/norm/RA_orig_seg_lognorm.rds")
celltype = "coarse_ids"

    
#sc
myles.skin = readRDS(input_data[1])$seur
print(table(myles.skin$coarse_ids))
myles.skin = myles.skin[, !grepl("dbl", myles.skin$coarse_ids)]
myles.skin = myles.skin[, !grepl("dPT_CUBN", myles.skin$coarse_ids)]
#xenium
sc = readRDS(input_data[2])

library(data.table)
library(tidyr)

CustomDotPlot <- function(seurat_obj, gene_cluster_map,
                          assay = "RNA", slot = "data",
                          cluster_col = "cluster",
                          scale_color = TRUE,
                          color_low = "lightgrey", color_high = "red") {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  if (!all(c("gene_label", "gene", "cluster") %in% colnames(gene_cluster_map))) {
    stop("gene_cluster_map must contain columns: gene_label, gene, cluster")
  }
    
  genes = gene_cluster_map$gene
  labels = gene_cluster_map$gene_label
    
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)[genes, , drop = FALSE]    
  expr_mat <- as.matrix(expr_mat > 0)
  cell_clusters <- seurat_obj@meta.data[[cluster_col]]    
  names(cell_clusters) <- rownames(seurat_obj@meta.data)   
  rownames(expr_mat) = labels

  # Melt to long format
  expr_dt <- as.data.table(as.table(expr_mat))  # columns: gene, cell, expr
  setnames(expr_dt, c("gene", "cell", "expr"))    
  print(dim(expr_dt))
    
  expr_dt[, cluster := cell_clusters[as.character(cell)]]
  pct_dt <- expr_dt[, .(pct_expr = mean(expr) * 100), by = .(gene, cluster)]
  print(dim(pct_dt))

  # Compute average expression (faster via Seurat)
    
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)[genes, , drop = FALSE]    
  expr_mat <- as.matrix(expr_mat)
  cell_clusters <- seurat_obj@meta.data[[cluster_col]]    
  names(cell_clusters) <- rownames(seurat_obj@meta.data)   
  rownames(expr_mat) = labels
  expr_dt <- as.data.table(as.table(expr_mat))  # columns: gene, cell, expr
  setnames(expr_dt, c("gene", "cell", "expr"))    
  print(dim(expr_dt))
    
  expr_dt[, cluster := cell_clusters[as.character(cell)]]
  avg_dt <- expr_dt[, .(avg_expr = mean(expr)), by = .(gene, cluster)]
  print(dim(avg_dt))
    
#   # Join both
   dot_dt <- merge(avg_dt, pct_dt, by = c("gene", "cluster"))
   dot_dt
    
}


overlap.genes = intersect(rownames(myles.skin), rownames(sc))
print(length(overlap.genes))


myles.skin.sub = myles.skin[overlap.genes, ]
sc.sub = subset(sc, features=overlap.genes)

rm(myles.skin)
rm(sc)

sc.sub@meta.data$tech = "xenium"
myles.skin.sub@meta.data$tech = "sc"

spatial <- sc.sub
    #CreateSeuratObject(sc.sub@assays$RNA, meta.data = sc.sub@meta.data, project = "spatial")

set.seed(99)
spatial <- FindClusters(spatial, graph.name='humap_fgraph', resolution=1.2)


library(dplyr)
de <- getDE2(myles.skin.sub, cluster='coarse_ids', top_n=5)

pdf("Rplots.pdf", width=38, height=12)
DotPlot_scCustom(spatial, features=de, group.by='seurat_clusters')
dev.off()

annot_labels = c(
    "0" = "F",
    "1" = "T",
    "2" = "PT",
    "3" = "M",
    "4" = "PT",
    "5" = "EC",
    "6" = "DLC_ALC",
    "7" = "PT",
    "8" = "M",
    "9" = "Principal",
    "10" = "PT",
    "11" = "ALC_DLC",
    "12" = "ALC_DLC",
    "13" = "Mural",
    "14" = "Intercalated",
    "15" = "EC",
    "16" = "Podocyte",
    "17" = "NOS",
    "18" = "T",
    "19" = "Plasma",
    "20" = "B",
    "21" = "Parietal",
    "22" = "M",
    "23" = "NOS",
    "24" = "EC",
    "25" = "T"
)

levels(spatial@meta.data$seurat_clusters) = annot_labels


saveRDS(spatial, 'spatial_sopa_newlineagelabels.rds')

svg("Rplots_i.svg", width=42, height=12)
DotPlot_scCustom(spatial, features=de, group.by='seurat_clusters') + get_theme(angle=45, size=25)
dev.off()

svg("Rplots_dimi.svg", width=13, height=12)
DimPlot_scCustom(spatial, group.by='seurat_clusters')
dev.off()

query.set <- as.SingleCellExperiment(spatial)
ref.set   <- as.SingleCellExperiment(myles.skin.sub)
col_data <- colData(ref.set)

print(table(spatial@meta.data$seurat_clusters))
print(sc.sub)
print(colnames(myles.skin.sub@meta.data))

print(table(spatial@meta.data$seurat_clusters))

##pred.grun2 <- SingleR(test = query.set, ref = ref.set, labels = col_data[, celltype], num.threads=16, de.method="t", method='cluster',  clusters=spatial$seurat_clusters)
###pred.grun2 <- SingleR(test = query.set, ref = ref.set, labels = col_data[, celltype], num.threads=16, de.method="wilcox", de.n=15, clusters=spatial$seurat_clusters, num.threads=10)
##
##
###saveRDS(pred.grun2, glue("sc_crossxen_singleR.rds"))
##
##
##png(glue("cross_xen_diagnostics.png"), width=1600, height=1000)
##print(plotScoreHeatmap(pred.grun2))
##dev.off()
##
##png(glue("cross_xen_delta.png"), width=1000, height=800)
##plotDeltaDistribution(pred.grun2)
##dev.off()
##
#### # fill the reference single cell annotation using the same column
#### merged.obj@meta.data = merged.obj@meta.data %>%
####     mutate(fine_ids = ifelse(is.na(fine_ids), !!sym(celltype), fine_ids))
##
##print('---')
##print(colnames(sc.sub@meta.data))
##
##pdf("seurat_cluster.pdf", width=13, height=6)
##DimPlot_scCustom(sc.sub, reduction='humap', group.by='seurat_clusters')
##dev.off()
##
##levels(sc.sub@meta.data$seurat_clusters) <- ifelse(is.na(pred.grun2$labels), 'NOS', pred.grun2$labels)
###sc.sub@meta.data$seurat_clusters = pred.grun2[match(rownames(sc.sub@meta.data), rownames(pred.grun2)),'pruned.labels']
##
##print('---')
##
##svg("singleR_cluster.svg", width=13, height=6)
##DimPlot_scCustom(sc.sub, reduction='humap', group.by='seurat_clusters')
##dev.off()
##
##
##sc.sub <- run_sctype(sc.sub, known_tissue_type="Kidney", slot="RNA")
##sc.sub@meta.data = sc.sub@meta.data %>% mutate(cell_type=ifelse(grepl("intercalated", sctype_classification), "Intercalated", sctype_classification))
##
##svg("sctype_cluster.svg", width=13, height=6)
##DimPlot_scCustom(sc.sub, reduction='humap', group.by='cell_type')
##dev.off()
##
