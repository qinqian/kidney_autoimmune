## https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(dplyr); library(condiments);
  library(RColorBrewer); library(scales)
  library(viridis); library(UpSetR)
  library(pheatmap); library(msigdbr)
  library(fgsea); library(knitr)
  library(ggplot2); library(gridExtra)
  library(tradeSeq); library(cowplot)
})
theme_set(theme_classic())
library(tidyverse)


import_TFGB <- function()
{
    raw <- GEOquery::getGEOSuppFiles("GSE114687", baseDir = tempdir(), 
        filter_regex = "GSE114687_pseudospace_cds")
    raw_loc <- rownames(raw)
    system(paste0("gunzip", " ", raw_loc))
    file <- stringr::str_remove(raw_loc, ".gz")
    cds <- readRDS(file)
    file.remove(file)
    counts <- cds@assayData$exprs
    phenoData <- pData(cds@phenoData)
    rd <- SimpleList(tSNEorig = cbind(cds@phenoData@data$TSNE.1, 
        cds@phenoData@data$TSNE.2))
    filt <- apply(counts, 1, function(x) {
        sum(x >= 2) >= 15
    })
    counts <- counts[filt, ]
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), 
        colData = phenoData, reducedDims = rd)
    return(sce)
}

library(GEOquery)
library(R.utils)
#####library(gunzip)

library(condiments)
library(slingshot)

imm.niche <- readRDS("250711_niches.rds")

data("tgfb", package = "bioc2021trajectories")

## df <- bind_cols(
##   as.data.frame(reducedDims(tgfb)$UMAP),
##   as.data.frame(colData(tgfb)[, !colnames(colData(tgfb)) == "slingshot"])
##   )
## ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = treatment_id)) +
##   geom_point(size = .7) +
##   scale_color_brewer(palette = "Accent") +
##   labs(col = "Treatment")

## ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = spatial_id)) +
##   geom_point(size = .7) +
##   scale_color_brewer(palette = "Dark2") +
##   labs(col = "Spatial ID")

library(Seurat)
library(harmony)

imm.niche.sub <- subset(imm.niche, subset = niche_label %in% c("Fibrosis & Interstitium", "Immune"))

imm.niche.sub <- NormalizeData(imm.niche.sub, scale.factor=median(colSums(imm.niche.sub[['RNA']]$counts)))

set.seed(9)
imm.niche.sub <- ScaleData(imm.niche.sub) %>% RunPCA(dims=1:50) %>% RunHarmony(c("condition", "sample_id")) %>% FindNeighbors() %>% RunUMAP(dims=1:50)

imm.niche.sce = as.SingleCellExperiment(imm.niche.sub)

tgfb <- condiments::imbalance_score(imm.niche.sce, conditions = imm.niche.sce$condition,
                                    k = 20, smooth = 40, dimred = "UMAP")

tgfb <- slingshot(tgfb, reducedDim = 'UMAP', clusterLabels = tgfb$niche_label_fine,
                  start.clus = 7, extend="n", reweight = FALSE, reassign = FALSE)

tgfb@int_metadata$slingshot <- tgfb$slingshot
mst <- slingMST(tgfb, as.df = TRUE)

df <- data.frame(Embeddings(imm.niche.sub, 'umap'), Cluster=imm.niche.sub@meta.data[, "niche_label_fine"])

df$scores <- tgfb$scores$scaled_scores

ggplot(df, aes(x = umap_1, y = umap_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Score")
ggplot(df, aes(x = umap_1, y = umap_2, col = Cluster)) +
  geom_point(size = .7, alpha = .3) +
  ## scale_color_brewer(palette = "Dark2") +
  labs(col = "Spatial ID") +
  geom_path(data = mst, col = "black", size = 1.5) +
  geom_point(data = mst, aes(col = Cluster), size = 5)
dev.off()

set.seed(100)
topologyTest(tgfb, conditions = tgfb$condition, rep = 50)

sdss <- slingshot_conditions(tgfb, tgfb$condition, approx_points = FALSE,
                             extend = "n", reweight = FALSE, reassign = FALSE)

msts <- lapply(sdss, slingMST, as.df = TRUE) %>%
  bind_rows(.id = "Batch") %>%
  arrange(Batch)

pdf("skeleton.pdf")
ggplot(df, aes(x = umap_1, y = umap_2, col = Cluster)) +
  geom_point(size = .7, alpha = .1) +
  ## scale_color_brewer(palette = "Accent") +
  geom_point(data = msts, size = 3) +
  geom_path(data = msts, aes(group = interaction(Lineage, Batch)),
            size = 2)
dev.off()

lineages <- lapply(sdss, slingCurves, as.df = TRUE) %>%
  bind_rows(.id = "Batch") %>%
  arrange(Order)


position <- data.frame(
  "tSNE1" = c(40, -30, 45),
  "tSNE2" = c(50, -50, -50),
  "Batch" = "H2122A",
  "text" = paste0("Lineage ", 1:3)
)

pdf("path.pdf")
ggplot() + 
  geom_point(data = df, aes(x = umap_1, y = umap_2, col = Cluster), size = .7, alpha = .2) +
  geom_path(data = lineages, size = 1.5, aes(x=umap_1, y=umap_2, group = interaction(Lineage, Batch)))
dev.off()


mapping <- matrix(rep(1:4, each = 2), nrow = 4, ncol = 2, byrow = TRUE)

sds <- merge_sds(sdss[[1]], sdss[[2]], 
                 condition_id = names(sdss), mapping = mapping)

df <- full_join(
  df %>% mutate(cells=rownames(df)) %>% 
    select(cells, umap_1, umap_2, Cluster),
  slingPseudotime(sds) %>%
    as.data.frame() %>%
    mutate(cells = rownames(.))
) %>%
    pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "pst")

pdf("differential_progression.pdf")
ggplot(df, aes(x = pst)) +
  geom_density(alpha = .4, aes(fill = Cluster), col = "transparent") +
  geom_density(aes(col = Cluster), fill = "transparent", size = 1.5) +
  guides(col = FALSE) +
  ## scale_fill_brewer(palette = "Accent") +
  ## scale_color_brewer(palette = "Accent") +
  labs(x = "Pseudotime", fill = "Type") +
  facet_wrap(~ Curve, scales = "free_x")
dev.off()

progressionTest(sds, conditions = imm.niche.sub$condition, lineages = TRUE)

df <- bioc2021trajectories::sling_reassign(sds) %>% 
  as.data.frame() %>%
  mutate(cells = rownames(.)) %>%
  dplyr::rename("Lineage1" = V1, "Lineage2" = V2, "Lineage3" = V3) %>%
  pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "weights") %>%
  full_join(df) %>%
  group_by(cells) %>%
  select(-pst) %>% 
  mutate(weights = weights / sum(weights)) %>% 
  ungroup()

pdf("differential_progression2.pdf")
ggplot(df %>% group_by(Cluster, Curve) %>% 
         summarise(weights = mean(weights), .groups = NULL),
       aes(x = Curve, fill = Cluster, y = weights)) +
  geom_col(position = "dodge") +
  ## scale_fill_brewer(palette = "Accent") +
  theme(legend.position = c(.7, .7)) +
  labs(x = "", y = "Mean weight")
dev.off()

## ggplot(df %>% pivot_wider(names_from = "Curve", values_from = "weights"),
##        aes(x = Lineage1, y = Lineage3)) +
##   geom_hex() +
##   scale_fill_viridis_c(direction = -1) +
##   facet_wrap(~Cluster, scales = "free") +
##   geom_abline(slope = -1, intercept = 1, linetype = "dotted") +
##   geom_abline(slope = -1, intercept = 2/3, linetype = "dotted") +
##   geom_abline(slope = -1, intercept = 1/3, linetype = "dotted") +
##   annotate("text", x = .53, y = .53, label = "w3 = 0", angle = -52) +
##   annotate("text", x = .62, y = .1, label = "w3 = 1/3", angle = -52) +
##   annotate("text", x = .14, y = .14, label = "w3 = 2/3", angle = -52) +
##   theme(legend.position = "bottom") +
##   labs(x = "Weights for Lineage 1 (w1)", y = "Weights for Lineage 2 (w2)",
##        fill = "counts per hexagon")
