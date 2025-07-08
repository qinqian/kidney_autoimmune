library(sccomp)
library(ggthemes)
library(Seurat)
library(tidyverse)
library(tidyr)
library(forcats)
library(dplyr)
library(spatula)
library(sf)


tesera_obj <- readRDS("../phaseB_niches/tessera/vignettes/ICI_tesera.rds")

table(tesera_obj@meta.data$seurat_cluster)
table(tesera_obj@meta.data$condition)

## 0 (blue) - fibrosis
## 1 (orange) - inflammation (lymphocytes + macrophages mostly; some eosinophils)
## 2 (green) - distal tubules
## 3 (red) - proximal tubules
## 4 (purple) - *
## 5 (brown) - *
## 6 (pink) - vessels/muscle
## 7 (gray-ish) - I am not able to discern what this is highlighting*
## 8 (green/yellow) - glomeruli
levels(tesera_obj@meta.data$seurat_cluster) <- c("fibrosis", "inflammation", "distal tubules",
                                                 "proximal tubules", "NA1", "NA2",
                                                 "vessels/muscle", "NA3", "glomeruli")

obj = tesera_obj
maps1 = purrr::map(1:4, function(i) {
    sample_id_test=unique(obj[, obj@meta.data$condition=='Control']@meta.data$sample_ids)[i]
    ggplot(obj[, obj@meta.data$condition=='Control']@meta.data %>% filter(sample_ids==sample_id_test)) + 
    geom_sf(aes(geometry = shape, 
                fill = seurat_clusters),
            color = NA) + 
    theme_void(base_size = 16) + ggtitle(sample_id_test) + 
    coord_sf(expand = FALSE) + 
    scale_fill_tableau('Classic 10') +
    theme(legend.position="bottom")
})

maps2 = purrr::map(1:4, function(i) {
    sample_id_test=unique(obj[, obj@meta.data$condition=='Case  ']@meta.data$sample_ids)[i]
    ggplot(obj[, obj@meta.data$condition=='Case  ']@meta.data %>% filter(sample_ids==sample_id_test)) + 
    geom_sf(aes(geometry = shape, 
                #colour=seurat_clusters,
                #linewidth = 0, 
                color = NA,
                fill = seurat_clusters), color = NA) + 
    theme_void(base_size = 16) + 
    coord_sf(expand = FALSE) + ggtitle(sample_id_test) + 
    scale_fill_tableau('Classic 10')+  theme(legend.position="bottom")
})


res =
    tesera_obj |>
    sccomp_estimate(
      formula_composition = ~ condition,
      .sample = sample_ids, .cell_group = seurat_cluster,
      cores = 1, verbose=FALSE
    )

res = res |>  sccomp_test()

pdf("sccomp_ici_tesera.pdf")
print(res |> sccomp_boxplot(factor = "condition"))
dev.off()

pdf("sccomp_ici_tesera_1d.pdf", width=9.6, height=5)
print((res |>  plot_1D_intervals()) + theme_bw(base_size=15) + ylab("Niche group"))
dev.off()

pdf("sccomp_ici_tesera_2d.pdf")
print(res |>   plot_2D_intervals())
dev.off()


## add confounding factor
meta = read.csv("shruti_meta_clean.csv")
meta = as_tibble(meta[, c("slide_id", "age", "bmi", "sex", "eGFR_base", "malignancy")])
meta = meta %>% rename(sample_ids = slide_id)

tesera_obj@meta.data = tesera_obj@meta.data %>% left_join(meta, by="sample_ids")

res.mixed =
    tesera_obj |>
    sccomp_estimate(
      formula_composition = ~ condition + age + sex,
      .sample = sample_ids, .cell_group = seurat_cluster,
      bimodal_mean_variability_association = TRUE,
      cores = 1, verbose=FALSE
    )

res.mixed = res.mixed |> sccomp_test()
pdf("sccomp_ici_tesera_1d_mixed.pdf", width=10, height=8)
print((res.mixed |>  plot_1D_intervals()))
dev.off()
