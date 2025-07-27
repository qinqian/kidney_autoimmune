## We load the required packages
library(Seurat)
library(decoupleR)
# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(UCell)


data = readRDS("250711_niches.rds")

net <- decoupleR::get_progeny(organism = 'human', top = 200)

# Extract the normalized log-transformed counts

mat <- as.matrix(data[['RNA']]$data)

# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)


# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
                         tidyr::pivot_wider(id_cols = 'source', 
                                            names_from = 'condition',
                                            values_from = 'score') %>%
                         tibble::column_to_rownames(var = 'source') %>%
                         Seurat::CreateAssayObject(.)

# Change assay
Seurat::DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- Seurat::ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster = Seurat::Idents(data)) %>% mutate(condition=data$condition) %>%
      unite("cluster_cond", c(cluster, condition), remove=F) %>% select(-cluster, -condition) %>% 
      pivot_longer(cols = -cluster_cond, 
                   names_to = "source", 
                   values_to = "score") %>%
      dplyr::group_by(cluster_cond, source) %>%
      dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
                tidyr::pivot_wider(id_cols = cluster_cond,
                                   names_from = 'source',
                                   values_from = 'mean') %>%
                tibble::column_to_rownames(var = 'cluster_cond') %>%
                as.matrix()

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15, cluster_rows=F,
#                   cluster_rows = FALSE,treeheight_row = 0)
                   treeheight_col = 20, filename='Fig3b2.pdf')


# Change assay
Seurat::DefaultAssay(object = data) <- "RNA"


genes <- net%>% filter(source=='JAK-STAT')%>% select(target) %>% pull()
genes2 <- net%>% filter(source=='NFkB')%>% select(target)%>%pull()

signatures <- list(
    JAK_STAT = intersect(genes, rownames(data)),
    NFkB = intersect(genes2, rownames(data))
)

library(UCell)
library(Seurat)
library(ggplot2)
library(patchwork)

seurat.object <- AddModuleScore_UCell(data, 
                                      features=signatures, name=NULL)

pdf("Fig3b3.pdf", width=8, height=8)
p2 <- FeaturePlot(seurat.object, reduction = "umap", features = names(signatures), split.by='condition')
dev.off()
