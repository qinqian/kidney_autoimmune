# Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Required packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview", "DOSE"))


library(Seurat)
library(tidyverse)
library(fgsea)
library(data.table)

imm.niche <- readRDS("250711_niches.rds")

plot.list <- list()
for (cl in unique(imm.niche$niche_label)) {
    imm.niche.cl = subset(imm.niche, subset = niche_label == cl)
    if (cl == 'Skeletal Muscle') {
        next
    }
    markers <- presto::wilcoxauc(imm.niche.cl, 'condition')%>% filter(group=='Case') %>% filter()
    markers$group = cl
    plot.list[[cl]] <- markers
}


plot.df <- plot.list %>% rbindlist()

library(ggrepel)
p1 = ggplot() +
    geom_point(data=plot.df, aes(x = logFC, y = -log10(padj), color=padj <= 1e-50 & auc >= 0.6 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(padj <= 1e-50 & auc >= 0.6 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 12) + facet_wrap(~group) + 
    ## labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme(legend.position = "none")


library(msigdbr)
library(fgsea)
library(tidyverse)
library(msigdbr)

genesets <- msigdbr() %>% filter(gs_collection=='C2')

pathway_list <- genesets %>% filter(grepl("KEGG", gs_name)) %>% 
  dplyr::group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

test.imm.niche <- fgsea(pathways=pathway_list,
      stats=setNames(plot.list[['Immune']]$logFC, plot.list[['Immune']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.fib.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Fibrosis & Interstitium']]$logFC, plot.list[['Fibrosis & Interstitium']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

pdf("test.pdf", width=12, height=6)
p2 = ggplot(test.imm.niche%>%filter(padj<=1e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 6)  + ggtitle("Immune GSEA")  
p3 = ggplot(test.fib.niche%>%filter(padj<=1e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 6)  + ggtitle("Fibrosis & Interstitium GSEA")  
print((p1+(p2/p3))+patchwork::plot_layout(width=c(3.2, 1.2)))
dev.off()
