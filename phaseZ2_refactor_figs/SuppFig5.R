library(Seurat)
library(tidyverse)
library(fgsea)
library(data.table)

imm.niche <- readRDS("../phaseZ_finalize_figs/250711_niches.rds")

plot.list <- list()
for (cl in unique(imm.niche$niche_label)) {
    imm.niche.cl = subset(imm.niche, subset = niche_label == cl)
    if (cl == 'Skeletal Muscle') {
        next
    }
    markers <- presto::wilcoxauc(imm.niche.cl, 'condition')%>% filter(group=='Case')
    markers$group = cl
    plot.list[[cl]] <- markers
}


plot.df <- plot.list %>% rbindlist()

library(ggrepel)

p1 = ggplot() +
    geom_point(data=plot.df %>% filter(!(group %in% c('Fibrosis & Interstitium', 'Immune', "Collecting Duct"))), aes(x = logFC, y = -log10(padj), color=padj <= 0.05 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(!(group %in% c('Fibrosis & Interstitium', 'Immune', "Collecting Duct"))) %>% filter(padj <= 0.05 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 13.5) + facet_wrap(group~., ncol=6) + 
    theme(legend.position = "none")

print(p1)
dev.off()

library(msigdbr)
library(fgsea)
library(tidyverse)
library(msigdbr)

H.genesets <- msigdbr() %>% filter(gs_collection=='H')
print(H.genesets)

H.pathway_list <- H.genesets %>% 
  dplyr::group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

genesets <- msigdbr() %>% filter(gs_collection=='C2')

pathway_list <- genesets %>% filter(grepl("KEGG", gs_name)) %>% 
  dplyr::group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

test.vessel.niche <- fgsea(pathways=pathway_list,
      stats=setNames(plot.list[['Vessel']]$logFC, plot.list[['Vessel']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.itub.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Injured Proximal Tubule']]$logFC, plot.list[['Injured Proximal Tubule']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.tub.niche <- fgsea(pathways=pathway_list,
      stats=setNames(plot.list[['Proximal Tubule']]$logFC, plot.list[['Proximal Tubule']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.Glomerulus.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Glomerulus']]$logFC, plot.list[['Glomerulus']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.limb.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Thick Ascending Limb']]$logFC, plot.list[['Thick Ascending Limb']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.distal.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Distal Convoluted Tubule']]$logFC, plot.list[['Distal Convoluted Tubule']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)


test.vessel.niche.H <- fgsea(pathways=H.pathway_list,
      stats=setNames(plot.list[['Vessel']]$logFC, plot.list[['Vessel']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.itub.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Injured Proximal Tubule']]$logFC, plot.list[['Injured Proximal Tubule']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.tub.niche.H <- fgsea(pathways=H.pathway_list,
      stats=setNames(plot.list[['Proximal Tubule']]$logFC, plot.list[['Proximal Tubule']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.Glomerulus.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Glomerulus']]$logFC, plot.list[['Glomerulus']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.limb.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Thick Ascending Limb']]$logFC, plot.list[['Thick Ascending Limb']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.distal.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Distal Convoluted Tubule']]$logFC, plot.list[['Distal Convoluted Tubule']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)


write_tsv(test.itub.niche, 'itub.niche.fgsea.tsv')
write_tsv(test.tub.niche, 'tub.niche.fgsea.tsv')
write_tsv(test.limb.niche, 'limb.niche.fgsea.tsv')
write_tsv(test.distal.niche, 'distal.niche.fgsea.tsv')
write_tsv(test.Glomerulus.niche, 'Glomerulus.niche.fgsea.tsv')


write_tsv(test.itub.niche.H, 'H.itub.niche.fgsea.tsv')
write_tsv(test.tub.niche.H, 'H.tub.niche.fgsea.tsv')
write_tsv(test.limb.niche.H, 'H.limb.niche.fgsea.tsv')
write_tsv(test.distal.niche.H, 'H.distal.niche.fgsea.tsv')
write_tsv(test.Glomerulus.niche.H, 'H.Glomerulus.niche.fgsea.tsv')


pdf("SuppFig5_v3.pdf", width=31, height=16)
p2 = ggplot(test.itub.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Injured Proximal Tubule KEGG GSEA")  +ylab("")
p3 = ggplot(test.tub.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Proximal Tubule KEGG GSEA")    +ylab("")
p4 = ggplot(test.limb.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Thick Ascending Limb KEGG GSEA")    +ylab("")
p5 = ggplot(test.distal.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Distal Convoluted Tubule KEGG GSEA")  
p6 = ggplot(test.Glomerulus.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Glomerulus KEGG GSEA")    +ylab("")
p7 = ggplot(test.vessel.niche%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Vessel KEGG GSEA")    +ylab("")
p20 = ggplot(test.itub.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Injured Proximal Tubule Hallmarks GSEA")    +ylab("")
p30 = ggplot(test.tub.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Proximal Tubule Hallmarks GSEA")    +ylab("")
p40 = ggplot(test.limb.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Thick Ascending Limb Hallmarks GSEA")    +ylab("")
p50 = ggplot(test.distal.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Distal Convoluted Tubule Hallmarks GSEA")    +ylab("")
p60 = ggplot(test.Glomerulus.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Glomerulus Hallmarks GSEA")    +ylab("")
p70 = ggplot(test.vessel.niche.H%>%filter(padj<=5e-2)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Vessel Hallmarks GSEA")    +ylab("")
print(p1/(p5+p6+p2+p3+p4+p7+p50+p60+p20+p30+p40+p70+patchwork::plot_layout(guides='collect', ncol=6, width=c(2, 2, 2, 2, 2, 2)))+patchwork::plot_layout(nrow=2, height=c(1.6, 3.5)))
dev.off()


