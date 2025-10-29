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
    geom_point(data=plot.df %>% filter(group %in% c('Fibrosis & Interstitium', 'Immune', "Collecting Duct", "Glomerulus")), aes(x = logFC, y = -log10(padj), color=padj <= 0.05 & auc >= 0.6 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(group %in% c('Fibrosis & Interstitium', 'Immune', "Collecting Duct", "Glomerulus")) %>% filter(padj <= 0.05 & auc >= 0.6 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 12) + facet_wrap(~group, ncol=1) + 
    theme(legend.position = "none")

p1a = ggplot() +
    geom_point(data=plot.df %>% filter(group %in% c("Collecting Duct")), aes(x = logFC, y = -log10(padj), color=padj <= 0.05 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(group %in% c("Collecting Duct")) %>% filter(padj <= 0.05 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 12) + #facet_wrap(~group, ncol=1) + 
    theme(legend.position = "none")

p1b = ggplot() +
    geom_point(data=plot.df %>% filter(group %in% c('Fibrosis & Interstitium')), aes(x = logFC, y = -log10(padj), color=padj <= 0.05 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(group %in% c('Fibrosis & Interstitium')) %>% filter(padj <= 0.05 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 12) + #facet_wrap(~group, ncol=1) + 
    theme(legend.position = "none")

p1c = ggplot() +
    geom_point(data=plot.df %>% filter(group %in% c('Immune')), aes(x = logFC, y = -log10(padj), color=padj <= 0.05 & abs(logFC) >= 1), alpha = 0.6, size = 1) +
    geom_text_repel(data=plot.df %>% filter(group %in% c('Immune')) %>% filter(padj <= 0.05 & abs(logFC) >= 1), aes(x=logFC, y = -log10(padj), label = feature), max.overlaps = 100, size = 3) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_minimal(base_size = 12) + #facet_wrap(~group, ncol=1) + 
    theme(legend.position = "none")


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

test.colduct.niche <- fgsea(pathways=pathway_list,
      stats=setNames(plot.list[['Collecting Duct']]$logFC, plot.list[['Collecting Duct']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)
test.Glomerulus.niche <- fgsea(pathways=pathway_list,
      stats=sort(setNames(plot.list[['Glomerulus']]$logFC, plot.list[['Glomerulus']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)


test.imm.niche.H <- fgsea(pathways=H.pathway_list,
      stats=setNames(plot.list[['Immune']]$logFC, plot.list[['Immune']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)
test.fib.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Fibrosis & Interstitium']]$logFC, plot.list[['Fibrosis & Interstitium']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.colduct.niche.H <- fgsea(pathways=H.pathway_list,
      stats=setNames(plot.list[['Collecting Duct']]$logFC, plot.list[['Collecting Duct']]$feature),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)

test.Glomerulus.niche.H <- fgsea(pathways=H.pathway_list,
      stats=sort(setNames(plot.list[['Glomerulus']]$logFC, plot.list[['Glomerulus']]$feature)),
      eps      = 0.0,
      minSize  = 15,
      maxSize  = 500) %>% arrange(padj)


write_tsv(test.imm.niche, 'imm.niche.fgsea.tsv')
write_tsv(test.fib.niche, 'fib.niche.fgsea.tsv')
write_tsv(test.imm.niche.H, 'H.imm.niche.fgsea.tsv')
write_tsv(test.fib.niche.H, 'H.fib.niche.fgsea.tsv')

write_tsv(test.colduct.niche, 'colduct.niche.fgsea.tsv')
write_tsv(test.Glomerulus.niche, 'Glomerulu.niche.fgsea.tsv')

write_tsv(test.colduct.niche.H, 'H.imm.niche.fgsea.tsv')
write_tsv(test.Glomerulus.niche.H, 'H.Glomerulus.niche.fgsea.tsv')



pdf("Fig4.pdf", width=28, height=12)
p2 = ggplot(test.imm.niche%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Immune GSEA")  
p3 = ggplot(test.fib.niche%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Fibrosis & Interstitium GSEA")  
p4 = ggplot(test.colduct.niche%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Collecting Duct GSEA")  
p5 = ggplot(test.Glomerulus.niche%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Glomerulus GSEA")  
p6 = ggplot(test.imm.niche.H%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Immune GSEA")  
p7 = ggplot(test.fib.niche.H%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Fibrosis & Interstitium GSEA")  
p8 = ggplot(test.colduct.niche.H%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Collecting Duct GSEA")  
p9 = ggplot(test.Glomerulus.niche.H%>%filter(padj<=5e-2), aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 9)  + ggtitle("Glomerulus GSEA")  
print(p1+(p2/p6)+(p3/p7)+(p4/p8)+patchwork::plot_layout(ncol=4, width=c(7, 2, 2, 2.3)))
dev.off()


net <- decoupleR::get_progeny(organism = 'human', top = 200)

mat <- as.matrix(imm.niche[['RNA']]$data)

# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', center=T,
                           minsize = 5)


# Extract mlm and store it in pathwaysmlm in data
imm.niche[['pathwaysmlm']] <- acts %>%
                         tidyr::pivot_wider(id_cols = 'source', 
                                            names_from = 'condition',
                                            values_from = 'score') %>%
                         tibble::column_to_rownames(var = 'source') %>%
                         Seurat::CreateAssayObject(.)

# Change assay
Seurat::DefaultAssay(object = imm.niche) <- "pathwaysmlm"

imm.niche@assays$pathwaysmlm@data <- imm.niche@assays$pathwaysmlm@data

# Extract activities from object as a long dataframe
df <- t(as.matrix(imm.niche@assays$pathwaysmlm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster = Seurat::Idents(imm.niche)) %>% mutate(condition=imm.niche$condition)


pathway_cols <- df %>%
  select(-condition, -cluster) %>%
  colnames()

# Pivot longer for tidyverse operations
df_long <- df %>%
  pivot_longer(cols = all_of(pathway_cols),
               names_to = "Pathway",
               values_to = "Activity")

# Compute log2 fold change per pathway
log2fc_df <- df_long %>%filter(!(Pathway %in% c("Hypoxia", "Androgen", "EGFR", "Estrogen"))) %>% filter(!(cluster%in%c('Skeletal Muscle'))) %>% 
  group_by(cluster, Pathway, condition) %>%
  summarise(Mean_Activity = mean(Activity) + 1e-6, .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = Mean_Activity) %>%
  mutate(Control = ifelse(is.na(Control), 1e-6, Control)) %>% 
  mutate(FC = `Case` - `Control`)%>% 
  select(cluster, Pathway, FC) %>% tidyr::pivot_wider(id_cols=c(cluster), names_from=Pathway, values_from="FC") %>%
  tibble::column_to_rownames(var = 'cluster') %>% 
  as.matrix()

# Compute Wilcoxon Rank-Sum Test per pathway
wilcox_pval_df <- df_long %>% filter(!(Pathway %in% c("Hypoxia", "Androgen", "EGFR", "Estrogen"))) %>% filter(!(cluster%in%c('Skeletal Muscle'))) %>% 
    unite("cluster_condition", c(Pathway, cluster), sep="____", remove=F) %>% 
    group_by(cluster_condition) %>%
    summarise(Wilcoxon_p = wilcox.test(Activity ~ condition)$p.value, .groups = 'drop')

# Optional: Adjust p-values (FDR)
wilcox_pval_df <- wilcox_pval_df %>%
  mutate(FDR = p.adjust(Wilcoxon_p, method = "fdr"))

wilcox_pval_df <- wilcox_pval_df %>% separate(cluster_condition, c("cluster", "condition"), sep="____")  %>% tidyr::pivot_wider(id_cols=c(cluster), names_from=condition, values_from="FDR") %>% tibble::column_to_rownames(var = 'cluster') %>% 
  as.matrix()


# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

# Convert p-value matrix into significance stars
signif_labels <- ifelse(wilcox_pval_df < 0.001, "***",
                  ifelse(wilcox_pval_df < 0.01, "**",
                  ifelse(wilcox_pval_df < 0.05, "*", "")))

# Make a list of annotations for each cell
annotation_text <- matrix(signif_labels, nrow = nrow(log2fc_df))

# Plot
## pheatmap::pheatmap(mat = top_acts_mat,
pheatmap::pheatmap(mat = log2fc_df,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15, display_numbers = annotation_text,
                   cellheight = 15, cluster_rows=F,
#                   cluster_rows = FALSE,treeheight_row = 0)
                   treeheight_col = 20, filename='Fig4_part2.pdf')


# Change assay
Seurat::DefaultAssay(object = imm.niche) <- "RNA"


genes <- net%>% filter(source=='JAK-STAT')%>% select(target) %>% pull()
genes2 <- net%>% filter(source=='NFkB')%>% select(target)%>%pull()

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

H.pathway_list[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]

signatures <- list(
    OXIDATIVE_PHOS = intersect(rownames(imm.niche), unlist(test.colduct.niche.H %>% filter(padj<=0.05) %>% filter(pathway=='HALLMARK_OXIDATIVE_PHOSPHORYLATION') %>% select(leadingEdge) %>% pull())),
    JAK_STAT = intersect(rownames(imm.niche), unlist(test.fib.niche %>% filter(padj<=0.05) %>% filter(pathway=='KEGG_JAK_STAT_SIGNALING_PATHWAY') %>% select(leadingEdge) %>% pull())),
    chemokine = intersect(rownames(imm.niche), unlist(test.imm.niche %>% filter(padj<=0.05) %>% filter(pathway=='KEGG_CHEMOKINE_SIGNALING_PATHWAY') %>% select(leadingEdge) %>% pull())),
    IFNG = intersect(rownames(imm.niche), unlist(test.imm.niche.H %>% filter(padj<=0.05) %>% filter(pathway=='HALLMARK_INTERFERON_GAMMA_RESPONSE') %>% select(leadingEdge) %>% pull()))
)

library(UCell)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sf)

seurat.object <- AddModuleScore_UCell(imm.niche, 
                                      features=signatures, name=NULL)

case1 <- "BS21-N65682A2"
case2 <- "BS23_52206A2"
cont1 <- "BS22_12012A1"
cont2 <- "BS2_61615A1"

library(patchwork)

feature = 'OXIDATIVE_PHOS'
feature2 = 'oxidative_phos'
df1 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2')
max_ox = max(df1$oxidative_phos, df2$oxidative_phos)
p21 = ggplot() + 
    geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')+ggtitle("OXIDATIVE PHOS") + coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE)
p22 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right') + coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE)
feature = 'JAK_STAT'
feature2 = 'jak_stat'
df1 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2')
max_jak = max(df1$jak_stat, df2$jak_stat)
p31 = ggplot() + 
    geom_sf(data = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) +  coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) +
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')+ggtitle("JAK STAT")
p32 = ggplot() + 
    geom_sf(data = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2'),
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) + coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + 
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
feature = 'chemokine'
feature2 = 'chemokine'
max_jak = max(df1$chemokine, df2$chemokine)
p41 = ggplot() + 
    geom_sf(data = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) + coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + 
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')+ggtitle("Chemokine")
p42 = ggplot() + 
    geom_sf(data = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2'),
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) + coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) +
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
feature = 'IFNG'
feature2 = 'ifng'
df1 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2')
max_ox = max(df1$ifng, df2$ifng)
p51 = ggplot() + 
    geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) + coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')+ggtitle("IFNG")
p52 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=15) + coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
part2<-(p21/p22)|(p31/p32)|(p41/p42)|(p51/p52)
ggsave("Fig4_part2.svg", plot = part2, width=25, height=8)


pdf("Fig4.pdf", width=15.5, height=12.5)
p2 = ggplot(test.imm.niche%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Immune KEGG GSEA")  +ylab("")
p3 = ggplot(test.fib.niche%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Fibrosis & Interstitium KEGG GSEA")   +ylab("")
p4 = ggplot(test.colduct.niche%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Collecting Duct KEGG GSEA")   +ylab("")
p5 = ggplot(test.Glomerulus.niche%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Glomerulus GSEA")  
p6 = ggplot(test.imm.niche.H%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Immune Hallmarks GSEA")   +ylab("")
p7 = ggplot(test.fib.niche.H%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Fibrosis & Interstitium GSEA")   +ylab("")
p8 = ggplot(test.colduct.niche.H%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Collecting Duct Hallmarks GSEA")   +ylab("")
p9 = ggplot(test.Glomerulus.niche.H%>%filter(padj<=5e-2)%>%slice_head(n=15)%>%mutate(pathway=factor(pathway, levels=as.character(pathway))), aes(x = abs(NES), y = gsub("HALLMARK_|KEGG_", "", pathway), fill = NES > 0)) +
    geom_col() + 
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    geom_vline(xintercept = 0, color = "black") +
    theme_minimal(base_size = 14)  + ggtitle("Glomerulus Hallmarks GSEA")   +ylab("")
##print((p1a|p4|p8)/(p1b|p3|p7)/(p1c|p2|p6)) #+patchwork::plot_layout(ncol=3, width=c(2, 2, 2)))) #/part2)
print((p4|p8)/(p3|p7)/(p2|p6)) #+patchwork::plot_layout(ncol=3, width=c(2, 2, 2)))) #/part2)
# plot_annotation(
#    tag_levels = list(c("A", "B", "C", "D", "E", 'F', 'G', 'H', 'I', 'J', 'K', 'L','M','N')), 
# ) &  theme(plot.tag = element_text(face = "bold", size = 25)))
dev.off()

pdf("Fig4b.pdf", width=15.5, height=12.5)
print(part2)
dev.off()
