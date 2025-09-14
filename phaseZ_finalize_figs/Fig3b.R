## We load the required packages
library(ggsci)
library(circlize)
library(scico)
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
                           .mor = 'weight', center=T,
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
## data <- Seurat::ScaleData(data)

data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@data

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
      as.data.frame() %>%
      dplyr::mutate(cluster = Seurat::Idents(data)) %>% mutate(condition=data$condition)


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
                   treeheight_col = 20, filename='Fig3b2.pdf')


# Change assay
Seurat::DefaultAssay(object = data) <- "RNA"


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
    JAK_STAT = intersect(genes, rownames(data)),
    NFkB = intersect(genes2, rownames(data)),
    OXIDATIVE_PHOS = H.pathway_list[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]
)

library(UCell)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sf)

seurat.object <- AddModuleScore_UCell(data, 
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
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
p22 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
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
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
p32 = ggplot() + 
    geom_sf(data = cbind(seurat.object@meta.data, FetchData(seurat.object, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS2'),
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_jak), low = 'white', high = '#832424')+ theme(aspect.ratio = 1, legend.position='right')
pall<-p21+p22+p31+p32+plot_layout(ncol=2, guides='collect')
ggsave("Fig3b3.svg", plot = pall, width=10, height=8)

## sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
## lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
## imm.niche <- readRDS("250711_niches.rds")

## sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(lennard_label=lennard.subtype@meta.data$lennard_label)
## sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(tile_label = imm.niche@meta.data[match(sc.niche$obj@meta.data$tile_id, rownames(imm.niche@meta.data)), 'niche_label'])
## sc.niche.heatmap<-sc.niche$obj@meta.data %>% count(tile_label, lennard_label) %>% pivot_wider(names_from=lennard_label, values_from=n, values_fill=0) %>% filter(!is.na(tile_label))
## sc.niche.heatmap = as.data.frame(sc.niche.heatmap)
## rownames(sc.niche.heatmap) = sc.niche.heatmap$tile_label
## sc.niche.heatmap = sc.niche.heatmap[,-1]
## sc.niche.heatmap <- scale(sc.niche.heatmap)
## col_fun <- colorRamp2(c(-2, 0, 2), scico(3, palette = "vik"))  # "vik" is diverging

## sc.niche.heatmap <- sc.niche.heatmap[rownames(sc.niche.heatmap)!="Skeletal Muscle",]

## niche_cols = pal_npg("nrc")(9)

## names(niche_cols) <- rownames(sc.niche.heatmap)
## patient <- read_tsv("my_sample_col.tsv")
## imm.niche@meta.data <- imm.niche@meta.data %>% mutate(patient_id=gsub("__2.*", "", sample_id))
## lennard.subtype@meta.data = lennard.subtype@meta.data %>% mutate(patient_id=gsub("output-XETG00150__0018462__|output-XETG00392__0045655__", "", gsub("__2.*", "", sample_ids)))

## p21 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
##             aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
##     geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==cont2),
##                aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
##     scale_fill_gradient(low = 'white', high = '#832424')+
##     guides(fill = guide_legend(nrow = 5),
##            color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='bottom',
##                                          legend.text = element_text(size = 3),      # legend labels
##                                          legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
##     coord_sf(xlim=c(6600, 7200), ylim=c(600, 900), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')
## dev.off()
