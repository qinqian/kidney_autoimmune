suppressPackageStartupMessages({
    library(Seurat)
    library(scCustomize)
    library(spatula)
    library(ComplexHeatmap)
    library(circlize)
    library(rcna)
    library(circlize)
    library(ggrastr)
    library(tidyverse)
    library(scico)
    library(circlize)    
    library(ggsci)
})
library(ggrastr)

get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_blank(), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal")
    defined_theme
}

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}


heatmap_col_fun <- colorRamp2(c(-2, 0, 2), scico(3, palette = "vik"))  # "vik" is diverging


celltype_colors <- c(
    # cell groups
  "T" = "#1b7837",
  "Regulatory T" = "#238b45",
  "Proliferating T cell" = "#006d2c",
  "Immune (LowQ)" = "#ccebc5",

  # B/Plasma (purples)
  "B" = "#d0d1e6",
  "Plasma" = "#54278f",

  # Myeloid (yellows & oranges)
  "Tissue Myeloid" = "#ffffcc",
  "Monocyte" = "#ffeda0",
  "Inflammatory Myeloid" = "#feb24c",
  "cDC1" = "#fee0d2",
  "pDC" = "#fee0d2",
  "Basophil" = "#fee0d2",    
    

  # Endothelial (blues)
  "Endothelial" = "#c6dbef",

  # Kidney cells (reds & browns)
  "Podocyte" = "#fbb4ae",
  "Parietal" = "#f768a1",
  "Proximal Tubule" = "#c51b8a",
  "Thick Ascending Limb" = "#7a0177",
  "Thin Ascending Limb" = "#fdae6b",
  "Thin Descending Limb" = "#e6550d",
  "Interstitial" = "#a63603",
  "Distal Convoluted Tubule" = "#6e016b",
  "Collecting Duct-IC" = "#6e016b",
  "Collecting Duct-PC" = "#6e016b",
  "Connecting Tubule" = "#6e016b"    
)

niche_cols = c(
  # Immune (purples)
  "Immune"              = "#6a51a3",  # medium-deep purple

  # Myeloid / Skeletal Muscle (yellows/oranges)
  #"Skeletal Muscle"      = "#fdcc8a",  # warm orange-yellow

  # Endothelial / Fibrosis (blues)
  "Fibrosis & Interstitium" = "#6baed6",  # light sky blue
  "Vessel"                 = "#08519c",  # deep blue

  # Kidney cells (reds → pinks → violets, smooth gradient)
  "Injured Proximal Tubule"  = "#fcbba1",  # soft pink
  "Proximal Tubule"          = "#fb6a4a",  # salmon
  "Thick Ascending Limb"     = "#ef3b2c",  # tomato red
  "Distal Convoluted Tubule" = "#a50f15",  # dark red
  "Collecting Duct"          = "#67000d",  # very dark red

  # Duct / Specialized (violet)
  "Glomerulus" = "#fdcc8a"  # violet (distinct from immune purple)
)


sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")

lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
imm.niche <- readRDS("250711_niches.rds")

niche.merge = imm.niche 
obj.merge = lennard.subtype
orig.merge = readRDS("all_KPMP_integrate_singlet_umap_umapnn_labels_umap.rds")

meta = read.csv("shruti_meta_clean (3).csv")
input_meta = meta[,c('slide_id', 'age', 'sex', 'case_ctrl', 'ICPi',  'malignancy', 'eGFR_base')] %>% arrange(case_ctrl)

cells_to_keep <- colnames(orig.merge)[orig.merge$tech=='xenium']
orig.merge.xen <- subset(orig.merge, cells = cells_to_keep)


# Main cells
main_cells <- colnames(orig.merge.xen)

# Assay cells
assay_cells <- colnames(orig.merge.xen@assays$RNA)  # or whichever assay you're using

# PCA cells
pca_cells <- rownames(orig.merge.xen@reductions$pca@cell.embeddings)

# Graph cells
graph_cells <- colnames(orig.merge.xen[['humap_fgraph']])  # adapt if using another graph

# Active identity names
ident_cells <- names(Idents(orig.merge.xen))

# Check mismatches
length(setdiff(assay_cells, main_cells))
length(setdiff(pca_cells, main_cells))
length(setdiff(graph_cells, main_cells))
length(setdiff(ident_cells, main_cells))

orig.merge.xen@meta.data  = orig.merge.xen@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1))
orig.merge.xen@meta.data  = orig.merge.xen@meta.data %>% mutate(case_ctrl=str_trim(input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl']))
orig.merge.xen@meta.data$case_ctrl_num = as.numeric(factor(str_trim(orig.merge.xen@meta.data$case_ctrl), levels=c("Control", "Case")))


obj.merge@meta.data  = obj.merge@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1))
obj.merge@meta.data  = obj.merge@meta.data %>% mutate(case_ctrl=str_trim(input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl']))
obj.merge@meta.data$case_ctrl_num = as.numeric(factor(str_trim(obj.merge@meta.data$case_ctrl), levels=c("Control", "Case")))
obj.merge@meta.data$cell_label = gsub(" Cell", "", obj.merge@meta.data$lennard_label)
obj.merge@meta.data = obj.merge@meta.data %>% mutate(cell_label = ifelse(cell_label=='Immune', 'Immune (LowQ)', cell_label))


#orig.baysor <- readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")
orig.baysor <- readRDS("kidney_orig_seg_merged.rds")
orig.baysor@meta.data <- orig.baysor@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)
lennard.subtype@meta.data <- lennard.subtype@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)
xy <- Embeddings(orig.baysor, 'spatial')[match(lennard.subtype@meta.data$uniq_id, orig.baysor@meta.data$uniq_id),]
rm(orig.baysor)

markers = c(            "CALB1", "HSD11B2", "SCNN1G",

    "TNXB", "COL5A1", "PDGFRA",
                        "PLA2R1", "PODXL", "WT1",
            "CXCL9", "SLAMF7", "CD38", # "IL2RG",
            "HAVCR1", "CDH6", "SOX9",            
            "HNF4A", "CUBN", "LRP2", #            
            "CASR", "UMOD", "SIM2",  # "CDH3"
            "NOTCH3", "CSPG4", "MCAM"
)

names(markers) <- c(
                        "C-Duct", "C-Duct", "C-Duct",
        "Fib", "Fib", "Fib",
                    "Glo", "Glo", "Glo",
    
                    "Immune", "Immune", "Immune",
                    "I-PT", "I-PT", "I-PT",
    
                    "PT", "PT", "PT",
    
                    "TAL", "TAL", "TAL",
                    "Vessel", "Vessel", "Vessel"
                    )


## axis <- ggh4x::guide_axis_truncated(
##   trunc_lower = unit(0, "npc"),
##   trunc_upper = unit(1, "cm")
## )

print('------')
p2 = DotPlot_scCustom(subset(imm.niche, subset=niche_label != 'Skeletal Muscle'), features=unname(markers), group.by='niche_label') + 
    theme(axis.text.x=element_text(face="bold", angle=90, hjust = 1, vjust=1), axis.text.y=element_text(face="bold")) + xlab("") + ylab("") + get_theme(angle=90, size=12) + theme(legend.box="vertical", legend.margin=margin(), plot.margin = unit(c(0,0,0,0), "cm"), ) + scale_size_continuous(range = c(0.1, 2)) +ggtitle("")

fig.size(8, 13)
fig2a = p2

markers = c("CD8A", "CD2", "ZAP70", "CTLA4", "FOXP3", "TIGIT", "CD3E", "MKI67", "TUBB", "XBP1", "CD38", "FCRL5", "MS4A1", "CD19", "CD79A", "CD163", "MRC1", "F13A1", "CD14", "CIITA", "FCN1", "CXCL9", "CXCL10", "MMP9", "WDFY4", "CLEC9A", "IRF8", "GZMB", "LILRA4", "IL3RA", "KIT", "HDC",  "PECAM1", "PLVAP", "PDGFRB", "PODXL", "FGF1", "NES", "SHANK3",  "ITGB3", "CFH", "BMP7", "TNC", "AEBP1", "LRP2", "CUBN", "PAH", "ITGB6", "CA12", "MUC1",  "EPCAM", "PROM1", "PAX8", "HSD11B2", "KCNJ10", "SERPINA5", "UMOD", "CASR", "SCNN1A", "SLC4A1", "DMRT2", "CLNK", "SCNN1G", "GATA3", "PFKFB3", "PKHD1", "CALB1", "KCNJ1")

library(glue)
niche.merge@meta.data$condition_num = as.numeric(factor(str_trim(niche.merge@meta.data$condition), levels=c("Control", "Case")))

niche.cna <- association.Seurat(
    seurat_object = niche.merge,
    test_var = 'condition_num',
    samplem_key = 'sample_id',
    graph_use = 'RNA_snn',
    verbose = TRUE,
    batches = NULL, ## no batch variables to include
    #covs = c("age", "sex", "ICPi") ## no covariates to include
)



p3=FeaturePlot_scCustom(niche.cna, features = c('cna_ncorrs_fdr10'), raster=T, raster.dpi=c(150, 150))[[1]] + #
    scale_color_gradient2(high = "#de2d26", mid = "white", low = "#2c7fb8", midpoint = 0,  guide = guide_colorbar(direction = "vertical"))+
    labs(title = 'ICI-AIN-associated niches', subtitle = 'Filtered for FDR<0.10', color = 'Correlation')+ ggplot2::theme(legend.position = "right")
fig2e = p3

sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(lennard_label=lennard.subtype@meta.data$lennard_label)
sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(tile_label = imm.niche@meta.data[match(sc.niche$obj@meta.data$tile_id, rownames(imm.niche@meta.data)), 'niche_label'])
sc.niche.heatmap<-sc.niche$obj@meta.data %>% count(tile_label, lennard_label) %>% pivot_wider(names_from=lennard_label, values_from=n, values_fill=0) %>% filter(!is.na(tile_label))
sc.niche.heatmap = as.data.frame(sc.niche.heatmap)
rownames(sc.niche.heatmap) = sc.niche.heatmap$tile_label
sc.niche.heatmap = sc.niche.heatmap[,-1]
sc.niche.heatmap <- scale(sc.niche.heatmap)

sc.niche.heatmap <- sc.niche.heatmap[rownames(sc.niche.heatmap)!="Skeletal Muscle",]

#names(niche_cols) <- rownames(sc.niche.heatmap)

fig.size(5,9)
## Use niche colors here
ha1 = rowAnnotation(samples = rownames(sc.niche.heatmap),
                    col=list(samples=niche_cols),
                    annotation_name_gp = gpar(fontsize = 12),    
                    show_legend = FALSE,
                    annotation_legend_param = list(
                    samples = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5), direction = "horizontal")))

## Use heatmap colors

ht <- Heatmap(
  sc.niche.heatmap,
  name = "Cell Composition Ratio",
  col = heatmap_col_fun,
  width=6.5,
  height=1,
  column_names_gp = gpar(fontsize = 13),
  heatmap_legend_param = list(
    legend_height = unit(5, "cm"),
    title_gp = gpar(fontsize = 15), # Increase legend title font size
    labels_gp = gpar(fontsize = 12),    
    direction = "horizontal"
  ), 
        cluster_columns = TRUE, 
        show_column_dend = FALSE, 
    left_annotation = ha1)

fig1f = ht


library(patchwork)
ht_plot <- wrap_elements(full=grid.grabExpr(draw(fig1f,
                              merge_legend = TRUE, heatmap_legend_side = "top",
                              annotation_legend_list = NULL)))

#fig2b1 = DimPlot_scCustom(subset(niche.merge, subset = condition == 'Case'), group.by='niche_label') + scale_color_manual(values=niche_cols) + ggtitle("ICI-AIN")
#fig2b2 = DimPlot_scCustom(subset(niche.merge, subset = condition == 'Control'), group.by='niche_label') + scale_color_manual(values=niche_cols) + ggtitle("ICI-ATN")
fig2b = DimPlot_scCustom(niche.merge, group.by='niche_label', raster.dpi=c(150, 150), raster=T) + scale_color_manual(values=niche_cols)

library(patchwork)
#pa_d <- p0 + p1 + fig2b + ht_plot + plot_layout(ncol=4, guides='collect') & theme(legend.position='bottom')

library(sf)
imm.niche@meta.data <- imm.niche@meta.data %>% mutate(patient_id=gsub("__2.*", "", sample_id))

case1 <- "BS21-N65682A2"
#case2 <- "BS23_52206A2"
cont1 <- "BS22_12012A1"
cont2 <- "BS2_61615A1"

library(tidyverse)
p20 = ggplot() +
    rasterise(geom_sf(data = as_tibble(imm.niche@meta.data, niche_label=imm.niche$niche_label_fine)%>% filter(patient_id==case1),
            aes(geometry = shape, fill = niche_label), color = NA), dpi=100) + theme_bw(base_size=12) + scale_fill_manual(values=niche_cols)+ theme(legend.position = 'none')+
    coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + ggtitle("ICI-AIN case patient 1 niche tiles")


fig2c = p20


fig.size(15, 20)
((fig2b + plot_layout(ncol=2, guides='collect') & theme(legend.position='none')) | ht_plot) / (fig2b |  fig2a) + 
 plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E")), 
 ) &  theme(plot.tag = element_text(face = "bold", size = 28))

ggsave('Fig2_part1.pdf', width=20, height=13)

library(tidyverse)
p10 = ggplot() +
    rasterise(geom_sf(data = as_tibble(imm.niche@meta.data, niche_label=imm.niche$niche_label_fine)%>% filter(patient_id==cont1),
              aes(geometry = shape, fill = niche_label), color = NA), dpi=100) + theme_bw(base_size=12) + scale_fill_manual(values=niche_cols) + theme(legend.position = 'none')+
    coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + ggtitle("ICI-ATN control patient 1 niche tiles")

fig2d = p10

feature = 'PODXL'
feature2 = 'podxl'
df1 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS22')
max_ox = max(df1$podxl, df2$podxl)

p51 = ggplot() + 
    rasterise(geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA), dpi=100) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + ggtitle("PODXL expression")

p52 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + ggtitle("PODXL expression")  

feature = 'CD79B'
feature2 = 'cd79b'
df1 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS22')
max_ox = max(df1$cd79b, df2$cd79b)

p51_2 = ggplot() + 
    rasterise(geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA), dpi=100) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + ggtitle("CD79B expression")

p52_2 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + ggtitle("CD79B expression")  

feature = 'UMOD'
feature2 = 'umod'
df1 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS22')
max_ox = max(df1$umod, df2$umod)

p51_3 = ggplot() + 
    rasterise(geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA), dpi=100) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + ggtitle("UMOD expression")

p52_3 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + ggtitle("UMOD expression")  

feature = 'CUBN'
feature2 = 'cubn'
df1 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id==case1)
df2 = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>%
                janitor::clean_names() %>% mutate(patient_id=`orig_ident`) %>% 
                filter(patient_id=='BS22')
max_ox = max(df1$cubn, df2$cubn)

p51_4 = ggplot() + 
    rasterise(geom_sf(data = df1,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA), dpi=100) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(1000, 4000), ylim=c(1000, 2500), expand = FALSE) + ggtitle("CUBN expression")

p52_4 = ggplot() + 
    geom_sf(data = df2,
            aes(geometry = shape, fill = !!sym(feature2)), color = NA) + theme_bw(base_size=12) +
    scale_fill_gradient(limits = c(0,max_ox), low = 'white', high = '#832424')+ theme(legend.position='right')+
    coord_sf(xlim=c(2000, 5000), ylim=c(500, 2000), expand = FALSE) + ggtitle("CUBN expression")  

fig2c + fig2d + p51 + p52 + plot_layout(ncol=2)
ggsave('Fig2_part2.pdf', width=20, height=20)

(fig2b | ht_plot) / (fig2c | fig2d) / (((p51 + p51_2) / (p51_3 + p51_4)) | ((p52+p52_2)/(p52_3+p52_4))) / (fig2e |  fig2a) + 
 plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E", 'F', 'G', 'H', 'I', 'J', 'K', 'L','M','N')), 
 ) &  theme(plot.tag = element_text(face = "bold", size = 28))

ggsave('Fig2_v3.pdf', width=18.5, height=21.5)

