library(Seurat)
library(dplyr)
library(sf)
library(tmap)
library(ggrastr)
library(tidyverse)
library(patchwork)


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


imm.niche <- readRDS("250711_niches.rds")

lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
orig.baysor <- readRDS("../phaseF_newpipeline/sopa_seg/comb_h5ad/kidney_orig_seg_merged.rds")

orig.baysor@meta.data <- orig.baysor@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)
lennard.subtype@meta.data <- lennard.subtype@meta.data%>%unite("uniq_id", c(sample, cell_id), remove=F)


xy <- Embeddings(orig.baysor, 'spatial')[match(lennard.subtype@meta.data$uniq_id, orig.baysor@meta.data$uniq_id),]
intensity <- Embeddings(orig.baysor, 'intensities')[match(lennard.subtype@meta.data$uniq_id, orig.baysor@meta.data$uniq_id),]

## comb_h5ad/kidney_orig_seg_merged.rds

patient <- read_tsv("my_sample_col.tsv")
imm.niche@meta.data <- imm.niche@meta.data %>% mutate(patient_id=gsub("__2.*", "", sample_id))


tile_sf <- imm.niche@meta.data %>% st_as_sf()
print(tile_sf)

print(st_crs(tile_sf))
st_set_crs(tile_sf, 4326)
print(st_crs(tile_sf))

common_xlim <- c(0, 5000)
common_ylim <- c(0, 5000)

case1 <- "BS21-N65682A2"
case2 <- "BS23_52206A2"
cont1 <- "BS22_12012A1"
cont2 <- "BS2_61615A1"


pdf("Fig2b.pdf", height=5, width=8)
tmap_options(component.autoscale = T)
## tmap_options(tmap.limits = c(facets.view=16, facets.plot=128))
#map1<-tm_shape(tile_sf %>% filter(patient_id %in% c(case1, case2))) +
map1<-tm_shape(tile_sf %>% filter(patient_id %in% c(case1))) +
  tm_polygons(fill='niche_label', col = NULL) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)
#map2<-tm_shape(tile_sf %>% filter(patient_id %in% c(cont1, cont2))) +
map2<-tm_shape(tile_sf %>% filter(patient_id %in% c(cont2))) +
  tm_polygons(fill='niche_label', col = NULL) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
#map1=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(case1, case2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
#map2=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(cont1, cont2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
## tmap_arrange(map1, map2)
feature = 'CXCL9'
#map1m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
map1m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)
#map2m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
map2m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont2))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
#map1=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(case1, case2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
#map2=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(cont1, cont2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
## tmap_arrange(map1, map2)
feature = 'COL5A1'
#map1m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
map1m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
#map2m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
map2m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont2))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
feature = 'PODXL'
#map1m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
map1m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1))) +
  tm_polygons(fill='PODXL', col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
#map2m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
map2m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont2))) +
  tm_polygons(fill='PODXL', col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
#map1=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(case1, case2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
#map2=tm_shape(tile_sf%>% st_as_sf(), # %>% filter(patient_id %in% c(cont1, cont2)),
#              tm_polygons(fill = "niche_label", col = NULL) + 
#              #tm_facets(c("sample_id"), ncol=2) + 
#              tm_layout(legend.title.size = 8, legend.text.size = 6.5,
#                        legend.frame = F, panel.label.size=6.5, 
#                        panel.label.height = 2,
#             )
#             )
# Immune | Glomerui | Fib 
tmap_arrange(map1, map1m1, map1m2, map1m0, map2, map2m1, map2m2, map2m0, nrow=2, ncol=4)
## bbox_union <- st_bbox(tile_sf)
## common_xlim <- bbox_union[c("xmin", "xmax")]
## common_xlim[0] = common_xlim[0] - 500
## common_xlim[1] = common_xlim[1] + 500
## common_ylim <- bbox_union[c("ymin", "ymax")]
## common_ylim[0] = common_ylim[0] - 180
## common_ylim[1] = common_ylim[1] + 200
## p1 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424')+theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p2 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case2),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424') +theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p3 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont1),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424')+theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p4 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424') +theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## print((p1/p2/p3/p4)+plot_layout(height=c(1,1,1,1), ncol=1))
## feature = 'PODXL'
## p1 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424')+theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p2 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case2),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424') +theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p3 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont1),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424')+theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## p4 = ggplot() + 
##     geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
##             aes(geometry = shape, fill = !!sym(feature)), color = "lightgray", size = 0.01) + 
##     scale_fill_gradient(low = 'white', high = '#832424') +theme_void()+coord_sf(xlim = common_xlim, ylim = common_ylim, expand = T) #+ theme(aspect.ratio = 1)
## print((p1/p2/p3/p4)+plot_layout(height=c(1,1,1,1), ncol=1))
dev.off()

library(ggsci)
library(circlize)
library(scico)

sc.niche <- readRDS("../phaseF_newpipeline/sopa_seg/sopa_baysor_tessera.rds")
lennard.subtype <- readRDS("250721_cells_annotated_lennard.rds")
imm.niche <- readRDS("250711_niches.rds")

sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(lennard_label=lennard.subtype@meta.data$lennard_label)
sc.niche$obj@meta.data = sc.niche$obj@meta.data %>% mutate(tile_label = imm.niche@meta.data[match(sc.niche$obj@meta.data$tile_id, rownames(imm.niche@meta.data)), 'niche_label'])
sc.niche.heatmap<-sc.niche$obj@meta.data %>% count(tile_label, lennard_label) %>% pivot_wider(names_from=lennard_label, values_from=n, values_fill=0) %>% filter(!is.na(tile_label))
sc.niche.heatmap = as.data.frame(sc.niche.heatmap)
rownames(sc.niche.heatmap) = sc.niche.heatmap$tile_label
sc.niche.heatmap = sc.niche.heatmap[,-1]
sc.niche.heatmap <- scale(sc.niche.heatmap)
col_fun <- colorRamp2(c(-2, 0, 2), scico(3, palette = "vik"))  # "vik" is diverging

sc.niche.heatmap <- sc.niche.heatmap[rownames(sc.niche.heatmap)!="Skeletal Muscle",]

niche_cols = pal_npg("nrc")(9)

names(niche_cols) <- rownames(sc.niche.heatmap)

library(scCustomize)
cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
names(cols) <- unique(lennard.subtype@meta.data$lennard_label)

patient <- read_tsv("my_sample_col.tsv")
imm.niche@meta.data <- imm.niche@meta.data %>% mutate(patient_id=gsub("__2.*", "", sample_id))
lennard.subtype@meta.data = lennard.subtype@meta.data %>% mutate(patient_id=gsub("output-XETG00150__0018462__|output-XETG00392__0045655__", "", gsub("__2.*", "", sample_ids)))

pdf("Rplots.pdf", width=24, height=6)
feature = 'PODXL'
p0 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, imm.niche$niche_label)%>% filter(patient_id==case1),
            aes(geometry = shape, fill = niche_label), color=NA) + theme_bw(base_size=12) +
    guides(fill = guide_legend(nrow = 3)) +
    coord_sf(expand = FALSE) + NULL+
    scale_fill_d3('category20c') + theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 6),      # legend labels
                                         legend.title = element_text(size = 7))+scale_fill_manual(values=niche_cols)
p1 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==case1),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(1700, 2100), ylim=c(900, 1200), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')
feature = 'CD38'
p2 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==case1),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==case1),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(3000, 3500), ylim=c(1500, 1800), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')
p01 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, imm.niche$niche_label)%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = niche_label), color=NA) + theme_bw(base_size=12) +
    guides(fill = guide_legend(nrow = 3)) +
    coord_sf(expand = FALSE) + NULL+
    scale_fill_d3('category20c') + theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 6),      # legend labels
                                         legend.title = element_text(size = 7))+scale_fill_manual(values=niche_cols)
p11 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==cont2),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='none',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(1600, 2300), ylim=c(1200, 1600), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')
feature = 'CD38'
p21 = ggplot() + 
    geom_sf(data = cbind(imm.niche@meta.data, FetchData(imm.niche, feature))%>% filter(patient_id==cont2),
            aes(geometry = shape, fill = !!sym(feature)), color = 'lightgray') + theme_bw(base_size=12) +
    geom_point_rast(data = lennard.subtype@meta.data %>% mutate(x=xy[,1], y=xy[,2]) %>% filter(patient_id==cont2),
               aes(x=x, y=y, color=lennard_label), size=0.25, alpha=0.9)+
    scale_fill_gradient(low = 'white', high = '#832424')+
    guides(fill = guide_legend(nrow = 5),
           color = guide_legend(override.aes = list(size = 2)))+ theme(aspect.ratio = 1, legend.position='bottom',
                                         legend.text = element_text(size = 3),      # legend labels
                                         legend.title = element_text(size = 4))+scale_color_manual(values=cols)+
    coord_sf(xlim=c(6600, 7200), ylim=c(600, 900), expand = FALSE) + NULL + theme(aspect.ratio = 1, legend.position='none')
print(p0+p1+p2+p01+p11+p21+plot_layout(ncol=6, widths=c(3,3,3,3,3,3)))
dev.off()
