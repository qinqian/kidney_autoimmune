library(Seurat)
library(dplyr)
library(sf)
library(tmap)
library(tidyverse)
library(patchwork)


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_text(size=size), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


imm.niche <- readRDS("250711_niches.rds")

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


pdf("Fig2b.pdf", height=8, width=8)
tmap_options(component.autoscale = T)
## tmap_options(tmap.limits = c(facets.view=16, facets.plot=128))
map1<-tm_shape(tile_sf %>% filter(patient_id %in% c(case1, case2))) +
  tm_polygons(fill='niche_label', col = NULL) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)
map2<-tm_shape(tile_sf %>% filter(patient_id %in% c(cont1, cont2))) +
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
map1m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)
map2m1<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
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
map1m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
map2m0<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
  tm_polygons(fill=feature, col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
feature = 'PODXL'
map1m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(case1, case2))) +
  tm_polygons(fill='PODXL', col = NULL, fill.scale = tm_scale_continuous()) +
  tm_facets(c("patient_id"), ncol=1)+
  tm_layout(legend.show = FALSE)    
map2m2<-tm_shape(cbind(tile_sf, FetchData(imm.niche, feature)) %>% filter(patient_id %in% c(cont1, cont2))) +
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
## tmap_arrange(map1, map2, map1m, map2m, map1m2, map2m2, nrow=4, ncol=4)
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
