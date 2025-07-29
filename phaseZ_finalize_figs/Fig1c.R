library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
suppressPackageStartupMessages({
    library(tessera)
    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(sf)
    ## Plotting functions 
    ## Not imported by Tessera
    library(ggplot2)
    library(ggthemes)
    library(viridis)
    library(patchwork)
    library(harmony)
    library(Seurat)
    library(dplyr)
})


fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}


get_theme <- function(size=12, angle=45) {
    defined_theme = theme_bw(base_size=size) + theme(legend.title=element_blank(), strip.text=element_text(size=size), legend.text=element_text(size=size), axis.title.x=element_text(size=size), axis.title.y=element_text(size=size), axis.text.y=element_text(size=size), axis.text.x=element_text(size=size, angle=angle, hjust = 1, vjust=1.05), legend.position="bottom", legend.box = "horizontal") 
    defined_theme
}


library(sccomp)

obj.merge = readRDS("250721_cells_annotated_lennard.rds")


meta = read.csv("~/shruti_meta_clean (3).csv")
input_meta = meta[,c('slide_id', 'age', 'sex', 'case_ctrl', 'ICPi',  'malignancy', 'eGFR_base')] %>% arrange(case_ctrl)
#tile_obj$condition = str_trim(meta[match(str_extract(as.character(obj.merge$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

obj.merge@meta.data$cell_label = gsub(" Cell", "", obj.merge@meta.data$lennard_label)
obj.merge@meta.data = obj.merge@meta.data %>% mutate(cell_label = ifelse(cell_label=='Immune', 'Immune (LowQ)', cell_label))

obj.merge@meta.data  = obj.merge@meta.data %>% mutate(sample_id=str_extract(sample, "__(BS\\d*[_-].*)__2024", group=1)) 
obj.merge@meta.data  = obj.merge@meta.data %>% mutate(case_ctrl=input_meta[match(sample_id, input_meta$slide_id), 'case_ctrl'])
obj.merge@meta.data$case_ctrl_num = as.numeric(factor(obj.merge@meta.data$case_ctrl))

set.seed(999)

obj.merge@meta.data$case_ctrl = factor(str_trim(obj.merge@meta.data$case_ctrl), levels=c("Control", "Case"))

res = obj.merge |>
    sccomp_estimate(
      formula_composition = ~ case_ctrl,
      sample = "sample_id", cell_group = "cell_label",
      cores = 2, verbose=FALSE
    )

res = res |>  sccomp_test()

pdf("Fig1e.pdf", height=5, width=4.6)
p <- (res |>  plot_1D_intervals(show_fdr_message=F)) + theme_bw(base_size=8) + theme(legend.position="bottom") + ylab("Cell type")
print(p)
dev.off()

meta = read.csv("~/shruti_meta_clean (3).csv")
input_meta = meta[,c('slide_id', 'age', 'sex', 'case_ctrl', 'ICPi',  'malignancy', 'eGFR_base')] %>% arrange(case_ctrl)
#tile_obj$condition = str_trim(meta[match(str_extract(as.character(obj.merge$sample_id), "(.*)__", group=1), meta$slide_id), 'case_ctrl'])

test = obj.merge@meta.data %>% filter(tech=='xenium') %>% janitor::clean_names()
test$condition = str_trim(meta[match(str_extract(as.character(subset(obj.merge, subset=tech=='xenium')$orig.ident), "__(BS.*A[1,2])__(2.+)", group=1), meta$slide_id), 'case_ctrl'])

sample_ids = test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
    mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(sample_ids) %>% pull()
sample_mat = as.matrix(test %>% filter(tech=='xenium') %>% count(sample_ids, lennard_label) %>% group_by(sample_ids) %>% 
    mutate(ratio = n/sum(n)) %>% select(-n) %>% pivot_wider(names_from=lennard_label, values_from=ratio, values_fill = 0) %>% ungroup() %>% select(-sample_ids))
rownames(sample_mat) = str_extract(as.character(sample_ids), "__(BS.*A[1,2])__(2.+)", group=1)

my_sample_col <- data.frame(sample = str_trim(meta$case_ctrl))
row.names(my_sample_col) <- meta$slide_id

my_sample_col <- my_sample_col %>% arrange(sample)

## my_sample_col = my_sample_col[rownames(sample_mat), ,drop=F]

my_sample_col = my_sample_col %>%
    mutate(group=case_when(
               sample=="Case" ~ "ICI-AIN",
               sample=="Control" ~ "ICI-ATN")) %>%
    mutate(group_simple = case_when(
               sample=="Case" ~ "AIN",
               sample=="Control" ~ "ATN")) %>% 
    mutate(patient=rep(seq(1, 4), 2))  %>%
    unite("patient_id", c(group_simple, patient), remove=F)

write_tsv(my_sample_col, "my_sample_col.tsv")


sample_mat = sample_mat[rownames(my_sample_col),]
rownames(sample_mat) = my_sample_col$patient_id

sample_mat <- scale(sample_mat)

my_sample_col = my_sample_col %>% rownames_to_column()
my_sample_col = my_sample_col %>% column_to_rownames(var='patient_id')

# # Define color function
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht_opt("heatmap_row_names_gp" = gpar(fontsize = 5))
ha1 = rowAnnotation(samples = my_sample_col$group, 
                    col=list(samples=c('ICI-AIN'='gray', 'ICI-ATN'='black')),
#                    annotation_label_location = "top",  # move label to top
                    annotation_name_gp = gpar(fontsize = 5),  
                    annotation_legend_param = list(
                    samples = list(title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5), direction = "horizontal")))


colnames(sample_mat)[colnames(sample_mat)=='Immune Cell'] = 'Immune (LowQ)'

cell_label = res %>% filter(!grepl("Intercep", parameter)) %>% arrange(desc(c_effect)) %>% select(cell_label) %>% pull()

colnames(sample_mat) = gsub(" Cell", "", colnames(sample_mat))

sample_mat = sample_mat[, cell_label]

# Desired cell size
cell_width <- unit(5, "mm")   # width per column
cell_height <- unit(5, "mm")   # height per row

# Main heatmap
ht <- Heatmap(
  sample_mat,
  name = "Expression",
  col = col_fun,
  ## width = cell_width * ncol(sample_mat),
  ## height = cell_height * nrow(sample_mat),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    direction = "horizontal"
  ),cluster_columns = FALSE,
  column_names_side = "top",
  right_annotation = ha1)

library(grid)

pdf("Fig1c.pdf", height=8, width=4.5)
## draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
##      annotation_legend_side = "bottom")
ht_grob <- grid.grabExpr(draw(ht, newpage = FALSE))
ht_plot <- wrap_elements(full = ht_grob)
combined <- ht_plot + (p+ggtitle("ICI-AIN vs ICI-ATN")) + plot_layout(ncol = 1, height=c(1.0, 2.3))
print(combined)
dev.off()
