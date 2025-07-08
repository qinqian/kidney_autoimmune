library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyplots)
library(stringr)
library(sf)

setwd("/Users/qian/Documents/active_projects/xenium/others_check_rds")
kidn_seu_obj = readRDS("../data/miles/shapes_seur_kidney_v5.rds")
metadata = read_csv("../data/kidney/shruti_meta_clean.csv")
kidn_seu_obj$seur@meta.data = kidn_seu_obj$seur@meta.data %>% mutate(slide_id=str_extract(sample_id, paste(metadata %>% select(slide_id) %>% pull(), collapse="|")))

kidn_seu_obj$seur@meta.data = kidn_seu_obj$seur@meta.data %>% mutate(group=kidn_seu_obj$seur@meta.data %>% left_join(metadata %>% select(slide_id, case_ctrl), by="slide_id") %>% select(case_ctrl))

## kidn_seu_obj$seur@meta.data = kidn_seu_obj$seur@meta.data %>% unite("group_sample", c(slide_id, group))

print(head(kidn_seu_obj$seur@meta.data))

pdf("kidney_groupby_condition.pdf", width=8, height=12)
print(FeaturePlot(kidn_seu_obj$seur, features=c("CXCL9", "CXCL10", "CXCL11", "CXCL12"), split.by="group"))
dev.off()


kidn_seu_obj$seur$sample_group <- paste(kidn_seu_obj$seur$slide_id, kidn_seu_obj$seur$group, sep = "_")
kidn_seu_obj$seur$sample_group <- factor(kidn_seu_obj$seur$sample_group, levels=unique(kidn_seu_obj$seur$sample_group)[order(str_extract(unique(kidn_seu_obj$seur$sample_group), "Control|Case"))])

pdf("kidney_groupby_condition2.pdf", width=25, height=12)
print(FeaturePlot(kidn_seu_obj$seur, features=c("CXCL9", "CXCL10", "CXCL11", "CXCL12"), split.by="slide_id"))
print(FeaturePlot(kidn_seu_obj$seur, features=c("CXCL9", "CXCL10", "CXCL11", "CXCL12"), split.by="sample_group"))
dev.off()

