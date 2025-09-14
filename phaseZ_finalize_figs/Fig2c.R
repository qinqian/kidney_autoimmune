suppressPackageStartupMessages({
    library(tessera)
    ## Downstream analysis in Seurat V5
    library(Seurat)
    library(scCustomize)
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
    library(cowplot)
    library(rcna)
    library(readr)
    library(stringr)
    library(spatula)
})

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

source("https://raw.githubusercontent.com/kris-nader/sp-type/main/sp-type.R")

tile_obj <- readRDS("250711_niches.rds")

tile_obj@meta.data = tile_obj@meta.data %>% mutate(niche_label_input = ifelse(niche_label=='Fibrosis & Interstitium', niche_label, niche_label_fine))

tile_df = tile_obj@meta.data

library(spatula)
library(purrr)
source("../bin/cell_interaction.R")

niche_interaction = list()
set.seed(999)

for (cont in unique(tile_df %>% filter(tile_df$condition=='Control') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                                  map.c.plain %>% select(X, Y),
                                  niches_ind, compartments=NULL, nsteps=3)
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'control') %>% mutate(sample_id=cont)
}


for (cont in unique(tile_df %>% filter(tile_df$condition=='Case') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                    map.c.plain %>% select(X, Y),
                    niches_ind, compartments=NULL, nsteps=3)    
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'case') %>% mutate(sample_id=cont)
}


niche_interaction.df = niche_interaction %>% rbindlist() %>% filter(grepl('Fib', index_type)) 
write_tsv(niche_interaction.df, "niche_interaction_step3.tsv")

    ## filter((cond=='case' & fdr <= 0.05)|(cond=='control')) %>% arrange(fdr)


pdf("test2.pdf", width=20, height=12)
p1 = niche_interaction.df %>% 
    ggplot(aes(y=freq, colour=cond, x=paste0(index_type, ":", type))) + 
    geom_point(position = position_dodge(width = 0.75)) +
    xlab("Fibrogenic niche interaction with other niches") + 
    ggpubr::stat_compare_means(
        ## comparisons = list(c("case", "control")),
        aes(group = cond),
        paired = F, 
        method = "wilcox.test",
        method.args = list(alternative = "greater"),        
        label = "p.signif", 
        label.y.npc = "top"
        )+
    ggpubr::theme_pubclean(base_size=9) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, size=15))
dev.off()

for (cont in unique(tile_df %>% filter(tile_df$condition=='Control') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                                  map.c.plain %>% select(X, Y),
                                  niches_ind, compartments=NULL, nsteps=1)
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'control') %>% mutate(sample_id=cont)
}


for (cont in unique(tile_df %>% filter(tile_df$condition=='Case') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                    map.c.plain %>% select(X, Y),
                    niches_ind, compartments=NULL, nsteps=1)    
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'case') %>% mutate(sample_id=cont)
}


niche_interaction.df = niche_interaction %>% rbindlist() %>% filter(grepl('Fib', index_type)) 
write_tsv(niche_interaction.df, "niche_interaction_step1.tsv")

    ## filter((cond=='case' & fdr <= 0.05)|(cond=='control')) %>% arrange(fdr)


pdf("test2_step1.pdf", width=20, height=12)
p2 = niche_interaction.df %>% 
    ggplot(aes(y=freq, colour=cond, x=paste0(index_type, ":", type))) + 
    geom_point(position = position_dodge(width = 0.75)) +
    xlab("Fibrogenic niche interaction with other niches") + 
    ggpubr::stat_compare_means(
        ## comparisons = list(c("case", "control")),
        aes(group = cond),
        paired = F, 
        method = "wilcox.test",
        method.args = list(alternative = "greater"),        
        label = "p.signif", 
        label.y.npc = "top"
        )+
    ggpubr::theme_pubclean(base_size=9) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, size=15))
dev.off()

for (cont in unique(tile_df %>% filter(tile_df$condition=='Control') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                                  map.c.plain %>% select(X, Y),
                                  niches_ind, compartments=NULL, nsteps=2)
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'control') %>% mutate(sample_id=cont)
}


for (cont in unique(tile_df %>% filter(tile_df$condition=='Case') %>% select(sample_id) %>% pull())) {
    map.c.plain <- tile_df %>% filter(sample_id==cont) %>% filter(niche_label_input!='Skeletal Muscle')
    niches = unique(map.c.plain$niche_label_input)
    niches_ind = as.character(map.c.plain$niche_label_input)
    coloc.mapc <- coloc_all_types(niches,
                    map.c.plain %>% select(X, Y),
                    niches_ind, compartments=NULL, nsteps=2)    
    niche_interaction[[cont]] = coloc.mapc %>% mutate(cond = 'case') %>% mutate(sample_id=cont)
}


niche_interaction.df = niche_interaction %>% rbindlist() %>% filter(grepl('Fib', index_type)) 
write_tsv(niche_interaction.df, "niche_interaction_step2.tsv")

    ## filter((cond=='case' & fdr <= 0.05)|(cond=='control')) %>% arrange(fdr)


pdf("test2_step2.pdf", width=20, height=36)
p3 = niche_interaction.df %>% 
    ggplot(aes(y=freq, colour=cond, x=paste0(index_type, ":", type))) + 
    geom_point(position = position_dodge(width = 0.75)) +
    xlab("Fibrogenic niche interaction with other niches") + 
    ggpubr::stat_compare_means(
        ## comparisons = list(c("case", "control")),
        aes(group = cond),
        paired = F, 
        method = "wilcox.test",
        method.args = list(alternative = "greater"),        
        label = "p.signif", 
        label.y.npc = "top"
        )+
    ggpubr::theme_pubclean(base_size=9) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, size=15))
print(p2/p3/p1)
dev.off()
