library(scCustomize)
library(glue)
library(Seurat)
library(ggplot2)
library(argparse)
library(patchwork)
library(tidyplots)
library(tidyverse)
library(ggplot2)
library(tessera)
library(plot1cell)


parser <- ArgumentParser(prog="subset_cells", description="")
parser$add_argument("data", help="input xenium directory or rds")
parser$add_argument("--cell", help="cell type")
parser$add_argument("--meta", help="output prefix")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()

match_celltype <- function(x) {    
    x@meta.data$cell_type = meta$coarse_ids[
                                    match(as.character(singlet@meta.data$cell_ident),
                                          rownames(meta))]
    x
}

meta = readRDS(args$meta)
print(head(meta))
singlet = readRDS(args$data)

singlet@meta.data = singlet@meta.data %>% mutate(sample_id=gsub("_orig_seg.h5ad", "", sample)) %>% 
   unite("cell_ident", c(sample_id, cell_id),  na.rm = TRUE, sep="_", remove = FALSE)

singlet = match_celltype(singlet)
print(head(singlet@meta.data, 3))
print(head(setdiff(singlet@meta.data$cell_ident, rownames(meta))))
print(length(setdiff(singlet@meta.data$cell_ident, rownames(meta))))

print(table(meta$coarse_ids))
print(table(singlet@meta.data$cell_type))

print(head(singlet@meta.data$cell_ident))
print(head(rownames(meta)))

saveRDS(subset(singlet, subset= coarse_ids == args$cell), file=glue("{args$output}.rds"))
