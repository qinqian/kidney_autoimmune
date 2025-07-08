library(scCustomize)
library(glue)
library(Seurat)
library(ggplot2)
library(argparse)
library(patchwork)
library(tidyplots)
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
    x@meta.data$cell_type = meta$coarse_ids[match(gsub("_orig_seg_adddbl.rds_", "",
                                                              gsub("(pre_lowQC|pre|post|post2|pre2|post_lowQC)_orig_seg_adddbl.rds__", "", rownames(x@meta.data))),
                                    gsub("_[12]_([a-z]*)", "_\\1", rownames(meta)))]
    x
}

meta = readRDS(args$meta)
print(head(meta))
singlet = readRDS(args$data)
print(head(singlet@meta.data, 3))
print(head(setdiff(rownames(singlet@meta.data), rownames(meta))))
print(length(setdiff(rownames(singlet@meta.data), rownames(meta))))

singlet = match_celltype(singlet)

print(table(meta$coarse_ids))
print(table(singlet@meta.data$cell_type))
saveRDS(subset(singlet, subset= coarse_ids == args$cell), file=glue("{args$output}.rds"))
