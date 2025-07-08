library(Seurat)
library(Matrix)
library(dplyr)
library(glue)
library(ggplot2)
library(scCustomize)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

parser <- ArgumentParser(prog="subset_cell_type.r", description="")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")

args = parser$parse_args()
output = args$output
print(args)

xen.l = list()
xen.c = list()

cell.ids = basename(args$data) %>% 
   str_replace("_orig_seg\\.rds", "")

for (idx in seq_along(args$data)) {
    r = args$data[idx]
    print(r)
    xen = readRDS(r)
    print(xen)
    if (idx == 1) {
        genes = rownames(xen)
    }
    xen.l[[r]] = xen@meta.data
    xen.c[[r]] = xen@assays[['RNA']]$counts
    rownames(xen.l[[r]]) = paste0(cell.ids[idx], "__", rownames(xen.l[[r]]))
    colnames(xen.c[[r]]) = rownames(xen.l[[r]])
    print(object.size(xen.l))
    print(object.size(xen.c[[r]]))
    print(dim(xen.c[[r]]))
    print(sum(rownames(xen.l[[r]]) == colnames(xen.c[[r]])))
    print(sum(rownames(xen.c[[r]]) == genes))
}


meta.data = bind_rows(xen.l)
print(dim(meta.data))
counts = do.call(cbind, xen.c)
print(dim(counts))
print(counts[1:5, 1:5])
print(colnames(counts)[1:5])
merged.xen = CreateSeuratObject(counts=counts, meta.data=meta.data, project='RA', assay='RNA')
print(object.size(merged.xen))
saveRDS(merged.xen, glue("comb/{args$output}_orig_seg.rds"))

