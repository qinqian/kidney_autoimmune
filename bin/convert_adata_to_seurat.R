library(Seurat)
library(sceasy)
library(glue)
library(ggplot2)
library(argparse)
library(stringr)

parser <- ArgumentParser(prog="convert_seurat_to_adata.R", description="convert seurat data to adata h5ad")
parser$add_argument("data", help="rds")
parser$add_argument("--output", help="output prefix")
args = parser$parse_args()

use_condaenv("scconverter")

sceasy::convertFormat(args$data, to="seurat", from="anndata", outFile=args$output)
print(packageVersion("SeuratObject"))
print(packageVersion("Seurat"))

print(sessionInfo())


