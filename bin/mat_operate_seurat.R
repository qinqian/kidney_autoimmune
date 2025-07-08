library(Seurat)
library(Matrix)
library(dplyr)
library(BPCells)
library(glue)
library(ggplot2)
library(scCustomize)
library(argparse)
library(stringr)
library(tidyverse)
library(patchwork)

# run uwot umap
Run_uwot_umap <- function(SeuratObj, min_dist = 0.3, spread = 0.8){
    HU <- uwot::umap(SeuratObj@reductions$harmony@cell.embeddings, min_dist = min_dist, 
                 spread = spread, ret_extra = 'fgraph', fast_sgd = FALSE)
    colnames(HU$embedding) = c('HUMAP1', 'HUMAP2')
    rownames(HU$fgraph) = colnames(HU$fgraph) = Cells(SeuratObj)
    SeuratObj[['humap']] <- Seurat::CreateDimReducObject(
        embeddings = HU$embedding,
        assay = 'RNA',
        key = 'HUMAP_',
        global = TRUE
    )
    HU_graph <- Seurat::as.Graph(HU$fgraph)
    DefaultAssay(HU_graph) <- DefaultAssay(SeuratObj)
    SeuratObj[['humap_fgraph']] <- HU_graph
    return(SeuratObj)
}


sc.normalize <- function(x, method="lognorm") {
    print(x)
    mat = LayerData(x, layer='counts')
    print(dim(mat))
    print(dim(mat))
    print("normalize")
    print("variable genes..")

    if (exists("RA_mat_objs/var_stats.rds")) {
        mat.norm = readRDS("RA_mat_objs/lognorm_mat.rds")
        stats = readRDS("RA_mat_objs/var_stats.rds")
    } else {
        mat.norm <- multiply_cols(mat, 1/Matrix::colSums(mat))
        if (c("nCount_Xenium") %in% colnames(x@meta.data)) {
        mat.norm <- log1p(mat.norm * median(x@meta.data$nCount_Xenium))
        } else {
        mat.norm <- log1p(mat.norm * median(x@meta.data$nCount_RNA))
        }
        stats <- matrix_stats(mat.norm, row_stats="variance")
        saveRDS(mat.norm, "RA_mat_objs/lognorm_mat.rds")
        saveRDS(stats, "RA_mat_objs/var_stats.rds")
    }

    if (exists("RA_mat_objs/variable_genes.rds")) {
       variable_genes <- readRDS("RA_mat_objs/variable_genes.rds")
    } else {
       variable_genes <- order(stats$row_stats["variance",], decreasing=TRUE) %>% 
         head(0.2*nrow(mat)) %>% 
         sort()
       saveRDS(variable_genes, "RA_mat_objs/variable_genes.rds")
    }

    print("scaling...")
    if (exists("RA_mat_objs/lognorm_scaledata.rds")) {
        mat_norm <- readRDS("RA_mat_objs/lognorm_scaledata.rds")
    } else {
        mat_norm <- mat[variable_genes,]
        mat_norm <- mat_norm %>% write_matrix_dir(tempfile("mat"))
        gene_means <- stats$row_stats["mean",variable_genes]
        gene_vars <- stats$row_stats["variance", variable_genes]
        mat_norm <- (mat_norm - gene_means) / gene_vars
        saveRDS(mat_norm, "RA_mat_objs/lognorm_scaledata.rds")
    }

    print("pca..")
    set.seed(9)
    if (exists("RA_mat_objs/pca_embeddings.rds")) {
        pca <- readRDS("RA_mat_objs/pca_embeddings.rds")
    } else {
        svd <- BPCells::svds(mat_norm, k=30)
        #svd <- irlba::irlba(mat_norm, nv=50)
        pca <- multiply_cols(svd$v, svd$d)
        saveRDS(pca, "RA_mat_objs/pca_embeddings.rds")
    }

    seu = CreateAssay5Object(counts=mat, data=mat.norm)
    seu = CreateSeuratObject(seu, meta.data=x@meta.data, assay="RNA")
    print(rownames(seu)[1:5])
    print(colnames(seu)[1:5])
    rownames(pca) <- colnames(seu)
    colnames(pca) <- paste0('PC_', seq(1, 30))

    seu[['pca']] <- Seurat::CreateDimReducObject(
        embeddings = pca,
        assay = 'RNA',
        key = 'PC_',
        global = TRUE
    )
    print(seu)

    # run harmony
    print("harmony..")
    set.seed(9)
    if (exists("RA_mat_objs/harmony_pca_embeddings.rds")) {
        harmonyObj <- readRDS("RA_mat_objs/harmony_pca_embeddings.rds")
    } else {
        harmonyObj <- harmony::RunHarmony(
                    data_mat = pca, ## PCA embedding matrix of cells
                    meta_data = x@meta.data, ## dataframe with cell labels
                    vars_use = 'orig.ident', ## variable to integrate out
                    nclust = 15, ## number of clusters in Harmony model
                    max_iter = 10, ## stop after initialization
                    return_object = F ## return the full Harmony model object
                    )
        saveRDS(harmonyObj, "RA_mat_objs/harmony_pca_embeddings.rds")
    }
    print(str(harmonyObj))
    rownames(harmonyObj) <- colnames(seu)
    colnames(harmonyObj) <- paste0('HPC_', seq(1, 30))
    seu[['harmony']] <- Seurat::CreateDimReducObject(
        embeddings = harmonyObj,
        assay = 'RNA',
        key = 'HPC_',
        global = TRUE
    )

    cat(sprintf("PCA dimensions: %s\n", toString(dim(pca))))
    print('umap...')
    set.seed(12341512)
    seu = Run_uwot_umap(seu)
    #failed here with large ef 500-2000
    #Error in names(x) <- value :
    #  'names' attribute [154015] must be the same length as the vector [1]
    #Calls: sc.normalize -> plot_embedding -> collect_features -> colnames<-
    #    In addition: Warning message:
    #In knn_hnsw(pca, ef = 500) :
    #  KNN search didn't find self-neighbor for 16 datapoints. Try higher ef value
    set.seed(999)
    seu <- FindClusters(seu, graph.name = 'humap_fgraph', resolution = 0.8, verbose = TRUE)

    #clusts <- knn_hnsw(pca, k=15, ef=200) %>% # Find approximate nearest neighbors
    #  knn_to_snn_graph() %>% # Convert to a SNN graph
    #  cluster_graph_louvain() # Perform graph-based clustering
    #saveRDS(clusts, "RA_mat_objs/clusts_ef200.rds")

    pdf("RA_mat_objs/umap.pdf", width=10, height=6)
    print(DimPlot(seu, reduction='humap'))
    dev.off()

    seu
}


parser <- ArgumentParser(prog="normalize_xenium.R", description="a wrapper for different normalization in single cells")
parser$add_argument("data", metavar="N", type="character", nargs="+", help="input xenium directory or rds")
parser$add_argument("--output", help="output prefix")
parser$add_argument("--method", help="normalization method")

args = parser$parse_args()

method = args$method
output = args$output
input_data = args$data

adata = readRDS(input_data)
print(dim(adata))
adata@meta.data$orig.ident = basename(as.character(adata@meta.data$sample))
print(head(adata@meta.data))

adata_dbl = adata[, adata@meta.data[, 'scDblFinder.class'] != 'singlet']
adata_dbl = sc.normalize(adata_dbl)
saveRDS(adata_dbl, glue("RA_mat_objs/{args$output}_orig_seg_lognorm_dbl.rds"))

adata_single = adata[, adata@meta.data[, 'scDblFinder.class'] == 'singlet']
adata_single = sc.normalize(adata_single)
saveRDS(adata_single, glue("RA_mat_objs/{args$output}_orig_seg_lognorm_singlet.rds"))

