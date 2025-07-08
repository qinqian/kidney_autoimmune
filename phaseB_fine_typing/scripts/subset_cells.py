import os
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from spatialdata_io import xenium
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation


def retype_T_withscimilarity(ad, plot=False, cell_threshold=10, gene_threshold=4500):

    # Instantiate the CellAnnotation object
    # Set model_path to the location of the uncompressed model
    model_path = "../phaseA_process/notebook/model_v1.1"
    ca = CellAnnotation(model_path=model_path)

    adams = align_dataset(ad, ca.gene_order, gene_overlap_threshold=gene_threshold)
    adams = lognorm_counts(adams)
    adams.obsm["X_scimilarity"] = ca.get_embeddings(adams.X)
    sc.pp.neighbors(adams, use_rep="X_scimilarity")
    sc.tl.umap(adams)

    predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_knn(
        adams.obsm["X_scimilarity"]
    )
    adams.obs["predictions_unconstrained"] = predictions.values

    cell_type_count = adams.obs["predictions_unconstrained"].value_counts()
    cell_type_count = cell_type_count.index[cell_type_count >= cell_threshold]
    ca.safelist_celltypes(cell_type_count)

    adams = ca.annotate_dataset(adams)
    print(adams.obs.columns)

    return adams


def run(data_path, out_path, threads, config):
    df = pd.read_csv(data_path.m_annot)

    adata = sc.read(data_path.h5ad)

    print(config.subtype)
    if config.subtype:
        adata = adata[adata.obs.cell_id.isin(df.loc[df.coarse_ids == "T", :].cell_id), :]
        print(adata.shape)
        adata.X = adata.layers['counts'].copy()
        print(adata.obs.columns)
        print('--------')
        adata = retype_T_withscimilarity(adata, cell_threshold=20)
        print('--------')
        print(adata.obs.columns)
    else:
        # borrow miles' annotation
        print(adata.shape)
        adata = adata[adata.obs.cell_id.isin(df.cell_id), :]
        adata.obs.loc[:, 'preannotation'] = pd.Categorical(adata.obs.cell_id.map(dict(zip(df.cell_id, df.coarse_ids))))
        print(adata.obs['preannotation'].isnull().sum())
        key = 'preannotation'
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
        sq.gr.nhood_enrichment(adata, cluster_key=key, numba_parallel=False, show_progress_bar=False, n_jobs=1)
        sq.gr.co_occurrence(adata, cluster_key=key, n_jobs=1, show_progress_bar=False)

    adata.write(out_path.h5ad)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.params)

