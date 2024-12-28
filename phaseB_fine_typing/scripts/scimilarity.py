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


def annot_xenium_adata_withscimilarity(ad, plot=False, cell_threshold=10, gene_threshold=4500):
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
    return adams


def sanity_check_plot(ad, out_fig):
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=[
            "leiden",
        ],
        wspace=0.4,
    )
    plt.savefig(out_fig)


def run(data_path, out_path, threads, config):
    print(config)
    adata = sc.read(data_path.h5ad)
    adata = annot_xenium_adata_withscimilarity(adata, cell_threshold=config.cell_count)
    sanity_check_plot(adata, out_path.fig)
    adata.write(out_path.h5ad)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.params)

