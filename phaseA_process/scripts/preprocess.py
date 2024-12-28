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


def process_xenium_adata(ad, plot=False):
    sc.pp.calculate_qc_metrics(ad, percent_top=(10, 20, 50, 150), inplace=True)
    cprobes = (
        ad.obs["control_probe_counts"].sum() / ad.obs["total_counts"].sum() * 100
    )
    cwords = (
        ad.obs["control_codeword_counts"].sum() / ad.obs["total_counts"].sum() * 100
    )
    print(f"Negative DNA probe count % : {cprobes}")
    print(f"Negative decoding count % : {cwords}")    
    sc.pp.filter_cells(ad, min_counts=20)
    sc.pp.filter_genes(ad, min_cells=5)
    sc.pp.normalize_total(ad, inplace=True)
    sc.pp.log1p(ad)
    sc.pp.pca(ad)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    sc.tl.leiden(ad)    

    if plot:
        fig, axs = plt.subplots(1, 4, figsize=(15, 4))
        axs[0].set_title("Total transcripts per cell")
        sns.histplot(
            ad.obs["total_counts"],
            kde=False,
            ax=axs[0],
        )
        
        axs[1].set_title("Unique transcripts per cell")
        sns.histplot(
            ad.obs["n_genes_by_counts"],
            kde=False,
            ax=axs[1],
        )
        
        axs[2].set_title("Area of segmented cells")
        sns.histplot(
            ad.obs["cell_area"],
            kde=False,
            ax=axs[2],
        )
        
        axs[3].set_title("Nucleus ratio")
        sns.histplot(
            ad.obs["nucleus_area"] / ad.obs["cell_area"],
            kde=False,
            ax=axs[3],
        )    
        sc.pl.umap(
            ad,
            color=[
                "total_counts",
                "n_genes_by_counts",
                "leiden",
                "nucleus_count"
            ],
            wspace=0.4,
        )    
    return ad


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


def load_data(data_path, out_path):
    sdata = xenium(data_path)
    adata = sdata.tables["table"]
    adata.layers["counts"] = adata.X.copy()
    adata = process_xenium_adata(adata)
    adata.write(out_path.h5ad)
    sanity_check_plot(adata, out_path.fig)


def do_something(data_path, out_path, threads, config):
    # python code
    #../data/kidney/20241025__200743__BWH_20241025_SHRUTI_RACHEL/output-XETG00392__0045655__BS21-N65682A2__20241025__201009
    #h5ad/output-XETG00392__0045655__BS21-N65682A2__20241025__201009.h5ad
    print(data_path.case)
    print(out_path)
    print(config)
    load_data(data_path.case, out_path)


do_something(snakemake.input, snakemake.output, snakemake.threads, snakemake.config)


