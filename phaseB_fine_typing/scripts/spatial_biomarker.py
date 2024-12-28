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


def visualize_spatial_marker(ads, out_fig, marker=['CXCL9', 'CXCL10', "FOXP3"]):
    fig, ax = plt.subplots(2, 4)
    fig.set_size_inches(25, 10)
    ad = sc.read(ads.case_h5ads[0])
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[0],
        wspace=0,
        ax=ax[0][0],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[1], wspace=0,
        ax=ax[0][1],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[2], wspace=0,
        ax=ax[0][2],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color="celltype_hint", wspace=0,
        ax=ax[0][3],
    )
    ad = sc.read(ads.cont_h5ads[0])
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[0],
        wspace=0,
        ax=ax[1][0],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[1],
        wspace=0,
        ax=ax[1][1],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=marker[2],
        wspace=0,
        ax=ax[1][2],
        legend_loc=None, colorbar=False
    )
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color="celltype_hint", wspace=0,
        ax=ax[1][3],
    )
    fig.subplots_adjust(right=0.5, left=0.02, top=0.95, bottom=0.1)
    plt.savefig(out_fig)


def run(data_path, out_path, threads, config):
    visualize_spatial_marker(data_path, out_path.fig)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.config)

