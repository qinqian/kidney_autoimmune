import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from spatialdata_io import xenium
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellAnnotation


def visualize_spatial_marker(ads, out_fig, marker=['CXCL9', 'CXCL10', "CXCL11"]):
    fig, ax = plt.subplots(8, 3)
    fig.set_size_inches(12, 10)
    groups = []
    for index, case in enumerate(ads.case_h5ads):
        group = case.split('/')[-3]
        groups.append(group)

    order = np.argsort(groups)
    groups = np.array(groups)[order]
    ads.case_h5ads = np.array(ads.case_h5ads)[order]

    for index, case in enumerate(ads.case_h5ads):
        group = case.split('/')[-3]
        ad = sc.read(case)
        for idx, m in enumerate(marker):
            sq.pl.spatial_scatter(
                ad,
                library_id="spatial",
                shape=None,
                color=m,
                wspace=0, hspace=0,
                ax=ax[index][idx],
                legend_loc=None, colorbar=False, frameon=False, fig=fig)
            ax[index][idx].set_axis_off()
            ax[index][idx].set_title(f"{group} {m}")

        #sq.pl.spatial_scatter(
        #    ad,
        #    library_id="spatial",
        #    shape=None,
        #    color="celltype_hint", wspace=0,
        #    ax=ax[index][3],
        #)

    fig.subplots_adjust(right=0.95, left=0.01, top=0.95, bottom=0.05, wspace=0, hspace=0.5)
    plt.savefig(out_fig)


def run(data_path, out_path, threads, config):
    visualize_spatial_marker(data_path, out_path.fig)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.config)

