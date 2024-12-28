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



def sanity_check_plot(ad, out_fig):
    sc.tl.diffmap(ad, n_comps=10)
    sc.tl.dpt(ad)
    sq.pl.spatial_scatter(
        ad,
        library_id="spatial",
        shape=None,
        color=[
            "leiden",
            "dpt_pseudotime",
            "celltype_hint",
        ],
        wspace=0.4,
    )
    plt.savefig(out_fig)


def run(data_path, out_path, threads, config):
    print(config)
    adata = sc.read(data_path.h5ad)
    sanity_check_plot(adata, out_path.fig)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.params)

