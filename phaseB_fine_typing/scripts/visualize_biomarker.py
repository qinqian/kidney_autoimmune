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


def colocalize_analysis(ad, key='celltype_hint2', target="TregFOXP3high"):
    sq.gr.spatial_neighbors(ad, coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(ad, cluster_key=key, numba_parallel=False, show_progress_bar=False, n_jobs=1)
    sq.gr.co_occurrence(ad, cluster_key=key, n_jobs=1, show_progress_bar=False)

    print('*************')
    print(ad.obs['celltype_hint2'])
    #### sq.gr.spatial_neighbors(ad, coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(ad, cluster_key=key, numba_parallel=False, show_progress_bar=False, n_jobs=1)
    sq.gr.co_occurrence(ad, cluster_key=key, n_jobs=1, show_progress_bar=False)
    print('-------------')

    out = ad.uns['celltype_hint2_co_occurrence']['occ']
    interval = ad.uns['celltype_hint2_co_occurrence']['interval'][1:]
    categories = ad.obs[key].cat.categories
    # find cluster neighborhood matrix
    assert target in categories, "f{target} not in [list(categories)]}"

    idx = np.where(categories == target)[0][0]
    df  = pd.DataFrame(out[idx, :, :].T, columns=categories).melt(var_name=key, value_name="probability")
    df["distance"] = np.tile(interval, len(categories))
    return df


def plot_codep(ad, out_fig, conn_key="TregFOXP3high"):
    celltype_hint = np.array(list(ad.obs["celltype_hint"].values))
    x = ad[:, 'FOXP3'].X.toarray().flatten()
    celltype_hint[(x >= 1) & (celltype_hint == 'regulatory T cell')] = 'TregFOXP3high'
    ad.obs.loc[:, 'celltype_hint2'] = pd.Categorical(celltype_hint)
    print(ad.obs.loc[:, 'celltype_hint2'])

    df = colocalize_analysis(ad, key='celltype_hint2')
    sq.pl.co_occurrence(
        ad,
        cluster_key="celltype_hint2",
        clusters=[conn_key, 'CD8-positive, alpha-beta T cell'],
        figsize=(18, 6)
    )
    plt.savefig(out_fig)
    return df


def visualize_spatial_marker(ads, out_fig, marker=["FOXP3"]):
    dfs = []

    fig, ax = plt.subplots(8, 2)
    fig.set_size_inches(18, 12)
    n = 0
    for case in ads.case_h5ads:
        group = case.split('/')[-3]
        ad = sc.read(case)
        celltype_hint = np.array(list(ad.obs["celltype_hint"].values))
        #x = ad[celltype_hint == 'regulatory T cell', 'FOXP3'].X.toarray()
        #ax[0][3].hist(x)
        try:
            df = plot_codep(ad, out_fig.fig2)
            df['sample'] = os.path.basename(case).replace('output−XETG00150__0018462__', '').replace('__20241025__201009_Tcells.h5ad', '').replace('output−XETG00392__0045655__', '').replace('__20240803__183643_Tcells.h5ad', '')
            df['group']  = group
            dfs.append(df)
            sc.pl.umap(ad, ax=ax[n][0], color=marker, show=False)
            sc.pl.umap(ad, ax=ax[n][1], color=['celltype_hint'], show=False)
            n += 1
        except:
            n += 1

    fig.subplots_adjust(right=0.5, left=0.02, top=0.95, bottom=0.1, wspace=0.4, hspace=0.4)
    fig.savefig(out_fig.fig1)

    dfs = pd.concat(dfs, axis=0)
    dfs.to_csv(out_fig.tabs, sep="\t")


def run(data_path, out_path, threads, config):
    visualize_spatial_marker(data_path, out_path)


run(snakemake.input, snakemake.output, snakemake.threads, snakemake.params)

