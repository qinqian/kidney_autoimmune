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


def colocalize_analysis(ad, key='celltype_hint2', target="T"):
    out = ad.uns[f'{key}_co_occurrence']['occ']
    interval = ad.uns[f'{key}_co_occurrence']['interval'][1:]
    categories = ad.obs[key].cat.categories
    # find cluster neighborhood matrix
    assert target in categories, "f{target} not in [list(categories)]}"

    idx = np.where(categories == target)[0][0]
    df  = pd.DataFrame(out[idx, :, :].T, columns=categories).melt(var_name=key, value_name="probability")
    df["distance"] = np.tile(interval, len(categories))
    return df


def plot_codep(ad, out_fig, conn_key="T"):
    celltype_hint = np.array(list(ad.obs["preannotation"].values))
    df = colocalize_analysis(ad, key='preannotation')
    #sq.pl.co_occurrence(
    #    ad,
    #    cluster_key="preannotation",
    #    clusters=[conn_key, 'Macrophage', 'Fibroblast'],
    #    figsize=(18, 6)
    #)
    fig, ax = plt.subplots(1, 2, figsize=(13, 7))
    sq.pl.nhood_enrichment(
        ad,
        cluster_key="preannotation",
        figsize=(8, 8),
        title="Neighborhood enrichment adata",
        ax=ax[0],
    )
    sq.pl.spatial_scatter(ad, color="preannotation", shape=None, size=2, ax=ax[1])
    ax[1].legend(title="cell types", bbox_to_anchor=(0.5, -0.1), loc='upper center', ncols=3) # bbox_to_anchor=(0.3, -1, 0.5, 0.5), 
    fig.savefig(out_fig)
    return df


def visualize_spatial_marker(ads, out_fig, marker=["FOXP3"]):
    dfs = []

    fig, ax = plt.subplots(8, 2)
    fig.set_size_inches(18, 12)
    n = 0
    for case in ads.case_h5ads:
        group = case.split('/')[-3]
        ad = sc.read(case)
        celltype_hint = np.array(list(ad.obs["preannotation"].values))
        #x = ad[celltype_hint == 'regulatory T cell', 'FOXP3'].X.toarray()
        #ax[0][3].hist(x)
        try:
            df = plot_codep(ad, out_fig.fig2[n])
            df['sample'] = os.path.basename(case).replace('output−XETG00150__0018462__', '').replace('__20241025__201009_Tcells.h5ad', '').replace('output−XETG00392__0045655__', '').replace('__20240803__183643_Tcells.h5ad', '')
            df['group']  = group
            dfs.append(df)
            sc.pl.umap(ad, ax=ax[n][0], color=marker, show=False)
            sc.pl.umap(ad, ax=ax[n][1], color=['preannotation'], show=False)
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

