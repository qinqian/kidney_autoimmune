from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scvelo as scv
from numpy import ndarray
from scipy.stats import mannwhitneyu, rankdata

##from scanorama import correct_scanpy
scv.logging.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")
np.random.seed(99)

from scipy.sparse import csr_matrix, issparse

__all__ = [
    "census_normalize",
    "remove_zero_mvg",
    "compute_similarity2",
    "compute_similarity1",
    "compute_gcs",
    "convert_to_markov",
    "any",
    "find_nonzero",
    "FNNLSa",
    "nu",
    "eps",
    "regressed",
    "diffused",
    "compare_cytotrace",
    "cytotrace",
    "batch_cytotrace",
    "compare_cytotrace_ncores",
    "cytotrace_ncore",
    "align_diffrate",
    "plot_multiadata",
    "cumulative_boxplot",
    "run_cytotrace",
]


# TODO: refactor to subpackage exposing a minimal public interface
# TODO: migrate all justifiably eager imports to the tops of new modules
# TODO: enable runtime type-checking
# TODO: remove unused comments and print statements
# TODO: add unit tests


def census_normalize(mat, count):
    "RNA-seq census normalization to correct cell lysis"
    if issparse(mat):
        # x = mat.copy()
        # x.data = 2**x.data - 1
        # gene x cell
        x = mat.multiply(count).multiply(1.0 / mat.sum(axis=0)).tocsr()
        return x.log1p()
    else:
        x = 2**mat - 1
        # gene x cell
        x = (x * count) / x.sum(axis=0)
        return np.log2(x + 1)


def remove_zero_mvg(mat):
    "remove cells not expressing any of the top 1000 variable genes"
    # gene x cell
    cell_count = (mat > 0).sum(axis=1)
    X_topcell = mat[cell_count >= 0.05 * mat.shape[1], :]
    var = X_topcell.var(axis=1)
    mean = X_topcell.mean(axis=1)
    disp = var / mean
    return X_topcell[
        disp >= (np.sort(disp)[::-1])[:1000][-1], :
    ]  # in case less than 1000 genes


def compute_similarity2(O: ndarray, P: ndarray) -> ndarray:
    "Compute pearson correlation between two matrices O and P using einstein summation"
    (n, t) = O.shape  # n traces of t samples
    (n_bis, m) = P.shape  # n predictions for each of m candidates

    DO = O - (
        np.einsum("nt->t", O, optimize="optimal") / np.double(n)
    )  # compute O - mean(O)
    DP = P - (
        np.einsum("nm->m", P, optimize="optimal") / np.double(n)
    )  # compute P - mean(P)

    cov = np.einsum("nm,nt->mt", DP, DO, optimize="optimal")

    varP = np.einsum("nm,nm->m", DP, DP, optimize="optimal")
    varO = np.einsum("nt,nt->t", DO, DO, optimize="optimal")
    tmp = np.einsum("m,t->mt", varP, varO, optimize="optimal")
    return cov / np.sqrt(tmp)


# adapt from https://github.com/ikizhvatov/efficient-columnwise-correlation/blob/master/columnwise_corrcoef_perf.py
def compute_similarity1(A):
    "Compute pairwise correlation of all columns in matrices A"
    n, t = A.shape  # genes x samples
    DO = A - (
        np.einsum("nt->t", A, optimize="optimal") / np.double(n)
    )  # compute O - mean(O)
    cov = np.einsum("nm,nt->mt", DO, DO, optimize="optimal")
    varP = np.einsum("nt,nt->t", DO, DO, optimize="optimal")
    return cov / np.sqrt(np.einsum("m,t->mt", varP, varP, optimize="optimal"))


def compute_gcs(mat, count, top_n_genes=200):
    "Compute gene set enrichment scores by correlating gene count and gene expression"
    corrs = compute_similarity2(mat.T, count.reshape(-1, 1))[0, :]
    # avoid nan, put all nan to -1
    corrs[np.isnan(corrs)] = -1
    gcs = mat[np.argsort(corrs)[::-1][:top_n_genes], :].mean(axis=0)
    return gcs


def convert_to_markov(sim):
    """Convert the Pearson correlation to Markov matrix

    TODO: use velocity graph to replace this markov matrix
    """
    Ds = np.copy(sim)

    cutoff = np.mean(Ds)
    np.fill_diagonal(Ds, 0)

    Ds[Ds <= 0] = 0
    Ds[Ds <= cutoff] = 0

    zero_rows = Ds.sum(axis=1) == 0
    zero_cols = Ds.sum(axis=0) == 0

    Ds[~zero_rows, :] = (
        Ds[~zero_rows, :].T / Ds[~zero_rows, :].sum(axis=1)
    ).T  # tips
    Ds[:, zero_cols] = 0
    return Ds


# http://xrm.phys.northwestern.edu/research/pdf_papers/1997/bro_chemometrics_1997.pdf
# matlab reference: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/3388/versions/1/previews/fnnls.m/index.html

from numpy import (
    abs,
    arange,
    argmax,
    finfo,
    float64,
    int64,
    min,
    newaxis,
    nonzero,
    sum,
    zeros,
)
from scipy.linalg import solve

nu = newaxis
import numpy as np

# machine epsilon
eps = finfo(float64).eps


def any(a):
    # assuming a vector, a
    larger_than_zero = sum(a > 0)
    if larger_than_zero:
        return True
    else:
        return False


def find_nonzero(a):
    # returns indices of nonzero elements in a
    return nonzero(a)[0]


def FNNLSa(XtX, Xty, tol=None):
    """Faster NNLS imported from https://github.com/delnatan/FNNLSa
    A fast non-negativity-constrained least squares algorithm. Journal of chemometrics
    """

    if tol is None:
        tol = eps

    M, N = XtX.shape
    # initialize passive set, P. Indices where coefficient is >0
    P = zeros(N, dtype=int64)
    # and active set. Indices where coefficient is <=0
    Z = arange(N) + 1
    # working active set
    ZZ = arange(N) + 1
    # initial solution vector, x
    x = zeros(N, dtype=float64)
    # weight vector
    w = Xty - XtX @ x
    # iteration counts and parameter
    it = 0
    itmax = 30 * N
    # MAIN LOOP
    # continue as long as there are indices within the active set Z
    # or elements in inner loop active set is larger than 'tolerance'
    piter = 0
    while any(Z) and any(w[ZZ - 1] > tol):
        piter += 1
        t = argmax(w[ZZ - 1]) + 1  # find largest weight
        t = ZZ[t - 1]
        P[t - 1] = t  # move to passive set
        Z[t - 1] = 0  # remove from active set
        PP = find_nonzero(P) + 1
        ZZ = find_nonzero(Z) + 1
        NZZ = ZZ.shape

        # compute trial solution, s
        s = zeros(N, dtype=float64)

        if len(PP) == 1:
            s[PP - 1] = Xty[PP - 1] / XtX[PP - 1, PP - 1]
        else:
            s[PP - 1] = solve(XtX[PP - 1, PP[:, nu] - 1], Xty[PP - 1])
        s[ZZ - 1] = 0.0  # set active coefficients to 0

        while any(s[PP - 1] <= tol) and it < itmax:
            it = it + 1
            QQ = find_nonzero((s <= tol) * P) + 1
            alpha = min(x[QQ - 1] / (x[QQ - 1] - s[QQ - 1]))
            x = x + alpha * (s - x)
            ij = find_nonzero((abs(x) < tol) * (P != 0)) + 1
            Z[ij - 1] = ij
            P[ij - 1] = 0
            PP = find_nonzero(P) + 1
            ZZ = find_nonzero(Z) + 1
            if len(PP) == 1:
                s[PP - 1] = Xty[PP - 1] / XtX[PP - 1, PP - 1]
            else:
                s[PP - 1] = solve(XtX[PP - 1, PP[:, nu] - 1], Xty[PP - 1])
            s[ZZ - 1] = 0.0
        # assign current solution, s, to x
        x = s
        # recompute weights
        w = Xty - XtX @ x
    return x, w


def regressed(markov, gcs, solver="fnnls"):
    """solve markov @ weight = gcs problems,

    solver: fnnls (default) is faster in larger dataset, e.g., above 20,000 cells
            nnls is faster in smaller dataset, e.g. less than 5,000 cells
    """
    print(solver)
    if solver == "fnnls":
        coef, err = FNNLSa(markov.T @ markov, markov.T @ gcs)
    elif (
        solver == "jfnnls"
    ):  # pip install fnnls # https://github.com/jvendrow/fnnls
        from fnnls import fnnls

        coef, err = fnnls(markov, gcs)
    elif solver == "nnls":
        from scipy.optimize import nnls

        print((markov > 0).sum())
        coef, err = nnls(markov, gcs)  # from 1987...
    elif solver == "lsq_linear":
        from scipy import sparse
        from scipy.optimize import lsq_linear

        # slow as well
        sol = lsq_linear(
            sparse.csr_matrix(markov),
            gcs,
            bounds=(0, np.inf),
            lsmr_tol="auto",
            verbose=1,
        )  # from 1997...
        coef = sol.x
    elif solver == "lasso":
        from scipy import sparse
        from sklearn.linear_model import Lasso

        lasso = Lasso(
            alpha=1e-7, max_iter=100000, fit_intercept=False, positive=True
        )
        den = (markov > 0).sum() / np.prod(markov.shape)
        print(den)
        if den >= 0.5:
            lasso.fit(markov, gcs)
        else:
            # lasso.fit(sparse.csr_matrix(markov), gcs)
            lasso.fit(sparse.coo_matrix(markov), gcs)
        coef = lasso.coef_
    elif solver == "nnlsgd":  # somehow not consistent with nnls and fnnls...
        coef = NNLS_GD(markov.T @ markov, markov.T @ gcs)
    elif solver == None:  # do not fit
        return gcs
    else:  # tnt-nn
        raise NotImplementedError
    # faster nnls
    return np.dot(markov, coef)


def diffused(markov, gcs, ALPHA=0.9):
    """Compute diffusion process"""
    v_prev = np.copy(gcs)
    v_curr = np.copy(gcs)

    for i in range(10000):
        v_prev = np.copy(v_curr)
        v_curr = ALPHA * (markov.dot(v_curr)) + (1 - ALPHA) * gcs

        diff = np.mean(np.abs(v_curr - v_prev))
        if diff <= 1e-6:
            break
    return v_curr


from scipy.sparse import issparse


def compare_cytotrace(
    adata,
    layer="all",
    cell_count=10,
    condition="age",
    solver="nnls",
    is_normalized=False,
    n_cores=4,
    top_n_genes=200,
):
    "Main interface of cytotrace reimplementation used for single dataset with multiple conditions"
    # condition common steps

    proc_time = time()
    if not is_normalized:
        if layer == "all":
            X = (adata.layers["spliced"] + adata.layers["unspliced"]).toarray()
        else:
            try:
                X = (
                    adata.layers[layer].toarray()
                    if issparse(adata.layers[layer])
                    else adata.layers[layer]
                )
            except:
                X = adata.X.toarray() if issparse(adata.X) else adata.X

        cells_selected = np.arange(X.shape[0])
        genes_selected = np.arange(X.shape[1])
        pgenes = (pd.isnull((X > 0).sum(axis=0))) | (X.var(axis=0) == 0)
        # cell x gene
        X = X[:, ~pgenes]
        X = (X.T / X.sum(axis=1)) * 1e6

        pqcells = (pd.isnull((X > 0).sum(axis=0))) | (
            (X > 0).sum(axis=0) <= cell_count
        )
        X = X[:, ~pqcells]

        genes_selected = genes_selected[~pgenes]
        cells_selected = cells_selected[~pqcells]

        X = np.log2(X + 1)
        counts = (X > 0).sum(axis=0)
        mat2 = census_normalize(X, counts)

        mvg = remove_zero_mvg(mat2)
        cells_selected = cells_selected[mvg.sum(axis=0) != 0]

        mat2 = mat2[:, mvg.sum(axis=0) != 0]
        counts = counts[mvg.sum(axis=0) != 0]

        # joint two conditions to determine gcs feature
        # using the same set of genes
        gcs = compute_gcs(mat2, counts, top_n_genes=top_n_genes)
        adata.uns["gcs"] = gcs
        adata.uns["cells_selected"] = cells_selected
        adata.uns["counts"] = counts
        adata.uns["mat2"] = mat2
        adata.uns["mvg"] = mvg
        adata.uns["genes_selected"] = genes_selected
    else:
        cells_selected = adata.uns["cells_selected"]
        counts = adata.uns["counts"]
        mat2 = adata.uns["mat2"]
        mvg = adata.uns["mvg"]
        genes_selected = adata.uns["genes_selected"]
        gcs = adata.uns["gcs"]

    # for indepent gcs feature
    # using different set of genes is very variable
    # gcs = np.zeros(len(cells_selected))

    # condition-specific steps
    for cond in np.unique(adata.obs.loc[:, condition].values[cells_selected]):
        selection = adata.obs.loc[:, condition].values[cells_selected] == cond

        gcs_cond = gcs[selection]
        counts_cond = counts[selection]
        mat2_cond = mat2[:, selection]

        # compute gcs feature independently for two conditions
        # this would select different gene set, and make the
        # cytotrace score incomparable
        # gcs_cond = compute_gcs(mat2_cond, counts_cond)

        mvg_cond = mvg[:, selection]

        trans_time = time()
        print("preprocessing: %s" % (trans_time - proc_time))
        # get transition matrix
        corr_cond = compute_similarity1(mvg_cond)
        markov_cond = convert_to_markov(corr_cond)

        gc_time = time()
        print("markov: %s" % (gc_time - trans_time))

        # get Gene count signature
        gcs_cond = regressed(markov_cond, gcs_cond, solver=solver)
        gcs_cond = diffused(markov_cond, gcs_cond)
        gcs[selection] = gcs_cond

        # condition specific cytotrace correlated gene sets
        scores_cond = rankdata(gcs_cond) / gcs_cond.shape[0]
        score_time = time()
        print("score time: %s" % (score_time - gc_time))

        corrs_cond = compute_similarity2(
            mat2_cond.T, scores_cond.reshape(-1, 1)
        )[0, :]
        final_time = time()
        print("final time: %s" % (final_time - score_time))
        adata.var["cytotrace_corrs_%s" % cond] = None
        adata.var.iloc[
            genes_selected,
            adata.var.columns.get_loc("cytotrace_corrs_%s" % cond),
        ] = corrs_cond

    adata.obs["gcs"] = None
    adata.obs["cytotrace"] = None
    adata.obs["counts"] = None
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("gcs")
    ] = gcs  # used for re-rank when integrating multiple-samples
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("cytotrace")] = (
        rankdata(gcs) / gcs.shape[0]
    )  # cytotrace scores
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("counts")
    ] = counts  # used for re-rank when integrating multiple-samples


def cytotrace_sparse(
    adata,
    layer="raw",
    cell_count=20,
    solver="nnls",
    top_n_features=200,
    skip_regress=False,
):
    "optimized version"
    proc_time = time()
    if not issparse(adata.layers[layer]):
        raise NotImplementedError
    else:
        X = adata.layers[layer]

    cells_selected = np.arange(X.shape[0])
    features_selected = np.arange(X.shape[1])

    n_cells = X.shape[0]
    feature_mean = (X.sum(0) / n_cells).A1
    feature_mean_sq = (X.multiply(X).sum(0) / n_cells).A1
    feature_var = (feature_mean_sq - feature_mean**2) * (
        n_cells / (n_cells - 1)
    )
    non_zero_features = (X > 0).sum(axis=0).A1
    pfeatures = np.isnan(non_zero_features) | (feature_var == 0)

    # cell x feature, filter feature
    X = X[:, ~pfeatures]
    features_selected = features_selected[~pfeatures]

    # gene x cell
    X = X.T
    X = X.multiply(1.0 / csr_matrix.sum(X, axis=0)) * 1e6
    X = X.tocsr()

    # filter cell, feature x cell
    feature_count_per_cell = (X > 0).sum(axis=0).A1
    pcells = np.isnan(feature_count_per_cell) | (
        feature_count_per_cell < cell_count
    )

    # feature x cell
    X = X[:, ~pcells]
    cells_selected = cells_selected[~pcells]
    feature_count_per_cell = feature_count_per_cell[~pcells]

    # census normalize, feature x cell
    census_X = (
        X.multiply(feature_count_per_cell)
        .multiply(1.0 / csr_matrix.sum(X, axis=0))
        .log1p()
        .tocsr()
    )

    # top variable features, feature x cell
    cell_count_per_feature = (census_X > 0).sum(axis=1).A1
    census_X_topcell = census_X[
        cell_count_per_feature >= 0.05 * census_X.shape[1], :
    ].tocsr()
    n_cells = census_X_topcell.shape[1]
    feature_mean = (census_X_topcell.sum(1) / n_cells).A1
    feature_mean_sq = (
        census_X_topcell.multiply(census_X_topcell).sum(1) / n_cells
    ).A1
    feature_var = (feature_mean_sq - feature_mean**2) * (
        n_cells / (n_cells - 1)
    )
    disp = feature_var / feature_mean
    mvg = census_X_topcell[
        disp >= (np.sort(disp)[::-1])[:1000][-1], :
    ]  # in case less than 1000 features

    # top 1000 variable features for markov matrix
    selection = (mvg.sum(axis=0) != 0).A1
    mvg = mvg[:, selection]
    corr = compute_similarity1(mvg.A)
    markov = convert_to_markov(corr)

    # filter census output with nonzero cells in top variable features matrix
    cells_selected = cells_selected[selection]
    census_X = census_X[:, selection]
    feature_count_per_cell = feature_count_per_cell[selection]

    # calculate GC scores
    corrs = compute_similarity2(
        census_X.T.A, feature_count_per_cell.reshape(-1, 1)
    )[0, :]
    # avoid nan, put all nan to -1
    corrs[np.isnan(corrs)] = -1
    gcs = census_X[np.argsort(corrs)[::-1][:top_n_features], :].mean(axis=0).A1

    if not skip_regress:
        # print("regress...")
        # from scipy.optimize import lsq_linear
        # # slow as well
        # sol = lsq_linear(csr_matrix(markov), gcs, bounds=(0, np.inf), lsmr_tol='auto', verbose=1)    # from 1997...
        # coef = sol.x

        from scipy.optimize import nnls

        coef, err = nnls(markov, gcs)  # from 1987...
        gcs = np.dot(markov, coef)

    # print(markov.shape, gcs.shape)
    gcs = diffused(markov, gcs)
    rank = rankdata(gcs)
    scores = rank / gcs.shape[0]

    # print(gcs)
    adata.obs["gcs"] = np.nan
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("gcs")] = gcs

    cytoGenes = compute_similarity2(census_X.T.A, scores.reshape(-1, 1))[0, :]

    adata.obs["cytotrace"] = np.nan
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("cytotrace")
    ] = scores

    adata.obs["counts"] = np.nan
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("counts")] = gcs

    adata.var["cytotrace"] = False
    adata.var.iloc[
        features_selected, adata.var.columns.get_loc("cytotrace")
    ] = True

    adata.var["cytotrace_corrs"] = np.nan
    adata.var.iloc[
        features_selected, adata.var.columns.get_loc("cytotrace_corrs")
    ] = np.array(corrs, dtype=np.float32)
    return {
        "CytoTRACE": scores,
        "CytoTRACErank": rank,
        "GCS": gcs,
        "Counts": feature_count_per_cell,
        "exprMatrix": census_X,
        "cytoGenes": cytoGenes,
        "gcsGenes": corrs,
        "filteredCells": np.setdiff1d(
            np.arange(adata.shape[0]), cells_selected
        ),
    }


def cytotrace(
    adata, layer="all", cell_count=10, solver="nnls", top_n_genes=200
):
    "Main interface of cytotrace reimplementation used for single dataset with one condition"
    proc_time = time()
    if layer == "all":
        X = (adata.layers["spliced"] + adata.layers["unspliced"]).toarray()
    #        X = (adata.layers['spliced'] + adata.layers['unspliced'])
    else:
        try:
            X = (
                adata.layers[layer].toarray()
                if issparse(adata.layers[layer])
                else adata.layers[layer]
            )
            # X = adata.layers[layer] if issparse(adata.layers[layer]) else adata.layers[layer]
        except:
            X = adata.X.toarray() if issparse(adata.X) else adata.X
            # X = adata.X if issparse(adata.X) else adata.X

    cells_selected = np.arange(X.shape[0])
    genes_selected = np.arange(X.shape[1])

    pgenes = (pd.isnull((X > 0).sum(axis=0))) | (X.var(axis=0) == 0)
    # cell x gene
    X = X[:, ~pgenes]
    # gene x cell
    X = (X.T / X.sum(axis=1)) * 1e6

    pqcells = (pd.isnull((X > 0).sum(axis=0))) | (
        (X > 0).sum(axis=0) <= cell_count
    )
    X = X[:, ~pqcells]

    genes_selected = genes_selected[~pgenes]
    cells_selected = cells_selected[~pqcells]

    X = np.log2(X + 1)
    counts = (X > 0).sum(axis=0)
    mat2 = census_normalize(X, counts)

    mvg = remove_zero_mvg(mat2)
    selection = mvg.sum(axis=0) != 0
    mvg = mvg[:, selection]
    mat2 = mat2[:, selection]
    counts = counts[selection]

    cells_selected = cells_selected[selection]

    trans_time = time()
    print("processing: %s" % (trans_time - proc_time))

    # get transition matrix
    corr = compute_similarity1(mvg)
    markov = convert_to_markov(corr)

    gc_time = time()
    print("markov: %s" % (gc_time - trans_time))

    # get Gene count signature
    # gcs = compute_gcs(mat2, counts)
    # gcs = compute_gcs(mat2, counts, top_n_genes=top_n_genes)

    corrs = compute_similarity2(mat2.T, counts.reshape(-1, 1))[0, :]
    # avoid nan, put all nan to -1
    corrs[np.isnan(corrs)] = -1
    gcs = mat2[np.argsort(corrs)[::-1][:top_n_genes], :].mean(axis=0)

    gcs_time = time()
    print("gcs: %s" % (gcs_time - gc_time))

    gcs = regressed(markov, gcs, solver=solver)
    regress_time = time()
    print("regression: %s" % (regress_time - gcs_time))

    gcs = diffused(markov, gcs)
    rank = rankdata(gcs)
    scores = rankdata(gcs) / gcs.shape[0]
    score_time = time()
    print("diffusion:%s " % (score_time - regress_time))

    cyto_corrs = compute_similarity2(mat2.T, scores.reshape(-1, 1))[0, :]
    gene_time = time()
    print("genes:%s" % (gene_time - score_time))

    adata.obs["gcs"] = None
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("gcs")] = gcs

    adata.obs["cytotrace"] = None
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("cytotrace")
    ] = scores

    adata.obs["counts"] = None
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("counts")] = counts

    adata.var["cytotrace"] = None
    adata.var.iloc[
        genes_selected, adata.var.columns.get_loc("cytotrace")
    ] = True

    adata.var["cytotrace_corrs"] = None
    adata.var.iloc[
        genes_selected, adata.var.columns.get_loc("cytotrace_corrs")
    ] = cyto_corrs
    return {
        "CytoTRACE": scores,
        "CytoTRACErank": rank,
        "GCS": gcs,
        "Counts": counts,
        "exprMatrix": mat2,
        "cytoGenes": cyto_corrs,
        "gcsGenes": corrs,
        "filteredCells": np.setdiff1d(
            np.arange(adata.shape[0]), cells_selected
        ),
    }


def visualize(adata, metrics, name="test"):
    import scanpy as sc
    import seaborn as sns

    fig, ax = plt.subplots(2)
    fig.set_size_inches(8, 6)
    sc.pl.umap(adata, color="cytotrace", show=False, ax=ax[0])
    order = (
        adata.obs.groupby("clusters").median().sort_values(["cytotrace"]).index
    )
    bplot = sns.boxplot(
        x="clusters",
        y="cytotrace",
        data=adata.obs,
        order=order,
        width=0.35,
        ax=ax[1],
    )
    ax[1].set_xticklabels(order, rotation=40, ha="right")
    fig.savefig(f"{name}_figure.pdf")


def batch_cytotrace(mvg_batch, gcs_batch, solver="jfnnls"):
    corr_batch = compute_similarity1(mvg_batch)
    markov_batch = convert_to_markov(corr_batch)
    gcs_time = time()
    gcs_batch = regressed(markov_batch, gcs_batch, solver=solver)
    regress_time = time()
    print("regression: %s" % (regress_time - gcs_time))
    print(markov_batch.shape, gcs_batch.shape)
    gcs_batch = diffused(markov_batch, gcs_batch)
    return gcs_batch


from scipy.sparse import issparse


def compare_cytotrace_ncores(
    adata,
    layer="all",
    cell_count=10,
    condition="age",
    solver="nnls",
    is_normalized=False,
    ncores=4,
    batch_cell=2000,
):
    "Main interface of cytotrace reimplementation used for single dataset with multiple conditions"
    # condition common steps

    proc_time = time()
    if not is_normalized:
        if layer == "all":
            X = (
                adata.layers["spliced"] + adata.layers["unspliced"]
            )  # .toarray()
        else:
            try:
                X = (
                    adata.layers[layer]
                    if issparse(adata.layers[layer])
                    else adata.layers[layer]
                )
            except:
                X = adata.X if issparse(adata.X) else adata.X

        cells_selected = np.arange(X.shape[0])
        genes_selected = np.arange(X.shape[1])

        from sklearn.utils import sparsefuncs_fast

        mean, var = sparsefuncs_fast.csr_mean_variance_axis0(X)

        pgenes = (np.isnan(np.array((X > 0).sum(axis=0))[0])) | (var == 0)

        # cell x gene
        X = X[:, ~pgenes]
        X = (X.T).multiply(1.0 / csr_matrix.sum(X.T, axis=0)).tocsr()
        pqcells = np.array(
            (np.isnan((X > 0).sum(axis=0)))
            | ((X > 0).sum(axis=0) <= cell_count)
        )[0]
        X = X[:, ~pqcells]

        genes_selected = genes_selected[~pgenes]
        cells_selected = cells_selected[~pqcells]

        counts = np.array((X > 0).sum(axis=0))[0]
        mat2 = census_normalize(X, counts).toarray()
        mvg = remove_zero_mvg(mat2)

        # sel = np.array(mvg.sum(axis=0))[0]!=0
        sel = mvg.sum(axis=0) != 0
        mat2 = mat2[:, sel]
        counts = counts[sel]
        cells_selected = cells_selected[sel]

        # joint two conditions to determine gcs feature
        # using the same set of genes
        gcs = compute_gcs(mat2, counts)
        adata.uns["gcs"] = gcs
        adata.uns["cells_selected"] = cells_selected
        adata.uns["counts"] = counts
        adata.uns["mat2"] = mat2
        adata.uns["mvg"] = mvg
        adata.uns["genes_selected"] = genes_selected
    else:
        cells_selected = adata.uns["cells_selected"]
        counts = adata.uns["counts"]
        mat2 = adata.uns["mat2"]
        mvg = adata.uns["mvg"]
        genes_selected = adata.uns["genes_selected"]
        gcs = adata.uns["gcs"]

    # for indepent gcs feature
    # using different set of genes is very variable
    # gcs = np.zeros(len(cells_selected))

    # parallel part
    import multiprocessing as mp

    # https://stackoverflow.com/questions/8533318/multiprocessing-pool-when-to-use-apply-apply-async-or-map
    # condition-specific steps
    pool = mp.Pool(ncores)

    for cond in np.unique(adata.obs.loc[:, condition].values[cells_selected]):
        print(cond)
        selection = adata.obs.loc[:, condition].values[cells_selected] == cond

        gcs_cond = gcs[selection]
        counts_cond = counts[selection]
        mat2_cond = mat2[:, selection]
        mvg_cond = mvg[:, selection]

        mvg_gcs_batches = []
        for i in range(0, gcs_cond.shape[0], batch_cell):
            mvg_batch = mvg_cond[:, i : (i + batch_cell)].copy()
            gcs_batch = gcs_cond[i : i + batch_cell].copy()
            print(mvg_batch.shape, gcs_batch.shape)
            mvg_gcs_batches.append([mvg_batch, gcs_batch, solver])

        # https://stackoverflow.com/questions/8533318/multiprocessing-pool-when-to-use-apply-apply-async-or-map
        # http://jennguyen1.github.io/nhuyhoa/software/Parallel-Processing.html
        result = pool.starmap_async(batch_cytotrace, mvg_gcs_batches)
        result.wait()

        gcs_batches = result.get()
        # gcs_batches = result.get()
        gcs_cond = np.hstack(gcs_batches)

        trans_time = time()
        print("preprocessing: %s" % (trans_time - proc_time))

        gc_time = time()
        print("markov: %s" % (gc_time - trans_time))
        gcs[selection] = gcs_cond

        # condition specific cytotrace correlated gene sets
        scores_cond = rankdata(gcs_cond) / gcs_cond.shape[0]
        score_time = time()
        print("score time: %s" % (score_time - gc_time))

        corrs_cond = compute_similarity2(
            mat2_cond.T, scores_cond.reshape(-1, 1)
        )[0, :]
        final_time = time()
        print("final time: %s" % (final_time - score_time))
        adata.var["cytotrace_corrs_%s" % cond] = None
        adata.var.iloc[
            genes_selected,
            adata.var.columns.get_loc("cytotrace_corrs_%s" % cond),
        ] = corrs_cond

    pool.close()
    pool.join()

    adata.obs["gcs"] = None
    adata.obs["cytotrace"] = None
    adata.obs["counts"] = None
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("gcs")
    ] = gcs  # used for re-rank when integrating multiple-samples
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("cytotrace")] = (
        rankdata(gcs) / gcs.shape[0]
    )  # cytotrace scores
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("counts")
    ] = counts  # used for re-rank when integrating multiple-samples


def cytotrace_ncore(
    adata,
    layer="all",
    cell_count=10,
    solver="nnls",
    ncores=4,
    batch_cell=3000,
    shuffle=3,
):
    "optimized version"
    proc_time = time()
    if layer == "all":
        X = adata.layers["spliced"] + adata.layers["unspliced"]
    else:
        try:
            X = (
                adata.layers[layer]
                if issparse(adata.layers[layer])
                else adata.layers[layer]
            )
        except:
            X = adata.X if issparse(adata.X) else adata.X

    cells_selected = np.arange(X.shape[0])
    genes_selected = np.arange(X.shape[1])

    from sklearn.utils import sparsefuncs_fast

    mean, var = sparsefuncs_fast.csr_mean_variance_axis0(X)

    pgenes = (np.isnan(np.array((X > 0).sum(axis=0))[0])) | (var == 0)
    # pgenes = np.isnan(mean * adata.shape[0]) | (var == 0)

    # cell x gene
    X = X[:, ~pgenes]
    # gene x cell
    # X = (X.T / np.array(X.sum(axis=1))[0]) * 1e6
    # sparse operation
    X = (X.T).multiply(1.0 / csr_matrix.sum(X.T, axis=0)).tocsr()

    pqcells = np.array(
        (np.isnan((X > 0).sum(axis=0))) | ((X > 0).sum(axis=0) <= cell_count)
    )[0]
    X = X[:, ~pqcells]

    genes_selected = genes_selected[~pgenes]
    cells_selected = cells_selected[~pqcells]

    counts = np.array((X > 0).sum(axis=0))[0]
    mat2 = census_normalize(X, counts).toarray()

    # cell_count = np.array((mat2 > 0).sum(axis=1))[:, 0]
    # X_topcell = mat2[cell_count >= 0.05 * mat2.shape[1], :]

    # mean, var = sparsefuncs_fast.csr_mean_variance_axis0(X_topcell)
    # disp = var / mean
    # mvg = X_topcell[disp >= (np.sort(disp)[::-1])[:1000][-1], :]
    mvg = remove_zero_mvg(mat2)

    # sel = np.array(mvg.sum(axis=0))[0]!=0
    sel = mvg.sum(axis=0) != 0
    mat2 = mat2[:, sel]
    counts = counts[sel]
    cells_selected = cells_selected[sel]

    # mat2 = mat2.toarray()
    # mvg = mvg.toarray()
    print(mvg.shape, mat2.shape)

    trans_time = time()
    print("processing: %s" % (trans_time - proc_time))

    gcs = compute_gcs(mat2, counts)

    # get transition matrix
    # parallel part
    import multiprocessing as mp

    pool = mp.Pool(ncores)

    # https://stackoverflow.com/questions/8533318/multiprocessing-pool-when-to-use-apply-apply-async-or-map
    mvg_gcs_batches = []
    for i in range(0, adata.shape[0], batch_cell):
        mvg_batch = mvg[:, i : (i + batch_cell)].copy()
        gcs_batch = gcs[i : i + batch_cell].copy()
        mvg_gcs_batches.append([mvg_batch, gcs_batch, solver])

    markov_diffusion_time = time()

    # https://stackoverflow.com/questions/8533318/multiprocessing-pool-when-to-use-apply-apply-async-or-map
    result = pool.starmap_async(batch_cytotrace, mvg_gcs_batches)
    pool.close()
    pool.join()
    score_time = time()
    print(score_time - markov_diffusion_time)
    gcs_batches = result.get()
    gcs_batches = np.hstack(gcs_batches)

    scores = rankdata(gcs_batches) / gcs_batches.shape[0]
    corrs = compute_similarity2(mat2.T, scores.reshape(-1, 1))[0, :]

    gene_time = time()
    print("genes:%s" % (gene_time - score_time))
    adata.obs["gcs"] = None
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("gcs")] = gcs

    adata.obs["cytotrace"] = None
    adata.obs.iloc[
        cells_selected, adata.obs.columns.get_loc("cytotrace")
    ] = scores

    adata.obs["counts"] = None
    adata.obs.iloc[cells_selected, adata.obs.columns.get_loc("counts")] = counts

    adata.var["cytotrace"] = None
    adata.var.iloc[
        genes_selected, adata.var.columns.get_loc("cytotrace")
    ] = True

    adata.var["cytotrace_corrs"] = None
    adata.var.iloc[
        genes_selected, adata.var.columns.get_loc("cytotrace_corrs")
    ] = corrs


def align_diffrate(
    adatas, labels, field="condition", type="A", outfield="cytotrace", ax=None
):
    "this is used for differentiation rate comparison across samples"
    scores = [
        adatas[0].obs.loc[:, outfield][adatas[0].obs.loc[:, field] == type],
        adatas[1].obs.loc[:, outfield][adatas[1].obs.loc[:, field] == type],
    ]
    if outfield in ["latent_time", "velocity_pseudotime"]:
        pvalg = mannwhitneyu(scores[0], scores[1], alternative="greater")[1]
        pvall = mannwhitneyu(scores[0], scores[1], alternative="less")[1]
    else:
        pvalg = mannwhitneyu(
            1 - scores[0], 1 - scores[1], alternative="greater"
        )[1]
        pvall = mannwhitneyu(1 - scores[0], 1 - scores[1], alternative="less")[
            1
        ]

    for s, l in zip(scores, labels):
        if outfield in ["latent_time", "velocity_pseudotime"]:
            s = np.sort(s)
        else:
            if outfield != "cytotrace":
                s = np.sort(-s)
            else:
                s = np.sort(1 - s)
        ax.step(
            np.concatenate([s, s[[-1]]]),
            np.arange(s.size + 1) / s.size,
            label=l,
        )
        ax.scatter(
            np.concatenate([s, s[[-1]]]), np.arange(s.size + 1) / s.size, s=5
        )
        if outfield == "cytotrace":
            ax.set_xlabel("Differentiation level")
        else:
            ax.set_xlabel(outfield)
        ax.set_ylabel("Cumulative proportion of cells")
        ax.set_title(type)
    ax.text(
        np.percentile(s, 0.6),
        0.4,
        f"{labels[0]} > {labels[1]}: {pvalg:.3E}",
    )
    ax.text(
        np.percentile(s, 0.6),
        0.48,
        f"{labels[0]} < {labels[1]}: {pvall:.3E}",
    )
    ax.legend()


def plot_multiadata(adatas):
    for a in adatas:
        print(a.obs.age[0])
        fig, ax = plt.subplots(1, 5)
        fig.set_size_inches(32, 9)
        scv.pl.velocity_embedding_grid(
            a,
            scale=0.2,
            color="louvain",
            show=False,
            basis="pca",
            legend_loc="on data",
            ax=ax[0],
            title=a.obs.age[0],
        )
        scv.pl.scatter(
            a,
            color="velocity_pseudotime",
            basis="pca",
            size=80,
            color_map="gnuplot",
            ax=ax[1],
            show=False,
            title=a.obs.age[0],
        )
        scv.pl.scatter(
            a,
            color="latent_time",
            basis="pca",
            size=80,
            color_map="gnuplot",
            ax=ax[2],
            show=False,
            title=a.obs.age[0],
        )
        scv.pl.scatter(
            a,
            color="end_points",
            basis="pca",
            size=80,
            color_map="gnuplot",
            ax=ax[3],
            show=False,
            title=a.obs.age[0],
        )
        scv.pl.scatter(
            a,
            color="root_cells",
            basis="pca",
            size=80,
            color_map="gnuplot",
            ax=ax[4],
            show=False,
            title=a.obs.age[0],
        )


# https://stackoverflow.com/questions/15408371/cumulative-distribution-plots-python
def cumulative_boxplot():
    data = adata_copy[adata_copy.obs.age == "A",].obs.n_counts
    # evaluate the histogram
    values, base = np.histogram(data, bins=10)
    # evaluate the cumulative
    cumulative = np.cumsum(values)
    # plot the cumulative function
    plt.scatter(base[:-1], cumulative / cumulative[-1], c="blue")
    # plot the survival function
    # plt.plot(base[:-1], len(data)-cumulative, c='green')
