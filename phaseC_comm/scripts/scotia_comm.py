import os
import scotia
from pathlib import Path
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import sys


def main():
    datap = '../phaseB_check_rds/notebook/kidney_T_miles.h5ad'
    ##adata_test = sc.read_h5ad('../../SCOTIA/example/merfish_liver_example.h5ad')
    adata_test = sc.read_h5ad(datap)
    known_lr_pairs = pd.read_csv("../../SCOTIA/example/lr_gene.list", sep = '\t', header = None, index_col = None)
    known_lr_pairs.columns = ['l_gene','r_gene']
    
    known_lr_pairs = known_lr_pairs.loc[known_lr_pairs.l_gene.str.upper().isin(adata_test.var.index), :]
    print(known_lr_pairs)
    #known_lr_pairs = known_lr_pairs.loc[:, known_lr_pairs.r_gene.str.upper().isin(adata_test.var.index)]
    print(adata_test.var.index)


    print(known_lr_pairs.head())
    print(adata_test.obs.head())
    print(adata_test.obs.columns)
    adata_test.obs['sample'] = 'A'
    adata_test.obs['fov'] = 9
    print(adata_test.obs['fov'].value_counts())
    print(adata_test.obs['sample'].value_counts())

    folder_path = Path("test_output/")
    folder_path.mkdir(parents=True, exist_ok=True)

    scotia.run_scotia.lr_score(adata = adata_test,
         lr_list = known_lr_pairs,
         sample_col = 'sample',
         fov_col = 'fov',
         #celltype_col = 'annotation',
         celltype_col = 'celltype_hint',
         output_path = 'test_output/')

if __name__ == '__main__':
    main()
