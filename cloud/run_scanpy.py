import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from datetime import timedelta

sc.settings.n_jobs = 64
sc.settings.verbosity = 4
rand_seed = 0

input_file = "/projects/MantonBM/MantonBM_nonmix_10x.h5"
output_file = "scanpy_analysis_result.h5ad"

print("Reading ICA (bone marrow) dataset")
adata = sc.read_10x_h5(input_file, genome = 'GRCh38')
adata.var_names_make_unique()

num_cells = adata.shape[0]
sc.pp.filter_cells(adata, min_counts = 100, min_genes = 500, max_counts = 6e5, max_genes = 6000)
sc.pp.filter_genes(adata, min_cells = num_cells * 0.0005)

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
adata.obs['n_counts'] = adata.X.sum(axis = 1).A1

adata = adata[adata.obs['percent_mito'] < 0.1, :]

sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e5)
sc.pp.log1p(adata)

adata.raw = adata.copy()

sc.pp.highly_variable_genes(adata, n_top_genes = 2000)
adata = adata[:, adata.var['highly_variable']]

sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata, n_comps = 50, random_state = rand_seed)

bbknn(adata, batch_key = 'Channel')

sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)

sc.tl.tsne(adata, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)

adata.uns['neighbors']['params']['n_neighbors'] = 15
sc.tl.umap(adata, random_state = rand_seed)

sc.tl.rank_genes_groups(adata, groupby = 'louvain', method = 't-test')

adata.write(output_file, compression = 'gzip')