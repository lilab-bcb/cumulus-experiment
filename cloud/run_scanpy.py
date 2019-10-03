import os, sys, time, math
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from datetime import timedelta

sc.settings.n_jobs = 32
sc.settings.verbosity = 4
rand_seed = 0

input_file = "/projects/benchmark/MantonBM/MantonBM_nonmix_10x.h5"
output_file = "scanpy_analysis_result.h5ad"

print("Reading ICA (bone marrow) dataset")
adata = sc.read_10x_h5(input_file, genome = 'GRCh38')
adata.var_names_make_unique()
adata.obs['Channel'] = adata.obs.index.map(lambda s: s.split('-')[0]).astype('category')

start_overall = time.time()

sc.pp.filter_cells(adata, min_genes = 500)
sc.pp.filter_cells(adata, max_genes = 6000)

num_cells = adata.shape[0]
sc.pp.filter_genes(adata, min_cells = math.floor(num_cells * 0.0005))

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
adata.obs['n_counts'] = adata.X.sum(axis = 1).A1

adata = adata[adata.obs['percent_mito'] < 0.1, :]

sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e5)
sc.pp.log1p(adata)

adata.raw = adata.copy()

start_hvg = time.time()
sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key = 'Channel')
end_hvg = time.time()
print("Time spent for HVG: {}.".format(timedelta(seconds = end_hvg - start_hvg)))

adata = adata[:, adata.var['highly_variable']]

start_pca = time.time()
sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata, n_comps = 50, random_state = rand_seed)
end_pca = time.time()
print("Time spent for PCA: {}.".format(timedelta(seconds = end_pca - start_pca)))

start_bbknn = time.time()
bbknn(adata, batch_key = 'Channel')
end_bbknn = time.time()
print("Time spent for BBKNN: {}.".format(timedelta(seconds = end_bbknn - start_bbknn)))

start_louvain = time.time()
sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)
end_louvain = time.time()
print("Time spent for Louvain: {}.".format(timedelta(seconds = end_louvain - start_louvain)))

start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
print("Time spent for TSNE: {}.".format(timedelta(seconds = end_tsne - start_tsne)))

start_umap = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_umap = time.time()
print("Time spent for UMAP: {}.".format(timedelta(seconds = end_umap - start_umap)))

start_de = time.time()
sc.tl.rank_genes_groups(adata, groupby = 'louvain', method = 't-test')
end_de = time.time()
print("Time spent for DE analysis: {}.".format(timedelta(seconds = end_de - start_de)))

end_overall = time.time()
print("Overall Time: {}.".format(timedelta(seconds = end_overall - start_overall)))

adata.write(output_file, compression = 'gzip')