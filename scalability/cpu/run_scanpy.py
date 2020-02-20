import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from scanpy.neighbors import compute_connectivities_umap
from umap.umap_ import nearest_neighbors
from numpy.random import RandomState
from datetime import timedelta

sc.settings.n_jobs = int(sys.argv[1])
sc.settings.verbosity = 4
rand_seed = 0
print("SCANPY experiment using {} cores:".format(sc.settings.n_jobs))

filter_norm_data = "/data/MantonBM_nonmix_filter_norm.h5ad"
hvf_file = "/data/MantonBM_nonmix_hvf.txt"
pca_data = "/data/MantonBM_nonmix_filter_norm_pca.h5ad"
corrected_data = "/data/MantonBM_nonmix_corrected.h5ad"

adata = sc.read_h5ad(filter_norm_data)
start_hvg = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel')
end_hvg = time.time()
logstr_hvg = "Time spent for HVG: {}.".format(timedelta(seconds = end_hvg - start_hvg))
print(logstr_hvg)
del adata

adata = sc.read_h5ad(corrected_data)
adata_hvf = adata[:, adata.var['highly_variable_features']].copy()
start_pca = time.time()
sc.pp.scale(adata_hvf, max_value = 10)
sc.tl.pca(adata_hvf, random_state = rand_seed)
end_pca = time.time()
logstr_pca = "Time spent for PCA: {}.".format(timedelta(seconds = end_pca - start_pca))
print(logstr_pca)
del adata

adata = sc.read_h5ad(corrected_data)
rs = RandomState(rand_seed)
start_knn = time.time()
knn_indices, knn_distances, _ = nearest_neighbors(adata.obsm['X_pca'], n_neighbors = 100, random_state = rs, metric = 'euclidean', metric_kwds = {}, angular = False, verbose = False)
end_knn = time.time()
logstr_knn = "Time spent for KNN: {}.".format(timedelta(seconds = end_knn - start_knn))
print(logstr_knn)
del adata

# Construct affinity matrix using KNN information from Pegasus results, but reduce to 15 neighbors.
adata = sc.read_h5ad(corrected_data)
n_cells = adata.shape[0]
n_neighbors = 15
knn_indices = np.concatenate((np.arange(n_cells).reshape(-1, 1), adata.uns['pca_knn_indices'][:, 0:(n_neighbors - 1)]), axis = 1)
knn_distances = np.concatenate((np.zeros(n_cells).reshape(-1, 1), adata.uns['pca_knn_distances'][:, 0:(n_neighbors - 1)]), axis = 1)
distances, connectivities = compute_connectivities_umap(knn_indices, knn_distances, n_obs = n_cells, n_neighbors = n_neighbors)
adata.uns['neighbors'] = {}
adata.uns['neighbors']['indices'] = knn_indices
adata.uns['neighbors']['distances'] = distances
adata.uns['neighbors']['connectivities'] = connectivities
adata.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': 'umap', 'metric': 'euclidean'}

start_umap = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_umap = time.time()
logstr_umap = "Time spent for UMAP: {}.".format(timedelta(seconds = end_umap - start_umap))
print(logstr_umap)
del adata

# Use KNN from Pegasus
adata = sc.read_h5ad(corrected_data)
n_cells = adata.shape[0]
knn_indices = np.concatenate((np.arange(n_cells).reshape(-1, 1), adata.uns['pca_knn_indices']), axis = 1)
knn_distances = np.concatenate((np.zeros(n_cells).reshape(-1, 1), adata.uns['pca_knn_distances']), axis = 1)
distances, connectivities = compute_connectivities_umap(knn_indices, knn_distances, n_obs = n_cells, n_neighbors = 100)
adata.uns['neighbors'] = {}
adata.uns['neighbors']['indices'] = knn_indices
adata.uns['neighbors']['distances'] = distances
adata.uns['neighbors']['connectivities'] = connectivities
adata.uns['neighbors']['params'] = {'n_neighbors': 100, 'method': 'umap', 'metric': 'euclidean'}

start_louvain = time.time()
sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)
end_louvain = time.time()
logstr_louvain = "Time spent for louvain: {}.".format(timedelta(seconds = end_louvain - start_louvain))
print(logstr_louvain)

start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
logstr_tsne = "Time spent for tSNE: {}.".format(timedelta(seconds = end_tsne - start_tsne))
print(logstr_tsne)

#adata.write("scanpy_result.h5ad", compression = 'gzip')