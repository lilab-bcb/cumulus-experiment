import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from scanpy.neighbors import compute_connectivities_umap
from umap.umap_ import nearest_neighbors
from numpy.random import RandomState
from datetime import timedelta, datetime

sc.settings.n_jobs = os.cpu_count()
sc.settings.verbosity = 4
rand_seed = 0

filter_norm_data = "/data/MantonBM_nonmix_filter_norm.h5ad"
hvf_file = "/data/MantonBM_nonmix_hvf.txt"
pca_data = "/data/MantonBM_nonmix_filter_norm_pca.h5ad"
corrected_data = "/data/MantonBM_nonmix_corrected.h5ad"

print("Reading ICA (bone marrow) dataset")
adata = sc.read_h5ad(filter_norm_data)
print("Finding highly variable genes")
start_hvg = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel')
end_hvg = time.time()
print("HVG starts at: {}".format(datetime.fromtimestamp(start_hvg)))
print("HVG ends at: {}".format(datetime.fromtimestamp(end_hvg)))
logstr_hvg = "Time spent for HVG: {}.".format(timedelta(seconds = end_hvg - start_hvg))
print(logstr_hvg)


print("Batch correction using BBKNN")
df_hvf = pd.read_csv(hvf_file)
df_hvf['highly_variable_features'] = [True] * df_hvf.shape[0]
df_hvf.set_index('index', inplace = True)
df_hvf.drop(columns = ['gene_ids'], inplace = True)
adata = sc.read_h5ad(pca_data)
adata.var.drop(columns = ['highly_variable_features'], inplace = True)
df_join = adata.var.join(df_hvf)
df_join['highly_variable_features'].fillna(False, inplace = True)
adata_hvf = adata[:, df_join['highly_variable_features']].copy()
start_bbknn = time.time()
bbknn(adata_hvf, batch_key = 'Channel')
end_bbknn = time.time()
print("BBKNN starts at: {}".format(datetime.fromtimestamp(start_bbknn)))
print("BBKNN ends at: {}".format(datetime.fromtimestamp(end_bbknn)))
logstr_batch = "Time spend for BBKNN batch correction: {}.".format(timedelta(seconds = end_bbknn - start_bbknn))
print(logstr_batch)


print("Running PCA")
adata = sc.read_h5ad(corrected_data)
adata_hvf = adata[:, adata.var['highly_variable_features']].copy()
start_pca = time.time()
sc.pp.scale(adata_hvf, max_value = 10)
sc.tl.pca(adata_hvf, random_state = rand_seed)
end_pca = time.time()
print("PCA starts at: {}".format(datetime.fromtimestamp(start_pca)))
print("PCA ends at: {}".format(datetime.fromtimestamp(end_pca)))
logstr_pca = "Time spent for PCA: {}.".format(timedelta(seconds = end_pca - start_pca))
print(logstr_pca)

print("Computing KNN")
adata = sc.read_h5ad(corrected_data)
rs = RandomState(rand_seed)
start_knn = time.time()
knn_indices, knn_distances, _ = nearest_neighbors(adata.obsm['X_pca'], n_neighbors = 100, random_state = rs, metric = 'euclidean', metric_kwds = {}, angular = False, verbose = False)
end_knn = time.time()
print("KNN starts at: {}".format(datetime.fromtimestamp(start_knn)))
print("KNN ends at: {}".format(datetime.fromtimestamp(end_knn)))
logstr_knn = "Time spent for KNN: {}.".format(timedelta(seconds = end_knn - start_knn))
print(logstr_knn)

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

print("Computing UMAP embedding")
start_umap = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_umap = time.time()
print("UMAP starts at: {}".format(datetime.fromtimestamp(start_umap)))
print("UMAP ends at: {}".format(datetime.fromtimestamp(end_umap)))
logstr_umap = "Time spent for UMAP: {}.".format(timedelta(seconds = end_umap - start_umap))
print(logstr_umap)

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

print("Finding Clusters using Louvain")
start_louvain = time.time()
sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)
end_louvain = time.time()
print("Louvain starts at: {}".format(datetime.fromtimestamp(start_louvain)))
print("Louvain ends at: {}".format(datetime.fromtimestamp(end_louvain)))
logstr_louvain = "Time spent for louvain: {}.".format(timedelta(seconds = end_louvain - start_louvain))
print(logstr_louvain)

print("Finding Clusters using Leiden")
start_leiden = time.time()
sc.tl.leiden(adata, resolution = 1.3, random_state = rand_seed)
end_leiden = time.time()
print("Leiden starts at: {}".format(datetime.fromtimestamp(start_leiden)))
print("Leiden ends at: {}".format(datetime.fromtimestamp(end_leiden)))
logstr_leiden = "Time spent for leiden: {}.".format(timedelta(seconds = end_leiden - start_leiden))
print(logstr_leiden)

print("Computing tSNE embedding")
start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
print("tSNE starts at: {}".format(datetime.fromtimestamp(start_tsne)))
print("tSNE ends at: {}".format(datetime.fromtimestamp(end_tsne)))
logstr_tsne = "Time spent for tSNE: {}.".format(timedelta(seconds = end_tsne - start_tsne))
print(logstr_tsne)


print("Computing diffmap")
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 100)
end_df = time.time()
print("Diffmap starts at: {}".format(datetime.fromtimestamp(start_df)))
print("Diffmap ends at: {}".format(datetime.fromtimestamp(end_df)))
logstr_df = "Time spent for diffmap: {}.".format(timedelta(seconds = end_df - start_df))
print(logstr_df)

adata.write("scanpy_mantonbm_result_no_fle.h5ad", compression = 'gzip')


adata = sc.read_h5ad("scanpy_mantonbm_result_no_fle.h5ad")
print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, iterations = 5000)
end_fle = time.time()
print("FLE starts at: {}".format(datetime.fromtimestamp(start_fle)))
print("FLE ends at: {}".format(datetime.fromtimestamp(end_fle)))
logstr_fle = "Time spent for FLE: {}.".format(timedelta(seconds = end_fle - start_fle))
print(logstr_fle)

adata.write("scanpy_result.h5ad", compression = 'gzip')