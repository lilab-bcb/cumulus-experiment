import time
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from scanpy.neighbors import compute_connectivities_umap
from umap.umap_ import nearest_neighbors
from numpy.random import RandomState
from datetime import timedelta

def record_time(duration_time, step_name, fp):
    fp.write("Time spent for {step}: {time}.\n\n".format(step = step_name, time = timedelta(seconds = duration_time)))

sc.settings.n_jobs = 28
sc.settings.verbosity = 4
rand_seed = 0

log_file = "1m_scanpy_cpu_{}.log".format(sc.settings.n_jobs)
fp = open(log_file, 'w')

filter_norm_data = "/data/1M_neurons_filter_norm.h5ad"
hvf_file = "/data/1M_neurons_hvf.txt"
pca_data = "/data/1M_neurons_filter_norm_pca.h5ad"
corrected_data = "/data/1M_neurons_corrected.h5ad"

print("Reading Mouse Neuron dataset")
adata = sc.read_h5ad(filter_norm_data)
print("Finding highly variable genes")
start_hvf = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel')
end_hvf = time.time()
record_time(end_hvf - start_hvf, "HVG", fp)

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
record_time(end_bbknn - start_bbknn, "BBKNN", fp)

print("Running PCA")
adata = sc.read_h5ad(corrected_data)
adata_hvf = adata[:, adata.var['highly_variable_features']].copy()
start_pca = time.time()
sc.pp.scale(adata_hvf, max_value = 10)
sc.tl.pca(adata_hvf, random_state = rand_seed)
end_pca = time.time()
record_time(end_pca - start_pca, "PCA", fp)

print("Computing neighborhood graph")
adata = sc.read_h5ad(corrected_data)
rs = RandomState(rand_seed)
start_knn = time.time()
knn_indices, knn_distances, _ = nearest_neighbors(adata.obsm['X_pca'], n_neighbors = 100, random_state = rs, metric = 'euclidean', metric_kwds = {}, angular = False, verbose = False)
end_knn = time.time()
record_time(end_knn - start_knn, "KNN", fp)

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
record_time(end_umap - start_umap, "UMAP", fp)

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
record_time(end_louvain - start_louvain, "Louvain", fp)

print("Finding Clusters using Leiden")
start_leiden = time.time()
sc.tl.leiden(adata, resolution = 1.3, random_state = rand_seed)
end_leiden = time.time()
record_time(end_leiden - start_leiden, "Leiden", fp)

print("Computing tSNE embedding")
start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
record_time(end_tsne - start_tsne, "tSNE", fp)

print("Computing diffmap")
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 100)
end_df = time.time()
record_time(end_df - start_df, "Diffmap", fp)

adata.write("1M_neurons_scanpy_result_no_fle.h5ad", compression = 'gzip')

fp.close()