import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
from scanpy.neighbors import compute_connectivities_umap
from datetime import timedelta
from termcolor import cprint

sc.settings.n_jobs = os.cpu_count()
sc.settings.verbosity = 4
rand_seed = 0

filter_norm_data = "../MantonBM_nonmix_10x_filter_norm.h5ad"
hvg_file = "../hvg.txt"
pca_data = "../MantonBM_nonmix_10x_filter_norm_pca.h5ad"
corrected_data = "../MantonBM_nonmix_corrected.h5ad"

f = open("scanpy.log", "w")

print("Reading ICA (bone marrow) dataset")
adata = sc.read_h5ad(filter_norm_data)
print("Finding highly variable genes")
start_hvg = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel')
end_hvg = time.time()
logstr_hvg = "Time spent for HVG: {}.".format(timedelta(seconds = end_hvg - start_hvg))
cprint(logstr_hvg, "yellow")
f.write(logstr_hvg + '\n')


print("Batch correction using BBKNN")
df_hvg = pd.read_csv(hvg_file)
df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
df_hvg.set_index('index', inplace = True)
df_hvg.drop(columns = ['gene_ids'], inplace = True)
adata = sc.read_h5ad(pca_data)
df_join = adata.var.join(df_hvg)
df_join['highly_variable_genes'].fillna(False, inplace = True)
adata_hvg = adata[:, df_join['highly_variable_genes']].copy()
start_bbknn = time.time()
bbknn(adata_hvg, batch_key = 'Channel')
end_bbknn = time.time()
logstr_batch = "Time spend for BBKNN batch correction: {}.".format(timedelta(seconds = end_bbknn - start_bbknn))
cprint(logstr_batch, "yellow")
f.write(logstr_batch + '\n')


print("Running PCA")
adata = sc.read_h5ad(corrected_data)
adata_hvg = adata[:, adata.var['highly_variable_genes']].copy()
start_pca = time.time()
sc.pp.scale(adata_hvg, max_value = 10)
sc.tl.pca(adata_hvg, random_state = rand_seed)
end_pca = time.time()
logstr_pca = "Time spent for PCA: {}.".format(timedelta(seconds = end_pca - start_pca))
cprint(logstr_pca, "yellow")
f.write(logstr_pca + '\n')


print("Computing neighborhood graph")
adata = sc.read_h5ad(corrected_data)
start_nn = time.time()
sc.pp.neighbors(adata, n_neighbors = 100, random_state = rand_seed)
end_nn = time.time()
logstr_nn = "Time spent for knn: {}.".format(timedelta(seconds = end_nn - start_nn))
cprint(logstr_nn, "yellow")
f.write(logstr_nn + '\n')


# Construct affinity matrix using KNN information from scCloud results.
adata = sc.read_h5ad(corrected_data)
n_cells = adata.shape[0]
knn_indices = np.concatenate((np.arange(n_cells).reshape(-1, 1), adata.uns['knn_indices']), axis = 1)
knn_distances = np.concatenate((np.zeros(n_cells).reshape(-1, 1), adata.uns['knn_distances']), axis = 1)
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
logstr_louvain = "Time spent for louvain: {}.".format(timedelta(seconds = end_louvain - start_louvain))
cprint(logstr_louvain, "yellow")
f.write(logstr_louvain + '\n')

print("Finding Clusters using Leiden")
start_leiden = time.time()
sc.tl.leiden(adata, resolution = 1.3, random_state = rand_seed)
end_leiden = time.time()
logstr_leiden = "Time spent for leiden: {}.".format(timedelta(seconds = end_leiden - start_leiden))
cprint(logstr_leiden, "yellow")
f.write(logstr_leiden + '\n')


print("Computing UMAP embedding")
start_umap = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_umap = time.time()
logstr_umap = "Time spent for UMAP: {}.".format(timedelta(seconds = end_umap - start_umap))
cprint(logstr_umap, "yellow")
f.write(logstr_umap + '\n')


print("Computing tSNE embedding")
start_tsne = time.time()
sc.tl.tsne(adata, n_pcs = 2, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
logstr_tsne = "Time spent for tSNE: {}.".format(timedelta(seconds = end_tsne - start_tsne))
cprint(logstr_tsne, "yellow")
f.write(logstr_tsne + '\n')


print("Computing diffmap")
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 50)
end_df = time.time()
logstr_df = "Time spent for diffmap: {}.".format(timedelta(seconds = end_df - start_df))
cprint(logstr_df, "yellow")
f.write(logstr_df + '\n')

f.close()

print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, n_jobs = sc.settings.n_jobs)
end_fle = time.time()
logstr_fle = "Time spent for FLE: {}.".format(timedelta(seconds = end_fle - start_fle))
cprint(logstr_fle, "yellow")