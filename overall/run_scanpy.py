import numpy as np
import pandas as pd
import scanpy as sc
from bbknn import bbknn
import os, sys, time

sc.settings.n_jobs = os.cpu_count()
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
print("Time spent for HVG = {:.2f} seconds.".format(end_hvg - start_hvg))
f.write("Time spent for HVG = {:.2f} seconds.\n".format(end_hvg - start_hvg))

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
print("Time spend for BBKNN batch correction = {:.2f} seconds.".format(end_bbknn - start_bbknn))
f.write("Time spend for BBKNN batch correction = {:.2f} seconds.\n".format(end_bbknn - start_bbknn))

print("Running PCA")
adata = sc.read_h5ad(corrected_data)
adata_hvg = adata[:, adata.var['highly_variable_genes']].copy()
start_pca = time.time()
sc.pp.scale(adata_hvg, max_value = 10)
sc.tl.pca(adata_hvg, random_state = rand_seed)
end_pca = time.time()
print("Time spent for PCA = {:.2f} seconds.".format(end_pca - start_pca))
f.write("Time spent for PCA = {:.2f} seconds.\n".format(end_pca - start_pca))



print("Computing neighborhood graph")
adata = sc.read_h5ad(corrected_data)
start_nn = time.time()
sc.pp.neighbors(adata, n_neighbors = 100, random_state = rand_seed)
end_nn = time.time()
print("Time spent for knn = {:.2f} seconds.".format(end_nn - start_nn))
f.write("Time spent for knn = {:.2f} seconds.\n".format(end_nn - start_nn))


print("Finding Clusters")
start_cluster = time.time()
sc.tl.leiden(adata, resolution = 1.3, random_state = rand_seed)
end_cluster = time.time()
print("Time spent for leiden = {:.2f} seconds.".format(end_cluster - start_cluster))
f.write("Time spent for leiden = {:.2f} seconds.\n".format(end_cluster - start_cluster))

print("Computing UMAP embedding")
start_umap = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_umap = time.time()
print("Time spent for UMAP = {:.2f} seconds.".format(end_umap - start_umap))
f.write("Time spent for UMAP = {:.2f} seconds.\n".format(end_umap - start_umap))

print("Computing tSNE embedding")
start_tsne = time.time()
sc.tl.tsne(adata, n_pcs = 2, use_rep = 'X_pca', random_state = rand_seed, use_fast_tsne = True)
end_tsne = time.time()
print("Time spent for tSNE = {:.2f} seconds.".format(end_tsne - start_tsne))
f.write("Time spent for tSNE = {:.2f} seconds.\n".format(end_tsne - start_tsne))

print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, n_jobs = sc.settings.n_jobs)
end_fle = time.time()
print("Time spent for FLE = {:.2f} seconds.".format(end_fle - start_fle))
f.write("Time spent for FLE = {:.2f} seconds.\n".format(end_fle - start_fle))

print("Computing diffmap")
adata = sc.read_h5ad(corrected_data)
sc.pp.neighbors(adata, n_neighbors = 100, random_state = rand_seed, method = 'gauss')
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 50)
end_df = time.time()
print("Time spent for diffmap = {:.2f} seconds.".format(end_df - start_df))
f.write("Time spent for diffmap = {:.2f} seconds.\n".format(end_df - start_df))


f.close()

#sc.pl.tsne(adata, color = ['leiden'], save = '.pdf')
#sc.pl.umap(adata, color = ['leiden'], save = '.pdf')
#sc.pl.draw_graph(adata, color = 'leiden', save = '.pdf')