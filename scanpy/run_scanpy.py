#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
from bbknn import bbknn
from utils import str_time

sc.settings.n_jobs = 8
rand_seed = 0
hvg_file = "../hvg.txt"
f = open("./scanpy_local.log", "w")

print("Reading ICA (bone marrow) dataset")
start = time.time()
adata = sc.read_10x_h5("../MantonBM_nonmix_tiny_10x.h5", genome='GRCh38')
adata.obs['channel'] = list(map(lambda s : s.split('-')[0], adata.obs.index.values))
adata.var_names_make_unique()

print("Filtering cells")
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, max_genes=5999)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata = adata[adata.obs['percent_mito'] < 0.1, :]

print("Filtering genes")
sc.pp.filter_genes(adata, min_cells=138)

print("Normalizing data")
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e5)
sc.pp.log1p(adata)


print("Getting Highly Variable Genes from scCloud")
df_hvg = pd.read_csv(hvg_file)
print("Highly variable gene set of size " + str(df_hvg.shape[0]))
df_hvg.rename(columns = {'index': 'name'}, inplace = True)
df_hvg.set_index('gene_ids', inplace = True)
df_hvg.drop(columns = ['name'], inplace = True)
df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
df_var = adata.var.reset_index().set_index('gene_ids').rename(columns = {'index': 'name'})
df_join = df_var.join(df_hvg)
df_join['highly_variable_genes'].fillna(False, inplace = True)
adata_variable_genes = adata[:, df_join['highly_variable_genes']]

#print("Finding variable genes")
#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, inplace=True)
#adata_variable_genes = adata[:, adata.var['highly_variable']]


print("Scale expression matrix of variable genes")
sc.pp.scale(adata_variable_genes, max_value=10)

print("Running PCA")
start_pca = time.time()
sc.tl.pca(adata_variable_genes)
end_pca = time.time()
adata.obsm = adata_variable_genes.obsm
print("Time spent for PCA = " + str_time(end_pca - start_pca) + ".")
f.write("Time spent for PCA = " + str_time(end_pca - start_pca) + ".\n")

print("Computing neighborhood graph")
start_nn = time.time()
sc.pp.neighbors(adata, n_neighbors=100)
end_nn = time.time()
print("Time spent for knn = " + str_time(end_nn - start_nn) + ".")
f.write("Time spent for knn = " + str_time(end_nn - start_nn) + ".\n")

print("Finding Clusters using Louvain algorithm")
start_lv = time.time()
sc.tl.louvain(adata, resolution = 1.3)
end_lv = time.time()
print("Time spent for louvain = " + str_time(end_lv - start_lv) + ".")
f.write("Time spent for louvain = " + str_time(end_lv - start_lv) + ".\n")

print("Finding Clusters using Leiden algorithm")
start_ld = time.time()
sc.tl.leiden(adata, resolution = 1.3)
end_ld = time.time()
print("Time spent for leiden = " + str_time(end_ld - start_ld) + ".")
f.write("Time spent for leiden = " + str_time(end_ld - start_ld) + ".\n")

print("Computing UMAP embedding")
start_up = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_up = time.time()
print("Time spent for UMAP = " + str_time(end_up - start_up) + ".")
f.write("Time spent for UMAP = " + str_time(end_up - start_up) + ".\n")

print("Computing tSNE embedding")
start_te = time.time()
sc.tl.tsne(adata, n_pcs = 2, random_state = rand_seed, use_fast_tsne = True)
end_te = time.time()
print("Time spent for tSNE = " + str_time(end_te - start_te) + ".")
f.write("Time spent for tSNE = " + str_time(end_te - start_te) + ".\n")

print("Computing diffmap")
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 50)
end_df = time.time()
print("Time spent for diffmap = " + str_time(end_df - start_df) + ".")
f.write("Time spent for diffmap = " + str_time(end_df - start_df) + ".\n")

adata.write("MantonBM_nonmix_tiny_scanpy.h5ad")
end = time.time()
print("Total time = " + str_time(end - start) + ".")
f.write("Total time = " + str_time(end - start) + ".\n")
f.close()

sc.pl.tsne(adata, color = ['leiden'], save = '.png')
sc.pl.umap(adata, color = ['leiden'], save = '.png')