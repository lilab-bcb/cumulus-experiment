#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import os, sys, re
import time
from termcolor import cprint
from mnnpy import mnn_correct
from utils import str_time

sc.settings.n_jobs = 12
rand_seed = 0

f = open("./scanpy_mnn.log", "w")

batch_folder = "/Users/yy939/Documents/manton-bm/raw_batch/tiny"
f_list = [f for f in os.listdir(batch_folder) if re.match('.+_HiSeq_1.h5', f)]
ad_list = []
cnt = 0
for f_name in f_list:
	cnt += 1
	cprint("{seq} Processing {file_name}...".format(seq = cnt, file_name = f_name), "yellow")
	
	adata = sc.read_10x_h5(batch_folder + '/' + f_name, genome='GRCh38')
	n_sample = adata.shape[0]
	batch_name = f_name.split('.')[0]
	adata.obs['batch'] = [batch_name] * n_sample
	adata.var_names_make_unique()
	ad_list.append(adata)


print("Getting Highly Variable Sets:")
hvg_list = list(pd.read_csv("../../hvg.txt", header = None)[0].values)
print(hvg_list[0:5])

print("Batch Correction Using MNN:")
start_mnn = time.time()
corrected = mnn_correct(*ad_list, var_subset = hvg_list, n_jobs = 8)
end_mnn = time.time()
print("Time spent = " + str_time(end_mnn - start_mnn))
f.write("Time spent for MNN = " + str_time(end_mnn - start_mnn) + '\n')

adata = corrected[0]
adata.write("MantonBM_nonmix_scanpy_mnn.h5ad")
print("Corrected data are saved.")

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

# Should be replaced by hvg
print("Finding variable genes")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, inplace=True)
adata_variable_genes = adata[:, adata.var['highly_variable']]


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
sc.tl.umap(adata, min_dist = 0.1, random_state = rand_seed)
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

adata.write("MantonBM_nonmix_scanpy_mnn_processed.h5ad")
end = time.time()
print("Total time = " + str_time(end - start) + ".")
f.write("Total time = " + str_time(end - start) + ".\n")

f.close()