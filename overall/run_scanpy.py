#!/usr/bin/env python

import numpy as np
import scanpy as sc
import os, sys, time

sc.settings.n_jobs = os.cpu_count()
rand_seed = 0
output_dir = "scanpy_output"

if output_dir not in os.listdir('.'):
	if os.system("mkdir {}".format(output_dir)):
		sys.exit(1)

f = open(output_dir + "/scanpy.log", "w")

print("Reading ICA (bone marrow) dataset")
start = time.time()
adata = sc.read_10x_h5("../MantonBM_nonmix_10x.h5", genome='GRCh38')
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

print("Finding variable genes")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5)
adata_variable_genes = adata[:, adata.var['highly_variable']]

print("Scale expression matrix of variable genes")
sc.pp.scale(adata_variable_genes, max_value=10)

print("Running PCA")
start_pca = time.time()
sc.tl.pca(adata_variable_genes, random_state = rand_seed)
end_pca = time.time()
adata.obsm = adata_variable_genes.obsm
print("Time spent for PCA = " + str(end_pca - start_pca) + " seconds.")
f.write("Time spent for PCA = " + str(end_pca - start_pca) + " seconds.\n")

print("Computing neighborhood graph")
start_nn = time.time()
sc.pp.neighbors(adata, n_neighbors=100)
end_nn = time.time()
print("Time spent for knn = " + str(end_nn - start_nn) + " seconds.")
f.write("Time spent for knn = " + str(end_nn - start_nn) + " seconds.\n")

print("Finding Clusters")
start_cluster = time.time()
sc.tl.leiden(adata, resolution = 1.3, random_state = rand_seed)
end_cluster = time.time()
print("Time spent for leiden = " + str(end_lv - start_lv) + " seconds.")
f.write("Time spent for leiden = " + str(end_lv - start_lv) + " seconds.\n")

print("Computing UMAP embedding")
start_up = time.time()
sc.tl.umap(adata, random_state = rand_seed)
end_up = time.time()
print("Time spent for UMAP = " + str(end_up - start_up) + " seconds.")
f.write("Time spent for UMAP = " + str(end_up - start_up) + "seconds.\n")

print("Computing tSNE embedding")
start_te = time.time()
sc.tl.tsne(adata, n_pcs = 2, random_state = rand_seed, use_fast_tsne = True, n_jobs = sc.settings.n_jobs)
end_te = time.time()
print("Time spent for tSNE = " + str(end_te - start_te) + " seconds.")
f.write("Time spent for tSNE = " + str(end_te - start_te) + " seconds.\n")

print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, n_jobs = sc.settings.n_jobs)
end_fle = time.time()
print("Time spent for FLE = " + str(end_fle - start_fle) + " seconds.")
f.write("Time spent for FLE = " + str(end_fle - start_fle) + " seconds.\n")

print("Computing diffmap")
start_df = time.time()
sc.tl.diffmap(adata, n_comps = 50)
end_df = time.time()
print("Time spent for diffmap = " + str(end_df - start_df) + " seconds.")
f.write("Time spent for diffmap = " + str(end_df - start_df) + " seconds.\n")

adata.write(output_dir + "/MantonBM_nonmix_scanpy.h5ad")
end = time.time()
print("Total time = " + str(end - start) + " seconds.")
f.write("Total time = " + str(end - start) + " seconds.")
f.close()

#sc.pl.tsne(adata, color = ['leiden'], save = '.pdf')
#sc.pl.umap(adata, color = ['leiden'], save = '.pdf')
#sc.pl.draw_graph(adata, color = 'leiden', save = '.pdf')