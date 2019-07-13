import numpy as np
import pandas as pd
import os, sys, time
import scanpy as sc
from bbknn import bbknn

sc.settings.n_jobs = os.cpu_count()
rand_seed = 0
filter_norm_data = "../MantonBM_nonmix_10x_filter_norm.h5ad"
hvg_file = "../hvg.txt"
pca_data = "../MantonBM_nonmix_10x_filter_norm_pca.h5ad"
corrected_data = "./sccloud_output/MantonBM_nonmix_sccloud_corrected.h5ad"

f = open(output_dir + "scanpy.log", "w")

print("Reading ICA (bone marrow) dataset")
start = time.time()
adata_filter_norm = sc.read_h5ad(filter_norm_data)
print("Finding variable genes")
start_hvg = time.time()
sc.pp.highly_variable_genes(adata_filter_norm, min_mean=0.0125, max_mean=7, min_disp=0.5)
end_hvg = time.time()
print("Time spent for HVG = {:.2f} seconds.".format(end_hvg - start_hvg))
f.write("Time spent for HVG = {:.2f} seconds.\n".format(end_hvg - start_hvg))

print("Getting Highly Variable Genes from scCloud")
df_hvg = pd.read_csv(hvg_file)
print("Highly variable gene set of size " + str(df_hvg.shape[0]))
df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
df_hvg.set_index('index', inplace = True)
df_hvg.drop(columns = ['gene_ids'], inplace = True)
adata_filter_norm.var.drop(columns = ['highly_variable_genes'], inplace = True)
df_join = adata_filter_norm.var.join(df_hvg)
df_join['highly_variable_genes'].fillna(False, inplace = True)
adata_variable_genes = adata_filter_norm[:, df_join['highly_variable_genes']].copy()


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
start_bbknn = time.time()
bbknn(adata, batch_key = 'Channel', metric = 'euclidean')
end_bbknn = time.time()
print("Time spent for bbknn = " + str(end_nn - start_nn) + " seconds.")
f.write("Time spent for bbknn = " + str(end_nn - start_nn) + " seconds.\n")

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