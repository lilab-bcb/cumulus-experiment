import os, time
import numpy as np
import scanpy as sc
from datetime import timedelta, datetime

sc.settings.n_jobs = os.cpu_count()
sc.settings.verbosity = 4
rand_seed = 0

src_data = "/data/MantonBM_nonmix_10x.h5"
log_file = "mantonbm_scanpy.log"

fp = open(log_file, 'w')

start_read = time.time()
adata = sc.read_10x_h5(src_data, genome = "GRCh38")
adata.var_names_make_unique()
adata.obs['Channel'] = adata.obs.index.map(lambda s: s.split('-')[0]).values
end_read = time.time()
fp.write("Read starts at: {}\n".format(datetime.fromtimestamp(start_read)))
fp.write("Read ends at: {}\n".format(datetime.fromtimestamp(end_read)))
fp.write("Time spent for Read: {}.\n\n".format(timedelta(seconds = end_read - start_read)))

start_filter = time.time()
n_cells = adata.shape[0]
sc.pp.filter_cells(adata, min_counts = 100)
sc.pp.filter_cells(adata, max_counts = 600000)
sc.pp.filter_cells(adata, min_genes = 500)
sc.pp.filter_cells(adata, max_genes = 6000)
sc.pp.filter_genes(adata, min_cells = n_cells * 0.0005)

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
adata.obs['n_counts'] = adata.X.sum(axis = 1).A1
adata = adata[adata.obs.percent_mito < 0.1, :]
end_filter = time.time()

fp.write("Filter starts at: {}\n".format(datetime.fromtimestamp(start_filter)))
fp.write("Filter ends at: {}\n".format(datetime.fromtimestamp(end_filter)))
fp.write("Time spent for Filter: {}.\n\n".format(timedelta(seconds = end_filter - start_filter)))

start_norm = time.time()
sc.pp.normalize_total(adata, target_sum = 1e5)
sc.pp.log1p(adata)
end_norm = time.time()
fp.write("LogNorm starts at: {}\n".format(datetime.fromtimestamp(start_norm)))
fp.write("LogNorm ends at: {}\n".format(datetime.fromtimestamp(end_norm)))
fp.write("Time spent for LogNorm: {}.\n\n".format(timedelta(seconds = end_norm - start_norm)))

#adata.raw = adata

start_hvg = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel')
end_hvg = time.time()

fp.write("HVG starts at: {}\n".format(datetime.fromtimestamp(start_hvg)))
fp.write("HVG ends at: {}\n".format(datetime.fromtimestamp(end_hvg)))
fp.write("Time spent for HVG: {}.\n\n".format(timedelta(seconds = end_hvg - start_hvg)))

adata = adata[:, adata.var.highly_variable]

start_pca = time.time()
sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata, n_comps = 50, random_state = rand_seed)
end_pca = time.time()
fp.write("PCA starts at: {}\n".format(datetime.fromtimestamp(start_pca)))
fp.write("PCA ends at: {}\n".format(datetime.fromtimestamp(end_pca)))
fp.write("Time spent for PCA: {}.\n\n".format(timedelta(seconds = end_pca - start_pca)))

start_knn = time.time()
sc.pp.neighbors(adata, n_neighbors = 100, n_pcs = 50, random_state = 0)
end_knn = time.time()
fp.write("KNN starts at: {}\n".format(datetime.fromtimestamp(start_knn)))
fp.write("KNN ends at: {}\n".format(datetime.fromtimestamp(end_knn)))
fp.write("Time spent for KNN: {}.\n\n".format(timedelta(seconds = end_knn - start_knn)))

start_louvain = time.time()
sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)
end_louvain = time.time()
fp.write("Louvain starts at: {}\n".format(datetime.fromtimestamp(start_louvain)))
fp.write("Louvain ends at: {}\n".format(datetime.fromtimestamp(end_louvain)))
fp.write("Time spent for Louvain: {}.\n\n".format(timedelta(seconds = end_louvain - start_louvain)))

start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', use_fast_tsne = True, random_state = 0)
end_tsne = time.time()
fp.write("tSNE starts at: {}\n".format(datetime.fromtimestamp(start_tsne)))
fp.write("tSNE ends at: {}\n".format(datetime.fromtimestamp(end_tsne)))
fp.write("Time spent for tSNE: {}.\n\n".format(timedelta(seconds = end_tsne - start_tsne)))

start_umap = time.time()
sc.tl.umap(adata, random_state = 0)
end_umap = time.time()
fp.write("UMAP starts at: {}\n".format(datetime.fromtimestamp(start_umap)))
fp.write("UMAP ends at: {}\n".format(datetime.fromtimestamp(end_umap)))
fp.write("Time spent for UMAP: {}.\n\n".format(timedelta(seconds = end_umap - start_umap)))

start_write = time.time()
adata.write("mantonbm_scanpy_result.h5ad")
end_write = time.time()
fp.write("Write starts at: {}\n".format(datetime.fromtimestamp(start_write)))
fp.write("Write ends at: {}\n".format(datetime.fromtimestamp(end_write)))
fp.write("Time spent for Write: {}.\n\n".format(timedelta(seconds = end_write - start_write)))

fp.close()