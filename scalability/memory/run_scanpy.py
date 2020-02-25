import os, sys, time
import numpy as np
import scanpy as sc
from datetime import timedelta, datetime

def record_time(start_time, end_time, step_name, fp):
    fp.write("{step} starts at: {time}\n".format(step = step_name, time = datetime.fromtimestamp(start_time)))
    fp.write("{step} ends at: {time}\n".format(step = step_name, time = datetime.fromtimestamp(end_time)))
    fp.write("Time spent for {step}: {time}.\n\n".format(step = step_name, time = timedelta(seconds = end_time - start_time)))

sc.settings.n_jobs = os.cpu_count()
sc.settings.verbosity = 3
rand_seed = 0

src_data = sys.argv[1]
out_name = sys.argv[2]
correct_batch = True if sys.argv[3] == '--batch' else False
log_file = "{}_scanpy.log".format(out_name)

fp = open(log_file, 'w')

start_read = time.time()
adata = sc.read_10x_h5(src_data, genome = "GRCh38")
adata.var_names_make_unique()
if correct_batch:
    adata.obs['Channel'] = adata.obs.index.map(lambda s: s.split('-')[0]).values
end_read = time.time()
record_time(start_read, end_read, "Read", fp)

start_filter = time.time()
n_cells = adata.shape[0]
sc.pp.filter_cells(adata, min_counts = 100)
sc.pp.filter_cells(adata, max_counts = 600000)
sc.pp.filter_cells(adata, min_genes = 500)
sc.pp.filter_cells(adata, max_genes = 6000)
sc.pp.filter_genes(adata, min_cells = n_cells * 0.0005)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
adata = adata[adata.obs.percent_mito < 0.1, :]
end_filter = time.time()
record_time(start_filter, end_filter, "Filter", fp)

start_norm = time.time()
sc.pp.normalize_total(adata, target_sum = 1e5)
sc.pp.log1p(adata)
end_norm = time.time()
record_time(start_norm, end_norm, "LogNorm", fp)

start_hvg = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=7, min_disp=0.5, n_top_genes = 2000, batch_key = 'Channel' if correct_batch else None)
end_hvg = time.time()
record_time(start_hvg, end_hvg, "HVG", fp)

adata = adata[:, adata.var.highly_variable]

start_pca = time.time()
sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata, n_comps = 50, random_state = rand_seed)
end_pca = time.time()
record_time(start_pca, end_pca, "PCA", fp)

start_knn = time.time()
sc.pp.neighbors(adata, n_neighbors = 100, n_pcs = 50, random_state = 0)
end_knn = time.time()
record_time(start_knn, end_knn, "KNN", fp)

start_louvain = time.time()
sc.tl.louvain(adata, resolution = 1.3, random_state = rand_seed)
end_louvain = time.time()
record_time(start_louvain, end_louvain, "Louvain", fp)

start_tsne = time.time()
sc.tl.tsne(adata, use_rep = 'X_pca', use_fast_tsne = True, random_state = 0)
end_tsne = time.time()
record_time(start_tsne, end_tsne, "tSNE", fp)

start_umap = time.time()
sc.tl.umap(adata, random_state = 0)
end_umap = time.time()
record_time(start_umap, end_umap, "UMAP", fp)

start_write = time.time()
adata.write("mantonbm_scanpy_result.h5ad")
end_write = time.time()
record_time(start_write, end_write, "Write", fp)

fp.close()