import scanpy as sc
import time

sc.settings.n_jobs = 8
sc.settings.verbosity = 4

adata = sc.read_h5ad("scanpy_knn.h5ad")

start_umap = time.time()
sc.tl.umap(adata, random_state = 0)
end_umap = time.time()
print("Time spent for UMAP = {:.2f} seconds.".format(end_umap - start_umap))

