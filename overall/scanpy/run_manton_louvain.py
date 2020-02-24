import numpy as np
import scanpy as sc

from scanpy.neighbors import compute_connectivities_umap
from datetime import timedelta, datetime

sc.settings.n_jobs = os.cpu_count()
sc.settings.verbosity = 4
rand_seed = 0

corrected_data = "/data/MantonBM_nonmix_corrected.h5ad"

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
print("Louvain starts at: {}".format(datetime.fromtimestamp(start_louvain)))
print("Louvain ends at: {}".format(datetime.fromtimestamp(end_louvain)))
logstr_louvain = "Time spent for louvain: {}.".format(timedelta(seconds = end_louvain - start_louvain))
print(logstr_louvain)

adata.write("manton_louvain.h5ad", compression = 'gzip')