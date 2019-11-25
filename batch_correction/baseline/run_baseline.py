import os, sys
import numpy as np
import pegasus as pg

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

adata = pg.read_input("/data/MantonBM_nonmix_tiny.h5sc")
pg.qc_metrics(adata)
pg.filter_data(adata)
pg.log_norm(adata)
pg.highly_variable_features(adata, consider_batch = False)

# Load precalculated PCA
cprint("Set precalculated PCA info...", "green")
adata.obsm['X_pca'] = np.load("/data/precalculated/batch_correction/pca.npy")
adata.uns['PCs'] = np.load("/data/precalculated/batch_correction/PCs.npy")
adata.uns['pca'] = {}
adata.uns['pca']['variance'] = np.load("/data/precalculated/batch_correction/pca_variance.npy")
adata.uns['pca']['variance_ratio'] = np.load("/data/precalculated/batch_correction/pca_variance_ratio.npy")

if os.system("pegasus cluster -p {cores} --umap /data/MantonBM_nonmix_tiny.h5sc baseline_result > baseline.log".format(cores = n_jobs)):
	sys.exit(1)