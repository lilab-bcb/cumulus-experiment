import os, sys
import pegasus as pg
from termcolor import cprint

n_cores = os.cpu_count()

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

pg.neighbors(adata)
pg.louvain(adata)
pg.umap(adata)
pg.write_output(adata, "ground")

del adata

if os.system("pegasus de_analysis -p {} --labels louvain_labels --t ground.h5ad ground.de.xlsx".format(n_cores)):
	sys.exit(1)

if os.system("pegasus annotate_cluster ground.h5ad ground.anno.txt"):
	sys.exit(1)

if os.system("pegasus plot scatter --basis umap --attributes louvain_labels,Channel ground.h5ad ground.umap.pdf"):
	sys.exit(1)