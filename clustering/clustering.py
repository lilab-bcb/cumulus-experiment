import os, sys
import scCloud
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

n_cores = os.cpu_count()
data_src = "MantonBM_nonmix_10x_corrected"
data_dst = "spectral_result"
spectral_label = "spectral_labels"

def run_spectral(data, rep_key, n_clusters = 20, K = 100, n_jobs = 1, random_state = 0, full_speed = False):

	X = data.obsm[rep_key].astype('float64')
	km = KMeans(n_clusters = n_clusters, n_jobs = n_jobs, random_state = random_state)
	km.fit(X)

	data.obs[spectral_label] = pd.Categorical(km.labels_)



if __name__ == '__main__':

	if os.system("scCloud cluster -p {jobs} --correct-batch-effect ../MantonBM_nonmix_10x.h5 {name}".format(jobs = n_cores, name = data_src)):
		sys.exit(1)

	adata = scCloud.tools.read_input(data_src + '.h5ad', mode = 'a')
	run_spectral(adata, 'X_diffmap', n_clusters = 20, n_jobs = n_cores)
	scCloud.tools.run_fitsne(adata, 'X_pca', n_jobs = n_cores)
	scCloud.tools.write_output(adata, data_dst)

	if os.system("scCloud plot scatter --basis fitsne --attributes {label},Individual {name}.h5ad {name}.fitsne.pdf".format(label = spectral_label, name = data_dst)):
		sys.exit(1)