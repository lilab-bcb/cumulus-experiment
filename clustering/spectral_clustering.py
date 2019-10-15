import os, sys
import pegasus as pg
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import time
from termcolor import cprint

n_cores = os.cpu_count()
data_src = "../MantonBM_nonmix_pegasus"
data_dst = "spectral_result"
spectral_label = "spectral_labels"

def run_spectral(data, rep_key, n_clusters, K = 100, n_jobs = 1, random_state = 0):

	X = data.obsm[rep_key].astype('float64')
	km = KMeans(n_clusters = n_clusters, n_jobs = n_jobs, random_state = random_state)
	km.fit(X)

	data.obs[spectral_label] = pd.Categorical(km.labels_)



if __name__ == '__main__':

	if (data_dst + '.h5ad') not in os.listdir('.'):
		adata = pg.read_input(data_src + '.h5ad')
		start_spec = time.time()
		run_spectral(adata, 'X_diffmap', n_clusters = 21, n_jobs = n_cores)
		end_spec = time.time()
		cprint("Time for spectral clustering = {} s.".format(end_spec - start_spec), "yellow")
		pg.fitsne(adata, rep = 'pca', n_jobs = n_cores)
		pg.write_output(adata, data_dst)

	if os.system("pegasus plot scatter --basis fitsne --attributes {label} {name}.h5ad {name}.fitsne.pdf".format(label = spectral_label, name = data_dst)):
		sys.exit(1)

	if os.system("pegasus de_analysis -p {jobs} --labels {label} --t {name}.h5ad {name}.de.xlsx".format(jobs = n_cores, label = spectral_label, name = data_dst)):
		sys.exit(1)

	if os.system("pegasus annotate_cluster {name}.h5ad {name}.anno.txt".format(name = data_dst)):
		sys.exit(1)