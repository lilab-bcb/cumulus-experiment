import sys, os, re, time
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from termcolor import cprint
from joblib import Parallel, delayed
from utils import str_time

folder = "../batch_correct_exp"


def get_NN_brute(data, num_threads = 8, K = 100, use_cache = False):
	if not use_cache:
		from sklearn.neighbors import NearestNeighbors
		start_brute = time.time()
		neighbors = NearestNeighbors(n_neighbors = K, algorithm = 'brute', metric = 'euclidean', n_jobs = num_threads)
		cprint("Calculating KNN with brute force algorithm...", "yellow")
		neighbors.fit(data.obsm['X_pca'])
		knn_distances, knn_indices = neighbors.kneighbors()
		end_brute = time.time()
		cprint("Finished with {alg} method!".format(alg = neighbors._fit_method), "yellow")
		cprint("Time used for brute force: " + str_time(end_brute - start_brute) + ".", "yellow")
		np.save('baseline_distances.npy', knn_distances)
		np.save('baseline_indices.npy', knn_indices)
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated brute force results...", "yellow")
		knn_indices = np.load('baseline_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_sccloud(data, num_threads = 8, K = 100, use_cache = False):
	if not use_cache:
		cprint("Calculating KNN with scCloud algorithm...", "yellow")
		from scCloud.tools.nearest_neighbors import calculate_nearest_neighbors
		start_sccloud = time.time()
		knn_indices, knn_distances = calculate_nearest_neighbors(data.obsm['X_pca'], num_threads, K = K, full_speed = True)
		end_sccloud = time.time()
		cprint("Finished!", "yellow")
		cprint("Time used for scCloud: " + str_time(end_sccloud - start_sccloud) + ".", "yellow")
		np.save('sccloud_indices.npy', knn_indices)
		np.save('sccloud_distances.npy', knn_distances)
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated scCloud results...", "yellow")
		knn_indices = np.load('sccloud_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_scanpy(data, num_threads = 8, K = 100, use_cache = False):
	"""
	Scanpy uses nearest_neighbors function from umap-learn. 
	Must activate scanpy conda environment to use the official version of umap-learn.
	"""
	if not use_cache:
		cprint("Calculating KNN with scanpy algorithm...", "yellow")
		from umap.umap_ import nearest_neighbors
		from numpy.random import RandomState
		rs = RandomState(0)
		start_scanpy = time.time()
		knn_indices, knn_distances, _ = nearest_neighbors(data.obsm['X_pca'], K, random_state = rs, metric = 'euclidean', metric_kwds = {}, angular = False, verbose = False)
		end_scanpy = time.time()
		cprint("Finished!", "yellow")
		cprint("Time used for scanpy: " + str_time(end_scanpy - start_scanpy) + ".")
		np.save('scanpy_indices.npy', knn_indices)
		np.save('scanpy_distances.npy', knn_distances)
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated scanpy results...", "yellow")
		knn_indices = np.load('scanpy_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_Seurat(data, use_cache = False):
	f_list = [f for f in os.listdir('.') if re.match('x_pca.txt', f)]
	if len(f_list) == 0:
		cprint("Saving PCA as txt file...", "yellow")
		np.savetxt('x_pca.txt', data.obsm['X_pca'])
		cprint("Finished!", "yellow")

	if not use_cache:
		cprint("Calculating KNN with seurat algorithm...", "yellow")
		if os.system("Rscript seurat_knn.R"):
			sys.exit(1)
		cprint("Finished!", "yellow")
	
	cprint("Loading calculated seurat results...", "yellow")
	knn_indices = np.loadtxt('seurat_indices.txt')
	knn_indices -= 1
	cprint("Finished!", "yellow")
	return knn_indices

def calc_accuracy_for_one_datapoint(idx, knn_indices, baseline_indices):
	# Add self point initially.
	num_correct = 1
	num_total = baseline_indices.shape[1] + 1

	for ptr in knn_indices[idx]:
		if np.sum(np.isin(ptr, baseline_indices[idx])) > 0:
			num_correct += 1

	return num_correct / num_total

def generate_knn_accuracy_dataframe(method, knn_indices, baseline_indices, K = 100, n_jobs = 8, temp_folder = None):

	cprint("Calculating {method} KNN accuracy...".format(method = method), "green")
	num_sample = baseline_indices.shape[0]

	# No need to search for self points.
	knn_indices = knn_indices[:, 0:K-1]
	baseline_indices = baseline_indices[:, 0:K-1]

	accuracy_arr = np.array(Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_accuracy_for_one_datapoint)(i, knn_indices, baseline_indices) for i in range(num_sample)))

	df = pd.DataFrame({'method': [method] * num_sample, 'accuracy': accuracy_arr})
	cprint("{method} result is generated.".format(method = method), "green")
	cprint("Mean = {mean:.2%}, Std = {std:.2%}".format(mean = df['accuracy'].mean(), std = df['accuracy'].std()), "yellow")

	return df

def plot_result(df_list):
	cprint("Plotting...", "yellow")
	df = pd.concat(df_list).reset_index()
	ax = sns.barplot(x = 'method', y = 'accuracy', data = df)
	ax.set(ylabel = 'Accuracy', xlabel = '')
	ax.set_ylim(0.95, 1.01)

	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])
	
	ax.get_figure().savefig("knn_accuracy.pdf")
	cprint("Finished!", "yellow")
	plt.close()

def main():

	# Check cache files.
	f_list = [f for f in os.listdir('.') if re.match('baseline_indices.npy', f)]
	cache_baseline = False if len(f_list) == 0 else True

	f_list = [f for f in os.listdir('.') if re.match('sccloud_indices.npy', f)]
	cache_sccloud = False if len(f_list) == 0 else True

	f_list = [f for f in os.listdir('.') if re.match('scanpy_indices.npy', f)]
	cache_scanpy = False if len(f_list) == 0 else True

	f_list = [f for f in os.listdir('.') if re.match('seurat_indices.txt', f)]
	cache_seurat = False if len(f_list) == 0 else True

	adata = None
	conda_env = os.environ['CONDA_DEFAULT_ENV']
	if (not cache_baseline) or (not cache_sccloud) or (not cache_scanpy):
		# Pick up appropriate h5ad reader function.
		if conda_env == 'scCloud':
			import scCloud
			adata = scCloud.tools.read_input('{folder}/MantonBM_nonmix.h5ad'.format(folder = folder), mode = 'a')
		elif conda_env == 'scanpy-benchmark':
			import scanpy as sc
			adata = sc.read_h5ad('{folder}/MantonBM_nonmix.h5ad'.format(folder = folder))
		else:
			raise ValueError("conda environment must be either \'scCloud\' or \'scanpy\'!")

	baseline_indices = get_NN_brute(adata, use_cache = cache_baseline)

	df_list = []
	data_ready = True

	# For scCloud
	if conda_env == 'scCloud' or cache_sccloud:
		knn_indices = get_NN_sccloud(adata, use_cache = cache_sccloud)
		df_sccloud = generate_knn_accuracy_dataframe('scCloud', knn_indices, baseline_indices)
		df_list.append(df_sccloud)
	else:
		data_ready = False

	# For scanpy
	if conda_env == 'scanpy' or cache_scanpy:
		knn_indices = get_NN_scanpy(adata, use_cache = cache_scanpy)
		knn_indices = knn_indices[:, 1:]
		df_scanpy = generate_knn_accuracy_dataframe('scanpy', knn_indices, baseline_indices)
		df_list.append(df_scanpy)
	else:
		data_ready = False

	# For seurat
	knn_indices = get_NN_Seurat(adata, use_cache = cache_seurat)
	knn_indices = knn_indices[:, 1:]
	df_seurat = generate_knn_accuracy_dataframe('Seurat', knn_indices, baseline_indices)
	df_list.append(df_seurat)

	if data_ready:
		plot_result(df_list)

if __name__ == '__main__':
	main()
