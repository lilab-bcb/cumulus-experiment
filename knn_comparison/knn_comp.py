import sys, os, re, time
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from termcolor import cprint
from joblib import Parallel, delayed
from utils import str_time

data_source = "../MantonBM_nonmix_corrected.h5ad"


def get_NN_brute(data, num_threads = 8, K = 100):
	f_list = [f for f in os.listdir('.') if f == 'baseline_indices.npy']

	if len(f_list) == 0:
		cprint("Use {} cores.".format(num_threads), "green")
		from sklearn.neighbors import NearestNeighbors
		start_brute = time.time()
		neighbors = NearestNeighbors(n_neighbors = K, algorithm = 'brute', metric = 'euclidean', n_jobs = num_threads)
		cprint("Calculating KNN with brute force algorithm...", "yellow")
		neighbors.fit(data.obsm['X_pca'])
		knn_distances, knn_indices = neighbors.kneighbors()
		end_brute = time.time()
		cprint("Finished with {alg} method!".format(alg = neighbors._fit_method), "yellow")
		cprint("Time used for brute force: " + str_time(end_brute - start_brute) + ".", "yellow")
		np.save('baseline_indices.npy', knn_indices)
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated brute force results...", "yellow")
		knn_indices = np.load('baseline_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_sccloud(data, num_threads = 8, K = 100):
	f_list = [f for f in os.listdir('.') if f == 'sccloud_indices.npy']

	if len(f_list) == 0:
		cprint("Use {} cores.".format(num_threads), "green")
		cprint("Calculating KNN with scCloud algorithm...", "yellow")
		from scCloud.tools.nearest_neighbors import calculate_nearest_neighbors
		start_sccloud = time.time()
		knn_indices, knn_distances = calculate_nearest_neighbors(data.obsm['X_pca'], num_threads, K = K, full_speed = True)
		end_sccloud = time.time()
		cprint("Finished!", "yellow")
		cprint("Time used for scCloud: " + str_time(end_sccloud - start_sccloud) + ".", "yellow")
		np.save('sccloud_indices.npy', knn_indices)
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated scCloud results...", "yellow")
		knn_indices = np.load('sccloud_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_scanpy(data, K = 100):
	"""
	Scanpy uses nearest_neighbors function from umap-learn. 
	Must activate scanpy conda environment to use the official version of umap-learn.
	"""
	f_list = [f for f in os.listdir('.') if f == 'scanpy_indices.npy']

	if len(f_list) == 0:
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
		cprint("Results are saved!", "yellow")
	else:
		cprint("Loading pre-calculated scanpy results...", "yellow")
		knn_indices = np.load('scanpy_indices.npy')
		cprint("Finished!", "yellow")

	return knn_indices

def get_NN_Seurat(data):
	f_list = [f for f in os.listdir('.') if f == 'seurat_indices.txt']

	if len(f_list) == 0:
		f_list = [f for f in os.listdir('.') if re.match('x_pca.txt', f)]
		if len(f_list) == 0:
			cprint("Saving PCA as txt file...", "yellow")
			np.savetxt('x_pca.txt', data.obsm['X_pca'])
			cprint("Finished!", "yellow")
		
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

	# violin plot
	ax = sns.violinplot(x = "method", y = "accuracy", data = df, order = ["scCloud", "scanpy", "Seurat V3"])
	ax.set(ylabel = 'Accuracy', xlabel = '')
	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])
	ax.get_figure().savefig("knn_accuracy_violinplot.pdf")
	cprint("Violin plot is generated!", "yellow")
	plt.close()

	# boxplot
	ax = sns.boxplot(x = "method", y = "accuracy", data = df, order = ["scCloud", "scanpy", "Seurat V3"])
	ax.set(ylabel = 'Accuracy', xlabel = '')
	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])
	ax.get_figure().savefig("knn_accuracy_boxplot.pdf")
	cprint("Boxplot is generated!", "yellow")
	plt.close()


if __name__ == '__main__':

	method = sys.argv[1]
	assert method in ["scCloud", "others", "plot"]

	n_cores = os.cpu_count()

	if method == 'scCloud':
		assert os.environ['CONDA_DEFAULT_ENV'] == 'scCloud'

		import scCloud
		adata = scCloud.tools.read_input(data_source, mode = 'a')

		get_NN_brute(adata, num_threads = n_cores)
		
		get_NN_sccloud(adata, num_threads = n_cores)


	elif method == 'others':
		assert os.environ['CONDA_DEFAULT_ENV'] == 'scanpy-benchmark'

		import scanpy as sc
		adata = sc.read_h5ad(data_source)

		# For scanpy
		get_NN_scanpy(adata)

		# For seurat
		get_NN_Seurat(adata)


	else:

		f_list = [f for f in os.listdir('.') if f in ["baseline_indices.npy", "sccloud_indices.npy", "scanpy_indices.npy", "seurat_indices.txt"]]

		if len(f_list) == 4:
			df_list = []
			baseline_indices = get_NN_brute(None)

			# For scCloud
			knn_indices = get_NN_sccloud(None)
			df_sccloud = generate_knn_accuracy_dataframe('scCloud', knn_indices, baseline_indices)
			df_list.append(df_sccloud)

			# For scanpy
			knn_indices = get_NN_scanpy(None)
			df_scanpy = generate_knn_accuracy_dataframe('scanpy', knn_indices, baseline_indices)
			df_list.append(df_scanpy)

			# For seurat
			knn_indices = get_NN_Seurat(None)
			df_seurat = generate_knn_accuracy_dataframe('Seurat V3', knn_indices, baseline_indices)
			df_list.append(df_seurat)
	
			plot_result(df_list)
