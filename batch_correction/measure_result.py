import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import read_mtx
from scipy import sparse
import scCloud
from scCloud.tools.diffusion_map import calculate_affinity_matrix
from termcolor import cprint
from sklearn.metrics import silhouette_score

method_list = ["baseline", "scCloud", "seurat", "mnn", "combat", "bbknn"]

measure_result = []

def process_baseline():
	cprint("For scCloud with no batch correction:", "red")
	f_list = [f for f in os.listdir("./ground") if f in ["ground.h5ad", "ground_cell_types.h5ad"]]
	if len(f_list) != 2:
		cprint("No processed data are found! Processing data using scCloud...", "green")
		if os.system("cd ./ground/ && ./run_sccloud.sh && cd .."):
			sys.exit(1)
		if os.system("cd ./ground/ && python gen_celltype.py && cd .."):
			sys.exit(1)

	cprint("Loading processed data...", "green")
	adata = scCloud.tools.read_input('./ground/ground.h5ad', mode = 'a')

	process_data(adata, './ground/ground_cell_types_result', method = 'baseline', processed = True)

def process_sccloud():
	cprint("For scCloud:", "red")
	f_list = [f for f in os.listdir("./scCloud") if f in ["tiny_sccloud_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using scCloud...", "green")
		if os.system("cd ./scCloud/ && python run_sccloud_batch_correct.py && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./scCloud/tiny_sccloud_corrected.h5ad', mode = 'a')

	process_data(adata, "./scCloud/tiny_sccloud_result", method = 'scCloud', processed = True)

def process_mnn():
	cprint("For MNN:", "red")
	f_list = [f for f in os.listdir("./mnn") if f in ["scanpy_mnn_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using MNN...", "green")
		if os.system("cd ./mnn/ && conda activate scanpy && python scanpy_mnn.py && conda deactivate && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./mnn/scanpy_mnn_corrected.h5ad', mode = 'a')

	cprint("Enforcing count matrix to be sparse...", "green")
	adata.X = sparse.csr_matrix(adata.X)

	cprint("Correcting cell names...", "green")
	adata.obs['cell_id'] = adata.obs.index.map(lambda s: s[:s.rfind('-')]).values
	adata.obs.set_index('cell_id', inplace = True)

	process_data(adata, "./mnn/scanpy_mnn_result", method = 'mnn')

def process_seurat():
	cprint("For Seurat:", "red")
	f_list = [f for f in os.listdir("./seurat") if f in ["gene_expression.mtx", "feature_names.txt", "sample_names.txt"]]
	if len(f_list) != 3:
		cprint("No corrected data are found! Correcting data using CCA...", "green")
		if os.system("cd ./seurat/ && Rscript seurat_cca.R && cd .."):
			sys.exit(1)

	cprint("Loading gene expression...", "green")
	adata = read_mtx("./seurat/gene_expression.mtx")

	cprint("Enforce count matrix to be sparse...", "green")
	adata.X = sparse.csr_matrix(adata.X)
	
	cprint("Loading feature names...", "green")
	df_features = pd.read_csv("./seurat/feature_names.txt", header = None)
	adata.var['index'] = df_features[0].values
	adata.var.set_index('index', inplace = True)

	cprint("Loading sample names...", "green")
	df_samples = pd.read_csv("./seurat/sample_names.txt", header = None)
	adata.obs['index'] = df_samples[0].values
	adata.obs.set_index('index', inplace = True)
	adata.obs['Channel'] = pd.Categorical(adata.obs.index.map(lambda s: s.split('-')[0]).values)
	adata.uns['genome'] = 'GRCh38'


	process_data(adata, "./seurat/seurat_result", method = 'seurat')

def process_combat():
	cprint("For ComBat:", "red")
	f_list = [f for f in os.listdir("./combat") if f in ["scanpy_combat_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using ComBat...", "green")
		if os.system("cd ./combat/ && conda activate scanpy && python scanpy_combat.py && conda deactivate && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./combat/scanpy_combat_corrected.h5ad', mode = 'a')

	cprint("Enforce count matrix to be sparse...", "green")
	adata.X = sparse.csr_matrix(adata.X)

	process_data(adata, "./combat/scanpy_combat_result", method = 'combat')

def process_bbknn():
	cprint("For BBKNN:", "red")
	f_list = [f for f in os.listdir("./bbknn") if f in ["scanpy_bbknn_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corredted data are found! Correcting data using BBKNN...", "green")
		if os.system("cd ./combat/ && conda activate scanpy && python scanpy_bbknn.py && conda deactivate && cd .."):
			sys.exit(1)

	cprint("loading corrected data...", "green")
	adata = scCloud.tools.read_input("./bbknn/scanpy_bbknn_corrected.h5ad", mode = 'a')

	cprint("Clustering...", "green")
	adata.uns['knn_indices'] = adata.uns['neighbors']['knn_indices'][:, 1:]
	adata.uns['knn_distances'] = adata.uns['neighbors']['knn_distances'][:, 1:]

	cprint("Computing UMAP...", "green")
	scCloud.tools.run_umap(adata, 'X_pca')

	cprint("For UMAP coordinates:", "yellow")
	process_data(adata, "./bbknn/scanpy_bbknn_result", method = 'bbknn', processed = True)

	scCloud.tools.write_output(adata, "./bbknn/scanpy_bbknn_result")


def process_data(data, output, method, processed = False):

	if not processed:

		cprint("Calculating PCA and KNN...", "green")
		data_c = data.copy()
		data_c.X = data_c.X.toarray()
		scCloud.tools.run_pca(data_c)
		data.obsm['X_pca'] = data_c.obsm['X_pca']
		scCloud.tools.get_kNN(data, 'X_pca', K = 100, n_jobs = 8)

		cprint("Computing UMAP...", "green")
		scCloud.tools.run_umap(data, 'X_pca')

		scCloud.tools.write_output(data, output)

	cprint("Calculating kBET measures on UMAP coordinates...", "green")
	kbet_stat, kbet_pvalue, kbet_ac_rate = scCloud.tools.calc_kBET(data, 'Channel', 'X_umap')
	cprint("Mean statistics is {stat:.4f}; Mean p-value is {pvalue:.4f}; Mean accept rate is {rate:.4f}.".format(stat = kbet_stat, pvalue = kbet_pvalue, rate = kbet_ac_rate), "yellow")

	cprint("Loading ground truth cell types...", "green")
	df_celltype = pd.read_csv("ground_cell_type.txt")
	assert np.sum(df_celltype['cell_id'] != data.obs.index.values) == 0
	data.obs['cell_types'] = pd.Categorical(df_celltype['cell_types'].values)

	# Set Individual
	if method == 'mnn':
		data.obs['Individual'] = pd.Categorical(data.obs.index.map(lambda s: s.split('_')[0][8:]))
	else:
		data.obs['Individual'] = pd.Categorical(data.obs['Channel'].apply(lambda s: s.split('_')[0][8:]))

	cprint("Calculating Mean Silhouette Score on UMAP coordinates...", "green")
	sil_score = silhouette_score(data.obsm['X_umap'], data.obs['cell_types'])
	cprint("Mean Silhouette Score on UMAP = {:.4f}.".format(sil_score), "yellow")

	cprint("Calculating kSIM on UMAP coordinates...", "green")
	ksim_mean, ksim_ac_rate = scCloud.tools.calc_kSIM(data, 'cell_types', rep_key = 'X_umap')

	cprint("Mean kSIM = {mean:.4f}, with accept rate {rate:.4f}.".format(mean = ksim_mean, rate = ksim_ac_rate), "yellow")

	measure_result.append((method, ksim_ac_rate, kbet_ac_rate))

	cprint("Plotting UMAP for cells with cell types...", "green")
	scCloud.tools.write_output(data, "{}_compare".format(method))
	if os.system("scCloud plot scatter --basis umap --attributes cell_types,Individual {name}_compare.h5ad {name}.celltypes.umap.pdf".format(name = method)):
		sys.exit(1)

def plot_scatter(precomputed = False):
	if precomputed:
		measure_result = [('Baseline', 0.8813, 0.1073), ('scCloud', 0.8137, 0.3540), ('seurat', 0.7412, 0.6070), ('MNN', 0.8553, 0.2678), ('ComBat', 0.8155, 0.3301), ('BBKNN', 0.8636, 0.2681)]

	method_l = []
	ksim_l = []
	kbet_l = []
	df_measure = pd.DataFrame({'method':[], 'kSIM':[], 'kBET':[]})
	for (method, ksim, kbet) in measure_result:
		method_l.append(method)
		ksim_l.append(ksim)
		kbet_l.append(kbet)
	
	df_measure['method'] = method_l
	df_measure['kSIM'] = ksim_l
	df_measure['kBET'] = kbet_l

	ax = sns.scatterplot(x = 'kSIM', y = 'kBET', hue = 'method', data = df_measure, legend = False)
	for line in range(0, df_measure.shape[0]):
		x_pos = df_measure.kSIM[line]
		y_pos = df_measure.kBET[line]
		if df_measure.method[line] == 'MNN':
			x_pos -= 0.018
			y_pos += 0.003
		elif df_measure.method[line] == 'BBKNN':
			x_pos += 0.003
			y_pos -= 0.003
		else:
			x_pos += 0.003
		ax.text(x_pos, y_pos, df_measure.method[line], horizontalalignment = 'left', size = 'medium', color = 'black')
	plt.xlabel('kSIM accept rate')
	plt.ylabel('kBET accept rate')
	plt.xlim(0.73, 0.91)


	plt.savefig("Figure_2C.pdf")
	plt.close()



if __name__ == "__main__":
	method = sys.argv[1]
	assert method in method_list or method == 'all' or method == 'plot'

	if method == 'baseline' or method == 'all':
		process_baseline()

	if method == 'scCloud' or method == 'all':
		process_sccloud()

	if method == 'seurat' or method == 'all':
		process_seurat()

	if method == 'mnn' or method == 'all':
		process_mnn()

	if method == 'combat' or method == 'all':
		process_combat()

	if method == 'bbknn' or method == 'all':
		process_bbknn()

	if method == 'all' or method == 'plot':
		precomputed = True if method == 'plot' else False
		plot_scatter(precomputed)
