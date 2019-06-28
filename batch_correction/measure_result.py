import os, sys
import numpy as np
import pandas as pd
from anndata import read_mtx
from scipy import sparse
import scCloud
from scCloud.tools.diffusion_map import calculate_affinity_matrix
from termcolor import cprint
from sklearn.metrics import silhouette_score

method_list = ["scCloud", "seurat", "mnn", "combat", "bbknn"]

def process_sccloud():
	cprint("For scCloud:", "red")
	f_list = [f for f in os.listdir("./scCloud") if f in ["tiny_sccloud_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using scCloud...", "green")
		if os.system("cd ./scCloud/ && ./run_sccloud_batch_correct.sh && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./scCloud/tiny_sccloud_corrected.h5ad', mode = 'a')

	process_data(adata, "./scCloud/tiny_sccloud_result", method = 'scCloud', processed = True)

	#cprint("For UMAP coordinates:", "green")
	#process_data(adata, "./scCloud/tiny_sccloud_result", method = 'scCloud', rep_key = 'X_umap', processed = True)

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

	cprint("Enforce count matrix to be sparse...", "green")
	adata.X = sparse.csr_matrix(adata.X)

	cprint("Clustering...", "green")
	adata.uns['knn_indices'] = adata.uns['neighbors']['knn_indices'][:, 1:]
	adata.uns['knn_distances'] = adata.uns['neighbors']['knn_distances'][:, 1:]

	W = calculate_affinity_matrix(adata.uns['knn_indices'], adata.uns['knn_distances'])
	adata.uns['W'] = W
	scCloud.tools.run_leiden(adata)

	cprint("Computing UMAP...", "green")
	scCloud.tools.run_umap(adata, 'X_pca')

	cprint("For UMAP coordinates:", "yellow")
	process_data(adata, "./bbknn/scanpy_bbknn_result", method = 'bbknn', rep_key = 'X_umap', processed = True)

	scCloud.tools.write_output(adata, "./bbknn/scanpy_bbknn_result")

def process_data(data, output, method, rep_key = 'X_pca', processed = False):

	if not processed:
		cprint("Calculating PCA and KNN...", "green")
		data.var['highly_variable_genes'] = [True] * data.shape[1]
		data_c = scCloud.tools.collect_highly_variable_gene_matrix(data)
		scCloud.tools.run_pca(data_c)
		data.obsm['X_pca'] = data_c.obsm['X_pca']
		scCloud.tools.run_diffmap(data, rep_key, n_jobs = 8, K = 100)
		scCloud.tools.get_kNN(data, 'X_diffmap', 100, n_jobs = 8)
		data.obsm['X_diffmap_pca'] = scCloud.tools.reduce_diffmap_to_3d(data.obsm['X_diffmap'])
		cprint("Clustering...", "green")
		data.uns['W'] = calculate_affinity_matrix(data.uns['knn_indices'], data.uns['knn_distances'])
		scCloud.tools.run_leiden(data)
		cprint("Computing FIt-SNE...", "green")
		scCloud.tools.run_fitsne(data, rep_key, n_jobs = 8)
		cprint("Computing UMAP...", "green")
		scCloud.tools.run_umap(data, rep_key)

		scCloud.tools.write_output(data, output)


	cprint("Calculating kBET measures on {} coordinates...".format(rep_key), "green")
	stat, pvalue, acc_rate = scCloud.tools.calc_kBET(data, 'Channel', rep_key)
	cprint("Mean statistics is {stat:.4f}; Mean p-value is {pvalue:.4f}; Mean accept rate is {rate:.4f}.".format(stat = stat, pvalue = pvalue, rate = acc_rate), "yellow")

	cprint("Loading ground truth cell types...", "green")
	df_celltype = pd.read_csv("ground_cell_type.txt")
	assert np.sum(df_celltype['cell_id'] != data.obs.index.values) == 0
	data.obs['cell_types'] = df_celltype['cell_types'].values
	subset_data = data[data.obs['cell_types'] != 10, :].copy()

	cprint("Calculating Mean Silhouette Score on UMAP coordinates...", "green")
	sil_score = silhouette_score(subset_data.obsm['X_umap'], subset_data.obs['cell_types'])
	cprint("Mean Silhouette Score on UMAP = {:.4f}.".format(sil_score), "yellow")

	cprint("Calculating kSIM on UMAP coordinates...", "green")
	ksim_mean, ksim_ac_rate = scCloud.tools.calc_kSIM(subset_data, 'cell_types', rep_key = 'X_umap')
	cprint("Mean kSIM = {mean:.4f}, with accept rate {rate:.4f}.".format(mean = ksim_mean, rate = ksim_ac_rate), "yellow")

	cprint("Plotting UMAP for cells with known cell types...", "green")
	scCloud.tools.write_output(subset_data, "temp")
	if os.system("scCloud plot scatter --basis umap --attributes leiden_labels,Channel temp.h5ad {}.celltypes.umap.pdf".format(method)):
		sys.exit(1)
	if os.system("rm temp.h5ad"):
		sys.exit(1)


if __name__ == "__main__":
	method = sys.argv[1]
	assert method in method_list or method == 'all'

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