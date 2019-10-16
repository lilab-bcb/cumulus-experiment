import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import read_mtx
import pegasus as pg
from termcolor import cprint
from sklearn.metrics import silhouette_score

method_list = ["baseline", "pegasus", "seurat", "mnn", "combat", "bbknn"]

figure_index = {'baseline' : 'B',
                'pegasus'  : 'C',
                'combat'   : 'D',
                'mnn'      : 'E',
                'bbknn'    : 'F',
                'seurat'   : 'G'
               }

method_print_name = {'baseline' : 'Baseline',
					 'pegasus'  : 'Pegasus',
					 'combat'   : 'ComBat',
					 'mnn'      : 'MNN',
					 'bbknn'    : 'BBKNN',
					 'seurat'   : 'Seurat v3'
}

measure_result = []

measure_precomputed_file = "correction_benchmark.txt"

def process_baseline():
	cprint("For pegasus with no batch correction:", "red")
	f_list = [f for f in os.listdir("./baseline") if f in ["baseline_result.h5ad"]]
	if len(f_list) != 1:
		cprint("No processed data are found! Processing data using pegasus without batch correction...", "green")
		if os.system("cd ./baseline/ && python run_baseline.py && cd .."):
			sys.exit(1)

	cprint("Loading processed data...", "green")
	adata = pg.read_input('./baseline/baseline_result.h5ad')

	process_data(adata, method = 'baseline', processed = True)

def process_pegasus():
	cprint("For pegasus:", "red")
	f_list = [f for f in os.listdir("./pegasus") if f in ["tiny_pegasus_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using pegasus...", "green")
		if os.system("cd ./pegasus/ && python run_pegasus_batch_correct.py && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = pg.read_input('./pegasus/tiny_pegasus_corrected.h5ad')

	process_data(adata, method = 'pegasus', processed = True)

def process_mnn():
	cprint("For MNN:", "red")
	f_list = [f for f in os.listdir("./mnn") if f in ["scanpy_mnn_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found!", "red")
		sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = pg.read_input('./mnn/scanpy_mnn_corrected.h5ad')

	cprint("Correcting cell names...", "green")
	adata.obs['cell_id'] = adata.obs.index.map(lambda s: s[:s.rfind('-')]).values
	adata.obs.set_index('cell_id', inplace = True)

	process_data(adata, method = 'mnn', output = "./mnn/scanpy_mnn_result")

def process_seurat():
	cprint("For Seurat:", "red")
	f_list = [f for f in os.listdir("./seurat") if f in ["matrix.mtx", "genes.txt", "barcodes.txt"]]
	if len(f_list) != 3:
		cprint("No corrected data are found!", "red")

	cprint("Loading gene expression...", "green")
	adata = read_mtx("./seurat/matrix.mtx")

	cprint("Loading gene names...", "green")
	df_genes = pd.read_csv("./seurat/genes.txt", header = None)
	adata.var['index'] = df_genes[0].values
	adata.var.set_index('index', inplace = True)

	cprint("Loading barcode names...", "green")
	df_barcodes = pd.read_csv("./seurat/barcodes.txt", header = None)
	adata.obs['index'] = df_barcodes[0].values
	adata.obs.set_index('index', inplace = True)
	adata.obs['Channel'] = pd.Categorical(adata.obs.index.map(lambda s: s.split('-')[0]).values)
	adata.uns['genome'] = 'GRCh38'

	process_data(adata, method = 'seurat', output = "./seurat/seurat_result")

def process_combat():
	cprint("For ComBat:", "red")
	f_list = [f for f in os.listdir("./combat") if f in ["scanpy_combat_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found!", "red")

	cprint("Loading corrected data...", "green")
	adata = pg.read_input('./combat/scanpy_combat_corrected.h5ad')

	process_data(adata, method = 'combat', output = "./combat/scanpy_combat_result")

def process_bbknn():
	cprint("For BBKNN:", "red")
	f_list = [f for f in os.listdir("./bbknn") if f in ["scanpy_bbknn_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corredted data are found!", "red")

	cprint("loading corrected data...", "green")
	adata = pg.read_input("./bbknn/scanpy_bbknn_corrected.h5ad")
	adata.uns['pca_knn_indices'] = adata.uns['neighbors']['knn_indices'][:, 1:]
	adata.uns['pca_knn_distances'] = adata.uns['neighbors']['knn_distances'][:, 1:]

	cprint("Computing UMAP...", "green")
	pg.umap(adata)

	cprint("For UMAP coordinates:", "yellow")
	process_data(adata, method = 'bbknn', processed = True)

	pg.write_output(adata, "./bbknn/scanpy_bbknn_result")


def process_data(data, method, output = None, processed = False):

	if not processed:

		cprint("Calculating PCA and KNN...", "green")
		pg.pca(data, features = "highly_variable_features" if "highly_variable_features" in data.var else None)
		pg.neighbors(data, n_jobs = 8)

		cprint("Computing UMAP...", "green")
		pg.umap(data)

		pg.write_output(data, output)

	cprint("Calculating kBET measures on UMAP coordinates...", "green")
	kbet_stat, kbet_pvalue, kbet_ac_rate = pg.calc_kBET(data, attr = 'Channel', rep = 'umap')
	cprint("Mean statistics is {stat:.4f}; Mean p-value is {pvalue:.4f}; Mean accept rate is {rate:.4f}.".format(stat = kbet_stat, pvalue = kbet_pvalue, rate = kbet_ac_rate), "yellow")

	cprint("Loading ground truth cell types...", "green")
	df_celltype = pd.read_csv("ground_cell_type.txt")
	assert np.sum(df_celltype['cell_id'] != data.obs.index.values) == 0
	data.obs['Cluster'] = pd.Categorical(df_celltype['cell_types'].values)

	# Set Individual
	if method == 'mnn':
		data.obs['Individual'] = pd.Categorical(data.obs.index.map(lambda s: s.split('_')[0][8:]))
	else:
		data.obs['Individual'] = pd.Categorical(data.obs['Channel'].apply(lambda s: s.split('_')[0][8:]))

	cprint("Calculating Mean Silhouette Score on UMAP coordinates...", "green")
	sil_score = silhouette_score(data.obsm['X_umap'], data.obs['Cluster'])
	cprint("Mean Silhouette Score on UMAP = {:.4f}.".format(sil_score), "yellow")

	cprint("Calculating kSIM on UMAP coordinates...", "green")
	ksim_mean, ksim_ac_rate = pg.calc_kSIM(data, attr = 'Cluster', rep = 'umap')
	cprint("Mean kSIM = {mean:.4f}, with accept rate {rate:.4f}.".format(mean = ksim_mean, rate = ksim_ac_rate), "yellow")

	measure_result.append((method_print_name[method], ksim_ac_rate, kbet_ac_rate))

	cprint("Plotting UMAP for cells with cell types...", "green")
	pg.write_output(data, "{}_compare".format(method))
	if os.system("pegasus plot scatter --basis umap --attributes Cluster,Individual {name}_compare.h5ad /output/Figure_S2{idx}.pdf".format(name = method, idx = figure_index[method])):
		sys.exit(1)

def plot_scatter(precomputed = False):
	if precomputed:
		df_measure = pd.read_csv(measure_precomputed_file)
	else:
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
		x_pos = df_measure.kSIM[line] + 0.003
		if df_measure.method[line] in ['Baseline', 'ComBat']:
			x_pos = df_measure.kSIM[line] - 0.003
		y_pos = df_measure.kBET[line]
		alignment = 'right' if df_measure.method[line] in ['Baseline', 'ComBat'] else 'left'
		ax.text(x_pos, y_pos, df_measure.method[line], horizontalalignment = alignment, size = 'medium', color = 'black')
	plt.xlabel('kSIM accept rate')
	plt.ylabel('kBET accept rate')


	plt.savefig("/output/Figure_2A.pdf")
	plt.close()



if __name__ == "__main__":

	method = sys.argv[1]
	assert method in method_list or method == 'all' or method == 'plot'

	if method == 'baseline' or method == 'all':
		process_baseline()

	if method == 'pegasus' or method == 'all':
		process_pegasus()

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
