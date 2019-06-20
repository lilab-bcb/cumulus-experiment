import os, sys
import numpy as np
import pandas as pd
from anndata import read_mtx
import scCloud
from termcolor import cprint

def process_sccloud():
	cprint("For scCloud:", "red")
	f_list = [f for f in os.listdir("./scCloud") if f in ["tiny_sccloud_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using scCloud...", "green")
		if os.system("cd ./scCloud/ && ./run_sccloud_batch_correct.sh && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./scCloud/tiny_sccloud_corrected.h5ad', mode = 'a')

	process_data(adata, processed = True)

def process_mnn():
	cprint("For MNN:", "red")
	f_list = [f for f in os.listdir("./mnn") if f in ["scanpy_mnn_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using MNN...", "green")
		if os.system("cd ./mnn/ && conda activate scanpy && python scanpy_mnn.py && conda deactivate && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./mnn/scanpy_mnn_corrected.h5ad', mode = 'a')

	process_data(adata)

def process_seurat():
	cprint("For Seurat:", "red")
	f_list = [f for f in os.listdir("./seurat") if f in ["gene_expression.mtx", "feature_names.txt", "sample_names.txt"]]
	if len(f_list) != 3:
		cprint("No corrected data are found! Correcting data using CCA...", "green")
		if os.system("cd ./seurat/ && Rscript seurat_cca.R && cd .."):
			sys.exit(1)

	cprint("Loading gene expression...", "green")
	adata = read_mtx("./seurat/gene_expression.mtx")
	
	cprint("Loading feature names...", "green")
	df_features = pd.read_csv("./seurat/feature_names.txt", header = None)
	adata.var['index'] = df_features[0].values
	adata.var.set_index('index', inplace = True)

	cprint("Loading sample names...", "green")
	df_samples = pd.read_csv("./seurat/sample_names.txt", header = None)
	adata.obs['index'] = df_samples[0].values
	adata.obs.set_index('index', inplace = True)
	adata.obs['Channel'] = pd.Categorical(adata.obs.index.map(lambda s: s.split('-')[0]).values)

	process_data(adata)

def process_combat():
	cprint("For ComBat:", "red")
	f_list = [f for f in os.listdir("./combat") if f in ["scanpy_combat_corrected.h5ad"]]
	if len(f_list) != 1:
		cprint("No corrected data are found! Correcting data using ComBat...", "green")
		if os.system("cd ./combat/ && conda activate scanpy && python scanpy_combat.py && conda deactivate && cd .."):
			sys.exit(1)

	cprint("Loading corrected data...", "green")
	adata = scCloud.tools.read_input('./combat/scanpy_combat_corrected.h5ad', mode = 'a')

	process_data(adata)

def process_data(data, processed = False):

	if not processed:
		cprint("Calculating PCA and KNN...", "green")
		scCloud.tools.run_pca(data)
		scCloud.tools.get_kNN(data, 'X_pca', K = 100)

	cprint("Calculating kBET measures...", "green")
	stat, pvalue, acc_rate = scCloud.tools.calc_kBET(data, 'Channel')
	cprint("Mean statistics is {stat:.4f}; Mean p-value is {pvalue:.4f}; Mean accept rate is {rate:.4f}.".format(stat = stat, pvalue = pvalue, rate = acc_rate), "yellow")

def main():
	#process_sccloud()
	#process_mnn()
	process_seurat()
	#process_combat()

if __name__ == "__main__":
	main()