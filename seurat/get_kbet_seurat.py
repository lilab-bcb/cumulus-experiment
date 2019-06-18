import numpy as np
import pandas as pd
from anndata import read_mtx
from termcolor import cprint
import scCloud

def main():
	cprint("Loading gene expression...", "green")
	adata = read_mtx("gene_expression.mtx")
	
	cprint("Loading feature names...", "green")
	df_features = pd.read_csv("feature_names.txt", header = None)
	adata.var['index'] = df_features[0].values
	adata.var.set_index('index', inplace = True)

	cprint("Loading sample names...", "green")
	df_samples = pd.read_csv("sample_names.txt", header = None)
	adata.obs['index'] = df_samples[0].values
	adata.obs.set_index('index', inplace = True)
	adata.obs['Channel'] = pd.Categorical(adata.obs.index.map(lambda s: s.split('-')[0]).values)

	cprint("Calculating PCA and KNN...", "green")
	scCloud.tools.run_pca(adata)
	scCloud.tools.get_kNN(adata, 'X_pca', K = 100)

	cprint("Calculating kBET measures...", "green")
	stat, pvalue, acc_rate = scCloud.tools.calc_kBET(adata, 'Channel')
	cprint("Mean statistics is {stat:.4f}; Mean p-value is {pvalue:.4f}; Mean accept rate is {rate:.4f}.".format(stat = stat, pvalue = pvalue, rate = acc_rate), "yellow")

if __name__ == "__main__":
	main()