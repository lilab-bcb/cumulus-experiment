import os, sys
import sccloud as scc
import pandas as pd
from termcolor import cprint

bm_full_data_source = "/projects/benchmark/MantonBM/MantonBM_nonmix.h5sc"
bm_tiny_data_source = "/projects/benchmark/MantonBM/MantonBM_nonmix_tiny.h5sc"
neuron_data_source = "/projects/benchmark/1M_Neurons/1M_neurons.h5"

def extract_hvf(input_file):

	fname_prefix = os.path.splitext(input_file)[0]

	# Run correction.
	cprint("Extracting HVGs from {}...".format(input_file), "green")
	adata = scc.read_input(input_file)
	scc.qc_metrics(adata)
	scc.filter_data(adata)
	scc.log_norm(adata)
	scc.highly_variable_features(adata, consider_batch = True, n_jobs = 8)
	scc.correct_batch(adata, features = "highly_variable_features")

	# Extract highly variable gene set to file.
	cprint("Extracting highly variable features to file...", "green")
	df_hvf = adata.var.loc[adata.var["highly_variable_features"]][["gene_ids"]]
	df_hvf.to_csv(fname_prefix + "_hvf.txt", index_label = "index")

	# PCA and KNN
	scc.pca(adata)
	scc.neighbors(adata, n_jobs = 20)
	scc.write_output(adata, fname_prefix + "_corrected")

	bdata = adata[:, adata.var['highly_variable_features']].copy()
	scc.write_output(bdata, fname_prefix + "_corrected_hvf")

def preprocess_data(in_file):

	cprint("Preprocessing {}...".format(in_file), "green")

	# Filtration and Normalization.
	suffix1 = "_filter_norm"
	output_name1 = os.path.splitext(in_file)[0] + suffix1
	output_hvf_name1 = output_name1 + "_hvf"

	adata = scc.read_input(in_file)
	scc.qc_metrics(adata)
	scc.filter_data(adata)
	scc.log_norm(adata)

	print("Getting Highly Variable Feature Set...")
	df_hvf = pd.read_csv(hvf_file)
	df_hvf['highly_variable_features'] = [True] * df_hvf.shape[0]
	df_hvf.set_index('index', inplace = True)
	df_hvf.drop(columns = ['gene_ids'], inplace = True)
	adata.var.drop(columns = ['highly_variable_features'], inplace = True)
	df_join = adata.var.join(df_hvf)
	df_join['highly_variable_features'].fillna(False, inplace = True)
	adata.var['highly_variable_features'] = df_join['highly_variable_features']
	adata_hvf = adata[:, adata.var['highly_variable_features']].copy()

	scc.write_output(adata, output_name1)
	scc.write_output(adata_hvf, output_hvf_name1)

	# Furthermore, do PCA.
	suffix2 = "_filter_norm_pca"
	output_name2 = os.path.splitext(in_file)[0] + suffix2
	output_hvf_name2 = output_name2 + "_hvf"

	scc.pca(adata)
	scc.write_output(adata, output_name2)

	scc.pca(adata_hvf)
	scc.write_output(adata_hvf, output_hvf_name2)

if __name__ == '__main__':
	extract_hvf(bm_full_data_source)
	extract_hvf(neuron_data_source)
	preprocess_data(full_data_source)
	preprocess_data(tiny_data_source)
	preprocess_data(neuron_data_source)