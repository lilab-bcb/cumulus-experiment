import os, sys
import scCloud
import pandas as pd
from termcolor import cprint

full_data_source = "MantonBM_nonmix_10x.h5"
tiny_data_source = "MantonBM_nonmix_tiny_10x.h5"

hvg_file = "hvg.txt"

def extract_hvg():

	# Run correction.
	cprint("Running scCloud batch correction on full data...", "green")
	adata = scCloud.tools.read_input(full_data_source)
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata)
	scCloud.tools.log_norm(adata, norm_count = 1e5)
	scCloud.tools.set_group_attribute(adata, attribute_string = None)
	scCloud.tools.estimate_adjustment_matrices(adata)
	scCloud.tools.select_highly_variable_genes(adata, consider_batch = True, n_jobs = 8)

	adata_c = scCloud.tools.collect_highly_variable_gene_matrix(adata)
	scCloud.tools.correct_batch_effects(adata_c)
	scCloud.tools.run_pca(adata_c)
	adata.obsm['X_pca'] = adata_c.obsm['X_pca']
	
	scCloud.tools.get_kNN(adata, 'X_pca', K = 100, n_jobs = 8)

	# Extract highly variable gene set to file.
	cprint("Extracting highly variable gene set to file...", "green")
	df_hvg = adata.var.loc[adata.var["highly_variable_genes"]][["gene_ids"]]
	df_hvg.to_csv(hvg_file)

	scCloud.tools.write_output(adata, "MantonBM_nonmix_corrected")

	bdata = adata[:, adata.var['highly_variable_genes']].copy()
	scCloud.tools.write_output(bdata, "MantonBM_nonmix_corrected_hvg")

def preprocess_data(in_file):

	cprint("Preprocessing {}...".format(in_file), "green")

	# Filtration and Normalization.
	suffix1 = "_filter_norm"
	output_name1 = in_file.split('.')[0] + suffix1
	output_hvg_name1 = output_name1 + "_hvg"

	adata = scCloud.tools.read_input(in_file)
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata)
	scCloud.tools.log_norm(adata, norm_count = 1e5)
	
	print("Getting Highly Variable Gene Set...")
	df_hvg = pd.read_csv(hvg_file)
	df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
	df_hvg.set_index('index', inplace = True)
	df_hvg.drop(columns = ['gene_ids'], inplace = True)
	adata.var.drop(columns = ['highly_variable_genes'], inplace = True)
	df_join = adata.var.join(df_hvg)
	df_join['highly_variable_genes'].fillna(False, inplace = True)
	adata_hvg = adata[:, df_join['highly_variable_genes']].copy()

	scCloud.tools.write_output(adata, output_name1)
	scCloud.tools.write_output(adata_hvg, output_hvg_name1)

	# Furthermore, do PCA.
	suffix2 = "_filter_norm_pca"
	output_name2 = in_file.split('.')[0] + suffix2
	output_hvg_name2 = output_name2 + "_hvg"
	adata_hvg_c = adata_hvg.copy()
	scCloud.tools.run_pca(adata_hvg_c)

	adata.obsm['X_pca'] = adata_hvg_c.obsm['X_pca']
	scCloud.tools.write_output(adata, output_name2)

	adata_hvg.obsm['X_pca'] = adata_hvg_c.obsm['X_pca']
	scCloud.tools.write_output(adata_hvg, output_hvg_name2)

if __name__ == '__main__':
	extract_hvg()
	#preprocess_data(full_data_source)
	#preprocess_data(tiny_data_source)