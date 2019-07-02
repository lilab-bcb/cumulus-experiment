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
	if os.system("scCloud cluster -p 8 --correct-batch-effect {src} MantonBM_nonmix_corrected".format(src = full_data_source)):
		sys.exit(1)

	# Extract highly variable gene set to file.
	cprint("Extracting highly variable gene set to file...", "green")
	adata = scCloud.tools.read_input("MantonBM_nonmix_corrected.h5ad", mode = 'a')
	df_hvg = adata.var.loc[adata.var["highly_variable_genes"]][["gene_ids"]]
	df_hvg.to_csv(hvg_file)

def preprocess_data(in_file):

	cprint("Preprocessing {}...".format(in_file), "green")

	# Filtration and Normalization.
	suffix1 = "_filter_norm"
	output_name1 = in_file + suffix1

	adata = scCloud.tools.read_input(in_file)
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata)
	scCloud.tools.log_norm(adata, norm_count = 1e5)
	
	scCloud.tools.write_output(adata, output_name1)
	
	# Furthermore, do PCA.
	suffix2 = "_filter_norm_pca"
	output_name2 = in_file + suffix2
	print("Getting Highly Variable Gene Set...")
	df_hvg = pd.read_csv(hvg_file)
	df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
	df_hvg.set_index('index', inplace = True)
	df_hvg.drop(columns = ['gene_ids'], inplace = True)
	adata.var.drop(columns = ['highly_variable_genes'], inplace = True)
	df_join = adata.var.join(df_hvg)
	df_join['highly_variable_genes'].fillna(False, inplace = True)
	adata_variable_genes = adata[:, df_join['highly_variable_genes']].copy()
	scCloud.tools.run_pca(adata_variable_genes)

	adata.obsm['X_pca'] = adata_variable_genes.obsm['X_pca']
	scCloud.tools.write_output(adata, output_name2)

if __name__ == '__main__':
	#extract_hvg()
	preprocess_data(full_data_source)
	preprocess_data(tiny_data_source)