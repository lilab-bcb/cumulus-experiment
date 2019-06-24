import scCloud
import sys
import pandas as pd

hvg_file = "./hvg.txt"

def process_data(in_file, calc_pca):
	suffix = "_filter_norm_pca" if calc_pca else "_filter_norm"
	output_name = in_file.split('/')[-1].split('.')[0] + suffix

	adata = scCloud.tools.read_input(in_file)
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata)
	scCloud.tools.log_norm(adata, norm_count = 1e5)
	
	if not calc_pca:
		scCloud.tools.write_output(adata, output_name)
	else:
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
		scCloud.tools.write_output(adata_variable_genes, output_name)

if __name__ == '__main__':
	in_file = sys.argv[1]
	calc_pca = False
	if len(sys.argv) > 2:
		if sys.argv[2] == "pca":
			calc_pca = True
	process_data(in_file, calc_pca)