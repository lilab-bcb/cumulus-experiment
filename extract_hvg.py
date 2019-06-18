import scCloud
import pandas as pd

def main():
	data_location = "MantonBM_nonmix_corrected.h5ad"
	output_name = "hvg.txt"

	adata = scCloud.tools.read_input(data_location, mode = 'a')
	df_hvg = adata.var.loc[adata.var["highly_variable_genes"]][["gene_ids"]]
	df_hvg.to_csv(output_name)

if __name__ == '__main__':
	main()