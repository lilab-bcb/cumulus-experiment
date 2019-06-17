import scCloud
import pandas as pd

data_location = "MantonBM_nonmix_corrected.h5ad"

adata = scCloud.tools.read_input(data_location, mode = 'a')
df_hvg = adata.var.loc[adata.var["highly_variable_genes"]][["gene_ids"]]
df_hvg.to_csv("hvg.txt")
