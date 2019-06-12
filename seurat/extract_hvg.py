import scCloud
import pandas as pd

data_location = "../MantonBM_nonmix_corrected.h5ad"

adata = scCloud.tools.read_input(data_location, mode = 'a')
hvg_arr = adata.var.loc[adata.var['selected']].index.values
df_hvg = pd.DataFrame({'hvg' : hvg_arr})
df_hvg.to_csv("hvg.txt", header = False, index = False)