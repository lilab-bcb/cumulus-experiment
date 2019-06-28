import scCloud
import pandas as pd

data_location = "ground_merged.h5ad"

adata = scCloud.tools.read_input(data_location, mode = 'a')
df_celltype = adata.obs[['cell_types']]
df_celltype.reset_index(inplace = True)
df_celltype.rename(columns = {'index': 'cell_id'}, inplace = True)

df_celltype.to_csv("../ground_cell_type.txt", index = False)