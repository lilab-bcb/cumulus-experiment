import sccloud as scc
import numpy as np
import pandas as pd
import os, sys
from sklearn.model_selection import train_test_split
from lightgbm import LGBMClassifier

data_location = "ground.h5ad"

adata = scc.read_input(data_location)

idx = np.isin(adata.obs['louvain_labels'], ['1', '4', '5', '6', '8'])
X_train, X_test, y_train, y_test = train_test_split(adata.X[idx], adata.obs.loc[idx, 'louvain_labels'], test_size = 0.1, random_state = 0, stratify = adata.obs.loc[idx, 'louvain_labels'])
lgb = LGBMClassifier(n_jobs = 8, metric = 'multi_error', importance_type = 'gain')
lgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], early_stopping_rounds = 1)

idx_unknown = np.isin(adata.obs['louvain_labels'], '10')
res = lgb.predict(adata.X[idx_unknown])

adata.obs['cell_types'] = adata.obs['louvain_labels'].astype('object')
adata.obs.loc[idx_unknown, 'cell_types'] = res
adata.obs['cell_types'] = pd.Categorical(adata.obs['cell_types'], categories = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '17', '18', '19'])
adata.obs['cell_types'] = adata.obs['cell_types'].apply(lambda s: '2' if s == '13' else s)
adata.obs['cell_types'] = pd.Categorical(adata.obs['cell_types']).rename_categories({'11':'10', '12':'11', '14':'12', '15':'13', '16':'14', '17':'15', '18':'16', '19':'17'})
scc.write_output(adata, "ground_cell_types")

df_celltype = adata.obs[['cell_types']]
df_celltype.reset_index(inplace = True)
df_celltype.rename(columns = {'index': 'cell_id'}, inplace = True)

df_celltype.to_csv("../ground_cell_type.txt", index = False)

if os.system("sccloud plot scatter --basis umap --attributes louvain_labels,cell_types {src} {outname}.celltype.umap.pdf".format(src = "ground_cell_types.h5ad", outname = os.path.splitext(data_location)[0])):
	sys.exit(1)