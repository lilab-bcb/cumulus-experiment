import scCloud
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from lightgbm import LGBMClassifier

data_location = "ground.h5ad"

adata = scCloud.tools.read_input(data_location, mode = 'a')

idx = np.isin(adata.obs['leiden_labels'], ['1', '4', '5', '6', '8'])
X_train, X_test, y_train, y_test = train_test_split(adata.X[idx], adata.obs.loc[idx, 'leiden_labels'], test_size = 0.1, random_state = 0, stratify = adata.obs.loc[idx, 'leiden_labels'])
lgb = LGBMClassifier(n_jobs = 8, metric = 'multi_error', importance_type = 'gain')
lgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], early_stopping_rounds = 1)

idx10 = np.isin(adata.obs['leiden_labels'], '10')
res = lgb.predict(adata.X[idx10])

adata.obs['new_labels'] = adata.obs['leiden_labels'].astype('object')
adata.obs.loc[idx10, 'new_labels'] = res
adata.obs['new_labels'] = pd.Categorical(adata.obs['new_labels'], categories = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '17', '18', '19'])

adata.obs['cell_types'] = adata.obs['new_labels'].apply(lambda s: '3' if s == '13' else s)
df_celltype = adata.obs[['cell_types']]
df_celltype.reset_index(inplace = True)
df_celltype.rename(columns = {'index': 'cell_id'}, inplace = True)

df_celltype.to_csv("../ground_cell_type.txt", index = False)