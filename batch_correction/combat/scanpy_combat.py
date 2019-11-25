import numpy as np
import pandas as pd
import scanpy as sc
import sys, time, os
from datetime import timedelta

sc.settings.n_jobs = os.cpu_count()
print("Using {} cores if possible.".format(sc.settings.n_jobs))
rand_seed = 0

src_file = "/data/MantonBM_nonmix_tiny_filter_norm.h5ad"
hvf_file = "/data/MantonBM_nonmix_tiny_hvf.txt"

print("Reading ICA (bone marrow) tiny dataset")
adata = sc.read_h5ad(src_file)

print("Getting Highly Variable Features from Pegasus")
df_hvf = pd.read_csv(hvf_file)
print("Highly variable gene set of size " + str(df_hvf.shape[0]))
df_hvf['highly_variable_features'] = [True] * df_hvf.shape[0]
df_hvf.set_index('index', inplace = True)
df_hvf.drop(columns = ['gene_ids'], inplace = True)
adata.var.drop(columns = ['highly_variable_features'], inplace = True)
df_join = adata.var.join(df_hvf)
df_join['highly_variable_features'].fillna(False, inplace = True)
adata_variable_genes = adata[:, df_join['highly_variable_features']].copy()
assert(adata_variable_genes.shape[1] == df_hvf.shape[0])


print("Combat to correct batch effects")
start_correction = time.time()
sc.pp.combat(adata_variable_genes, key = 'Channel')
end_correction = time.time()
print("ComBat Time = {}.".format(timedelta(seconds = end_correction - start_correction)))

print("Save corrected data...")
adata_variable_genes.write("scanpy_combat_corrected.h5ad")