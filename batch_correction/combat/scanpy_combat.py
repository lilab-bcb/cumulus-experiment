#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
from bbknn import bbknn
from utils import str_time

sc.settings.n_jobs = 8
rand_seed = 0

src_file = "../MantonBM_nonmix_tiny_10x_filter_norm.h5ad"
hvg_file = "../hvg.txt"
log_file = "./scanpy_mnn.log"

f = open(log_file, "w")

print("Reading ICA (bone marrow) tiny dataset")
start = time.time()
adata = sc.read_h5ad(src_file)

print("Getting Highly Variable Genes from scCloud")
df_hvg = pd.read_csv(hvg_file)
print("Highly variable gene set of size " + str(df_hvg.shape[0]))
df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]
df_hvg.set_index('index', inplace = True)
df_hvg.drop(columns = ['gene_ids'], inplace = True)
adata.var.drop(columns = ['highly_variable_genes'], inplace = True)
df_join = adata.var.join(df_hvg)
df_join['highly_variable_genes'].fillna(False, inplace = True)
adata_variable_genes = adata[:, df_join['highly_variable_genes']].copy()


print("Combat to correct batch effects")
start_combat = time.time()
sc.pp.combat(adata_variable_genes, key = 'Channel')
end_combat = time.time()
print("Time spent for combat = " + str_time(end_combat - start_combat) + ".")
f.write("Time spent for combat = " + str_time(end_combat - start_combat) + ".\n")

print("Save corrected data...")
adata_variable_genes.write("scanpy_combat_corrected.h5ad")