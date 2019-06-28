import numpy as np
import pandas as pd
import os, sys
import scanpy as sc
from bbknn import bbknn

src_file = "../MantonBM_nonmix_tiny_10x_filter_norm_pca.h5ad"

print("Reading ICA (bone marrow) dataset")
adata = sc.read_h5ad(src_file)

print("Computing neighborhood graph using BBKNN...")
bbknn(adata, batch_key = 'Channel', metric = 'euclidean')

adata_variable_genes.write("scanpy_bbknn_corrected.h5ad")