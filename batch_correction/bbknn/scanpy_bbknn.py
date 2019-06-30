import numpy as np
import pandas as pd
import os, sys, time
import scanpy as sc
from bbknn import bbknn

sc.settings.n_jobs = os.cpu_count()
print("Using {} cores if possible.".format(sc.settings.n_jobs))

src_file = "../MantonBM_nonmix_tiny_10x_filter_norm_pca.h5ad"

print("Reading ICA (bone marrow) dataset")
adata = sc.read_h5ad(src_file)

print("Computing neighborhood graph using BBKNN...")
start_correction = time.time()
bbknn(adata, batch_key = 'Channel', metric = 'euclidean')
end_correction = time.time()
print("BBKNN Time = {:.4f} s.".format(end_correction - start_correction))

adata.write("scanpy_bbknn_corrected.h5ad")