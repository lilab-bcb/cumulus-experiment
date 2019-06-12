#!/usr/bin/env python

import numpy as np
import scanpy as sc
import os, sys
import time
from termcolor import cprint
from mnnpy import mnn_correct
from utils import str_time

f = open("./scanpy_mnn.log", "w")

batch_folder = "/Users/yy939/Documents/manton-bm/raw_batch"

f_list = os.listdir(batch_folder)
ad_list = []
cnt = 0
for f_name in f_list:
	cnt += 1
	cprint("{seq} Processing {file_name}...".format(seq = cnt, file_name = f_name), "yellow")
	
	adata = sc.read_10x_h5(batch_folder + '/' + f_name, genome='GRCh38')
	n_sample = adata.shape[0]
	batch_name = f_name.split('.')[0]
	adata.obs['batch'] = [batch_name] * n_sample
	adata.var_names_make_unique()
	ad_list.append(adata)

	if cnt == 2: break


print("Getting Highly Variable Sets:")
bdata = sc.read_h5ad('../MantonBM_nonmix_corrected.h5ad')
hvg_list = list(bdata.var[bdata.var['selected'] == True].index.values)
print(hvg_list[0:5])

print("Batch Correction Using MNN:")
start_mnn = time.time()
corrected = mnn_correct(*ad_list, var_subset = hvg_list, n_jobs = 8)
end_mnn = time.time()
print("Time spent = " + str_time(end_mnn - start_mnn))
f.write("Time spent for MNN = " + str_time(end_mnn - start_mnn) + '\n')

adata = corrected[0]
adata.write("MantonBM_nonmix_scanpy_mnn.h5ad")

f.close()