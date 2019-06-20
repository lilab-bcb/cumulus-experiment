#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
import os, sys, re
import time
from termcolor import cprint
from mnnpy import mnn_correct
from utils import str_time

sc.settings.n_jobs = 8
rand_seed = 0

src_file = "../MantonBM_nonmix_tiny_10x_filter_norm.h5ad"
hvg_file = "../hvg.txt"
log_file = "./scanpy_mnn.log"


def calc_mnn():

	f = open(log_file, "w")	

	adata = sc.read_h5ad(src_file)

	# Split raw data into matrices due to Channel.
	print("Splitting raw data into Channels...")
	idx_dict = {}
	for i in range(adata.shape[0]):
		chan = adata.obs.iloc[i].name.split('-')[0]
		if chan in idx_dict:
			idx_dict[chan].append(i)
		else:
			idx_dict[chan] = [i]

	ad_dict = {}
	for (chan, indices) in idx_dict.items():
		ad_dict[chan] = adata[indices, :]


	# Get highly variable genes.
	df_hvg = pd.read_csv(hvg_file)
	cprint("Highly variable gene set of size " + str(df_hvg.shape[0]), "yellow")
	df_hvg.set_index('index', inplace = True)
	df_hvg.drop(columns = ['gene_ids'], inplace = True)
	df_hvg['highly_variable_genes'] = [True] * df_hvg.shape[0]

	cnt = 0

	for (chan, adata) in ad_dict.items():
		cnt += 1
		cprint("{seq}. Processing {chan_name}...".format(seq = cnt, chan_name = chan), "yellow")
	
		print("Filtering genes using highly variable gene set")
		adata.var.drop(columns = ['highly_variable_genes'], inplace = True)
		df_join = adata.var.join(df_hvg)
		df_join['highly_variable_genes'].fillna(False, inplace = True)
		adata_variable_genes = adata[:, df_join['highly_variable_genes']].copy()
		adata_variable_genes.var.drop(columns = ['n_cells', 'percent_cells', 'robust'], inplace = True)

		ad_dict[chan] = adata_variable_genes
		print("Dataset of shape: " + str(ad_dict[chan].shape))


	print("Batch Correction Using MNN:")
	start_mnn = time.time()
	corrected = mnn_correct(*list(ad_dict.values()), n_jobs = sc.settings.n_jobs, batch_key = 'Channel')
	end_mnn = time.time()
	print("Time spent = " + str_time(end_mnn - start_mnn))
	f.write("Time spent for MNN = " + str_time(end_mnn - start_mnn) + '\n')	

	adata = corrected[0]
	adata.write("scanpy_mnn_corrected.h5ad")

	print("Corrected data are saved.")
	f.write("Corrected data are saved." + "\n")
	f.close()


def main():
	calc_mnn()

if __name__ == "__main__":
	main()