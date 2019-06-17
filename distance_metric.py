from joblib import Parallel, delayed
from scipy.stats import entropy
from termcolor import cprint
from scCloud.tools.nearest_neighbors import calc_kBET
import numpy as np
import pandas as pd
import sys, os
import seaborn as sns
import matplotlib.pyplot as plt


def calc_JSD(P, Q):
	M = (P + Q) / 2
	return (entropy(P, M, base = 2) + entropy(Q, M, base = 2)) / 2.0

def calc_kBJSD_for_one_datapoint(pos, attr_values, knn_indices, ideal_dist):
	idx = np.append(knn_indices[pos], [pos])
	empirical_dist = pd.Series(attr_values[idx]).value_counts(normalize = True, sort = False).values
	return calc_JSD(ideal_dist, empirical_dist)

def calc_kBJSD(data, attr, knn_keyword = 'knn', K = 25, n_jobs = 1, temp_folder = None):
	assert attr in data.obs and data.obs[attr].dtype.name == 'category'

	ideal_dist = data.obs[attr].value_counts(normalize = True, sort = False).values # ideal no batch effect distribution
	nsample = data.shape[0]	
	nbatch = ideal_dist.size

	attr_values = data.obs[attr].values.copy()
	attr_values.categories = range(nbatch)

	knn_indices = data.uns[knn_keyword + '_indices'][:, 0 : K - 1]
	kBJSD_arr = np.array(Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_kBJSD_for_one_datapoint)(i, attr_values, knn_indices, ideal_dist) for i in range(nsample)))

	return (kBJSD_arr.mean(), kBJSD_arr)

def get_kDistances(adata, fn_prefix, attr, n_jobs = 1):
	# kBET
	kbet, pvalue, ac_rate = calc_kBET(adata, attr, n_jobs = n_jobs)
	cprint("kBET statistic for {filename} is {stat:.4f}, with p-value {pvalue:.4f} and accept rate {rate:.4f}.".format(filename = fn_prefix, stat = kbet, pvalue = pvalue, rate = ac_rate), "yellow")

	# kBJSD
	#kbjsd, kbjsd_arr = calc_kBJSD(adata, attr, n_jobs = n_jobs)
	#cprint("kBJSD for {filename} is {res:.4f}.".format(filename = fn_prefix, res = kbjsd), "yellow")
	#ax = sns.distplot(kbjsd_arr, hist = False, rug = True)
	#ax.get_figure().savefig("{prefix}_{attr}_kbjsd_dist.pdf".format(prefix = fn_prefix, attr = attr))
	#plt.close()

	return (kbet, pvalue, ac_rate)