import numpy as np
import pandas as pd
from scipy.stats import chi2

def calc_kBET_on_dataframe_for_one_datapoint(pos, attr_values, knn_indices, ideal_dist, K):
	indices = knn_indices[pos]
	dof = ideal_dist.size - 1

	observed_counts = pd.Series(attr_values[indices]).value_counts(sort = False).values
	expected_counts = ideal_dist * K
	static = np.sum(np.divide(np.square(np.subtract(observed_counts, expected_counts)), expected_counts))
	p_value = 1 - chi2.cdf(static, dof)
	return (static, p_value)


def calc_kBET_on_dataframe(data, nn_indices, df_cells, attr, K = 25, n_jobs = 1, temp_folder = None):
	"""
	This kBET metric is based on paper "A test metric for assessing single-cell RNA-seq batch correction" [M. Büttner, et al.] in Nature Methods, 2018.

	:return:
		static_mean: average chi-square static over all the data points.
		pvalue_mean: average p-value over all the data points.
	"""
	assert df_cells[attr].dtype.name == 'category'

	from joblib import Parallel, delayed

	ideal_dist = df_cells[attr].value_counts(normalize = True, sort = False).values # ideal no batch effect distribution
	nsample = data.shape[0]	
	nbatch = ideal_dist.size

	attr_values = df_cells[attr].values.copy()
	attr_values.categories = range(nbatch)

	knn_indices = nn_indices[:, 0 : K] - 1

	kBET_arr = np.array(Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_kBET_on_dataframe_for_one_datapoint)(i, attr_values, knn_indices, ideal_dist, K) for i in range(nsample)))

	res = kBET_arr.mean(axis = 0)
	stat_mean = res[0]
	pvalue_mean = res[1]

	return (stat_mean, pvalue_mean)

def process_data(title, pca_file, knn_file, cells_file):
	df_pca = pd.read_csv(pca_file)

	df_knn = pd.read_csv(knn_file)

	df_cells = pd.read_csv(cells_file, header = None)
	df_cells['Channel'] = pd.Series(df_cells[0].apply(lambda s: s.split('-')[0]), dtype = 'category')
	df_cells.set_index(0, inplace = True)

	stat, pvalue = calc_kBET_on_dataframe(df_pca.values, df_knn.values, df_cells, 'Channel')
	print("For {name} dataset, Mean statistic is {stat:.4f}, and Mean p-value is {pvalue:.4f}.".format(name = title, stat = stat, pvalue = pvalue))

	return (stat, pvalue)

def main():
	stat_raw, p_raw = process_data("Raw", pca_file = "seurat_raw_pca.txt", knn_file = "seurat_raw_knn.txt", cells_file = "seurat_raw_cells.txt")
	stat_cca, p_cca = process_data("CCA-Corrected", pca_file = "seurat_cca_pca.txt", knn_file = "seurat_cca_knn.txt", cells_file = "seurat_cca_cells.txt")

	print("Mean statistic is improved by {:.2%}.".format(np.abs(stat_raw - stat_cca) / stat_raw))
	print("Mean p-value is improved by {:.2%}.".format(np.abs(p_raw - p_cca) / p_raw))

if __name__ == "__main__":
	main()