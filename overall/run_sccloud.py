import os, sys

n_cores = os.cpu_count()
out_name = "MantonBM_nonmix_sccloud_benchmark"
data_src = "../MantonBM_nonmix.h5sc"

if os.system("sccloud cluster -p {jobs} --correct-batch-effect --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle {src_name} {out_name} > sccloud.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

