import os, sys

n_cores = os.cpu_count()
out_name = "pbmc_benchmark"
data_src = "/data/5k_pbmc_v3.h5"

if os.system("pegasus cluster -p {jobs} --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle {src_name} {out_name} > pegasus_pbmc.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

