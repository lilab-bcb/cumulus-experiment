import os, sys

n_cores = 8
out_name = "pbmc_pegasus_benchmark".format(n_cores)
data_src = "/data/5k_pbmc_v3.h5"

if os.system("pegasus cluster -p {jobs} --percent-mito 20 --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle {src_name} {out_name} > pbmc_pegasus_cpu_{jobs}.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)