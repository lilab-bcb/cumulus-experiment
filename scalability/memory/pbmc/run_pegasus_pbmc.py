import os, sys

n_cores = 8
data_src = "/data/5k_pbmc_v3.h5"

if os.system("pegasus cluster -p {jobs} --percent-mito 20 --knn-full-speed --spectral-louvain --fitsne --net-umap {src} pbmc_pegasus > pbmc_pegasus.log".format(jobs = n_cores, src = data_src)):
    sys.exit(1)