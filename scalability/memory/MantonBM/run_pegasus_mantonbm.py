import os, sys

data_src = "/data/MantonBM_nonmix.h5sc"
n_cores = 28

if os.system("pegasus cluster -p {jobs} --knn-full-speed --spectral-louvain --fitsne --net-umap {src} mantonbm_pegasus > mantonbm_pegasus.log".format(jobs = n_cores, src = data_src)):
    sys.exit(1)