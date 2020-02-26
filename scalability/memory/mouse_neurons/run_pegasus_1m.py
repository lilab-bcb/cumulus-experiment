import os, sys

data_src = "/data/1M_neurons.h5"
n_cores = os.cpu_count()

if os.system("pegasus cluster -p {jobs} --knn-full-speed --spectral-louvain --fitsne --net-umap {src} 1M_pegasus > 1M_pegasus.log".format(jobs = n_cores, src = data_src)):
    sys.exit(1)