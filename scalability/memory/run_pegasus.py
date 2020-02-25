import os, sys

n_cores = os.cpu_count()
data_src = sys.argv[1]
out_name = sys.argv[2]

if os.system("pegasus cluster -p {jobs} --knn-full-speed --spectral-louvain --fitsne --net-umap {src} {dst} > {dst}_pegasus.log".format(jobs = n_cores, src = data_src, dst = out_name)):
    sys.exit(1)