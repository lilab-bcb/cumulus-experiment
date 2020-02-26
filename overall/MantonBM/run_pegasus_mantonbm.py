import os, sys

n_cores = 28
out_name = "mantonbm_pegasus_benchmark"
data_src = "/data/MantonBM_nonmix.h5sc"

if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle {src_name} {out_name} > mantonbm_pegasus_cpu_{jobs}.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

