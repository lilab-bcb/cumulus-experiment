import os, sys

n_cores = sys.argv[1]
print("Pegasus experiment using {} cores:".format(n_cores))
out_name = "MantonBM_nonmix_{}_cores".format(n_cores)
data_src = "/data/MantonBM_nonmix.h5sc"

if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --spectral-louvain --fitsne --net-umap {src_name} {out_name} > pegasus_{jobs}_cores.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)