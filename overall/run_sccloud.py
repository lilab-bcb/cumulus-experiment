import os, sys

n_cores = os.cpu_count()
out_name = "MantonBM_nonmix_sccloud_benchmark"
data_src = "../MantonBM_nonmix_10x.h5"

if os.system("scCloud cluster -p {jobs} --correct-batch-effect --benchmark-time --diffmap-full-speed --run-approximated-louvain --run-approximated-leiden --run-fitsne --net-ds-full-speed --run-net-umap --run-net-fle {src_name} {out_name}".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

