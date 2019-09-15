import os, sys

n_cores = os.cpu_count()
out_name = "1M_neurons_benchmark"
data_src = "/projects/benchmark/1M_Neurons/1M_neurons.h5"

if os.system("sccloud cluster -p {jobs} --correct-batch-effect --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle --fle-memory 128 {src_name} {out_name} > 1m_sccloud.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

