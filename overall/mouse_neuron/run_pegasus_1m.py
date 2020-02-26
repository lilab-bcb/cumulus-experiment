import os, sys

n_cores = 28
out_name = "1M_neurons_pegasus_benchmark"
data_src = "/data/1M_neurons.h5"

if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle --fle-memory 128 {src_name} {out_name} > 1m_pegasus_cpu_{jobs}.log".format(jobs = n_cores, src_name = data_src, out_name = out_name)):
	sys.exit(1)

