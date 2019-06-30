import os, sys

n_jobs = os.cpu_count()
print("Use {} cores.".format(n_jobs))

if os.system("scCloud cluster -p {cores} --correct-batch-effect --run-leiden --run-approximated-leiden --run-fitsne --run-umap ../../MantonBM_nonmix_tiny_10x.h5 tiny_sccloud_corrected".format(cores = n_jobs)):
	sys.exit(1)