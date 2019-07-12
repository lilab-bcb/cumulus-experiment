import os, sys, time
from termcolor import cprint

n_cores = os.cpu_count()
output_dir = "sccloud_output"
data_src = "../MantonBM_nonmix_10x.h5"

if output_dir not in os.listdir('.'):
	if os.system("mkdir {}".format(output_dir)):
		sys.exit(1)

start_sccloud = time.time()
if os.system("scCloud cluster -p {jobs} --correct-batch-effect --diffmap-full-speed --run-louvain --run-leiden --run-approximated-louvain --run-approximated-leiden --run-tsne --run-fitsne --run-umap --run-fle --net-ds-full-speed --run-net-tsne --run-net-fitsne --run-net-umap --run-net-fle {src_name} {out_name}".format(jobs = n_cores, src_name = data_src, out_name = output_dir + '/MantonBM_nonmix_sccloud')):
	sys.exit(1)
end_sccloud = time.time()
cprint("scCloud finishes in {:.2f} seconds.".format(end_sccloud - start_sccloud), "yellow")