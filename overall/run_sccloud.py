import os, sys, time
from termcolor import cprint

n_cores = os.cpu_count()
output_dir = "sccloud_output"
out_name = "MantonBM_nonmix_sccloud"
data_src = "../MantonBM_nonmix_10x.h5"

if output_dir not in os.listdir('.'):
	if os.system("mkdir {}".format(output_dir)):
		sys.exit(1)

start_sccloud = time.time()
if os.system("scCloud cluster -p {jobs} --correct-batch-effect --diffmap-full-speed --run-approximated-leiden --run-fitsne --net-ds-full-speed --run-net-umap --run-net-fle {src_name} {out_name}".format(jobs = n_cores, src_name = data_src, out_name = output_dir + '/' + out_name)):
	sys.exit(1)
end_sccloud = time.time()
cprint("scCloud finishes in {:.2f} seconds.".format(end_sccloud - start_sccloud), "yellow")

if os.system("mv {src} ..".format(src = output_dir + '/' + out_name + '.h5ad')):
	sys.exit(1)