import time
import scanpy as sc
from datetime import timedelta

def record_time(duration_time, step_name, fp):
    fp.write("Time spent for {step}: {time}.\n\n".format(step = step_name, time = timedelta(seconds = duration_time)))

sc.settings.verbosity = 4
rand_seed = 0

log_file = "1M_neurons_scanpy.log".format(sc.settings.n_jobs)
fp = open(log_file, 'a')

adata = sc.read_h5ad("scanpy_neuron_result_no_fle.h5ad")
print("Computing FLE embedding")
start_fle = time.time()
sc.tl.draw_graph(adata, random_state = rand_seed, iterations = 5000)
end_fle = time.time()
record_time(end_fle - start_fle, "FLE", fp)

adata.write("1M_neurons_scanpy_result.h5ad")

fp.close()