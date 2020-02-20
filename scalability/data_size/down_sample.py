import os, re, sys
import pegasus as pg
import numpy as np

n_cpu = os.cpu_count()
#data_src = "/data/MantonBM_nonmix.h5sc"
data_src = "MantonBM_nonmix.h5sc"

size_list = [5000, 10000, 25000, 50000, 100000, 200000]

def down_sample(data, size):
    n_obs = data.shape[0]
    indices = np.random.choice(n_obs, size, replace = False)

    ds_data = data[indices, :].copy()
    return ds_data


if __name__ == '__main__':
    f_list = [f for f in os.listdir('.') if re.match('MantonBM_.*_ds.h5ad', f)]
    if len(f_list) < 6:

        adata = pg.read_input(data_src)
        pg.qc_metrics(adata)
        pg.filter_data(adata)

        for size in size_list:
            print("Generating a random sample of size {}:".format(size))
            bdata = down_sample(adata, size)
            pg.write_output(bdata, "MantonBM_{}_ds.h5ad".format(size))

    cnt = 1
    for size in size_list:
        print("{idx}. Experiment on subsample of MantonBM dataset of size {size}:".format(idx = cnt, size = size))

        if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --diffmap --spectral-louvain --spectral-leiden --fitsne --net-umap --net-fle {src_name} {out_name} > ds_{size}.log".format(jobs = n_cpu, src_name = "MantonBM_{}_ds.h5ad".format(size), out_name = "MantonBM_{}_ds_result.h5ad".format(size), size = size)):
            sys.exit(1)
        
        cnt += 1