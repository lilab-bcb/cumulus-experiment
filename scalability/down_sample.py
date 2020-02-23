import os, re, sys
import numpy as np
import pegasus as pg

data_src = "/data/MantonBM_nonmix.h5sc"

size_list = [5000, 10000, 25000, 50000, 100000, 200000]
cpu_list = [8, os.cpu_count()]

def down_sample(data, size):
    n_obs = data.shape[0]
    indices = np.random.choice(n_obs, size, replace = False)

    ds_data = data[indices, :].copy()
    return ds_data

def num_to_k(n):
    num_part = n / 1000
    if num_part == int(num_part):
        return str(int(num_part)) + 'k'
    else:
        return str(num_part) + 'k'


if __name__ == '__main__':
    cmd = sys.argv[1]
    assert cmd in ['sampling', 'benchmark']

    if cmd == 'sampling':
        adata = pg.read_input(data_src)
        pg.qc_metrics(adata)
        pg.filter_data(adata)

        for size in size_list:
            size_str = num_to_k(size)
            print("Generating a random sample of size {}:".format(size_str))
            bdata = down_sample(adata, size)
            out_name = "MantonBM_{}_ds".format(size_str)
            pg.write_output(bdata, "{}.h5ad".format(out_name))

    else:
        cnt = 1
        for size in size_list:
            size_str = num_to_k(size)

            for cpu in cpu_list:
                print("{idx}. Experiment on subsample of MantonBM dataset of size {size} with {cpu} cores:".format(idx = cnt, size = size_str, cpu = cpu))

                if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --louvain --spectral-louvain --leiden --spectral-leiden --fitsne --tsne --umap --net-umap --fle --net-fle {src_name} {out_name} > ds_{size}_cpu_{jobs}.log".format(jobs = cpu, src_name = "MantonBM_{}_ds.h5ad".format(size_str), out_name = "pegasus_{size}_ds_cpu_{cpu}_result.h5ad".format(size = size_str, cpu = cpu), size = size_str)):
                    sys.exit(1)

                cnt += 1

        for cpu in cpu_list:
            print("{idx}. Experiment on the full MantonBM dataset with {cpu} cores:".format(idx = cnt, cpu = cpu))

            if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --louvain --spectral-louvain --leiden --spectral-leiden --fitsne --tsne --umap --net-umap --fle --net-fle {src_name} {out_name} > ds_fullsize_cpu_{jobs}.log".format(jobs = cpu, src_name = "/data/MantonBM_nonmix.h5sc", out_name = "pegasus_fullsize_cpu_{}_result.h5ad".format(cpu))):
                sys.exit(1)

            cnt += 1
