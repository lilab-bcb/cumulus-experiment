import os, re, sys
import numpy as np

data_src = "/data/MantonBM_nonmix.h5sc"

size_list = [5000, 10000, 25000, 50000, 100000, 200000]
cpu_list = [8, 28]

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

def k_to_num(s):
    return int(float(s[:-1]) * 1000)


if __name__ == '__main__':
    cmd = sys.argv[1]
    assert cmd in ['sampling', 'benchmark']

    if cmd == 'sampling':
        import pegasus as pg
        adata = pg.read_input(data_src)
        pg.qc_metrics(adata)
        pg.filter_data(adata)

        for size in size_list:
            size_str = num_to_k(size)
            print("Generating a random sample of size {}:".format(size_str))
            bdata = down_sample(adata, size)
            out_name = "MantonBM_{}_ds".format(size_str)
            pg.write_output(bdata, "{}.h5ad".format(out_name))

    elif cmd == 'benchmark':
        cnt = 1
        for size in size_list:
            size_str = num_to_k(size)

            for cpu in cpu_list:
                print("{idx}. Experiment on subsample of MantonBM dataset of size {size} with {cpu} cores:".format(idx = cnt, size = size_str, cpu = cpu))

                if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --louvain --spectral-louvain --leiden --spectral-leiden --fitsne --tsne --umap --net-umap --fle --net-fle {src_name} {out_name} > ds_{size}_cpu_{jobs}.log".format(jobs = cpu, src_name = "MantonBM_{}_ds.h5ad".format(size_str), out_name = "pegasus_{size}_ds_cpu_{cpu}_result".format(size = size_str, cpu = cpu), size = size_str)):
                    sys.exit(1)

                cnt += 1

        for cpu in cpu_list:
            print("{idx}. Experiment on the full MantonBM dataset with {cpu} cores:".format(idx = cnt, cpu = cpu))

            if os.system("pegasus cluster -p {jobs} --correct-batch-effect --knn-full-speed --louvain --spectral-louvain --leiden --spectral-leiden --fitsne --tsne --umap --net-umap --fle --net-fle {src_name} {out_name} > ds_fullsize_cpu_{jobs}.log".format(jobs = cpu, src_name = "/data/MantonBM_nonmix.h5sc", out_name = "pegasus_fullsize_cpu_{}_result".format(cpu))):
                sys.exit(1)

            cnt += 1

    else:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        df = pd.read_csv("stats.csv")

        df['size_num'] = df['size'].apply(lambda s: k_to_num(s))
        df['Setting'] = df['cpu'].apply(lambda n: 'laptop' if n == 8 else 'server')

        ax = sns.lineplot(x = 'size_num', y = 'tsne', hue = 'Setting', style = 'Setting', markers = True, dashes = False, data = df)
        ax.set(xlabel = 'Number of cells in dataset', ylabel = 'Time')
        fig = ax.get_figure()
        fig.suptitle("Runtime of calculating tSNE embeddings")
        fig.savefig("test.pdf")
        plt.close()