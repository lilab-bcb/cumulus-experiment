import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def hms_to_hour(s):
    time_list = list(map(lambda a: a.strip(), s.split(' ')))
    assert len(time_list) == 3

    secs = 0.0
    for t in time_list:
        secs = 60 * secs + float(t[:-1])

    return secs / 3600

def get_amortized_cost(n_channels, unit_cost):
    cost_l = np.arange(n_channels) + 1
    return cost_l * unit_cost

def get_runtime(runtime_stat):
    n_channels = runtime_stat.size
    runtime_l = []

    for n in np.arange(n_channels) + 1:
        res = np.max(runtime_stat[0:n])
        runtime_l.append(res)

    return runtime_l

def generate_plots(df):
    fig, axes = plt.subplots(1, 2, figsize = (10, 6))

    # Plot runtime against number of channels.
    sns.lineplot(x = 'Channels', y = 'Runtime', data = df, ax = axes[0])
    axes[0].set_ylim(1, 9)
    axes[0].set(xlabel = "Number of channels", ylabel = "Runtime in hours")

    # Plot cost against number of channels.
    sns.lineplot(x = 'Channels', y = 'Cost', data = df, ax = axes[1])
    axes[1].set(xlabel = "Number of channels", ylabel = "Amortized Cost in US Dollars")

    fig.tight_layout(pad = 5.0)
    fig.suptitle("Total time and cost of Cumulus Cell Ranger Workflow regarding Number of 10X Channels")
    fig.savefig("channel_stats.pdf")

if __name__ == '__main__':
    df = pd.read_csv("channel_stats.csv")
    n_channels = df.shape[0]
    runtime_each_channel = df['Time'].apply(lambda s: hms_to_hour(s)).values

    cost_list = get_amortized_cost(n_channels = n_channels, unit_cost = 1.61)
    runtime_list = get_runtime(runtime_each_channel)

    df_plot = pd.DataFrame({'Channels': np.arange(n_channels) + 1, 'Runtime': runtime_list, 'Cost': cost_list})

    generate_plots(df_plot)