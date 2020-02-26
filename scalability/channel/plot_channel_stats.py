import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def hms_to_time(s, unit):
    assert unit in ['hour', 'second']

    time_list = list(map(lambda a: a.strip(), s.split(' ')))
    assert len(time_list) == 2

    minutes = 0.0
    for t in time_list:
        minutes = 60 * minutes + float(t[:-1])

    if unit == 'hour':
        return minutes / 60
    else:
        return minutes * 60

def get_runtime(runtime_stat):
    n_channels = runtime_stat.size
    runtime_l = []

    for n in np.arange(n_channels) + 1:
        res = np.max(runtime_stat[0:n])
        runtime_l.append(res)

    return runtime_l

def get_amortized_cost(runtime_arr, total_cost):
    unit_cost = total_cost / np.sum(runtime_arr)
    cost_per_channel = list(map(lambda n: unit_cost * n, runtime_arr))

    return cost_per_channel

def get_total_cost(cost_arr):
    cost_l = []
    acc = 0.0

    for cost in cost_arr:
        acc += cost
        cost_l.append(acc)

    return cost_l

def generate_plots(df):
    fig, axes = plt.subplots(1, 2, figsize = (10, 6))

    # Plot runtime against number of channels.
    axes[0].plot('Channels', 'Runtime', 'k.', data = df, markersize = 5)
    axes[0].set_ylim(0.5, 12.5)
    axes[0].set(xlabel = "10X channel", ylabel = "Runtime in hours")

    # Plot cost against number of channels.
    axes[1].plot('Channels', 'Cost', 'b', data = df)
    axes[1].plot('Channels', 'Cost', 'k.', data = df, markersize = 5)
    axes[1].set(xlabel = "Number of channels", ylabel = "Amortized total cost in US dollars")
    
    fig.tight_layout(pad = 5.0)
    fig.savefig("channel_stats.pdf")

if __name__ == '__main__':
    df = pd.read_csv("channel_stats.csv")
    n_channels = df.shape[0]
    runtime_list = df['Time'].apply(lambda s: hms_to_time(s, unit = 'hour')).values

    seconds_per_channel = df['Time'].apply(lambda s: hms_to_time(s, unit = 'second')).values
    df['Cost'] = get_amortized_cost(seconds_per_channel, total_cost = 101.43)
    cost_list = get_total_cost(df['Cost'].values)

    df['Cost'] = np.round(df['Cost'].values, decimals = 2)
    df.to_csv("channel_stats_updated.csv", index = False)

    df_plot = pd.DataFrame({'Channels': np.arange(n_channels) + 1, 'Runtime': runtime_list, 'Cost': cost_list})

    generate_plots(df_plot)