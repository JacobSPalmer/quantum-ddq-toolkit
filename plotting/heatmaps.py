# Check the STIM Getting Started guide for specific package version dependencies (as there are some specific versions of these that STIM requires)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns 

from generators import rsc_z_generator as scg

def create_heatmap(distance, i2e, title = None, def_coords = None, p_def = None, heat_range = [0.0, 0.05], fig_size = (10, 8), dpi = 120):
    _, c2i = scg.coord_circuit(distance)
    if def_coords:
        if not p_def:
            print(f"Must provide an physical error rate for specified defect coords at {def_coords}")
        i2e = scg.map_defects(c2i, i2e, def_coords, p_def)
    fig, ax = plt.subplots(figsize=fig_size)
    circuit_error_heatmap(distance ,i2e, ax, heat_range=heat_range)
    if title:
        ax.set_title(title) 
    fig.set_dpi(dpi)

def plot_lattice_errors_ax(i2e, distance, mu, sigma, ax):
    samples = i2e.values()
    x = np.linspace(min(samples), max(samples), num=len(samples))
    pdf = stats.norm.pdf(x, mu, sigma)
    ax.plot(x, pdf, 'r-', lw=2)
    ax.hist(samples, bins=15, alpha=0.7, density=True)
    ax.set_xlabel(r"Physical Error (p)")
    ax.set_ylabel("Frequency")

def c2e_to_df(c2i, i2e):
    df_c = pd.DataFrame({'index': list(c2i.values()),'coord': list(c2i.keys()) })
    df_e = pd.DataFrame({'index': list(i2e.keys()),'error': list(i2e.values()) })

    df = pd.merge(df_e, df_c, on='index')
    df['x'] = df['coord'].apply(lambda x: x[0])
    df['y'] = df['coord'].apply(lambda x: x[1])
    return df

def complete_range_for_ce_df(d, df):
    points = [i/10 for i in range(5, (d+1)*10, 5)]
    all_coords = pd.DataFrame(columns=['coord'])
    
    complete = []
    for x in points:
        for y in points:
            complete.append((x, y))
    all_coords['coord'] = complete
    all_coords['x'] = all_coords['coord'].apply(lambda x: x[0])
    all_coords['y'] = all_coords['coord'].apply(lambda x: x[1])

    merged_df = pd.merge(all_coords, df, on=['x', 'y', 'coord'], how='left')
    return merged_df

def circuit_error_heatmap(d, i2e, ax, heat_range = None, colormap='rocket_r'):
    _, c2i= scg.coord_circuit(d)
    df_qubits = c2e_to_df(c2i, i2e)
    df_comp_coords = complete_range_for_ce_df(d, c2e_to_df(c2i, i2e))
    table = df_comp_coords.pivot(index='x', columns='y', values='error')

    if heat_range is None:
        heat_range = [df_comp_coords['error'].min(),df_comp_coords['error'].max()]
        # print(f"Heat Range:{heat_range}")
    g = sns.heatmap(table, cmap=colormap, ax = ax, vmin=heat_range[0], vmax=heat_range[1])

    g.xaxis.set_ticks_position("top")
    g.xaxis.set_label_position("top")
    g.yaxis.set_label_text('Y')
    g.xaxis.set_label_text('X')

    ax.collections[0].colorbar.set_label("Physical Error (p)")

    data_qubits_locs = df_qubits[df_qubits['x'].apply(lambda x: x == int(x))]
    measure_qubit_locs = df_qubits[df_qubits['x'].apply(lambda x: x != int(x))]
    ax.scatter((data_qubits_locs['y']*2)-.5, (data_qubits_locs['x']*2)-.5, 
            c='gray',
            marker='+',
            s=200)
    ax.scatter((measure_qubit_locs['y']*2)-.5, (measure_qubit_locs['x']*2)-.5, 
            c='gray',
            marker='x',
            s=120)
    ax.scatter((data_qubits_locs['y']*2)-.5, (data_qubits_locs['x']*2)-.5, 
               edgecolor='gray',
               facecolor='w',
               s=50)
    ax.scatter((measure_qubit_locs['y']*2)-.5, (measure_qubit_locs['x']*2)-.5, 
               edgecolor='gray',
               facecolor='w',
               s=50)
    
    ax.grid(which='both',linewidth=1, alpha=0.1)