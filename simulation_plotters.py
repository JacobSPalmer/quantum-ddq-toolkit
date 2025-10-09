# You will need to install the sinter, stim, matplotlib, pymatching, and numpy packages!
# Check the STIM Getting Started guide for specific package version dependencies (as there are some specific versions of these that STIM requires)
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from adjustText import adjust_text

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'font.size': 12})

# These should be default packages!
import os
from typing import List
from pprint import pprint as pp
from collections import OrderedDict

from helpers import *
import ddq_circuit_generators as scg


def sort_legend(lables_handles, convert_to_percent = False, key=str):
    handles = lables_handles[0]
    labels = lables_handles[1]

    by_label = OrderedDict(zip(labels, handles))
    sorted_labels = sorted(by_label, key=key)
    sorted_handles = [by_label[label] for label in sorted_labels]
    
    #Turn into percentage
    if convert_to_percent:
        sorted_labels = [f'{float(label)*100:.0f}%' for label in sorted_labels]

    return (sorted_handles, sorted_labels)

def format_ax_data(df, ax, ler_by_round = True, group_samples_by=None, x_axis_col=None):
    if not group_samples_by or not x_axis_col:
        match df.iloc[0]['noise_model']:
            case 'heterogeneous':
                group_samples_by = 'p_sigma'
                x_axis_col = 'p_mu'
            case 'homogeneous':
                group_samples_by = 'p_defect'
                x_axis_col = 'p'
            case 'homogeneous-uniform':
                group_samples_by = 'distance'
                x_axis_col = 'p'
            case 'heterogeneous-defect':
                group_samples_by = 'p_defect'
                x_axis_col = 'p_mu'
    
    groups = df[group_samples_by].unique()
    g_label = group_samples_by
    row_lambda = lambda g: df[df[g_label] == g].sort_values(by=x_axis_col)
    x_lambda = lambda rows: rows[x_axis_col]

    y_lambda = lambda rows: rows['ler_round'] if ler_by_round else rows['ler_shot']
    
    groups = df[g_label].unique()
    colors = cm.tab10(np.linspace(0, 1, 10))
    shapes = iter(['o','^','s','p', 'h','o','d','+', 'X'])

    x_data_arr = []
    y_data_arr = []
    group_data_arr = []
    for g, c, s in zip(np.sort(groups), colors, shapes):
        rows = row_lambda(g)
        x_data = x_lambda(rows)
        y_data = y_lambda(rows)

        x_data_arr.append(x_data.tolist())
        y_data_arr.append(y_data.tolist())
        group_data_arr.append(g)
        if g_label == 'p_defect' and g == 0.0:
            ax.plot(x_data,y_data,c=[*c[:3], .6], marker=f"{s}", label=r'No Defect')
        else:
            ax.plot(x_data,y_data,c=[*c[:3], .6], marker=f"{s}", label=g)
    return x_data_arr, y_data_arr, group_data_arr

## Homogeneous Noise Model Plot Generators ##
def plot_simulation_data(df, use_loglog = True, ler_threshold = 0.0057, y_range = [10**-5, 0.1], scope_ler_by_rounds = True, ax = None, show_intersect=False, intersect_df = None, dpi=120):
    is_multiplot = True
    if not ax:
        fig, ax = plt.subplots(1,1)
        fig.set_dpi(dpi)
        is_multiplot = False

    x_data, y_data, group_data = format_ax_data(df, ax)
    noise_model = df.iloc[0]['noise_model']
    distance = df.loc[0]['distance']
    
    if ler_threshold:
        ax.hlines([ler_threshold], 0.0, 1.0, linestyles='dashed', label=f'{ler_threshold} logical threshold', color='gray')
    if noise_model == 'heterogeneous':
        x_variable = 'p_mu'
        legend_title = r'p$_{σ}$'
        xlabel = r'Mean Physical Error Rate (p$_{μ}$)'
        ax_title = f"Heterogeneous Noise; (d={distance})"
    elif noise_model == 'homogeneous':
        x_variable = 'p'
        legend_title = r'p$_{defect}$'
        xlabel = 'Physical Error Rate (p)'
        ax_title = f"Homogeneous Noise w/ {df.iloc[0]['defect_type'].title()} Qubit Defect; (d={distance})"
    elif noise_model == 'homogeneous-uniform':
        x_variable = 'p'
        legend_title = 'd'
        xlabel = 'Physical Error Rate (p)'
        ax_title = f"Uniform Homogeneous Noise"
    elif noise_model == 'heterogeneous-defect':
        x_variable = 'p_mu'
        legend_title = r'p$_{defect}$'
        xlabel = r'Mean Physical Error Rate (p$_{μ}$)'
        ax_title = f"Heterogeneous Noise w/ {df.iloc[0]['defect_type'].title()} Qubit Defect; (d={distance}, {r'p$_{σ}$'}={df.iloc[0]['p_sigma']})"

    if use_loglog: 
        ax.loglog()
        ax.set_xlim(min(df[x_variable]), max(df[x_variable]))
        ax.set_ylim(0.0001, top=.5) 


    ylabel = f"Logical Error Rate per {'round' if scope_ler_by_rounds else 'shot'} (ε)"

    ax.grid(which='major')
    ax.grid(which='minor')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(ax_title)
    
    if show_intersect:
        intersection_x_coords = []
        texts = []
        all_intersects_df = pd.DataFrame(columns=['distance', 'line', 'xy coord'])
        for x, y, g in zip(x_data, y_data, group_data):
            for i in range(len(x) - 1):
                x1, y1 = x[i], y[i]
                x2, y2 = x[i+1], y[i+1]

                if (y1 <= ler_threshold < y2) or (y2 <= ler_threshold < y1):
                    if y2 - y1 != 0:  # Avoid division by zero for horizontal segments
                        intersection_x = x1 + (ler_threshold - y1) * (x2 - x1) / (y2 - y1)
                        intersection_x_coords.append(intersection_x)
                        tb = ax.text(intersection_x, ler_threshold, f'{intersection_x:.5f}',
                                verticalalignment='top', horizontalalignment='left',
                                color='black', fontsize=10, alpha=0.8,
                                bbox=dict(facecolor='lightgray', edgecolor='black', linewidth=1, boxstyle='round,pad=0.5', alpha=0.5))
                        texts.append(tb)
                        ax.plot(intersection_x, ler_threshold, marker='x', markersize=6) # Mark the intersection
                        intersect_df.loc[len(intersect_df)] = {'distance': df.iloc[0]['distance'], 'line': g, 'xy coord': round_to_sig_figs(intersection_x, 4)}
                                

        if texts:
            adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5), # Add arrows to indicate original position
                expand_points=(1.2, 1.5), # Expand bounding box for point repulsion
                force_points=(0.2, 0.5), # Repel force from points
                force_text=(1.0, 1.6),
                force_static=(0.5, 0.8),
                explode_radius=(1.5),
                ax=ax, # Repel force from other texts
                lim=100)
        
    handles, labels = sort_legend(plt.gca().get_legend_handles_labels())

    ax.legend(title=legend_title, labels=labels, handles=handles, fontsize='small', loc='lower right' if use_loglog else 'best')

    if not is_multiplot:
        ax.plot()

def plot_multi_distance_simulation_data(arr_df, grid_dims = None, use_loglog = True, ler_threshold = 0.0057, x_range = [0.001, 0.005], y_range = [10**-5, 0.1], scope_ler_by_rounds = True, show_intersect=False, intersect_df = None, dpi=240):
    if grid_dims is None:
        grid_dims = [math.floor(len(arr_df)/2) + len(arr_df)%2 , 2]
    fig = plt.figure(figsize=(10*grid_dims[1], 8*grid_dims[0]))
    gs = fig.add_gridspec(grid_dims[0],grid_dims[1], hspace=0.09, wspace=0.04)
    gs.subplots().flatten()
    for df, ax in zip(arr_df, fig.get_axes()):
        plot_simulation_data(df, use_loglog, ler_threshold, y_range, scope_ler_by_rounds, ax, show_intersect, intersect_df, dpi)
    for ax in fig.get_axes():
        ax.label_outer()
    fig.set_dpi(dpi)