# You will need to install the sinter, stim, matplotlib, pymatching, and numpy packages!
# Check the STIM Getting Started guide for specific package version dependencies (as there are some specific versions of these that STIM requires)
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from adjustText import adjust_text
import re

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'font.size': 12})

# These should be default packages!
import os
from typing import List
from pprint import pprint as pp
from collections import OrderedDict

from helpers.helpers import *
import generators.rsc_z_generator as scg

def save_legend(figure_path, handles, labels, legend_title, dpi, loc='center', frameon=True):
    fig_legend = plt.figure(dpi=dpi)
    fig_legend.legend(handles=handles, labels=labels, loc=loc, frameon=frameon, title = legend_title)
    save_pyplot_as_image(path=figure_path+'-legend', fig=fig_legend, dpi=dpi, include_timestamps=False)
    plt.close(fig_legend) 

def format_defect_title(defect_types):
    if isinstance(defect_types, list):
        return ', '.join([i.title() for i in defect_types])
    else:
        return defect_types.title()

def sort_legend(lables_handles, convert_to_percent = False, key_type="float"):
    handles = lables_handles[0]
    labels = lables_handles[1]

    by_label = OrderedDict(zip(labels, handles))
    if key_type == 'str':
        sorted_labels = by_label.keys()
    else: 
        sorted_labels = sorted(by_label, key=lambda key: re.findall(r"[-+]?\d*\.?\d+", key))
    
    sorted_handles = [by_label[label] for label in sorted_labels]
    
    #Turn into percentage
    if convert_to_percent:
        sorted_labels = [f'{float(label)*100:.0f}%' for label in sorted_labels]

    return (sorted_handles, sorted_labels)
                
#TODO - remove calculation and use generate_BADs for intersection determination and use this as the method to plot the intersections
def generate_intersections(df, x_data, y_data, group_data, ler_target, intersect_df, ax):
    intersection_x_coords = []
    texts = []
    all_intersects_df = pd.DataFrame(columns=['distance', 'line', 'xy coord'])
    for x, y, g in zip(x_data, y_data, group_data):
        for i in range(len(x) - 1):
            x1, y1 = x[i], y[i]
            x2, y2 = x[i+1], y[i+1]

            if (y1 <= ler_target < y2) or (y2 <= ler_target < y1):
                if y2 - y1 != 0:  # Avoid division by zero for horizontal segments
                    intersection_x = x1 + (ler_target - y1) * (x2 - x1) / (y2 - y1)
                    intersection_x_coords.append(intersection_x)
                    tb = ax.text(intersection_x, ler_target, f'{intersection_x:.5f}',
                            verticalalignment='top', horizontalalignment='left',
                            color='black', fontsize=10, alpha=0.8,
                            bbox=dict(facecolor='lightgray', edgecolor='black', linewidth=1, boxstyle='round,pad=0.5', alpha=0.5))
                    texts.append(tb)
                    ax.plot(intersection_x, ler_target, marker='x', markersize=6) # Mark the intersection
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

def format_ax_data(df, ax, ler_by_round = True, group_samples_by_override=None, x_axis_col_overide = None):
    match df.iloc[0]['noise_model']:
        # If creating a special case (i.e. using an alpha scalar versus static deviations), you can create a special case here to specify what the plotter uses as the grouping columns (group_samples_by) and the column for the x data (x_axis_col)
        # Note - If you do create a special case that requires additional columns in the stats DF, ensure that the DF is properly setup in the helpers.py and reflected in the experiment function in the simulation_experiments.py 
        case 'homogeneous-uniform':
            group_samples_by = 'distance'
            x_axis_col = 'p'
        case 'homogeneous-defect':
            group_samples_by = 'p_defect'
            x_axis_col = 'p'
        case 'heterogeneous':
            group_samples_by = 'p_sigma'
            x_axis_col = 'mu_out'
        case 'heterogeneous-sig-scalar':
            group_samples_by = 'p_sig_scalar'
            x_axis_col = 'mu_out'
        case 'heterogeneous-defect' | 'heterogeneous-defect-sig-scalar':
            group_samples_by = 'p_defect'
            x_axis_col = 'mu_out'

    # In some cases, it's helpful to override the default plot setup for a specific noise model experiment. In this case, you can call the plotter function 
    if x_axis_col_overide:
        x_axis_col = x_axis_col_overide
    if group_samples_by_override:
        group_samples_by = group_samples_by_override
    
    groups = df[group_samples_by].unique()
    row_lambda = lambda g: df[df[group_samples_by] == g].sort_values(by=x_axis_col)
    # row_lambda = lambda g: df[df[group_samples_by] == g]
    x_lambda = lambda rows: rows[x_axis_col]
    y_lambda = lambda rows: rows['ler_round'] if ler_by_round else rows['ler_shot']

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
        if group_samples_by == 'p_defect' and g == 0.0:
            ax.plot(x_data,y_data,c=[*c[:3], .6], marker=f"{s}", label=r'No Defect')
        else:
            ax.plot(x_data,y_data,c=[*c[:3], .6], marker=f"{s}", label=g)
    return x_data_arr, y_data_arr, group_data_arr

def default_axis_by_simulation_data(df, x_col_override):
    match df.iloc[0]['noise_model']:
        case 'homogeneous-uniform':
            x_variable = 'p'
            legend_title = 'd'
            x_label = 'Physical Error Rate (p)'
            ax_title = f"Uniform Homogeneous Noise"
        case 'homogeneous-defect':
            x_variable = 'p'
            legend_title = r'p$_{defect}$'
            x_label = 'Physical Error Rate (p)'
            ax_title = f"Homogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{σ}$'
            ax_title = f"Heterogeneous Noise; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous-sig-scalar':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'α = p$_{σ}$/p$_{μ}$'
            ax_title = f"Heterogeneous Noise; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous-defect':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{defect}$'
            ax_title = f"Heterogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']}, {r'p$_{σ}$'}={df.iloc[0]['p_sigma']})"
        case 'heterogeneous-defect-sig-scalar':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{defect}$'
            ax_title = f"Heterogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']}, {r'α'}={df.iloc[0]['p_sig_scalar']})"
    
    return {'x_variable':x_variable,
            'x_label':x_label,
            'legend_title':legend_title,
            'ax_title':ax_title}


## Noise Model Plot Generators ##
def plot_simulation_data(df, 
                         use_loglog = True, 
                         ler_target = 0.005, 
                         y_range = None,
                         x_range = None,
                         scope_ler_by_rounds = True, 
                         ax = None, 
                         show_intersect=False, 
                         intersect_df = None, 
                         dpi=120, 
                         x_col_override = None,
                         group_col_override = None,
                         title_override = None,
                         y_label_override = None,
                         x_label_override = None,
                         legend_title_override = None,
                         save_figure_path = None,
                         remove_title = False,
                         remove_legend = False,
                         save_figure_legend_seperately = False,
                         include_timestamps_in_figure_file = False,
                         plot_dictionary = {}
                         ):
    
    is_multiplot = True
    if not ax:
        fig, ax = plt.subplots(1,1)
        fig.set_dpi(dpi)
        is_multiplot = False

    x_data, y_data, group_data = format_ax_data(df, ax, ler_by_round=scope_ler_by_rounds, x_axis_col_overide=x_col_override, group_samples_by_override=group_col_override)

    noise_model = df.iloc[0]['noise_model']
    distance = df.loc[0]['distance']
    
    if ler_target:
        ax.hlines([ler_target], 0.0, 1.0, linestyles='dashed', label=f'{ler_target} {r'LER$_{target}$'}', color='gray')

    y_label = f"Logical Error Rate per {'round' if scope_ler_by_rounds else 'shot'} (ε)"

    # These cases match the listed noise model for the experiment function in simulation_experiments.py
    # If you create a new noise model that does not use an existing noise model below, you must specify the x_variable, legend_title, x_label, and ax_title for this new experiment. These are used to label the graph.
    match noise_model:
        case 'homogeneous-uniform':
            x_variable = 'p'
            legend_title = 'd'
            x_label = 'Physical Error Rate (p)'
            ax_title = f"Uniform Homogeneous Noise"
        case 'homogeneous-defect':
            x_variable = 'p'
            legend_title = r'p$_{defect}$'
            x_label = 'Physical Error Rate (p)'
            ax_title = f"Homogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{σ}$'
            ax_title = f"Heterogeneous Noise; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous-sig-scalar':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'α = p$_{σ}$/p$_{μ}$'
            ax_title = f"Heterogeneous Noise; (d={distance}, r={df.iloc[0]['rounds']})"
        case 'heterogeneous-defect':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{defect}$'
            ax_title = f"Heterogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']}, {r'p$_{σ}$'}={df.iloc[0]['p_sigma']})"
        case 'heterogeneous-defect-sig-scalar':
            if x_col_override == 'p_mu':
                x_variable = 'p_mu'
                x_label = r'Mean Physical Error Rate (p$_{μ}$)'
            else:
                x_variable = 'mu_out'
                x_label = r'Actual Mean Physical Error Rate (p$_{μ-out}$)'
            legend_title = r'p$_{defect}$'
            ax_title = f"Heterogeneous Noise w/ {format_defect_title(df.iloc[0]['defect_type'])} Qubit Defect; (d={distance}, r={df.iloc[0]['rounds']}, {r'α'}={df.iloc[0]['p_sig_scalar']})"


    # TODO - clean this up (minimal working example)
    if x_col_override: x_variable = x_col_override

    if title_override: ax_title = title_override
    if y_label_override: y_label = y_label_override
    if x_label_override: x_label = x_label_override
    if legend_title_override: legend_title = legend_title_override


    if use_loglog: 
        ax.loglog()
        ax.set_xlim(min(df[x_variable]), max(df[x_variable]))
        ax.set_ylim(0.0001, top=.5) 
    
    # For scenarios where you want to specifically zoom in perhaps!
    if x_range: ax.set_xlim(x_range)
    if y_range: ax.set_ylim(y_range)
    
    ax.grid(which='major')
    ax.grid(which='minor')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if not remove_title: ax.set_title(ax_title)
    
    if show_intersect:
        generate_intersections(df, x_data, y_data, group_data, ler_target, intersect_df, ax)
        
    #TODO - this is custom handled since the keys in the homogeneous-uniform get ordered incorrectly since using ints not floats
    if noise_model == 'homogeneous-uniform' or group_col_override == 'distance':
        handles, labels = sort_legend(plt.gca().get_legend_handles_labels(), key_type='str')
    else:
        handles, labels = sort_legend(plt.gca().get_legend_handles_labels())

    #TODO simply out of hand number of parameters these plotter elements have accured but alas we persist. clean up later by doing a dictionary approach (like how matplotlib handles alot of their arguement parameters).
    ax.legend(title=legend_title, fancybox=True, framealpha=1.0, labels=labels, handles=handles, fontsize='small', loc='lower right' if use_loglog else 'best')
    if remove_legend:
        ax.get_legend().set_visible(False)

    if not is_multiplot:
        ax.plot()
        if save_figure_path:
            save_pyplot_as_image(save_figure_path, dpi, include_timestamps=include_timestamps_in_figure_file)
            if save_figure_legend_seperately:
                save_legend(save_figure_path, handles, labels, legend_title, dpi, frameon=True)


# Wrapper for plot to easily go over multiple distances
def plot_multi_distance_simulation_data(arr_df, 
                                        grid_dims = None, 
                                        fig_size_scalar = [10, 8],
                                        fig_hspace = 0.09,
                                        fig_wspace = 0.04,
                                        use_loglog = True, 
                                        ler_target = 0.005, 
                                        x_range = None, 
                                        y_range = None, 
                                        x_axis_override = None, 
                                        scope_ler_by_rounds = True, 
                                        show_intersect=False, 
                                        intersect_df = None, 
                                        dpi=120,
                                        save_figure_path = None,
                                        remove_title = False,
                                        remove_legend = False,
                                        save_figure_legend_seperately = False,
                                        include_alphanumeric_subplot_labels = True,
                                        include_timestamps_in_figure_file = False):
    if grid_dims is None:
        grid_dims = [math.floor(len(arr_df)/2) + len(arr_df)%2 , 2]

    fig = plt.figure(figsize=(fig_size_scalar[0]*grid_dims[1], fig_size_scalar[1]*grid_dims[0]))
    gs = fig.add_gridspec(grid_dims[0],grid_dims[1], hspace=fig_hspace, wspace=fig_wspace)
    gs.subplots().flatten()
    for df, ax in zip(arr_df, fig.get_axes()):
        plot_simulation_data(df, 
                             use_loglog=use_loglog, 
                             ler_target=ler_target, 
                             y_range=y_range, 
                             x_range=x_range, 
                             scope_ler_by_rounds=scope_ler_by_rounds, 
                             ax=ax, show_intersect=show_intersect, 
                             intersect_df=intersect_df, 
                             dpi=dpi,
                             remove_title=remove_title,
                             x_col_override=x_axis_override,
                             remove_legend = remove_legend,
                             save_figure_path=save_figure_path,
                             save_figure_legend_seperately = False,
                             include_timestamps_in_figure_file = include_timestamps_in_figure_file)

    labels = iter(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'])
    for ax in fig.get_axes():
        ax.label_outer()
        if include_alphanumeric_subplot_labels:
            ax.text(
                0.0, 1.0, next(labels), transform=(
                    ax.transAxes + ScaledTranslation(-20/120, +1/120, fig.dpi_scale_trans)),
                fontsize='x-large', fontweight='bold', va='bottom')
        

    fig.set_dpi(dpi)

    # If you would like to easily save the generated figure
    if save_figure_path:
        save_pyplot_as_image(save_figure_path, dpi, fig=fig, include_timestamps=include_timestamps_in_figure_file)
    if save_figure_legend_seperately:
        handles, labels = plt.gca().get_legend_handles_labels()
        title = plt.gca().get_legend().get_title().get_text()
        save_legend(figure_path=save_figure_path, dpi=dpi, handles=handles, labels=labels, legend_title=title)
