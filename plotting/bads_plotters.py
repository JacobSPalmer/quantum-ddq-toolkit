import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'font.size': 12})

from helpers.helpers import *

def plot_BADs(title, 
              df: pd.DataFrame, 
              x_col = 'BAD', 
              y_col='distance', 
              group_col = None, 
              bad_title = 'p',
              y_axis_title = None, 
              legend_title ='group',
              y_range = None,
              x_range = None, 
              show_shaded_region = False,
              figsize=(8, 6),
              dpi = 120,
              reverse=False):
    fig, ax = plt.subplots(1,1, figsize=figsize)
    fig.set_dpi(dpi)

    if reverse:
        y = y_col
        y_col = x_col
        x_col = y

    if group_col:
        groups = df[group_col].unique()

    colors = iter(cm.tab10(np.linspace(0, 1, 10)))
    shapes = iter(['o','^','s','p', 'h','o','d','+', 'X'])

    if group_col:
        groups = []
        for name, group in df.groupby('line'):
            c = next(colors)
            groups.append(group)
            ax.plot(group[x_col], group[y_col], label=name,color=c, marker= f'{next(shapes)}')
            if show_shaded_region: ax.fill_between(group[x_col], group[y_col], y2=0, alpha=0.3, color=c)
    else:
        ax.plot(df[x_col], df[y_col], marker= f'{next(shapes)}')
        if show_shaded_region: ax.fill_between(df[x_col],df[y_col], y2=0, alpha=0.2)


    ax.set_title(title)

    y_min = df[y_col].min() - df[y_col].min()*.1 if not y_range else y_range[0]
    y_max = df[y_col].max() + df[y_col].max()*.1 if not y_range else y_range[1]
    x_min = df[x_col].min() - df[x_col].max()*.1 if not x_range else x_range[0]
    x_max = df[x_col].max() + df[x_col].max()*.1 if not x_range else x_range[1]

    if reverse:
        ax.grid(which='major', axis='both')
        ax.set_title(title)

        ax.set_ylabel('BAD(' + bad_title +')')
        ax.set_xlabel(y_col.title())

        ax.set_xticks(df[x_col].unique().tolist())
        ax.set_ylim(y_min)
    else:
        ax.grid(which='major', axis='y')
        ax.set_title(title)

        ax.set_xlabel('BAD(' + bad_title +')')
        ax.set_ylabel(y_axis_title if y_axis_title else y_col.title())

        ax.set_yticks(df[y_col].unique().tolist())
        ax.set_ylim(y_min)

    if group_col:
        ax.legend(title=legend_title, loc='best', fontsize='small')
    ax.plot()


