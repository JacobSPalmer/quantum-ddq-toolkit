# You will need to install the sinter, stim, matplotlib, pymatching, and numpy packages!
# Check the STIM Getting Started guide for specific package version dependencies (as there are some specific versions of these that STIM requires)
import stim
import sinter
import pymatching
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import re

from adjustText import adjust_text
from datetime import datetime

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'font.size': 12})

# These should be default packages!
import os
from typing import List
from pprint import pprint as pp
from collections import OrderedDict

# This python script contains some modified IBM Repition Surface Code that allows you to specify a (at least for the moment) a single qubit and individually set it's defect rate
import ddq_circuit_generators as scg

## Helper Functions ##
def round_to_sig_figs(num, sig_figs):
    if num == 0: 
        return 0
    factor = 10**(sig_figs - 1 - int(math.floor(math.log10(abs(num)))))
    return round(num * factor) / factor

def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)
    
    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors

def get_true_stats_from_i2e(i2e):        
    mu = np.mean(list(i2e.values()))
    sigma = np.std(list(i2e.values()))
    med = np.median(list(i2e.values()))
    return mu, sigma, med

def get_coord_for_qubit_type(distance, type = 'center data'):
    #TODO - add type options for various spaces (i.e. edge data qubit or edge measure qubit)
    match type.lower():
        case 'edge data': 
            return (1, distance)
        case 'edge measure':
            return (distance-0.5, .5)
        case 'center data':
            return (int((distance + 1)/2), int((distance+1)/2))
        case 'center measure':
            return (int((distance + 1)/2) + 0.5, int((distance+1)/2) + 0.5)
        case _:
            #TODO - throw error here probably
            return None



## STIM Simulation Data Decomposers ##
def stats_to_df(stats):
    """
    Harvests the STIM simulation data into a easy-to-use DataFrame for additional analysis of simulation results and performance.
    Additional noise model DF can be created or extended as needed to reflect the setup specified in the experiment.
    """
    noise_model = stats[0].json_metadata['noise_model']
    match noise_model:
        case 'heterogeneous' | 'heterogeneous-sig-scalar':
            df = pd.DataFrame(columns=['p_sigma','p_mu', 'p_sig_scalar', 'distance', 'rounds', 'shots',  'errors', 'ler_shot','mu_out', 'sigma_out', 'ler_round', 'noise_model'])
            if noise_model == 'heterogeneous-sig-scalar':
                group_col = 'p_sig_scalar'
            else:
                group_col = 'p_sigma'
            x_col = 'mu_out'
            for s in stats:
                df.loc[len(df)] = {
                    'p_mu': s.json_metadata['p_mu'],
                    'p_sigma': s.json_metadata['p_sigma'],
                    'distance': s.json_metadata['distance'],
                    'rounds': s.json_metadata['rounds'],
                    'shots': s.shots,
                    'errors': s.errors,
                    'ler_shot': s.errors/s.shots,
                    'mu_out': s.json_metadata['mu_out'],
                    'sigma_out': s.json_metadata['sigma_out'],
                    'p_sig_scalar': s.json_metadata['p_sig_scalar'],
                    'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                    'noise_model': s.json_metadata['noise_model']
                }
        case 'homogeneous-defect':
            df = pd.DataFrame(columns=['p','p_defect','defect_coordinate','defect_type', 'distance', 'rounds', 'shots',  'errors', 'ler_shot', 'ler_round', 'noise_model'])
            group_col = 'p_defect'
            x_col = 'p'
            for s in stats:
                df.loc[len(df)] = {
                    'p': s.json_metadata['p'],
                    'p_defect': s.json_metadata['p_def'],
                    'defect_coordinate': s.json_metadata['defect_coord'],
                    'defect_type': s.json_metadata['defect_type'],
                    'distance': s.json_metadata['distance'],
                    'rounds': s.json_metadata['rounds'],
                    'shots': s.shots,
                    'errors': s.errors,
                    'ler_shot': s.errors/s.shots,
                    'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                    'noise_model': s.json_metadata['noise_model']
                }
        case 'heterogeneous-defect' | 'heterogeneous-defect-sig-scalar':
            df = pd.DataFrame(columns=['p_sigma','p_mu','p_sig_scalar','p_defect','defect_coordinate','defect_type', 'distance', 'rounds', 'shots',  'errors', 'ler_shot','mu_out', 'sigma_out', 'ler_round', 'noise_model'])
            group_col = 'p_defect'
            x_col = 'mu_out'
            for s in stats:
                df.loc[len(df)] = {
                    'p_mu': s.json_metadata['p_mu'],
                    'p_sigma': s.json_metadata['p_sigma'],
                    'p_defect': s.json_metadata['p_def'],
                    'defect_coordinate': s.json_metadata['defect_coord'],
                    'defect_type': s.json_metadata['defect_type'],
                    'distance': s.json_metadata['distance'],
                    'rounds': s.json_metadata['rounds'],
                    'shots': s.shots,
                    'errors': s.errors,
                    'ler_shot': s.errors/s.shots,
                    'mu_out': s.json_metadata['mu_out'],
                    'sigma_out': s.json_metadata['sigma_out'],
                    'p_sig_scalar': s.json_metadata['p_sig_scalar'],
                    'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                    'noise_model': s.json_metadata['noise_model']
                }
        # The default case is the homogenous uniform noise model, and contains the basic elements needed to plot performance data
        case _:
            group_col = 'distance'
            x_col = 'p'
            df = pd.DataFrame(columns=['p','distance', 'rounds', 'shots',  'errors', 'ler_shot', 'ler_round', 'noise_model', 'mu_out', 'sigma_out'])
            for s in stats:
                df.loc[len(df)] = {
                    'p': s.json_metadata['p'],
                    'distance': s.json_metadata['distance'],
                    'rounds': s.json_metadata['rounds'],
                    'shots': s.shots,
                    'errors': s.errors,
                    'ler_shot': s.errors/s.shots,
                    'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                    'noise_model': s.json_metadata['noise_model']
                }
                # df.sort_values(by=['distance','p'])
    
    df.sort_values(by=[group_col,x_col])
    return df

def concat_stats_df(df_l):
    df_concat = df_l[0].copy()
    for df in df_l[1:]:
        df_concat = pd.concat([df_concat, df], axis=0, ignore_index=True)
    return df_concat

def stats_df_to_csv(df, filename, subfolder):
    path = 'data'

    now = datetime.now()
    filename = f"{filename}-{now.strftime("%Y%m%d_%H%M")}.csv"

    subfolder_path = os.path.join(path, subfolder)
    filepath = os.path.join(path, subfolder, filename)
    

    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)

    
    if isinstance(df, list):
        concat_stats_df(df).to_csv(filepath, index=False)
    else:
        df.to_csv(filepath, index=False)
    print(f'Saved plot as {filepath}')
    
def save_pyplot_as_image(path, dpi):
    subfolder = re.search(r"^(.*)/", path).group(0)
    filename = re.search(r'[^/]+$', path).group(0)


    subfolder_path = os.path.join('data', subfolder)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)

    now = datetime.now()
    filepath = os.path.join(subfolder_path, f"{filename}-{now.strftime("%Y%m%d_%H%M")}.png")

    plt.savefig(fname=filepath, dpi=dpi, bbox_inches = 'tight')
    print(f'Saved plot as {filepath}')

def simulation_stats_from_csv(path_to_data, csv_file_name, split_df_on = None):
    df = pd.read_csv(f"{path_to_data}/{csv_file_name}")
    if split_df_on:
        df_list = []
        for s in df[split_df_on].unique():
            df_t = df.loc[df[split_df_on] == s].copy()
            df_list.append(df_t.reset_index(drop=True))
        return df_list
    else:
        return df



    