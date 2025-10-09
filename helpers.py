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
from adjustText import adjust_text

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
    Harvests the STIM simulation data into a easy-to-use DataFrame for additional analysis of simulation results and performance
    """
    noise_model = stats[0].json_metadata['noise_model']
    if noise_model == 'heterogeneous':
        df = pd.DataFrame(columns=['p_sigma','p_mu', 'distance', 'rounds', 'shots',  'errors', 'ler_shot', 'ler_round', 'noise_model'])
        group_col = 'p_sigma'
        x_col = 'p_mu'
        for s in stats:
            df.loc[len(df)] = {
                'p_mu': s.json_metadata['p_mu'],
                'p_sigma': s.json_metadata['p_sigma'],
                'distance': s.json_metadata['distance'],
                'rounds': s.json_metadata['rounds'],
                'shots': s.shots,
                'errors': s.errors,
                'ler_shot': s.errors/s.shots,
                'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                'noise_model': s.json_metadata['noise_model']
            }
        df.sort_values(by=['p_sigma', 'p_mu'])
    elif noise_model == 'homogeneous':
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
        df.sort_values(by=['p_defect','p'])
    elif noise_model == 'heterogeneous-defect':
        df = pd.DataFrame(columns=['p_sigma','p_mu','p_defect','defect_coordinate','defect_type', 'distance', 'rounds', 'shots',  'errors', 'ler_shot', 'ler_round', 'noise_model'])
        group_col = 'p_defect'
        x_col = 'p_mu'
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
                'ler_round': round_to_sig_figs(s.errors/(s.shots * s.json_metadata['rounds']), 10),
                'noise_model': s.json_metadata['noise_model']
            }
        df.sort_values(by=[group_col, x_col])
    else:
        group_col = 'distance'
        x_col = 'p'
        df = pd.DataFrame(columns=['p','distance', 'rounds', 'shots',  'errors', 'ler_shot', 'ler_round', 'noise_model'])
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
            df.sort_values(by=['distance','p'])
    return df

def concat_stats_df(df_l):
    df_concat = df_l[0].copy()
    for df in df_l[1:]:
        df_concat = pd.concat([df_concat, df], axis=0, ignore_index=True)
    return df_concat

def stats_df_to_csv(df, filename, subfolder):
    path = 'data'
    filename = f"{filename}.csv"
    subfolder_path = os.path.join(path, subfolder)
    filepath = os.path.join(path, subfolder, filename)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)
    
    if isinstance(df, list):
        return concat_stats_df(df).to_csv(filepath, index=False)
    else:
        return df.to_csv(filepath, index=False)