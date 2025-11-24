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

## DDQ Scripts ##
from helpers import *
import ddq_circuit_generators as scg

## These can be used as is, expanded, or used as a framework to develop your custom experiments that focus on a different aspect.
## If there is any interest in using this project for any case, please feel free to reach out with any questions! :)

## Case 1: Uniform Homogeneous (No Defect) ##
def uniform_homogeneous_noise_model_simulation_distance(distances: list, 
                                                        p_values = None,
                                                        rounds = 3,
                                                        use_d_for_rounds = False,
                                                        max_shots = 50_000):
    """
    Runs STIM simulations of the repetition code under a homogeneous noise model given the specified parameters with the option to specify defect qubit parameters.

    Args:
        distance (integer): Code distance for the generated circuits
        p_values (array): Array of physical qubit error values to generate samples over
        rounds (integer): # of complete cycles in the generated repitition code circuits
        max_shots (integer): # of times each circuit will be sampled
    """

    if not p_values:
        p_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    ## Converts the given defect type string to a coordinate appropriate to the specified code distance (i.e. d=5, defect_type = 'center-data' => def_coord = (3,3) ## 
    tasks = [
        sinter.Task(
            circuit = stim.Circuit(scg.surface_code_circuit_string(distance = d, 
                                                                    rounds= d if use_d_for_rounds else rounds, 
                                                                    p=p)),
            json_metadata={'p': p, 
                            'rounds': d if use_d_for_rounds else rounds,
                            'distance': d,
                            'noise_model': 'homogeneous-uniform'}
        )
        for p in p_values
        for d in distances
    ]

    stats: List[sinter.TaskStats] = sinter.collect(
        num_workers=os.cpu_count(),
        tasks=tasks,
        decoders=['pymatching'],
        max_shots=max_shots
    )
    return stats_to_df(stats)

## Case 2: Uniform Homogeneous with Defects ##
def homogeneous_noise_model_simulation(distance: int | list, 
                                       p_values = None, 
                                       p_defect_values = None, 
                                       defect_type = None, 
                                       rounds = 3, 
                                       use_d_for_rounds = False,
                                       max_shots = 50_000):
    """
    Runs STIM simulations of the repetition code under a homogeneous noise model given the specified parameters with the option to specify defect qubit parameters.

    Args:
        distance (integer): Code distance for the generated circuits
        p_values (array): Array of physical qubit error values to generate samples over
        p_defect_values (array): Array of physical qubit errors to generate samples over for the given defective qubit.
        defect_type (string): String representing the location on the lattice of the defective qubit. Valid strings are either 'center-data', 'center-measure', 'edge-data', 'edge-measure'. Given the defect type, the exact coordinate of qubit's changes depending on code distance and this allows for selecting the same qubit in simulations with any given distance.
        rounds (integer): # of complete cycles in the generated repitition code circuits
        max_shots (integer): # of times each circuit will be sampled
    """
    
    if not p_values:
        p_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    ## Converts the given defect type string to a coordinate appropriate to the specified code distance (i.e. d=5, defect_type = 'center-data' => def_coord = (3,3) ## 
   
    def task(d):
        def_coord = get_coord_for_qubit_type(d, defect_type)
        tasks = [
            sinter.Task(
                circuit = stim.Circuit(scg.surface_code_circuit_string(distance = d, 
                                                                        rounds= d if use_d_for_rounds else rounds, 
                                                                        p=p,
                                                                        p_def=[p_def],
                                                                        def_coord=[def_coord])),
                json_metadata={'p': p, 
                                'p_def': p_def, 
                                'rounds': d if use_d_for_rounds else rounds,
                                'distance': d,
                                'defect_coord': def_coord,
                                'defect_type': defect_type,
                                'noise_model': 'homogeneous-defect'}
            )
            for p in p_values
            for p_def in p_defect_values
        ]

        stats: List[sinter.TaskStats] = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=max_shots,
        )
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return task(distance)

## Case 3: Heterogeneous Noise With No Defects ##
def heterogeneous_noise_model_simulation(distance: int | list, 
                                         p_sigma_values, 
                                         use_p_sigma_values_as_scalars = False,
                                         p_mu_values = None, 
                                         rounds = 3,
                                         use_d_for_rounds = False,
                                         max_shots = 50_000,
                                         dist_func = None):
    """
    Runs STIM simulations of the parameterized repetition code with a heterogeneous noise model given the specified parameters.

    Args:
        distance (integer | array): Code distance for the generated circuits
        p_sigma_values (array): Array of std deviations (p_sigma) to generate circuits for. For each combination of sigma and mu values provided, a unique circuit will be generated where each qubit's physical error will be randomly sampled from a distribution created from the given mu and sigma values. 
        p_mu_values (array): Array of mean physical errors (p_mu) to create the distribution around.  If not specified, the mean physical errors will default to 20 points representing noise averages evenly distributed from 10^-3 to 10^-1 
        rounds (integer): # of complete cycles in the generated repitition code circuits
        max_shots (integer): # of times each circuit will be sampled
    """
    if not p_mu_values:
        p_mu_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    def task(d):
        tasks = []
        for p_mu in p_mu_values:
            for p_sigma in p_sigma_values:
                if use_p_sigma_values_as_scalars:
                    scalar = p_sigma
                    p_sigma = p_sigma * p_mu
                else:
                    scalar = p_sigma / p_mu
                circ_str, i2e, _ = scg.surface_code_circuit_string_with_sigma_noise(distance = d, 
                                                                                        rounds=d if use_d_for_rounds else rounds, 
                                                                                        p_mu=p_mu, 
                                                                                        p_sigma=p_sigma, 
                                                                                        return_qubit_mapping=True)
                mu_out, sigma_out, med = get_true_stats_from_i2e(i2e)
                task = sinter.Task(
                    circuit = stim.Circuit(circ_str),
                    json_metadata={'p_mu': p_mu, 
                                'p_sigma': p_sigma, 
                                'p_sig_scalar': scalar,
                                'rounds':d if use_d_for_rounds else rounds,
                                'distance':d,
                                'noise_model': 'heterogeneous' if not use_p_sigma_values_as_scalars else 'heterogeneous-sig-scalar',
                                'mu_out': mu_out,
                                'sigma_out': sigma_out
                    })
                tasks.append(task)

        stats: List[sinter.TaskStats] = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=max_shots
        )
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return task(distance)



## Case 4: Heterogeneous Noise With Defect ##
def heterogeneous_noise_model_with_defect_simulation(distance, 
                                                     scalar, 
                                                     use_p_sigma_values_as_scalars = False,
                                                     p_defect_values = None, 
                                                     defect_type = None, 
                                                     p_mu_values = None, 
                                                     rounds = 3, 
                                                     use_d_for_rounds = False,
                                                     max_shots = 50_000):
    """
    Runs STIM simulations of the parameterized repetition code with a heterogeneous noise model with defects given the specified parameters.

    Args:
        distance (integer | array): Code distance for the generated circuits
        p_sigma (int): The std deviations (p_sigma) to generate circuits for.
        p_defect_values (array): Array of physical qubit errors to generate samples over for the given defective qubit.
        defect_type (string): String representing the location on the lattice of the defective qubit. Valid strings are either 'center-data', 'center-measure', 'edge-data', 'edge-measure'. Given the defect type, the exact coordinate of qubit's changes depending on code distance and this allows for selecting the same qubit in simulations with any given distance.
        p_mu_values (array): Array of mean physical errors (p_mu) to create the distribution around.  If not specified, the mean physical errors will default to 20 points representing noise averages evenly distributed from 10^-3 to 10^-1 
        rounds (integer): # of complete cycles in the generated repitition code circuits
        max_shots (integer): # of times each circuit will be sampled
    """
    if not p_mu_values:
        p_mu_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    def task(d, rounds):
        def_coord = get_coord_for_qubit_type(d, defect_type)
        if def_coord is None:
            print(f'{defect_type} is not a valid special string location. Please use either "center data", "center measure", "edge data", "edge measure"')
            return
        
        tasks = []
        for p_def in p_defect_values:
            for p_mu in p_mu_values:
                p_sigma = scalar * p_mu
                circ_str, i2e, c2i = scg.surface_code_circuit_string_with_sigma_noise(distance = d, 
                                                                                rounds=d if use_d_for_rounds else rounds, 
                                                                                p_mu=p_mu, 
                                                                                p_sigma=p_sigma,
                                                                                p_def=[p_def] if p_def != 0.0 else None,
                                                                                def_coord=[def_coord] if p_def != 0.0 else None, 
                                                                                return_qubit_mapping=True)

                # Removes the defect from impacting the average of the sampled distribution, since this is not part of the distribution but instead a manually specified error rate 
                if p_def != 0.0: i2e.pop(c2i.get(def_coord))

                mu_out, sigma_out, med = get_true_stats_from_i2e(i2e)

                task = sinter.Task(
                    circuit = stim.Circuit(circ_str),
                    json_metadata={'p_mu': p_mu, 
                                'p_sigma': p_sigma, 
                                'p_sig_scalar': scalar,
                                'rounds':d if use_d_for_rounds else rounds,
                                'distance':d,
                                'p_def': p_def,
                                'defect_coord': def_coord,
                                'defect_type': defect_type,
                                'noise_model': 'heterogeneous-defect' if not use_p_sigma_values_as_scalars else 'heterogeneous-defect-sig-scalar',
                                'mu_out': mu_out,
                                'sigma_out': sigma_out}
                    )
                tasks.append(task)

        stats: List[sinter.TaskStats] = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=max_shots
        )
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d, rounds=d) for d in distance]
    else:
        return task(distance)