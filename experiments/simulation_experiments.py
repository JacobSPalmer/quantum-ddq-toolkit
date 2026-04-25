import stim
import sinter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import time
from adjustText import adjust_text

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'font.size': 12})

# These should be default packages!
import os
from typing import List
from pprint import pprint as pp
from collections import OrderedDict

## DDQ Scripts ##
from helpers.helpers import *
import generators.rsc_z_generator as scg

## These can be used as is, expanded, or used as a framework to develop your custom experiments that focus on a different aspect.
## If there is any interest in using this project for any case, please feel free to reach out with any questions! :)

## Case 1: Uniform Homogeneous (No Defect) ##
def uniform_homogeneous_noise_model_simulation_distance(distance: int | list, 
                                                        p_values = [],
                                                        rounds = None,
                                                        rounds_lambda = None,
                                                        max_shots = 50_000,
                                                        log_progress = True):
    """
    Runs STIM simulations of the repetition code under a homogeneous noise model given the specified parameters with the option to specify defect qubit parameters.

    Args:
        distance (integer|list): Code distance for the generated circuits
        rounds (integer|none): # of complete cycles in the generated repitition code circuits. If none, each circuit will have rounds equal to the distance of circuit.
        p_values (array): Array of physical qubit error values to generate samples over
        max_shots (integer): # of times each circuit will be sampled
    """

    if not p_values:
        p_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    if rounds:
        rounds_lambda = lambda _d: rounds
    elif rounds_lambda is None:
        rounds_lambda = lambda _d: _d

    def task(d):
        t = time.time()

    ## Converts the given defect type string to a coordinate appropriate to the specified code distance (i.e. d=5, defect_type = 'center-data' => def_coord = (3,3) ## 
        tasks = [
            sinter.Task(
                circuit = stim.Circuit(scg.surface_code_circuit_string(distance = d, 
                                                                        rounds = rounds_lambda(d), 
                                                                        p = p)),
                json_metadata={'p': p, 
                                'rounds': rounds_lambda(d),
                                'distance': d,
                                'noise_model': 'homogeneous-uniform'}
            )
            for p in p_values
            for d in distance
        ]

        stats: List[sinter.TaskStats] = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=max_shots
        )
        if log_progress: print(f'Completed d={d} - Runtime: {(time.time() - t):.4f} seconds\n')
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return task(distance)
    
## Case 2: Uniform Homogeneous with Defects ##
def homogeneous_noise_model_simulation(distance: int | list, 
                                       p_values = [], 
                                       p_defect_values = [],
                                       defect_type = [],
                                       max_shots = 50_000,
                                       rounds = None,
                                       rounds_lambda = None,
                                       log_progress = True):
    """
    Runs STIM simulations of the repetition code under a homogeneous noise model given the specified parameters with the option to specify defect qubit parameters.

    Args:
        distance (integer | list): Code distance(s) for the generated circuits.
        rounds (integer | none): # of complete cycles in the generated repitition code circuits. If none, each circuit will have rounds equal to the distance of circuit.
        rounds_lambda (lambda | none): Can be used to create a lambda function for specifying rounds based on some output
        defect_type (array): String representing the location on the lattice of the defective qubit. Valid strings are either 'center-data', 'center-measure', 'edge-data', 'edge-measure'. Given the defect type, the exact coordinate of qubit's changes depending on code distance and this allows for selecting the same qubit in simulations with any given distance.
        p_defect_values (array): Array of physical qubit errors to generate samples over for the given defective qubit.
        p_values (array): Array of physical qubit error values to generate samples over
        max_shots (integer): # of times each circuit will be sampled
        log_progress (bool): If true, progress is printed after each completed simulation for a given distance
    """
    
    if not p_values:
        p_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    ## Converts the given defect type string to a coordinate appropriate to the specified code distance (i.e. d=5, defect_type = 'center-data' => def_coord = (3,3) ## 
   
    if rounds:
        rounds_lambda = lambda _d: rounds
    elif rounds_lambda is None:
        rounds_lambda = lambda _d: _d

    def task(d):
        t = time.time()
        def_coord = []

        for loc in defect_type:
            try:
                def_coord.append(get_coord_for_qubit_type(d, loc))
            except:
                print(f'{loc} is not a valid special string location. Please use either "center data", "center measure", "edge data", "edge measure"')

        tasks = [
            sinter.Task(
                circuit = stim.Circuit(scg.surface_code_circuit_string(distance = d, 
                                                                        rounds= rounds_lambda(d), 
                                                                        p=p,
                                                                        p_def=p_def,
                                                                        def_coord=def_coord)),
                json_metadata={'p': p, 
                                'p_def': p_def, 
                                'rounds': rounds_lambda(d),
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
        if log_progress: print(f'Completed d={d} - Runtime: {(time.time() - t):.4f} seconds\n')
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return [task(distance)]




## Case 3: Heterogeneous Noise With No Defects ##
def heterogeneous_noise_model_simulation(distance: int | list, 
                                         p_sigma_values: list, 
                                         p_mu_values = [], 
                                         use_p_sigma_values_as_scalars = False,
                                         rounds = None,
                                         rounds_lambda = None,
                                         max_shots = 50_000,
                                         log_progress = True,
                                         distibution_func = None):
    """
    Runs STIM simulations of the parameterized repetition code with a heterogeneous noise model given the specified parameters.

    Args:
        distance (integer | array): Code distance for the generated circuits
        rounds (integer | none): # of complete cycles in the generated repitition code circuits. If none, each circuit will have rounds equal to the distance of circuit.
        rounds_lambda (lambda | none): Can be used to create a lambda function for specifying rounds based on some output
        p_sigma_values (array): Array of values that determine the std deviation to generate circuits with. For each combination of sigma and mu values provided, a unique circuit will be generated where each qubit's physical error will be randomly sampled from a distribution created from the given mu and sigma values. If use_p_sigma_values_as_scalars is set to true, then the values provided here will instead be used as scalars that represent a percentage of of the mean noise value to use as deviation.
        p_mu_values (array): Array of mean physical errors (p_mu) to create the distribution around.  If not specified, the mean physical errors will default to 20 points representing noise averages evenly distributed from 10^-3 to 10^-1 
        use_p_sigma_values_as_scalars (bool): If set to true, the scalar will be used to set the deviation as a percentage of the mean value. Otherwise, will be used as a static std deviation across all sampled values.
        max_shots (integer): # of times each circuit will be sampled
        log_progress (bool): If true, progress is printed after each completed simulation for a given distance
        distibution_func (lambda): Optional lambda parameter to specify how the distrubtion should be created
    """
    if not p_mu_values:
        p_mu_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]

    if rounds:
        rounds_lambda = lambda _d: rounds
    elif rounds_lambda is None:
        rounds_lambda = lambda _d: _d

    def task(d):
        t = time.time()
        tasks = []

        for p_mu in p_mu_values:
            for p_sigma in p_sigma_values:
                if use_p_sigma_values_as_scalars:
                    scalar = p_sigma
                    p_sigma = p_sigma * p_mu
                else:
                    scalar = p_sigma / p_mu
                circ_str, i2e, _ = scg.surface_code_circuit_string_with_sigma_noise(distance = d, 
                                                                                        rounds=rounds_lambda(d), 
                                                                                        p_mu=p_mu, 
                                                                                        p_sigma=p_sigma, 
                                                                                        return_qubit_mapping=True)
                mu_out, sigma_out, med = get_true_stats_from_i2e(i2e)
                task = sinter.Task(
                    circuit = stim.Circuit(circ_str),
                    json_metadata={'p_mu': p_mu, 
                                'p_sigma': p_sigma, 
                                'p_sig_scalar': scalar,
                                'rounds':rounds_lambda(d),
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
        if log_progress: print(f'Completed d={d} - Runtime: {(time.time() - t):.4f} seconds\n')

        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return [task(distance)]



## Case 4: Heterogeneous Noise With Defect ##
def heterogeneous_noise_model_with_defect_simulation(distance: int, 
                                                     scalar: float, 
                                                     defect_type = [],
                                                     p_defect_values = [], 
                                                     p_mu_values = [], 
                                                     rounds = None,
                                                     rounds_lambda = None,
                                                     use_p_sigma_values_as_scalars = True,
                                                     max_shots = 50_000,
                                                     log_progress = True,
                                                     distribution_func = None):
    """
    Runs STIM simulations of the parameterized repetition code with a heterogeneous noise model with defects given the specified parameters.

    Args:
        distance (integer | array): Code distance for the generated circuits
        rounds (integer | none): # of complete cycles in the generated repitition code circuits. If none, each circuit will have rounds equal to the distance of circuit.
        rounds_lambda (lambda | none): Can be used to create a lambda function for specifying rounds based on some output
        scalar (float): Scalar determines the std deviation to generate circuits for or the percentage of the mean noise value given for p_mu_values to use as std deviation. If use_p_sigma_values_as_scalars is set to true, then this will be used as a scalar value.
        defect_type (array): Array representing the locations on the lattice of the defective qubit. Valid strings are either 'center-data', 'center-measure', 'edge-data', 'edge-measure'. Given the defect type, the exact coordinate of qubit's changes depending on code distance and this allows for selecting the same qubit in simulations with any given distance.        p_defect_values (array): Array of physical qubit errors to generate samples over for the given defective qubit.
        p_defect_values (array): Array of physical qubit errors to generate samples over for the given defective qubit.
        p_mu_values (array): Array of mean physical errors (p_mu) to create the distribution around.  If not specified, the mean physical errors will default to 20 points representing noise averages evenly distributed from 10^-3 to 10^-1 
        use_p_sigma_values_as_scalars (bool): If set to true, the scalar will be used to set the deviation as a percentage of the mean value. Otherwise, will be used as a static std deviation across all sampled values.
        max_shots (integer): # of times each circuit will be sampled
        log_progress (bool): If true, progress is printed after each completed simulation for a given distance
        distribution_func (lambda): Optional lambda parameter to specify how the distrubtion should be created. If `None`, then the distribution will default to a Gaussian distribution.
    """
    if not p_mu_values:
        p_mu_values = [round_to_sig_figs(i, 15) for i in np.logspace(-3, -1, 20).tolist()]
    if rounds:
        rounds_lambda = lambda _d: rounds
    elif rounds_lambda is None:
        rounds_lambda = lambda _d: _d

    def task(d):
        t = time.time()
        def_coord = []

        if isinstance(defect_type, str):
            def_coord.append(get_coord_for_qubit_type(d, defect_type))
        else:
            for loc in defect_type:
                try:
                    def_coord.append(get_coord_for_qubit_type(d, loc))
                except:
                    print(f'{loc} is not a valid special string location. Please use either "center data", "center measure", "edge data", "edge measure"')
        
        tasks = []
        for p_def in p_defect_values:
            for p_mu in p_mu_values:
                if use_p_sigma_values_as_scalars:
                    p_sigma = scalar * p_mu
                else:
                    p_sigma = scalar
                circ_str, i2e, c2i = scg.surface_code_circuit_string_with_sigma_noise(distance = d, 
                                                                                rounds=rounds_lambda(d), 
                                                                                p_mu=p_mu, 
                                                                                p_sigma=p_sigma,
                                                                                p_def=p_def if p_def != 0.0 else None,
                                                                                def_coord=def_coord if p_def != 0.0 else None, 
                                                                                return_qubit_mapping=True)

                # Removes the defect from impacting the average of the sampled distribution, since this is not part of the distribution but instead a manually specified error rate 
                if p_def != 0.0: 
                    for loc in def_coord: i2e.pop(c2i.get(loc))

                mu_out, sigma_out, med = get_true_stats_from_i2e(i2e)

                task = sinter.Task(
                    circuit = stim.Circuit(circ_str),
                    json_metadata={'p_mu': p_mu, 
                                'p_sigma': p_sigma, 
                                'p_sig_scalar': scalar,
                                'rounds':rounds_lambda(d),
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
        if log_progress: print(f'Completed d={d} - Runtime: {(time.time() - t):.4f} seconds\n')
        return stats_to_df(stats)

    if isinstance(distance, list):
        return [task(d) for d in distance]
    else:
        return [task(distance)]