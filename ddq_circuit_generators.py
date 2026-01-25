import numpy as np
import math
from decimal import Decimal
import matplotlib.pyplot as plt
import scipy.stats as stats
# rng = np.random.default_rng(seed=7318151)

## Coordinates ##
def data_coords(distance):
    coords = []
    for row in range(1, distance+1):
        for col in range(1, distance+1):
            coords.append((col, row))
    return coords

def z_measure_coords(distance):
    coords = []
    for row in range(1, distance): # don't include the last row
        for col in range(1, distance+1, 2): # only take every other qubit
            if row%2: 
                coords.append((col-0.5, row+0.5))
            else:
                coords.append((col+0.5, row+0.5))
    return coords

def x_measure_coords(distance):
    coords = []
    for row in range(1, distance+2): # include extra for last row measures
        for col in range(2, distance, 2): # start from second column, ignore last
            if row%2:
                coords.append((col+0.5, row-0.5))
            else:
                coords.append((col-0.5, row-0.5))
    return coords

def coords_to_index(coords):
    return {tuple(c):i for i,c in dict(enumerate(coords)).items()}

def adjacent_coords(coord):
    col, row = coord
    adjacents = [(col-0.5, row-0.5), (col+0.5, row-0.5),
               (col-0.5, row+0.5), (col+0.5, row+0.5),
              ]
    return adjacents

## Utility/Helper Functions ##
def index_string(coord_list, c2i):
    # Returns the indicies for each coord in a list as space-delimited string.
    return ' '.join(str(c2i[coord]) for coord in coord_list)

def index_array(coord_list, c2i):
    return [c2i[coord] for coord in coord_list]

def index_arr_to_str(coord_list):
    return ' '.join(str(x) for x in coord_list)

def split_array_on_def(coord_list, c2i, def_coord):
    if def_coord not in coord_list:
        a = [index_array(coord_list, c2i)]
        return a
    defec = c2i[def_coord]
    a = index_array(coord_list, c2i)
    def_ind = a.index(defec)
    slice1 = a[:def_ind]
    slice2 = a[(def_ind+1):]

    return [slice1, [defec], slice2]

def chunk_array(coord_list, c2i, def_coord):
    def_ind = c2i[def_coord]
    a = index_array(coord_list, c2i)
    chunks = [a[i:i + 2] for i in range(0, len(a), 2)]
    index = 0
    for i, subarray in enumerate(chunks):
        if def_ind in subarray:
            index = i

    slice_1_c = sum(chunks[:index], [])
    def_slice = chunks[index]
    slice_2_c = sum(chunks[(index + 1):], [])
    return [slice_1_c, def_slice, slice_2_c]

def chunk_arr(a, size=2):
    return [a[i:i + 2] for i in range(0, len(a), 2)]

def count_sig_figs(n):
    d = Decimal(f"{n}")
    return len(d.normalize().as_tuple().digits)

def round_to_sig_figs(num, sig_figs):
    if num == 0: 
        return 0
    factor = 10**(sig_figs - 1 - int(math.floor(math.log10(abs(num)))))
    return round(num * factor) / factor

def validate_def_coords(def_coord, p_def, distance):
    if len(def_coord) != len(p_def):
        print("Number of defective coordinates does not match number of defective error rates!\n")
        return False
    for coord in def_coord:
        if coord[0] > .5+distance or coord[1] > .5+distance:
            print(f"Defective qubit coordinates ({coord}) invalid for d={distance} code!\n")
            print(f"Range is (1, 1) <= (coord x, coord y) <= ({distance + .5}, {distance + .5})\n")
            return False
        elif coord[0]%.5 != 0 or coord[1]%.5 != 0:
            print(f"Defective qubit coordinates ({coord})!\n")
            print(f"Coord x & y must be multiple of .5 (i.e., x%.5 == 0 & y%.5 == 0)\n")
            return False
    return True

## Gate String Formatters ##
def format_defec_string(operation_string, split_array, p, p_defec, is_one_qubit = True):
    # TODO - clean up
    if len(split_array) == 1:
        return f"""
        {operation_string}({p}) {index_arr_to_str(split_array[0])}
        """
    else:
        return f"""
        {operation_string}({p}) {index_arr_to_str(split_array[0])}
        {operation_string}({p_defec if is_one_qubit else p_defec*1.333}) {index_arr_to_str(split_array[1])}
        {operation_string}({p}) {index_arr_to_str(split_array[2])}
        """

def single_op_string(operation_string, coords, c2i, i2e):
    stim_string = ""
    for coord in coords:
        index = c2i[coord]
        p = i2e[index]
        stim_string += f"""
        {operation_string}({p}) {index}"""
    return stim_string

def controlled_op_string(operation_string, coords, c2i, i2e):
    stim_string = ""
    for cont_coords in chunk_arr(coords):
        index_1 = c2i[cont_coords[0]]
        index_2 = c2i[cont_coords[1]]
        p_1 = i2e[index_1]
        p_2 = i2e[index_2]
        #TODO - check if there is any frame of reference for the error rate on a 2 qubit operation???
        p_comb = round_to_sig_figs(((p_1 + p_2)/2) * 1.2, count_sig_figs(p_1))
        
        stim_string += f"""
        {operation_string}({p_comb}) {index_1} {index_2}"""
    return stim_string

## Mapping funcs for Coords, Indices, and Errors ##
def prepare_coords(distance):
    datas = data_coords(distance)
    x_measures = x_measure_coords(distance)
    z_measures = z_measure_coords(distance)
    c2i = coords_to_index(datas+x_measures+z_measures)
    return datas, x_measures, z_measures, c2i

## Prebuilt Distribution Generators ##
def prepare_error_map(distance, 
                      noise_dist_func = None, 
                      **kwargs):
    if noise_dist_func:
        return noise_dist_func(distance = distance, **kwargs)
    else: 
        return prepare_norm_dist_errors(distance = distance, **kwargs)
# Norm #
# def prepare_norm_dist_errors(distance, p_mu, p_sigma, precision=5):
#     rng = np.random.default_rng()
#     print(p_sigma)
#     return {key: round_to_sig_figs(rng.normal(loc=p_mu, scale=p_sigma), precision) for key in range(0, 2*distance**2 - 1)}

# Homogeneous Uniform Constant Error #
def prepare_constant_errors(distance, error):
    return {key: error for key in range(0, 2*distance**2 - 1)}

## Norm ##
# def prepare_norm_dist_errors(distance, mu, sigma, precision=5):
#     rng = np.random.default_rng()
#     return {key: round(rng.normal(loc=mu, scale=sigma),4) for key in range(0, 2*distance**2 - 1)}

## Truncated Norm ##
# TODO - review if this impacts any simulation in a major way. necessary to prevent negative error rates but there might be a better dist. to use entirely
def prepare_norm_dist_errors(distance, p_mu, p_sigma, precision=5):
    truncated_dist = stats.truncnorm((0.000001-p_mu)/p_sigma, (1-p_mu)/p_sigma, loc=p_mu, scale=p_sigma)
    return {key: round_to_sig_figs(truncated_dist.rvs(1)[0], precision) for key in range(0, 2*distance**2 - 1)}

def map_defects(c2i, i2e, def_coord, p_def):
    # Returns the index associated with the coordinate (simple helper)
    for c, p in zip(def_coord, p_def):
        i2e[c2i[c]] = p
    return i2e

def coord_circuit(distance):
    # Returns a Stim circuit string that adds a QUBIT_COORDS instruction for each
    #  qubit, based on the coordinate-to-index mapping.
    _, _, _, c2i = prepare_coords(distance)
    stim_circuit = ""
    for coord, index in c2i.items():
        stim_circuit += f"QUBIT_COORDS({','.join(map(str, coord))}) {index}\n"
    return stim_circuit, c2i

def label_indices(distance):
    datas, x_measures, z_measures, c2i = prepare_coords(distance)
    all_qubits = datas + x_measures + z_measures
    i = 0
    stim_string = ""
    for coord in datas:
        stim_string += f"Y_ERROR(0.{i:>02}) {c2i[coord]}\n"
        i += 1
    stim_string += "TICK\n"
    for coord in x_measures:
        stim_string += f"X_ERROR(0.{i:>02}) {c2i[coord]}\n"
        i += 1
    stim_string += "TICK\n"
    for coord in z_measures:
        stim_string += f"Z_ERROR(0.{i:>02}) {c2i[coord]}\n"
        i += 1

    return stim_string

def split_list(input_list, delimiter):
    """Splits a list into sublists based on a delimiter."""
    
    result = []
    current_sublist = []
    
    for item in input_list:
        if item == delimiter:
            result.append(current_sublist)
            current_sublist = []
        else:
            current_sublist.append(item)
    
    result.append(current_sublist)
    return result

## Format Stim String Functions ##

def lattice_with_noise(distance, i2e):
    _, x_measures, z_measures, c2i = prepare_coords(distance)

    stim_string = ""
    for i in range(4):
        cx_qubits = []
        for measure in z_measures:
            z_controls = adjacent_coords(measure)
            control = z_controls[i]
            if control in c2i:
                cx_qubits.extend([control, measure])
        
        for measure in x_measures:
            x_targets = adjacent_coords(measure)
            index_reorder = [0, 2, 1, 3]
            target = x_targets[index_reorder[i]]
            if target in c2i:
                cx_qubits.extend([measure, target])

        idle_qubits = [coord for coord in c2i.keys() if coord not in cx_qubits]
        
        stim_string += f"""
        CX {index_string(cx_qubits, c2i)}

        ## DEPOLARIZE2->cx_qubits ##
        {controlled_op_string("DEPOLARIZE2", cx_qubits, c2i, i2e)}

        ## DEPOLARIZE1->idle_qubits ##
        {single_op_string("DEPOLARIZE1", idle_qubits, c2i, i2e)}
        TICK
        """

    return stim_string

def stabilizers_with_noise(distance, i2e):
    datas, x_measures, z_measures, c2i = prepare_coords(distance)
    all_measures = x_measures + z_measures
    all_qubits = datas + all_measures
    
    stim_string = f"""
    R {index_string(all_measures, c2i)}

    ## X_ERROR->all_measures ##
    {single_op_string("X_ERROR", all_measures, c2i, i2e)}

    ## DEPOLARIZE1->datas ##
    {single_op_string("DEPOLARIZE1", datas, c2i, i2e)}

    TICK
    H {index_string(x_measures, c2i)}
    ## DEPOLARIZE1->all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}

    TICK
    """

    stim_string += lattice_with_noise(distance, i2e)
    
    stim_string += f"""
    H {index_string(x_measures, c2i)}
    
    ## DEPOLARIZE1->all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}

    TICK

    ## X_ERROR->all_measures ##
    {single_op_string("X_ERROR", all_measures, c2i, i2e)}

    ## DEPOLARIZE1->datas ##
    {single_op_string("DEPOLARIZE1", datas, c2i, i2e)}

    M {index_string(all_measures, c2i)}
    TICK
    """
    return stim_string

def initialization_step(distance, i2e):
    datas, x_measures, z_measures, c2i = prepare_coords(distance)
    all_measures = x_measures + z_measures
    all_qubits = datas + all_measures
    
    ##TODO - Double check that the depolarization should go on the all_qubits rather (as opposed to being on the x_measures)
    stim_string = f"""
    R {index_string(all_qubits, c2i)}
    ## X_ERROR -> all_qubits ##
    {single_op_string("X_ERROR", all_qubits, c2i, i2e)}
    TICK
    H {index_string(x_measures, c2i)}
    ## DEPOLARIZE1 -> all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}
    TICK
    """

    stim_string += lattice_with_noise(distance, i2e)

    
    stim_string += f"""
    H {index_string(x_measures, c2i)}

    ## DEPOLARIZE1 -> all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}
    TICK

    ## X_ERROR -> all_measures ##
    {single_op_string("X_ERROR", all_measures, c2i, i2e)}
    M {index_string(all_measures, c2i)}

    ## DEPOLARIZE1 -> datas ##
    {single_op_string("DEPOLARIZE1", datas, c2i, i2e)}
    
    TICK
    """

    for i in range(1, len(z_measures)+1):
        stim_string += f"DETECTOR({i}, 0) rec[{-i}]\n"
    return stim_string

def rounds_step(distance, rounds, i2e):
    if rounds <= 2:
        return "\n"
    datas, x_measures, z_measures, c2i = prepare_coords(distance)

    stim_string = f"REPEAT {rounds-2} {{\n"
    stim_string += stabilizers_with_noise(distance, i2e)

    num_measures_per_type = len(z_measures) # number of measures per type per round
    for i in range(1, num_measures_per_type+1): # offset to the previous round
        stim_string += f"DETECTOR({i}, 0) rec[{-i}] rec[{-(i+2*num_measures_per_type)}]\n"
    for i in range(1, num_measures_per_type+1): # offset to the other type and to the previous round
        stim_string += f"DETECTOR({i}, 0) rec[{-(i+num_measures_per_type)}] rec[{-(i+3*num_measures_per_type)}]\n"
        
    stim_string += """
    }
    """

    return stim_string

def final_step(distance, i2e):
    datas, x_measures, z_measures, c2i = prepare_coords(distance)
    all_measures = x_measures + z_measures
    all_qubits = datas + all_measures

    stim_string = f"""
    R {index_string(all_measures, c2i)}

    ## X_ERROR->all_measures ##
    {single_op_string("X_ERROR", all_measures, c2i, i2e)}

    ## DEPOLARIZE1->datas ##
    {single_op_string("DEPOLARIZE1", datas, c2i, i2e)}

    TICK
    H {index_string(x_measures, c2i)}

    ## DEPOLARIZE1->all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}
    
    TICK
    """

    stim_string += lattice_with_noise(distance, i2e)

    stim_string += f"""
    H {index_string(x_measures, c2i)}

    ## DEPOLARIZE1->all_qubits ##
    {single_op_string("DEPOLARIZE1", all_qubits, c2i, i2e)}
    
    TICK
    ## X_ERROR->all_qubits ##
    {single_op_string("X_ERROR", all_qubits, c2i, i2e)}

    M {index_string(all_qubits, c2i)}
    """

    num_measures_per_type = len(z_measures) 
    num_datas = len(datas)
    for i in range(1, num_measures_per_type+1): 
        stim_string += f"DETECTOR({i}, 0) rec[{-i}] rec[{-(i+2*num_measures_per_type+num_datas)}]\n"
    for i in range(1, num_measures_per_type+1): 
        stim_string += f"DETECTOR({i}, 0) rec[{-(i+num_measures_per_type)}] rec[{-(i+3*num_measures_per_type+num_datas)}]\n"
    
    coord_to_record_index = {coord: i-len(all_qubits) for i, coord in enumerate(all_qubits)}
    for i, measure in enumerate(z_measures):
        record_indices = []
        record_indices.append(coord_to_record_index[measure])
        adjacent_datas = adjacent_coords(measure)
        
        for data in adjacent_datas:
            if data in all_qubits:
                record_indices.append(coord_to_record_index[data])
        recs = [f"rec[{j}]" for j in record_indices]
        stim_string += f"DETECTOR({i}, 0) {' '.join(recs)}\n"
    
    obs_recs = [f"rec[{-(i+2*num_measures_per_type)}]" for i in range(1, distance+1)]
    stim_string += f"OBSERVABLE_INCLUDE(0) {' '.join(obs_recs)}"

    return stim_string

# Plot Generators
def plot_lattice_errors(i2e, distance, mu, sigma):
    samples = i2e.values()
    x = np.linspace(min(samples), max(samples), num=len(samples))
    # x = list(i2e.values())
    pdf = stats.norm.pdf(x, mu, sigma)

    plt.plot(x, pdf, 'r-', lw=2)

    plt.hist(samples, bins=15, alpha=0.7, density=True)
    plt.title("Normal Distribution for d={distance}, σ={sigma}, μ={mu}".format(distance=distance, sigma=sigma, mu=mu))
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.show()

## Generate Stim String Function ##

# p -> physical error rate for all qubits
# defects -> boolean; if true, the qubits at the specified def_coords with be overrided with the values in p_def
# p_def -> array of physical error rates for defective qubits
# defect_coords -> array of coordinates for defective qubits [(x1, y1), (x2, y2), ...], p_def[i] corresponds to defect_coords[i]
def surface_code_circuit_string(distance: int, 
                                rounds: int, 
                                p: int,
                                p_def = [],
                                def_coord = []):
    i2e = prepare_constant_errors(distance, p)

    string, c2i= coord_circuit(distance)

    if def_coord and p_def:
        # if only one p_def is given, give all defective coords the same error
        if type(p_def) is not list:
            p_def = [p_def] * len(def_coord)
        if len(p_def) == 1 and p_def[0] == 0:
            p_def[0] = p
        # Validate the specified defective coordinates (i.e., within range on the lattice for the given distance)
        if not validate_def_coords(def_coord, p_def, distance):
            return None
        else:
            map_defects(c2i, i2e, def_coord, p_def)

    #Circuit builder steps:
    string += initialization_step(distance, i2e)
    string += rounds_step(distance, rounds, i2e)
    string += final_step(distance, i2e)
    return string

# p_mu -> mean for normal dist. of physical error values
# p_sigma -> std. deviation for normal dist. of physical error values
# defects -> boolean; if true, the qubits at the specified def_coords with be overrided with the values in p_def
# p_def -> array of physical error rates for defective qubits
# defect_coords -> array of coordinates for defective qubits [(x1, y1), (x2, y2), ...], p_def[i] corresponds to defect_coords[i]
def surface_code_circuit_string_with_sigma_noise(distance: int, 
                                                 rounds: int, 
                                                 p_mu: float, 
                                                 p_sigma: float, 
                                                 p_def = [], 
                                                 def_coord = [],
                                                 noise_func = None,
                                                 show_norm_plot = False, 
                                                 return_qubit_mapping = False):
    string, c2i = coord_circuit(distance)

    # If deviation is 0, treat it like a homogeneous noise model. Helpful for gathering baselines to no deviation.
    if p_sigma == 0.0:
        if return_qubit_mapping:
            i2e = prepare_constant_errors(distance, p_mu)
            return surface_code_circuit_string(distance, rounds, p = p_mu, p_def= p_def, def_coord = def_coord), i2e, c2i
        else:
            return surface_code_circuit_string(distance, rounds, p = p_mu, p_def= p_def, def_coord = def_coord)

    i2e = prepare_error_map(distance=distance, noise_dist_func=noise_func, p_mu = p_mu, p_sigma=p_sigma)

    if def_coord and p_def:
        # if only one p_def is given, give all defective coords the same error
        if type(p_def) is not list:
            p_def = [p_def] * len(def_coord)
        # Validate the specified defective coordinates (i.e., within range on the lattice for the given distance)
        if not validate_def_coords(def_coord, p_def, distance):
            return None
        else:
            map_defects(c2i, i2e, def_coord, p_def)

    if show_norm_plot:
        plot_lattice_errors(i2e, distance, p_mu, p_sigma)
    
    #Circuit builder steps:
    string += initialization_step(distance, i2e)
    string += rounds_step(distance, rounds, i2e)
    string += final_step(distance, i2e)

    if return_qubit_mapping:
        return string, i2e, c2i
    else:
        return string