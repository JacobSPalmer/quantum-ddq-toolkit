# quantum-ddq-toolkit #
[![DOI](https://zenodo.org/badge/952059039.svg)](https://doi.org/10.5281/zenodo.21299538)

This toolkit aims to provide methodologies for creating quantum circuits & codes with granular control over specific qubit parameters, enabling simulation, diagnosis, and testing of the effects of defective/dead and heterogeneous qubits on quantum error correction codes. This tool is a living project and will be used and updated as this research develops.

The paper associated with this work, published within npj Unconventional Computing, can be found [here](https://www.nature.com/articles/s44335-026-00084-2).

Additionally, this package serves as the spiritual predecessor to our subsequent work which expands upon the work presented here to support chip-wide, multi-qubit simulations under granular, heterogeneous noise entitled [qSNOW](https://github.com/JacobSPalmer/qSNOW), which is actively under development.

## Installing Packages ##
All required packages (and their versions for the current release of this repo) are included in a Conda environment for ease of use.

### Install & Update ###
To install the Conda env:

'''
conda env create -f environment.yml
'''

To update the conda env
'''
conda env update -f environment.yml
'''

### Run ###
To run the Conda env:
'''
conda activate DDQ
'''
