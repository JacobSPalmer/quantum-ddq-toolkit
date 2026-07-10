# quantum-ddq-toolkit #
[![DOI](https://zenodo.org/badge/952059039.svg)](https://doi.org/10.5281/zenodo.21299538)

This toolkit aims to provide methodologies for creating quantum circuits & codes with granular control over specific qubit parameters, enabling simulation, diagnosis, and testing of the effects of defective/dead and heterogeneous qubits on quantum error correction codes. This tool is a living project and will be used and updated as this research develops. Please feel free to contact me if you have any questions on using the tool! :)

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
