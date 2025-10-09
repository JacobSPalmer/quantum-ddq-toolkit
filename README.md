# quantum-ddq-toolkit #
This toolkit aims to provide methodologies for the creation of quantum circuits &amp; codes with granular control of specific qubit parameters to allow for simulating, diagnosing, and testing the effects of defective/dead qubits on quantum error correction codes.

## Installing Packages ##
All the packages required (and the versions of each for the current release of this repo) are included as a Conda environment for ease of use.

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