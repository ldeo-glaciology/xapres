"""
xapres: A Python package for loading and processing ApRES radar data using xarray.

xapres provides tools for loading and processing data from the Autonomous phase-sensitive 
Radio-Echo Sounder (ApRES), an ice-penetrating radar system. The package leverages xarray 
for efficient handling of multidimensional radar data and provides methods for:

- Loading ApRES .dat files from local storage or cloud buckets
- Converting raw radar data to xarray datasets with proper coordinates and metadata
- Computing radar profiles from chirp data using FFT processing
- Analyzing displacement timeseries for ice motion studies
- Computing strain rates and velocity profiles
- Handling both attended and unattended measurement modes
- Processing polarimetric radar data

The package integrates with the BAS ApRES library for low-level file reading and extends
it with xarray-based data structures for scientific analysis workflows.

Main modules:
    load : Functions and classes for loading ApRES data files
    utils : Utility functions for data processing and analysis

Key classes:
    load.from_dats : Main class for loading and processing .dat files

Key functions:
    load.generate_xarray : Simple wrapper for loading data into xarray
    load.load_zarr : Load processed data from Zarr storage format

See the online documentation for detailed examples and API reference.
"""
import os
import sys

# add git submodule to path to allow imports to work
submodule_name = 'bas-apres'
#(parent_folder_path, current_dir) = os.path.split(os.path.dirname(__file__))
#sys.path.append(os.path.join(parent_folder_path, submodule_name))
#print(f"adding {(os.path.join(parent_folder_path, submodule_name))} to path in init")

#print(os.path.dirname(__file__))
#print( f"adding {os.path.join(os.path.dirname(__file__), submodule_name)} to path in init")
sys.path.append(os.path.join(os.path.dirname(__file__), submodule_name))


from . import load
from . import utils

utils.add_methods_to_xarrays()



