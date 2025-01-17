"""
xapres is a python package for loading and processing data from ApRES, an ice-penetrating radar, using xarray. 
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



