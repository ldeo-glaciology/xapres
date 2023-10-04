import pytest
from xapres_package import load 
import numpy

# Test the loading of a single dat file from the google bucket

def test_dat_file_loading():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    xa = load.load_from_dat(max_range=1400)
    xa.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0])

    
    assert numpy.isclose(xa.data.profile.mean().values.real, 8.093929591292018e-07)   