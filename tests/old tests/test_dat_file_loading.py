import pytest
from xapres import load 
import numpy

# Test the loading of a single dat file from the google bucket

def test_dat_file_loading():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    fs = load.from_dats(max_range=1400)
    fs.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0],
                bursts_to_process=[0])
    assert numpy.isclose(fs.data.chirp.mean().values, 0.02611298) 