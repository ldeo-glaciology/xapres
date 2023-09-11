import pytest
from xapres_package import ApRESDefs 
import numpy

# Test the loading of a single dat file from the google bucket

def test_dat_file_loading():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    xa1 = ApRESDefs.xapres(max_range=1400)
    xa1.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0])

    
    assert numpy.isclose(xa1.data.profile.mean().values.real,8.093929591292018e-07)   