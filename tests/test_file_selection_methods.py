import pytest
from xapres_package import load

## Use two different methods for selecting the same ApRES DAT file from a google bucket. 
#  In each case load it and then check that we have loaded the correct file. 

def test_file_selection_methods():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    xa1 = load.load_from_dat(max_range=1400)
    xa1.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0,1])

    xa2 = load.load_from_dat(max_range=1400)
    xa2.load_all(directory, 
                remote_load = True,
                file_names_to_process = xa1.dat_filenames_to_process)

    assert xa1.data.identical(xa2.data)