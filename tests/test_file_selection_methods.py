import pytest
from xapres_package import load

## Use two different methods for selecting the same ApRES DAT file from a google bucket. 
#  In each case load it and then check that we have loaded the correct file. 

def test_file_selection_methods():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    fs1 = load.from_dats(max_range=1400)
    fs1.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0,1])

    fs2 = load.load_from_dat(max_range=1400)
    fs2.load_all(directory, 
                remote_load = True,
                file_names_to_process = fs1.dat_filenames_to_process)

    assert fs1.data.identical(fs2.data)