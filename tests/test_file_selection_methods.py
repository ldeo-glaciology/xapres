import pytest
from xapres_package.utils import *
from xapres_package import ApRESDefs

## Use two different methods from selecting the same ApRES DAT file from a google bucket. 
#  In each case load it and then check that we have loaded the correct file. 

def test_file_selection():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    xa1 = ApRESDefs.xapres(max_range=1400)
    xa1.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0,1])

    xa2 = ApRESDefs.xapres(max_range=1400)
    xa2.load_all(directory, 
                remote_load = True,
                file_names_to_process = xa1.dat_filenames_to_process)

    assert xa1.data.identical(xa2.data)