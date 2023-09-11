import pytest
from xapres_package import ApRESDefs 
import numpy
import os
def test_file_search_methods():
    xa = ApRESDefs.xapres()
    print(os.getcwd())    
    data_directory = '../data'
    
    higher_level_list_of_dats = xa.list_files(data_directory + "/sample")

    print(higher_level_list_of_dats)
    print(xa.list_files(data_directory))
    
    # this checks that the list of files is not empty
    assert higher_level_list_of_dats  

    lower_level_list_of_dats = xa.list_files(data_directory + "/sample/polarmetric")
    # test that all the files found in a lower level directory were also found when searching in a higher level directory
    assert all(item in higher_level_list_of_dats for item in lower_level_list_of_dats)

    # test that the case of the extension (DAT vs dat) doesnt matter
    assert len(xa.list_files(data_directory + "/sample/different_case_examples")) == 2

    # test the search_suffix option is working
    assert len(xa.list_files(data_directory + "/sample/polarmetric", search_suffix='HH')) == 1