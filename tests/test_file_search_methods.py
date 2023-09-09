import pytest
from xapres_package import ApRESDefs 
import numpy

def test_file_selection():
    xa = ApRESDefs.xapres()

    # test that the recursive searching works. 
    #    i.e. when you search in a higher-level directory, you recover the contents of a lower-level directory
    higher_level_list_of_dats = xa.list_files("../../data/sample")
    lower_level_list_of_dats = xa.list_files("../../data/sample/polarmetric")
    assert all(item in higher_level_list_of_dats for item in lower_level_list_of_dats)

    # test that the case of the extension doesnt matter
    assert len(xa.list_files("../../data/sample/different_case_examples")) == 2

    # test the search_suffix option is working
    assert len(xa.list_files("../../data/sample/polarmetric", search_suffix='HH')) == 1