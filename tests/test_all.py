import pytest
import xapres
import xarray as xr
from xapres import load, utils
import numpy


def test_bound_methods_are_added_correctly():
    assert xr.DataArray.dB
    assert xr.DataArray.sonify
    assert xr.DataArray.displacement_timeseries
    

# Test the loading of a single dat file from the google bucket
def test_dat_file_loading():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    fs = load.from_dats()
    fs.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0],
                bursts_to_process=[0])
    assert numpy.isclose(fs.data.chirp.mean().values, 0.02611298) 
    
    
# test the displacement calculation
def test_displacement_calculation():
    from_zarr = load.load_zarr()  # lazily load a large APRES dataset from Greenland
    p1 = from_zarr.isel(time=2000).profile_stacked # select a profile 
    p2 = from_zarr.isel(time=2100).profile_stacked # select a different profile 
       
    utils.compute_displacement(p1, p2)    # calculate the displacement between the two profiles

    t = from_zarr.sel(time='2022-07-17').profile_stacked   # select all the profiles on a specfic date
    results = t.displacement_timeseries(bin_size = 30, offset = 3) # compute a time series of displacement from these data. Use non-default values for offset and bin_size 

    assert (abs(results)>0).all().load()



def test_file_search_methods():
    fs = load.from_dats()   
    data_directory = 'data'
    
    higher_level_list_of_dats = fs.list_files(data_directory + "/sample")
    
    # this checks that the list of files is not empty
    assert higher_level_list_of_dats  

    lower_level_list_of_dats = fs.list_files(data_directory + "/sample/polarmetric")
    # test that all the files found in a lower level directory were also found when searching in a higher level directory
    assert all(item in higher_level_list_of_dats for item in lower_level_list_of_dats)

    # test that the case of the extension (DAT vs dat) doesnt matter
    assert len(fs.list_files(data_directory + "/sample/different_case_examples")) == 2

    # test the search_suffix option is working
    assert len(fs.list_files(data_directory + "/sample/polarmetric", search_suffix='HH')) == 1



## Use two different methods for selecting the same ApRES DAT file from a google bucket. 
#  In each case load it and then check that we have loaded the correct file. 

def test_file_selection_methods():
    directory='gs://ldeo-glaciology/GL_apres_2022/A101'
    fs1 = load.from_dats()
    fs1.load_all(directory, 
                remote_load = True,
                file_numbers_to_process=[0,1])

    fs2 = load.from_dats()
    fs2.load_all(directory, 
                remote_load = True,
                file_names_to_process = fs1.dat_filenames_to_process)

    assert fs1.data.equals(fs2.data)

## test polarmetric local loading by loading the same waypoint twice as if it is two different ones and chaeking 
# that you get the same thing twice.

# tests the attended option and the polarmetric option from locally stored dat files (as opposed to cloud stored dat files)

def test_polarmetric_load():
    
    fs = load.from_dats()
    fs.load_all(attended=True, 
                directory=["data/sample/polarmetric", "data/sample/polarmetric"], 
                polarmetric=True)
    
    assert len(fs.data.waypoint) == 2
    assert all(fs.data.isel(waypoint=0).filename.values == fs.data.isel(waypoint=1).filename.values)



#  Test `generate_xarray` and `load_zarr` wrappers
def test_wrappers():
    
    from_DAT_unattended = load.generate_xarray(directory='data/sample/single_dat_file/', 
                file_numbers_to_process = [0], 
                bursts_to_process=[0],
                )

    from_zarr = load.load_zarr()
