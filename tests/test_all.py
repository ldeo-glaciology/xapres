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
    fd = load.from_dats()

    from_local = fd.load_all(directory='data/sample/multi-burst-dat-file/', legacy_fft=False) # load the data from a local directory
    p1 = from_local.isel(time=2, attenuator_setting_pair=0).profile # select a profile 
    p2 = from_local.isel(time=5, attenuator_setting_pair=0).profile # select a different profile 

    utils.compute_displacement(p1, p2, bin_size = 25)    # calculate the displacement between the two profiles. Use a non-default value for bin_size


    profiles = from_local.isel(attenuator_setting_pair=0).profile  # select all the profiles on a specfic date
    results = profiles.displacement_timeseries(bin_size = 30, offset = 3) # compute a time series of displacement from these data. Use non-default values for offset and bin_size. 
    assert results



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
    directory='data/sample/polarmetric'
    fs1 = load.from_dats()
    fs1.load_all(directory, legacy_fft=False, file_numbers_to_process=[0,1])

    fs2 = load.from_dats()
    fs2.load_all(directory, legacy_fft=False, file_names_to_process = fs1.dat_filenames_to_process)

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

from numpy import allclose as npc

def test_fft_calculations():
    # initialize
    fd = load.from_dats()
    directory='data/sample/single_dat_file/'
    
    # load the data from dat files
    ## load data from a local directory, compute the fft with the legacy method and dont correct the padding error
    load_oldfft_uncorrectedPad = fd.load_all(directory, legacy_fft=True, corrected_pad=False).profile
    ## load data from a local directory, compute the fft with the legacy method, but this time correct the padding error
    load_oldfft_correctedPad = fd.load_all(directory, legacy_fft=True, corrected_pad=True).profile

    ## load from a local directory, compute the fft with the new method
    load_newfft_full = fd.load_all(directory, legacy_fft=False)
    load_newfft = load_newfft_full.profile
    ## load from a local directory, compute the fft with the new method while  setting the crop-limits on the chirp to be their default values (this shouldnt effect the answer from the line above)
    load_newfft_defaultLimits = fd.load_all(directory, legacy_fft=False, addProfileToDs_kwargs={'crop_chirp_start': 0,'crop_chirp_end': 1}).profile
    ## load from a local directory, compute the fft with the new method while  setting the crop-limits on the chirp to some other values (this will effect the answer)
    load_newfft_nonDefaultLimits = fd.load_all(directory, legacy_fft=False, addProfileToDs_kwargs={'crop_chirp_start': 0,'crop_chirp_end': 0.5}).profile

    # Compute the ffts on pre-loaded data
    ## use the method .addProfileToDs() to compute the fft on a pre-loaded dataset
    afterLoad_newfft_ds  = load_newfft_full.addProfileToDs()
    ## use the method .computeProfile() to compute the fft on a pre-loaded chirp dataarray
    afterLoad_newfft_da  = load_newfft_full.chirp.computeProfile()

    # Change a constant used in the calculation of the range, this dosesnt effect the profiles, just the profile_range
    constants = load_newfft_full.attrs['constants']
    constants['c'] = 2e8
    afterLoad_newfft_da_differentConstants = load_newfft_full.chirp.computeProfile(constants=constants)

    assert not npc(load_oldfft_uncorrectedPad.values, load_oldfft_correctedPad.values)
    d = load_newfft.dims  #needed to transpose the dataarrays that use legacy_fft=True to be the same as those which use legacy_fft=False
    assert npc(load_oldfft_correctedPad.transpose(*d).values, load_newfft.values)
    assert npc(load_oldfft_correctedPad.transpose(*d).values, load_newfft_defaultLimits.values)
    assert npc(afterLoad_newfft_ds.profile.values, load_newfft_full.profile.values)
    assert npc(afterLoad_newfft_da.values, load_newfft_full.profile.values)
    assert npc(afterLoad_newfft_da_differentConstants.values, load_newfft_full.profile.values)
    assert not npc(afterLoad_newfft_da_differentConstants.profile_range.values, load_newfft_full.profile_range.values)