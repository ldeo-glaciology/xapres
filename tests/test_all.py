import pytest
import xarray as xr
from xapres import load, utils
import numpy as np
from numpy import allclose as npc
import os

def test_bound_methods_are_added_correctly():
    assert xr.DataArray.dB
    assert xr.DataArray.sonify
    assert xr.DataArray.displacement_timeseries

# Test the loading of a single dat file from the google bucket
def test_dat_file_loading():
    print('new version')
    directory='gs://ldeo-glaciology/apres/thwaites/continuous/ApRES_Lake1/SD2/DIR2023-01-16-0336'
    fs = load.from_dats()
    fs.list_files(directory)
    fs.load_all(directory)
    fs.load_all(directory,
            file_numbers_to_process=[0])
    load.generate_xarray(directory=directory)
    
    
# test the displacement calculation
def test_displacement_calculation():
    fd = load.from_dats()

    from_local = fd.load_all(directory='data/sample/multi-burst-dat-file/') # load the data from a local directory
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
    fs1.load_all(directory, file_numbers_to_process=[0,1])

    fs2 = load.from_dats()
    fs2.load_all(directory, file_names_to_process = fs1.dat_filenames_to_process)

    assert fs1.data.equals(fs2.data)

    fs3 = load.from_dats()
    fs3.load_all(directory, file_numbers_to_process="All")
    fs4 = load.from_dats()
    fs4.load_all(directory, file_names_to_process="All")

    with pytest.raises(ValueError):
        fs4 = load.from_dats()
        fs4.load_all(directory, file_numbers_to_process="All", file_names_to_process = fs1.dat_filenames_to_process) 

## test polarmetric local loading by loading the same waypoint twice as if it is two different ones and chaeking 
# that you get the same thing twice.

# tests the attended option and the polarmetric option from locally stored dat files (as opposed to cloud stored dat files)

def test_polarmetric_load_directory_list():
    
    fs = load.from_dats()
    fs.load_all(attended=True, 
                directory=["data/sample/polarmetric", "data/sample/polarmetric"], 
                polarmetric=True)
    
    assert len(fs.data.waypoint) == 2
    assert all(fs.data.isel(waypoint=0).filename.values == fs.data.isel(waypoint=1).filename.values)

def test_load_single():
    fs = load.from_dats()
    fs.load(dat_filename = 'data/sample/multi-burst-dat-file/DATA2022-05-22-1939.DAT')
    fs.load(dat_filename = 'data/sample/multi-burst-dat-file/DATA2022-05-22-1939.DAT',
        bursts_to_process=[0])
    fs.load(dat_filename = 'data/sample/multi-burst-dat-file/DATA2022-05-22-1939.DAT',
        bursts_to_process=[0, 2])
    fs.load(dat_filename = 'data/sample/single_dat_file/DATA2023-01-05-0315.DAT',
            attended = True) 
    fs.load(dat_filename = 'data/sample/single_dat_file/DATA2023-01-05-0315.DAT',
        bursts_to_process=[0, 2], attended = True)

def test_attended_load_all_directory_str():
    fd = load.from_dats()
    fd.load_all(attended=True, directory='data/sample/attended/')
    
def test_attended_load_all_directory_list():
    fd = load.from_dats()
    fd.load_all(attended=True, directory=['data/sample/attended/'])
  
def test_unattended_load_all_options():
    fd = load.from_dats()
    from_local = fd.load_all(directory='data/sample/multi-burst-dat-file/', disable_progress_bar=True) 
    from_local = fd.load_all(directory='data/sample/multi-burst-dat-file/', disable_progress_bar=False) 
    from_local = fd.load_all(directory='data/sample/multi-burst-dat-file/', max_range = 1500) 

def test_load_single_dat_file_with_multiple_bursts():
    dat_file = 'data/sample/multi-burst-dat-file/DATA2022-05-22-1939.DAT'
    fd = load.from_dats()
    ds = fd.load(dat_file, computeProfiles=False)
    ds = fd.load(dat_file, computeProfiles=True)

def test_load_single_dat_file_with_single_bursts():
    dat_file = 'data/sample/single_dat_file/DATA2023-01-05-0315.DAT'
    fd = load.from_dats()
    ds = fd.load(dat_file, computeProfiles=False)
    ds = fd.load(dat_file, computeProfiles=True)
    ds = fd.load(dat_file, computeProfiles=False, attended = True)
    ds = fd.load(dat_file, computeProfiles=True, attended = True)
    ds = fd.load(dat_file, computeProfiles=True, attended = True, max_range = 100)
    ds = fd.load(dat_file, computeProfiles=True, addProfileToDs_kwargs = {'crop_chirp_start': 0,'crop_chirp_end': 0.3, 'demean': False, 'scale_for_window': False, 'drop_noisy_chirps': True, 'clip_threshold': 0.8})

def test_chirp_time_dtype():
    fd = load.from_dats()
    directory='data/sample/single_dat_file/'
    ds = fd.load_all(directory)
    assert ds.chirp_time.dtype == 'float64'
    dat_file = 'data/sample/single_dat_file/DATA2023-01-05-0315.DAT'
    ds = fd.load(dat_file)
    assert ds.chirp_time.dtype == 'float64'

def test_not_supplying_a_directory():
    cwd = os.getcwd()
    os.chdir('data/sample/single_dat_file/')
    fd = load.from_dats(loglevel='DEBUG')
    fd.list_files()
    fd.load_all()
    os.chdir(cwd)

def test_when_no_files_are_found():
    fd = load.from_dats()
    fd.list_files('data/sample/empty')
    fd.load_all('data/sample/empty')

def test_selecting_bursts_to_process():
    directory='data/sample/multi-burst-dat-file'
    fs = load.from_dats()
    fs.load_all(directory, bursts_to_process='All')
    fs.load_all(directory, bursts_to_process=1)
    fs.load_all(directory, bursts_to_process=[1,2])
    fs.load_all(directory, bursts_to_process=[1,3, 20000])

#  Test `generate_xarray` and `load_zarr` wrappers
def test_wrappers_on_local():
    
    from_DAT_unattended = load.generate_xarray(directory='data/sample/single_dat_file/', 
                file_numbers_to_process = [0], 
                bursts_to_process=[0],
                )

    from_zarr = load.load_zarr()

def test_fft_calculations():
    # initialize
    fd = load.from_dats()
    directory='data/sample/single_dat_file/'

    ## define some options for the new fft to use the older defaults, i.e. before we corrected things to agree exactly with fmcw_load (https://github.com/ldeo-glaciology/xapres/pull/62)
    ops = {'demean': False, 'scale_for_window': False}
    
    ## load from a local directory, compute the fft with the new method
    load_newfft_full = fd.load_all(directory, addProfileToDs_kwargs=ops)
    ## load from a local directory, compute the fft with the new method while  setting the crop-limits on the chirp to be their default values (this shouldnt effect the answer from the line above)
    load_newfft_defaultLimits = fd.load_all(directory, addProfileToDs_kwargs={'crop_chirp_start': 0,'crop_chirp_end': 1} | ops).profile
    ## load from a local directory, compute the fft with the new method while  setting the crop-limits on the chirp to some other values (this will effect the answer)
    load_newfft_nonDefaultLimits = fd.load_all(directory, addProfileToDs_kwargs={'crop_chirp_start': 0,'crop_chirp_end': 0.5} | ops).profile

    # Compute the ffts on pre-loaded data
    ## use the method .addProfileToDs() to compute the fft on a pre-loaded dataset
    afterLoad_newfft_ds  = load_newfft_full.addProfileToDs(**ops)
    ## use the method .computeProfile() to compute the fft on a pre-loaded chirp dataarray
    afterLoad_newfft_da  = load_newfft_full.chirp.computeProfile(**ops)

    # Change a constant used in the calculation of the range, this dosesnt effect the profiles, just the profile_range
    constants = {'c': 1e8}
    afterLoad_newfft_da_differentConstants = load_newfft_full.chirp.computeProfile(constants=constants, **ops)

     ## choose some other options
    load_newfft_full.chirp.computeProfile(drop_noisy_chirps=True)
    load_newfft_full.chirp.computeProfile(detrend=True)
    load_newfft_full.chirp.computeProfile(stack=True)
    load_newfft_full.chunk({'time': 1, 'profile_range': 1}).chirp.computeProfile()

    assert npc(afterLoad_newfft_ds.profile.values, load_newfft_full.profile.values)
    assert npc(afterLoad_newfft_da.values, load_newfft_full.profile.values)
    assert npc(afterLoad_newfft_da_differentConstants.values, load_newfft_full.profile.values)
    assert not npc(afterLoad_newfft_da_differentConstants.profile_range.values, load_newfft_full.profile_range.values)

def test_comparison_with_matlab_code():
    """ Test that we get the same results as apres data loaded using fmcw_load, 
    a script typically used to load apres data:

    ApRES data is traditionally processed with code written in matlab supplied by bas (fmcw_plot, etc., https://github.com/ldeo-glaciology/phase-sensitive-radar-processing/tree/main/code/ApRES_Utils_Nicholls_250221). 

    Here we compare the results of running xpares and the results of running scripts from this collection of matlab code. 

    """
    # load the chirps, perform an fft, and put them all in an xarray
    directory = 'data/sample/thwaites/'
    p_data = load.generate_xarray(directory=directory, addProfileToDs_kwargs={'max_range': 2500})
    p_chirps = p_data['chirp']

    """
    Load output of matlab code
    The mat file loaded below was created using this version of the fmcw code: https://github.com/ldeo-glaciology/phase-sensitive-radar-processing/tree/main/code/ApRES_Utils_Nicholls_250221
    and these commands
    ```
    vdat = fmcw_load('../../data/sample/thwaites/DATA2023-02-12-0437.DAT'); %defaults to just the first burst
    [rc,~,spec_cor,spec] = fmcw_range(vdat, 2, 2500, @blackman);
    save('../data/sample/thwaites/DATA2023-02-12-0437_p4.mat')
    ```
    """
    # import mat file from the data directory
    import scipy.io as sio
    m_data = sio.loadmat('data/sample/thwaites/DATA2023-02-12-0437_p2.mat')

    m_profiles = m_data['spec_cor']

    vdat = m_data['vdat']

    m_chirps = vdat[0,0]['vif']
    m_chirp = m_chirps[0,]

    # extract just one chirp fro the xapres-loaded data
    p_chirp = p_chirps.isel(time=0, chirp_num=0).values.squeeze()

    # ensure that they are the same except for a 1.25 offset (which is not important)
    assert npc(m_chirp, p_chirp[:len(m_chirp)] + 1.25)

    # compute profiles using the same constants that are generated by the matlab code to get an exact comparison
    ## define a constants dict
    constants = {}
    constants['T'] = 1               # chirp duration [s]
    constants['f_1'] = 200e6         # starting frequency [Hz]
    constants['f_2'] = 400e6         # ending frequency [Hz]
    constants['B'] = vdat[0][0]['B'][0][0]         # bandwidth [Hz]
    constants['K'] = vdat[0][0]['K'][0][0]/2/np.pi            # rate of change of frequency [Hz/s]
    constants['c'] = 300000000.0     # speed of light in a vacuum [m/s]
    constants['ep'] = 3.18           # permittivity of ice
    constants['f_c'] = vdat[0][0]['fc'][0][0]# center frequency [Hz]
    constants['dt'] = 1/40000        # time step [s]

    # compute the profiles using these constants
    temp_new_constants = p_data.chirp.isel().computeProfile(constants=constants, max_range = 2500).isel(time=0).values.squeeze() 
   
    # compare all the profiles to the atlab-loaded ones.
    assert np.allclose(m_profiles, temp_new_constants)

def test_invalid_log_level():
    with pytest.raises(ValueError):
        load.from_dats(loglevel='INVALID')

def test_dB_method():
    directory='data/sample/single_dat_file/'
    fs = load.from_dats()
    ds = fs.load_all(directory)
    ds.profile.dB()

def test_sonify_method():
    directory='data/sample/single_dat_file/'
    fs = load.from_dats()
    ds = fs.load_all(directory)
    ds.chirp.isel(chirp_num=0).sonify(save=True)
    with pytest.raises(BaseException):
        ds.chirp.sonify() # rasies an error because we are asking for more than one chirp to be sonified
    
    # test the case when chirp_time is in timedelta format
    from_zarr = load.load_zarr() 
    from_zarr.isel(chirp_num=0, time = 100, attenuator_setting_pair = 0 ).chirp.sonify()

def test_fft_with_no_constants_supplied():
    directory='data/sample/single_dat_file/'
    fs = load.from_dats()
    ds = fs.load_all(directory)
    p1 = ds._replace(attrs={}).addProfileToDs().profile
    p2 = ds.addProfileToDs().profile
    assert p1.equals(p2)

def test_attended_fft():
    fd = load.from_dats()
    fd.load_all(attended=True, directory='data/sample/attended/').addProfileToDs()
    fd.load_all(attended=True, directory='data/sample/attended/').chirp.computeProfile()