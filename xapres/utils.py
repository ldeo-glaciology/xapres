"""
Functions for performing common operations on ApRES data.
"""
import numpy as np
import xarray as xr
import datetime
import dask.array as da

def displacement_timeseries(self: xr.DataArray, 
                            offset: int=1,  
                            bin_size: int=20, 
                            lower_limit_on_fit: float=800.0): 
    """
    Compute displacement, phase, coherence and associated uncertainties, as functions of depth and time, given a time series of complex ApRES profiles. 

    Profiles offset by a user-defined number of time steps are compared to compute the time series of displacement. 
    
    The time interval between outputted profiles is offset*dt where dt is the time between measurements. 

    The profiles of displacement etc. are binned in depth, with bin_size samples in each bin.

    Parameters:
    - self (xr.DataArray): The input data array containing a time series of complex profiles.
    - offset (int, optional): The time offset between the two time series. Default is 1.
    - bin_size (int, optional): Size of the vertical bins. Default is 20.

    Returns:
    xr.Dataset: Timeseries of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth.

    """
    # extract two time series
    profile1_unaligned = self.isel(time=slice(0,-offset))
    profile2_unaligned= self.isel(time=slice(offset,None))
    
    # compute the binned coherence, phase, and displacement, with uncertainties and put them into an xarray dataset  
    ds = compute_displacement(profile1_unaligned, 
                              profile2_unaligned, 
                              bin_size = bin_size,
                              lower_limit_on_fit = lower_limit_on_fit)

    # add attributes related to the this processing
    ds.attrs["offset"] = offset
    ds.attrs["processing"] = f"Created by the displacement_timeseries function in xapres using an offset of {offset} and bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def compute_displacement(profile1_unaligned: xr.DataArray, 
                        profile2_unaligned: xr.DataArray, 
                        bin_size: int=20, 
                        lower_limit_on_fit: float=800.0):
    """
    Compute displacement, coherence, and related uncertainties for binned time series data.

    Parameters:
    - profile1_unaligned (xr.DataArray): Unaligned time series data for the first measurement.
    - profile2_unaligned (xr.DataArray): Unaligned time series data for the second measurement.
    - bin_size (int, optional): Size of the vertical bins. Default is 20.

    Returns:
    xr.Dataset: Timeseries of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth.

    """
    if not isinstance(profile1_unaligned, xr.DataArray) or not isinstance(profile2_unaligned, xr.DataArray):
        raise TypeError("profile1_unaligned and profile2_unaligned must be xarray DataArrays")
    
    profiles = combine_profiles(profile1_unaligned, profile2_unaligned)

    profiles_binned = bin_profiles(profiles, bin_size)

    coherence = compute_coherence(profiles_binned.isel(shot_number=0), profiles_binned.isel(shot_number=1))
   
    # add attributes related to the bin_depth
    coherence.bin_depth.attrs["units"] = "m"
    coherence.bin_depth.attrs["long_name"] = "depth to the center of each bin"
    coherence.bin_depth.attrs["standard_name"] = "bin depth"         

    # add attributes related to the coherence
    coherence.attrs["units"] = "unitless"
    coherence.attrs["long_name"] = "complex coherence between measurements"

    # compute the phase and add attributes
    phase = -xr.apply_ufunc(np.angle, coherence, dask="allowed").rename("phase")
    phase.attrs["units"] = "radians"
    phase.attrs["long_name"] = "coherence phase"

    # compute phase uncertainties
    phase_uncertainty = ((1./abs(coherence))*np.sqrt((1.-abs(coherence)**2.)/(2.*bin_size))).rename('phase_uncertainty')
    phase_uncertainty.attrs["units"] = "radians"
    phase_uncertainty.attrs["long_name"] = "uncertainty in coherence phase"

    # compute the displacement
    displacement = phase2range(phase).rename("displacement")
    displacement.attrs["units"] = "m"
    displacement.attrs["long_name"] = "displacement since previous measurement"

    # compute the displacement uncertainty
    disp_uncertainty = phase2range(phase_uncertainty).rename('disp_uncertainty')
    disp_uncertainty.attrs["units"] = "m"
    disp_uncertainty.attrs["long_name"] = "uncertainty in displacement since previous measurement"

    dt_years = ((profiles.profile_time.sel(shot_number=2) - profiles.profile_time.sel(shot_number=1)) / np.timedelta64(1,'D') / 365.25).rename('dt_years')
    dt_years.attrs['units'] = 'years'
    dt_years.attrs['long_name'] = 'Time between shots'
    dt_years.attrs['description'] = 'Time in years between shots used in each measurement of displacement, vertical velocity, etc. dt_years[i] is the time between shot [j] and shot [j-1]'

    # vertical velocity
    velocity = (displacement / dt_years).rename('velocity')
    velocity.attrs['units'] = 'meters/year'
    velocity.attrs['long_name'] = 'Vertical velocity'

    # strain rates
    strain_rates = velocity.computeStrainRates(lower_limit_on_fit = lower_limit_on_fit)

    # combine to an xarray dataset
    da_list = [profiles, coherence, phase, phase_uncertainty, displacement, disp_uncertainty, velocity, strain_rates]
    ds = xr.merge(da_list)

    # add attributes related to this processing
    ds.attrs["bin_size"] = bin_size
    ds.attrs["description"] = "Time series of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth."
    ds.attrs["processing"] = f"Created by the compute_displacement function in xapres using a bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def combine_profiles(profile1_unaligned, profile2_unaligned):
    """Combine two timeseries of profiles. In the case of unattended data, record the midpoint time and the time of each computed profile"""
    
    if 'waypoint' in profile1_unaligned.coords:
        # data is taken in attended mode and we dont need to get the midpoint time and align
        # rename time as profile_time
        profile1_unaligned = profile1_unaligned.rename({'time':'profile_time'})
        profile2_unaligned = profile2_unaligned.rename({'time':'profile_time'})
        profiles = xr.concat([profile1_unaligned, profile2_unaligned], dim='shot_number')

    else:
        # in the case when we selected the time step with .isel(time=N), where N is an integer, we dont have time as a dimension. The following accounts for this scenario.
        if 'time' not in profile1_unaligned.dims:
            profile1_unaligned = profile1_unaligned.expand_dims(dim="time")
        if 'time' not in profile2_unaligned.dims:
            profile2_unaligned = profile2_unaligned.expand_dims(dim="time")

        #if 'time' not in profile1_unaligned.dims and 'time' in profile1_unaligned.coords:
            
        #else:
        
        # record the time interval between measurements
        t1 = profile1_unaligned.time.data
        t2 = profile2_unaligned.time.data
        dt = t2-t1

        # change the name of the time coordinates so that they can be retained when the profiles are concatenated, then drop the original 'time' coordinates (the latter is achieved with drop_vars)
        profile1_unaligned = profile1_unaligned.assign_coords(profile_time=("time", profile1_unaligned.time.data)).drop_vars("time")
        profile2_unaligned = profile2_unaligned.assign_coords(profile_time=("time", profile2_unaligned.time.data)).drop_vars("time")

        # concatenate the two profiles
        profiles = xr.concat([profile1_unaligned, profile2_unaligned], dim='shot_number', coords=['profile_time', 'burst_number', 'filename'])

        # add the midpoint time 
        profiles = profiles.assign_coords(time=("time", t1+dt/2))
        profiles.time.attrs["description"] = "mid-point time of two profiles used in the computation"

    profiles = profiles.assign_coords(shot_number=("shot_number", 1+profiles.shot_number.data))
    profiles.shot_number.attrs['long_name'] = 'shot number'
    profiles.shot_number.attrs['description'] = 'number of the shot used in each measurement'

    return profiles

def bin_profiles(profiles, bin_size):
    """Put the time series into vertical bins"""
    
    # Bin in depth
    profiles_binned = profiles.coarsen(profile_range=bin_size, boundary='trim').construct(profile_range=("bin_depth", "sample_in_bin"))

    # Compute the bin depth and add it to the DataArray
    bin_depth = profiles_binned.profile_range.mean(dim='sample_in_bin').data
    profiles_binned = profiles_binned.assign_coords(bin_depth=("bin_depth", bin_depth))

    return profiles_binned

def compute_coherence(b1_binned, b2_binned):
    # compute the coherence
    top = (b1_binned * np.conj(b2_binned)).sum(dim="sample_in_bin")
    bottom = np.sqrt( (np.abs(b1_binned)**2).sum(dim="sample_in_bin") * (np.abs(b2_binned)**2).sum(dim="sample_in_bin"))

    return (top/bottom).rename("coherence")

def phase2range(phi, 
                lambdac=0.5608, 
                rc=None, 
                K=2e8, 
                ci=1.6823e8):
        """
        Convert phase difference to range for FMCW radar
        Parameters
        ---------
        lambdac: float
            wavelength (m) at center frequency
        rc: float; optional
            coarse range of bin center (m)
        K:  float; optional
            chirp gradient (rad/s/s)
        ci: float; optional
            propagation velocity (m/s)
        ### Original Matlab File Notes ###
        Craig Stewart
        2014/6/10
        

        if not all([K,ci]) or rc is None:
            # First order method
            # Brennan et al. (2014) eq 15
            print('not precise')
            r = lambdac*phi/(4.*np.pi)
        else:
            print('Precise')
            # Appears to be from Stewart (2018) eqn 4.8, with tau = 2*R/ci and omega_c = 2 pi /lambdac, where R is the range
            r = phi/((4.*np.pi/lambdac) - (4.*rc[None,:]*K/ci**2.))
        """
        return lambdac*phi/(4.*np.pi)

def computeStrainRates(self, lower_limit_on_fit = 800):
    """Compute strain rates from a dataset of ApRES data. For use by the function `compute_displacement`"""
    velocity_cropped = self\
            .squeeze()\
            .where(self.bin_depth < lower_limit_on_fit)
    
    fit_ds = velocity_cropped.polyfit('bin_depth', 1, full = True)
            
    strain_rate = fit_ds.sel(degree = 1, drop =True).polyfit_coefficients.rename('strain_rate')

    surface_intercept =  fit_ds.sel(degree = 0, drop =True).polyfit_coefficients.rename('surface_intercept') 

    # R^2
    y_mean = velocity_cropped.mean(dim = 'bin_depth')
    SS_tot = ((velocity_cropped - y_mean)**2).sum(dim = 'bin_depth')
    R2 = (1 - (fit_ds.polyfit_residuals/SS_tot)).rename('r_squared')

    # add attrs
    strain_rate.attrs['units'] = '1/year'
    strain_rate.attrs['long_name'] = f"vertical strain rate in upper {lower_limit_on_fit} m"
    strain_rate.attrs['lower_limit_on_fit_meters'] = lower_limit_on_fit

    surface_intercept.attrs['units'] = 'meters/year'
    surface_intercept.attrs['long_name'] = 'vertical velocity at the surface from the linear fit'
    surface_intercept.attrs['lower_limit_on_fit_meters'] = lower_limit_on_fit
    
    R2.attrs['long_name'] = 'r-squared value for the linear fit'
    R2.attrs['units'] = '-' 
    
    return xr.merge([strain_rate, surface_intercept, R2])

def dB(self):
    """
    A function to convert profile data to decibels.
    
    The function is added to xarray dataarrays as a bound method in two functions. 
    
    """
    return 20*np.log10(np.abs(self)) 

def sonify(self, 
           play=True, 
           save=False, 
           wav_filename="chirp"):
    """
    A function to sonify a chirp - play the signal as a sound.
     
    The function is added to xarray dataarrays as a bound method in two functions. 
    
    It requires soundfile and sounddevice to be installed.
    """

    try:     
        import soundfile as sf
        import sounddevice as sd 
    except ImportError:
        print("sounddevice and soundfile are required to sonify the chirps. pip install them if you need this feature") 

    # make sure the input is just one chirp    
    if self.size != self.chirp_time.size: 
        raise BaseException('sonify only works for single chirps.')    

    # subset the chirp to remove popping
    chirp = self.isel(chirp_time =slice(5000,-500))

    # convert to a numpy array and remove singleton dimensions
    chirp_values = chirp.values.squeeze()

    t = chirp.chirp_time.values

    # convert the start and end time of the subsetted chirp to seconds, dealing with the cases when chirp.chirp_time is numpy array and when it is a timeselta64
    if isinstance(t[0], float):
        startTimeInSeconds = t[0] 
        endTimeInSeconds = t[-1] 
    elif isinstance(t[0], np.timedelta64):
        startTimeInSeconds = t[0] / np.timedelta64(1, 's')
        endTimeInSeconds = t[-1] / np.timedelta64(1, 's')
        print(startTimeInSeconds, endTimeInSeconds)

    # calculate the sample rate 
    samplerate = chirp.chirp_time.size / (endTimeInSeconds - startTimeInSeconds)
    samplerate = samplerate.astype(int)

    #
    if play:
        sd.play(chirp_values, samplerate=samplerate)

    if save:
        sf.write(f"{wav_filename} .wav", chirp_values, samplerate=samplerate)

def addProfileToDs(self: xr.Dataset, **kwargs):

    
    if 'constants' in self.attrs:
        profile = self.chirp.computeProfile(constants = self.attrs['constants'], **kwargs)
    else:
        profile = self.chirp.computeProfile(**kwargs)

    # remove profile variable and profile range, if they exist 
    if 'profile' in self.data_vars:
        out = self.drop_dims('profile_range')
    else: 
        out = self
        
    return xr.merge([out, profile], combine_attrs='override')


def computeProfile(self: xr.DataArray,
                   pad_factor=2, 
                   drop_noisy_chirps=False,
                   clip_threshold=1.2,
                   min_chirps=0,
                   demean=True,
                   detrend=False,
                   stack=False,
                   scale_for_window=True,
                   crop_chirp_start=0,
                   crop_chirp_end=1,
                   max_range=None,
                   constants={}):
    """
    Compute profiles from chirp data.
    -----------
    self : xr.DataArray
        The input chirp data array.
    pad_factor : int, optional
        Factor by which to pad the chirp data (default is 2).
    drop_noisy_chirps : bool, optional
        Whether to drop noisy chirps (default is False).
    clip_threshold : float, optional
        Threshold for clipping noisy chirps (default is 1.2).
    min_chirps : int, optional
        Minimum number of chirps required to keep a burst (default is 20).
    demean : bool, optional
        Whether to demean the chirp data (default is False).
    detrend : bool, optional
        Whether to detrend the chirp data (default is False).
    stack : bool, optional
        Whether to stack the chirp data (default is False).
    scale_for_window : bool, optional
        Whether to scale the chirp data for the window (default is True).
        This is only an option to allow tests to effectively compare the output of 
        this function with the output of the legacy fft method.
    crop_chirp_start : float, optional
        Start time for cropping chirps in seconds (default is 0, the start of the chirp).
    crop_chirp_end : float, optional
        End time for cropping chirps in seconds (default is 1, the end f the chirp).
    max_range : float, optional
        Maximum range for the profile in meters (default is None).
    constants : dict, optional
        Dictionary of user-defined constants for the radar system. 
        Any constants that are not supplied in constants are defined in default_constants()
    
    Returns:
    --------
    xr.DataArray
        The computed radar profile with range as the coordinate.
    """

    constants = default_constants() | constants

    B = constants['B']       # bandwidth [Hz]
    K = constants['K']       # rate of chnge of frequency [Hz/s]
    c = constants['c']       # speed of light in a vacuum [m/s]
    ep = constants['ep']     # permittivity of ice
    f_c = constants['f_c']   # center frequency [Hz]
    dt = constants['dt']     # time step [s]

    def rdei(x):
        """round down to the nearest even integer and return an integer"""
        return int(np.floor(x/2) * 2)
    
    def freq2range(frequencies):
        """"return the range for a given frequency"""
        return c * frequencies / (2*np.sqrt(ep)*K)

    Nt = rdei(self.chirp_time.size)   
    chirps = self.isel(chirp_time = slice(0, Nt))
    
    sampling_frequency = 1/dt 

    # if crop_chirp_start is not None:
    chirps = chirps.sel(chirp_time = slice(crop_chirp_start, crop_chirp_end))

    if drop_noisy_chirps:
        bad_chirps =  chirps.where(abs(chirps) > clip_threshold)
        good_bursts = bad_chirps.max(dim='chirp_time').count(dim='chirp_num') > min_chirps
        chirps = chirps.where(good_bursts)
        chirps = chirps.where(abs(chirps).max(dim='chirp_time')<clip_threshold)

    if demean:
        chirps = chirps - chirps.mean(dim='chirp_time')

    if detrend:
        p = chirps.polyfit('chirp_time', 1)
        fit = xr.polyval(chirps.chirp_time, p.polyfit_coefficients)
        chirps = chirps - fit

    if stack:   
        chirps = chirps.mean(dim='chirp_num', skipna=True)

    # note on variable naming below: s stands for the pre-fft signal, following others' notation, w stands for windowed, p stands for padding, and so on
   
    # window
    window = xr.DataArray(np.blackman(chirps.chirp_time.size), dims = 'chirp_time')
    s_w = chirps * window  
    
    # pad
    s_wp = s_w.pad(chirp_time=int((Nt*pad_factor-Nt)/2), constant_values=0)  

    # roll
    s_wpr = s_wp.roll(chirp_time=int(Nt*pad_factor/2))  

    if contains_dask_array(s_wpr):
        s_wpr = s_wpr.chunk({'chirp_time':-1}).load()

    # fft
    S_wpr = xr.apply_ufunc(np.fft.fft, 
                        s_wpr,
                        input_core_dims=[["chirp_time"]],
                        output_core_dims=[["chirp_time"]])
    
    # scale for padding
    S_wpr = S_wpr/s_wpr.chirp_time.size * np.sqrt(2*pad_factor) 
    
    if scale_for_window:
        S_wpr = S_wpr/ np.sqrt(np.mean((np.blackman(Nt))**2)) # scale for window
    
    # compute range
    indexes      = np.arange(s_wpr.chirp_time.size) 
    frequencies  = indexes * sampling_frequency/s_wpr.chirp_time.size
    profile_range = freq2range(frequencies)

    # reference array
    m = np.arange(len(S_wpr.chirp_time))/pad_factor
    phiref = 2*np.pi*f_c*m/B -  m * m * 2*np.pi * K/2/B**2
    S_wprr = S_wpr * np.exp(phiref*(-1j))
    S_wprr = S_wprr.rename('profile')
    S_wprr = S_wprr.rename({'chirp_time': 'profile_range'})
    S_wprr['profile_range'] = profile_range

    # crop to max_range
    if max_range is None:
        max_range = S_wprr.profile_range[-1]/2
        S_wprr = S_wprr.isel(profile_range = slice(0, S_wprr.profile_range.size//2-1))
    else:
        S_wprr = S_wprr.where(S_wprr.profile_range <= max_range, drop=True)
    
    # add attributes
    S_wprr.attrs['long_name'] = 'profile' 
    S_wprr.attrs['units'] = '-'
    S_wprr.attrs['description'] = 'complex profile computed from the fourier transform of the de-ramped chirp'
    S_wprr.profile_range.attrs['long_name'] = 'depth'
    S_wprr.profile_range.attrs['units'] = 'meters'
    S_wprr.attrs['constants'] = constants

    return S_wprr

def default_constants():
    constants = {}
    constants['T'] = 1               # chirp duration [s]
    constants['f_1'] = 200e6         # starting frequency [Hz]
    constants['f_2'] = 400e6         # ending frequency [Hz]
    constants['B'] = constants['f_2']-constants['f_1']          # bandwidth [Hz]
    constants['K'] = constants['B']/constants['T']            # rate of chnge of frequency [Hz/s]
    constants['c'] = 300000000.0     # speed of light in a vacuum [m/s]
    constants['ep'] = 3.18           # permittivity of ice
    constants['f_c'] = (constants['f_2']+constants['f_1'])/2   # center frequency [Hz]
    constants['dt'] = 1/40000        # time step [s]

    return constants
    
def add_methods_to_xarrays():
    
    da_methods = [dB, sonify, displacement_timeseries, computeProfile, computeStrainRates]
    for method in da_methods:
        setattr(xr.DataArray, method.__name__, method)

    ds_methods = [addProfileToDs]
    for method in ds_methods:
        setattr(xr.Dataset, method.__name__, method)

def contains_dask_array(dataarray):
    return isinstance(dataarray.data, da.Array)