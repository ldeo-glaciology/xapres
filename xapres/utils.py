# -*- coding: utf-8 -*-
"""
Created on Fri Jan 06, 2023 by George Lu

Utility functions that are called by multiple different classes

"""
import numpy as np
import xarray as xr
import datetime



def displacement_timeseries(self: xr.DataArray, 
                            offset: int=1,  
                            bin_size: int=20): 
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
    ds = compute_displacement(profile1_unaligned, profile2_unaligned, bin_size = bin_size)

    # add attributes related to the this processing
    ds.attrs["offset"] = offset
    ds.attrs["processing"] = f"Created by the displacement_timeseries function in xapres using an offset of {offset} and bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def compute_displacement(profile1_unaligned: xr.DataArray, 
                        profile2_unaligned: xr.DataArray, 
                        bin_size: int=20):
    """
    Compute displacement, coherence, and related uncertainties for binned time series data.

    Parameters:
    - profile1_unaligned (xr.DataArray): Unaligned time series data for the first measurement.
    - profile2_unaligned (xr.DataArray): Unaligned time series data for the second measurement.
    - bin_size (int, optional): Size of the vertical bins. Default is 20.

    Returns:
    xr.Dataset: Timeseries of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth.

    """

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

    # combine to an xarray dataset
    da_list = [profiles, coherence, phase, phase_uncertainty, displacement, disp_uncertainty]
    ds = xr.merge(da_list)

    # add attributes related to this processing
    ds.attrs["bin_size"] = bin_size
    ds.attrs["description"] = "Time series of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth."
    ds.attrs["processing"] = f"Created by the compute_displacement function in xapres using a bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def combine_profiles(profile1_unaligned, profile2_unaligned):
    """Combine two timeseries of profiles. In the case of unattended data, record the midpoint time and the time of each computed profile"""
    
    if 'time' not in profile1_unaligned.dims and 'time' in profile1_unaligned.coords:
        # data is taken in attended mode and we dont need to get the midpoint time and align
        profiles = xr.concat([profile1_unaligned, profile2_unaligned], dim='shot_number')
    else:

        # in the case when we selected the time step with .isel(time=N), where N is an integer, we dont have time as a dimension. THe following accounts for this scenario
        if 'time' not in profile1_unaligned.dims:
            profile1_unaligned = profile1_unaligned.expand_dims(dim="time")
        if 'time' not in profile2_unaligned.dims:
            profile2_unaligned = profile2_unaligned.expand_dims(dim="time")

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
    profiles_binned = profiles.coarsen(profile_range=bin_size, boundary='trim').construct(profile_range=("bin", "sample_in_bin"))

    # Compute the bin depth and add it to the DataArray
    bin_depth = profiles_binned.profile_range.mean(dim='sample_in_bin').data
    profiles_binned = profiles_binned.assign_coords(bin_depth=("bin", bin_depth))

    return profiles_binned

def compute_coherence(b1_binned, b2_binned):
    # compute the coherence
    top = (b1_binned * np.conj(b2_binned)).sum(dim="sample_in_bin")
    bottom = np.sqrt( (np.abs(b1_binned)**2).sum(dim="sample_in_bin") * (np.abs(b2_binned)**2).sum(dim="sample_in_bin"))
    coherence = (top/bottom).rename("coherence")

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
        """

        if not all([K,ci]) or rc is None:
            # First order method
            # Brennan et al. (2014) eq 15
            r = lambdac*phi/(4.*np.pi)
        else:
            # Precise
            # Appears to be from Stewart (2018) eqn 4.8, with tau = 2*R/ci and omega_c = 2 pi /lambdac, where R is the range
            r = phi/((4.*np.pi/lambdac) - (4.*rc[None,:]*K/ci**2.))
        return r

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

    # calculate the sample rate 
    samplerate = chirp.chirp_time.size / (endTimeInSeconds - startTimeInSeconds)
    samplerate = samplerate.astype(int)

    #
    if play:
        sd.play(chirp_values, samplerate=samplerate)

    if save:
        sf.write(f"{wav_filename} .wav", chirp_values, samplerate=samplerate)

def add_methods_to_xarrays():
    
    methods = [dB, sonify, displacement_timeseries]
    
    for method in methods:
        setattr(xr.DataArray, method.__name__, method)