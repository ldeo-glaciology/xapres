# -*- coding: utf-8 -*-
"""
Created on Fri Jan 06, 2023 by George Lu

Utility functions that are called by multiple different classes

"""
import numpy as np
import xarray as xr

def coherence(s1, s2):
        """
        Phase correlation between two elements of the scattering matrix
        Jordan et al. (2019) eq. 13
        Parameters
        ---------
        s1: array
            first acquisition
        s2:
            second acquisition
        Output
        ---------
        c:  array
            phase coherence
        """
        top = np.einsum('ij,ij->i', s1, np.conj(s2))
        bottom = np.sqrt(np.sum(np.abs(s1)**2,axis=1)*np.sum(np.abs(s2)**2,axis=1))
        c = top/bottom

        return c

def generate_range_diff(data1, data2, win_cor, step, range_ext=None, win_wrap=10, thresh=0.9, uncertainty='CR'):
    """
    Input data:
    data1, data2: xarray.DataArrays describing the "profile" variable
    range
    """

    # Get times of burst:
    t1 = data1.time.data
    t2 = data2.time.data
    dt = (t2-t1)/ np.timedelta64(1, 's')
    #self.logger.info(f"Time between bursts : {dt}s")
    # Get phase difference
    
    if range_ext is not None:
        # Fill a depth array which will be more sparse than the full Range vector
        idxs = np.arange(win_cor//2, range_ext.shape[0]-win_cor//2, step).astype(int)
        ds = range_ext[idxs]
    else:
        idxs = np.arange(win_cor//2, data1.shape[1]-win_cor//2, step).astype(int)
        ds = data1.profile_range[idxs]

    # Create data and coherence vectors
    acq1 = data1
    acq2 = data2
    co = np.empty_like(np.stack([ds.data]*data1.shape[0])).astype(np.cdouble)
    for i, idx in enumerate(idxs):
        # index two sub_arrays to compare
        arr1 = acq1[:,idx-win_cor//2:idx+win_cor//2]
        arr2 = acq2[:,idx-win_cor//2:idx+win_cor//2]
        # correlation coefficient between acquisitions
        # amplitude is coherence between acquisitions and phase is the offset
        co[:,i] = self.coherence(arr1.data, arr2.data)
    
    # Phase unwrapping
    phi = -np.angle(co).astype(float)
    for i in range(co.shape[1]-1):
        for t in range(co.shape[0]):
            idx = i+1
            if np.all(abs(co[t,idx-win_wrap:idx+win_wrap]) < thresh):
                continue
            if phi[t,idx]-phi[t,idx-1] > np.pi:
                phi[t,idx:] -= 2.*np.pi
            elif phi[t,idx]-phi[t,idx-1] < -np.pi:
                phi[t,idx:] += 2.*np.pi
    # Range difference calculation
    w = phase2range(phi,
                        0.5608,
                        ds.data,
                        2e8,
                        1.6823e8)

    # If the individual acquisitions have had uncertainty calculations
    
    if uncertainty == 'CR':
        # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
        sigma = (1./abs(co))*np.sqrt((1.-abs(co)**2.)/(2.*win_cor))
        # convert the phase offset to a distance vector
        w_err = phase2range(sigma,
                                0.5608,
                                ds.data,
                                2e8,
                                1.6823e8)

    '''elif uncertainty == 'noise_phasor':
        # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
        # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
        r_uncertainty = phase2range(self, self.unc1, self.header.lambdac) +\
            phase2range(self, self.unc2, self.header.lambdac)
        idxs = np.arange(win//2, len(self.data)-win//2, step)
        w_err = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])'''
    
    coords = {'time':(['time'],t2,{'units': 'seconds','long_name':'Time of second burst'}),'profile_range':(['profile_range'],ds.profile_range.data,{'units': 'm','long_name':'Depth'})}
    data_vars = {'time_diff':(['time'],np.cumsum(dt),{'units': 'seconds','long_name':'Time since first burst'}),
                    'range_diff':(['time','profile_range'], w, #np.cumsum(w,axis=0), 
                        {'units': 'm', 
                        'long_name':'Range difference'}),
        'err':(['time','profile_range'], w_err, 
                        {'units': 'm', 
                        'long_name':'Error'})}
    ds_xr = xr.Dataset(data_vars=data_vars, coords=coords)
    return ds_xr, co, phi # returning velocities in mm/day


def phase2range(phi, lambdac=0.5608, rc=None, K=2e8, ci=1.6823e8):
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