# -*- coding: utf-8 -*-
"""
Created on Fri Jan 06, 2023 by George Lu

Utility functions that are called by multiple different classes

"""
import numpy as np

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

    # cut out the start and end to remove popping
    chirp = self.isel(chirp_time =slice(5000,-500))
    chirp_values = chirp.values.squeeze()
    
    samplerate = chirp.chirp_time.size / ((chirp.chirp_time[-1].values - chirp.chirp_time[0].values))
    samplerate = samplerate.astype(int)

    if play:
        sd.play(chirp_values, samplerate=samplerate)

    if save:
        sf.write(f"{wav_filename} .wav", chirp_values, samplerate=samplerate)