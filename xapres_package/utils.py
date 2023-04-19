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
