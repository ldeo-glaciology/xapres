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
                            lower_limit_on_fit: float=None,
                            min_depth_for_ezz_fit: float=None,
                            max_depth_for_ezz_fit: float=None): 
    """
    Compute displacement, phase, coherence and associated uncertainties, as functions of depth and time, given a time series of complex ApRES profiles. 

    Profiles offset by a user-defined number of time steps are compared to compute the time series of displacement. 
    
    The time interval between outputted profiles is offset*dt where dt is the time between measurements. 

    The profiles of displacement etc. are binned in depth, with bin_size samples in each bin.

    Parameters:
    - self (xr.DataArray): The input data array containing a time series of complex profiles.
    - offset (int, optional): The time offset between the two time series. Default is 1.
    - bin_size (int, optional): Size of the vertical bins. Default is 20.
    - min_depth_for_ezz_fit (float, optional): Upper limit on the fit for computing strain rates. Default is 0.0.
    - max_depth_for_ezz_fit (float, optional): Lower limit on the fit for computing strain rates. Default is 800.0.
    - lower_limit_on_fit (float, optional): depreciated: use max_depth_for_ezz_fit instead.  
    
    Returns:

    xr.Dataset: Timeseries of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth.

    Notes:
    Computes the variance of the phase and displacement, and the vertical velocity using the Cramer-Rao bound (Rosen et al., 2000; eqn 68).


    References:
    P. A. Rosen et al., "Synthetic aperture radar interferometry," in Proceedings of the IEEE, vol. 88, no. 3, pp. 333-382, March 2000, doi: 10.1109/5.838084.
    
    """
    # extract two time series
    profile1_unaligned = self.isel(time=slice(0,-offset))
    profile2_unaligned= self.isel(time=slice(offset,None))
    
    # compute the binned coherence, phase, and displacement, with uncertainties and put them into an xarray dataset  
    ds = compute_displacement(profile1_unaligned, 
                              profile2_unaligned, 
                              bin_size = bin_size,
                              lower_limit_on_fit = lower_limit_on_fit,
                              min_depth_for_ezz_fit = min_depth_for_ezz_fit,
                              max_depth_for_ezz_fit = max_depth_for_ezz_fit)

    # add attributes related to the this processing
    ds.attrs["offset"] = offset
    ds.attrs["processing"] = f"Created by the displacement_timeseries function in xapres using an offset of {offset} and bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def compute_displacement(profile1_unaligned: xr.DataArray, 
                        profile2_unaligned: xr.DataArray, 
                        bin_size: int=20, 
                        lower_limit_on_fit: float=None,
                        min_depth_for_ezz_fit: float=None,
                        max_depth_for_ezz_fit: float=None):
    """
    Compute displacement, coherence, velocity, strain rates, and related uncertainties from ApRES profiles.

    Parameters:
    - profile1_unaligned (xr.DataArray): Complex ApRES profiles from the first measurement (or measurements in each pair). This can be a time series and can have multiple attenuator settings. 
    - profile2_unaligned (xr.DataArray): Complex ApRES profiles from the second measurement (or measurements in each pair). This can be a time series and can have multiple attenuator settings. 
    - bin_size (int, optional): Size of the bins to use when binned the profile data vertically before comparing profiles. Default is 20.
    - min_depth_for_ezz_fit (float, optional): Upper limit on the fit for computing strain rates. Default is 0.0.
    - max_depth_for_ezz_fit (float, optional): Lower limit on the fit for computing strain rates. Default is 800.0.
    - lower_limit_on_fit (float, optional): depreciated: use max_depth_for_ezz_fit instead.  
    
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

    # compute phase variance
    phase_variance = ((1/(abs(coherence)**2))*np.sqrt((1 - abs(coherence)**2)/(2*bin_size))).rename('phase_variance')
    phase_variance.attrs["units"] = "radians^2"
    phase_variance.attrs["long_name"] = "variance in coherence phase"

    # compute the displacement
    displacement = phase2range(phase).rename("displacement")
    displacement.attrs["units"] = "m"
    displacement.attrs["long_name"] = "displacement since previous measurement"

    # compute the displacement variance
    disp_variance = (phase2range(np.sqrt(phase_variance))**2).rename('disp_variance')
    disp_variance.attrs["units"] = "m^2"
    disp_variance.attrs["long_name"] = "variance in displacement since previous measurement"

    dt_years = ((profiles.profile_time.sel(shot_number=2) - profiles.profile_time.sel(shot_number=1)) / np.timedelta64(1,'D') / 365.25).rename('dt_years')
    dt_years.attrs['units'] = 'years'
    dt_years.attrs['long_name'] = 'Time between shots'
    dt_years.attrs['description'] = 'Time in years between shots used in each measurement of displacement, vertical velocity, etc. dt_years[i] is the time between shot [j] and shot [j-1]'

    # vertical velocity
    velocity = (displacement / dt_years).rename('velocity')
    velocity.attrs['units'] = 'meters/year'
    velocity.attrs['long_name'] = 'Vertical velocity'

    # vertical velocity variance
    velocity_variance = ((np.sqrt(disp_variance) / dt_years)**2).rename('velocity_variance')
    velocity_variance.attrs["units"] = "meters^2/year^2"
    velocity_variance.attrs["long_name"] = "variance in vertical velocity"

    # strain rates
    strain_rates = computeStrainRates(xr.merge([velocity, velocity_variance]),
                                      max_depth_for_ezz_fit = max_depth_for_ezz_fit, 
                                      min_depth_for_ezz_fit = min_depth_for_ezz_fit, 
                                      lower_limit_on_fit = lower_limit_on_fit)

    # combine to an xarray dataset
    da_list = [profiles, coherence, phase, phase_variance, displacement, disp_variance, velocity, velocity_variance, strain_rates]
    ds = xr.merge(da_list)

    # add attributes related to this processing
    ds.attrs["bin_size"] = bin_size
    ds.attrs["description"] = "Time series of profiles of coherence, phase, displacement, and associated uncertainties, binned in depth."
    ds.attrs["processing"] = f"Created by the compute_displacement function in xapres using a bin size of {bin_size} on {datetime.datetime.now() }"

    return ds

def combine_profiles(profile1_unaligned, profile2_unaligned):
    """Align two timeseries of ApRES profiles, (in the case of unattended data) recording the midpoint time and the time of each computed profile"""
    
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
    """
    Bin ApRES profiles vertically, while computing the depth of the center of each bin.
    
    Parameters:
    - profiles (xr.DataArray): The input data array containing a time series of complex profiles.
    - bin_size (int): Size of the bins measured in number of data points.

    Returns:
    xr.DataArray: The binned profiles.
    
    """
    # Bin in depth
    profiles_binned = profiles.coarsen(profile_range=bin_size, boundary='trim').construct(profile_range=("bin_depth", "sample_in_bin"))

    # Compute the bin depth and add it to the DataArray
    bin_depth = profiles_binned.profile_range.mean(dim='sample_in_bin').data
    profiles_binned = profiles_binned.assign_coords(bin_depth=("bin_depth", bin_depth))

    return profiles_binned

def compute_coherence(p1, p2):
    """
    Compute the complex coherence between two binned time series of complex profiles. 

    The time series must have previously been binned using `bin_profiles`. 

    Parameters:
    - p1 (xr.DataArray): The first time series of complex profiles.
    - p2 (xr.DataArray): The second time series of complex profiles.

    Returns:
    xr.DataArray: The complex coherence between the two time series of profiles.
    
    """  

    top = (p1 * np.conj(p2)).sum(dim="sample_in_bin")
    bottom = np.sqrt( (np.abs(p1)**2).sum(dim="sample_in_bin") * (np.abs(p2)**2).sum(dim="sample_in_bin"))

    return (top/bottom).rename("coherence")

def phase2range(phi, 
                lambdac=0.5608):
    """
    Convert phase difference to range.

    Parameters:
    - phi (xr.DataArray): The phase difference computed in vertical bins between two profiles.
    - lambdac (float, optional): The wavelength at the center frequency. Default is 0.5608 m.

    Returns:
    xr.DataArray: The displacement corresponding to the supplied phase differences.

    Notes: 
    An alternative version of this code, which comes from the original matlab code, offers the choice of computing a 'precise' range. 
    It does this by taking account (approximately) of both constant terms in the phase difference equation (eqn 13 in Brennan et al. 2014).
    The third term in that equation is always negligible (see https://github.com/ldeo-glaciology/xapres/blob/master/notebooks/guides/coherence_range.ipynb). 
    So we neglect it here. For clarity the original code, adapted from Craifg Stewart, 2014/6/10, is included as follows:

    lambdac: float
        wavelength (m) at center frequency
    rc: float; optional, default is None
        coarse range of bin center (m)
    K:  float; optional, default is 2e8
        chirp gradient (rad/s/s)
    ci: float; optional, default is 1.6823e8
        propagation velocity (m/s)

    ```
    if not all([K,ci]) or rc is None:
        # First order method
        # Brennan et al. (2014) eq 15
        print('not precise')
        r = lambdac*phi/(4.*np.pi)
    else:
        print('Precise')
        # From Stewart (2018) eqn 4.8, with tau = 2*R/ci and omega_c = 2 pi /lambdac, where R is the range
        r = phi/((4.*np.pi/lambdac) - (4.*rc[None,:]*K/ci**2.))

    ```
    """
    return lambdac*phi/(4.*np.pi)

def computeStrainRates(ds: xr.Dataset,
                       lower_limit_on_fit: float=None,
                       min_depth_for_ezz_fit: float=None,
                       max_depth_for_ezz_fit: float=None):
    """Compute strain rates from a dataset of ApRES data along with the variance in the strain-rate estimates. For use by the function `compute_displacement`"""
    
    # parse inputs
    if lower_limit_on_fit is not None:
        print("lower_limit_on_fit is depreciated. Use max_depth_for_ezz_fit instead.")
        if max_depth_for_ezz_fit is None:
            max_depth_for_ezz_fit = lower_limit_on_fit           
        else:
            print(f"Because you also set the value of max_depth_for_ezz_fit (= {max_depth_for_ezz_fit}), this value will be used instead of lower_limit_on_fit.")

    if min_depth_for_ezz_fit is None:
        min_depth_for_ezz_fit = 0.0
    if max_depth_for_ezz_fit is None:
        max_depth_for_ezz_fit = 800.0

    if min_depth_for_ezz_fit > max_depth_for_ezz_fit:
        print("min_depth_for_ezz_fit is greater than max_depth_for_ezz_fit. Swapping the two.")
        min_depth_for_ezz_fit, max_depth_for_ezz_fit = max_depth_for_ezz_fit, min_depth_for_ezz_fit
    
    # crop profiles
    ds_cropped = ds.squeeze().sel(bin_depth = slice(min_depth_for_ezz_fit, max_depth_for_ezz_fit))

    # perform weighted least squares fit
    fit_ds = weighted_least_squares(ds_cropped)

    # add attributes
    fit_ds.strain_rate.attrs['units'] = '1/year'
    fit_ds.strain_rate.attrs['long_name'] = f"vertical strain rate computed between {min_depth_for_ezz_fit} and {max_depth_for_ezz_fit} m"
    fit_ds.strain_rate.attrs['max_depth_for_ezz_fit_meters'] = max_depth_for_ezz_fit
    fit_ds.strain_rate.attrs['min_depth_for_ezz_fit_meters'] = min_depth_for_ezz_fit

    fit_ds.strain_rate_variance.attrs['units'] = '1/year'
    fit_ds.strain_rate_variance.attrs['long_name'] = f"variance in vertical strain rate computed between {min_depth_for_ezz_fit} and {max_depth_for_ezz_fit} m"
    fit_ds.strain_rate_variance.attrs['max_depth_for_ezz_fit_meters'] = max_depth_for_ezz_fit
    fit_ds.strain_rate_variance.attrs['min_depth_for_ezz_fit_meters'] = min_depth_for_ezz_fit

    fit_ds.surface_intercept.attrs['units'] = 'meters/year'
    fit_ds.surface_intercept.attrs['long_name'] = 'vertical velocity at the surface from the linear fit'
    fit_ds.surface_intercept.attrs['max_depth_for_ezz_fit_meters'] = max_depth_for_ezz_fit
    fit_ds.surface_intercept.attrs['min_depth_for_ezz_fit_meters'] = min_depth_for_ezz_fit

    fit_ds.surface_intercept_variance.attrs['units'] = 'meters/year'
    fit_ds.surface_intercept_variance.attrs['long_name'] = 'variance in vertical velocity at the surface from the linear fit'
    fit_ds.surface_intercept_variance.attrs['max_depth_for_ezz_fit_meters'] = max_depth_for_ezz_fit
    fit_ds.surface_intercept_variance.attrs['min_depth_for_ezz_fit_meters'] = min_depth_for_ezz_fit

    fit_ds.r_squared.attrs['units'] = '-'
    fit_ds.r_squared.attrs['long_name'] = 'r-squared value for the weighted linear fit'

    fit_ds.sum_squared_residuals.attrs['units'] = 'm^2/yr^2'
    fit_ds.sum_squared_residuals.attrs['long_name'] = 'sum of squared residuals between the weighted linear fit and the velocities'


    return fit_ds
    #return xr.merge([strain_rate, strain_rate_uncertainty, surface_intercept, surface_intercept_uncertainty, R2, fit_ds.polyfit_residuals])

def compute_residuals_weighted(ds):
    """ Compute the weighted sum of squared residuals of the least squares fit."""
    residuals = ds.velocity - (ds.strain_rate * ds.bin_depth + ds.surface_intercept)
    return (residuals**2).weighted(1/ds.velocity_variance).sum(dim='bin_depth').rename('sum_squared_residuals')

def compute_r_squared(ds):
    """ Compute the r-squared value for the weighted least squares fit."""
    sum_of_square_residuals = compute_residuals_weighted(ds)
    y_mean = ds.velocity.weighted(1/ds.velocity_variance).mean(dim = 'bin_depth')
    SS_tot = ((ds.velocity - y_mean)**2).weighted(1/ds.velocity_variance).sum(dim = 'bin_depth')
    return  xr.merge([(1 - (sum_of_square_residuals/SS_tot)).rename('r_squared'), sum_of_square_residuals])

def weighted_least_squares(ds, deg = 1, cov = 'unscaled'):
    """
    Perform weighted least squares fits on velocity profiles. 

    Parameters:
    - ds (xr.Dataset): The input dataset containing the velocity profiles and associated uncertainties.
    - deg (int, optional): The degree of the polynomial to fit. Default is 1.
    - cov (str, optional): The type of covariance matrix to use. Default is 'unscaled'. This produces uncertainties that are true representations of the uncertainty. The other option (cov = "full") produces relative uncertainties. 

    Returns:
    xr.Dataset: The parameters of the polynomial fit, their variances, and the residuals.

    Notes: 
    The weighting assumes that the uncertainty computed by compute_displacement is the variance of the measurement. 
    As described in the numpy.polyfit documentation for the case when we want to use inverse-variance weighting, 
    the weights supplied to polyfit in this function are reciprocal of the **standard deviation**, not the variance itself.
    So, given that we have the variance, we need to supply 1/sqrt(variance) to polyfit. 
    """

    def my_polyfit(x, y, sigma, cov, deg=1):
        #p, residuals,  rank, singular_values, rcond = np.polyfit(x, y,  w=1/sigma, deg=deg, full=True)
        p, V = np.polyfit(x, y,  w=1/sigma, deg=deg, cov=cov)
        p_variances = np.diag(V)
        return p, p_variances

    res = xr.apply_ufunc(
        my_polyfit,
        ds.bin_depth,
        ds.velocity,
        (ds.velocity_variance)**0.5,
        input_core_dims=[["bin_depth"], ["bin_depth"], ["bin_depth"]],
        output_core_dims=[["degree"], ["degree"]],
        kwargs = {'deg': deg, 'cov': cov},
        vectorize=True,   
        dask='parallelized',
        dask_gufunc_kwargs = {'output_sizes':{'degree': 2}},
        output_dtypes = [np.dtype(np.float64), np.dtype(np.float64)]
    )

    strain_rate = res[0].sel(degree = 0, drop =True).rename('strain_rate')
    strain_rate_variance = res[1].sel(degree = 0, drop =True).rename('strain_rate_variance')

    surface_intercept =  res[0].sel(degree = 1, drop =True).rename('surface_intercept') 
    surface_intercept_variance = res[1].sel(degree = 1, drop =True).rename('surface_intercept_variance')

    #sum_squared_residuals = res[2].squeeze().rename("sum_squared_residuals_from_polyfit")

    R2 = compute_r_squared(xr.merge([ds.velocity, ds.velocity_variance, strain_rate, surface_intercept]))

    #out = [res[0].rename("parameters"), res[1].rename("parameter_uncertainty") , res[2].squeeze().rename("polyfit_residuals")]
    return xr.merge([strain_rate, strain_rate_variance, surface_intercept, surface_intercept_variance, R2])

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