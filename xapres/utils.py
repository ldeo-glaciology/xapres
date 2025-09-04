"""
Utility functions for performing common operations on ApRES radar data.

This module provides a comprehensive set of functions for analyzing and processing
Autonomous phase-sensitive Radio-Echo Sounder (ApRES) data stored in xarray datasets.
The functions focus on displacement analysis, strain rate calculations, and signal
processing operations commonly used in glaciological radar studies.

Key functionality includes:
- Displacement timeseries computation from profile pairs
- Coherence analysis between measurement epochs
- Strain rate estimation using weighted least squares fitting
- Phase-to-range conversion for displacement measurements
- Signal processing utilities (windowing, FFT, filtering)
- Data visualization helpers (dB conversion, sonification)
- Method binding to extend xarray DataArray and Dataset functionality

The module follows a functional programming approach with functions designed to work
seamlessly with xarray data structures. Many functions are automatically bound as
methods to xarray objects for convenient access (e.g., `dataarray.displacement_timeseries()`).

Mathematical foundations are based on radar interferometry principles and the
Cramer-Rao bound for uncertainty estimation, following established practices in
ice-penetrating radar analysis.
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
    Compute displacement timeseries and strain rates from complex ApRES profiles.

    This function performs radar interferometry analysis on a time series of complex
    ApRES profiles to compute ice displacement, vertical velocities, and strain rates.
    It uses phase differences between profile pairs separated by a user-defined time
    offset to estimate ice motion.

    The analysis follows established radar interferometry principles, computing coherence,
    phase differences, and their associated uncertainties. Displacement measurements are
    binned vertically and used to estimate strain rates through weighted least squares
    fitting.

    Parameters
    ----------
    self : xarray.DataArray
        Time series of complex ApRES profiles with dimensions including 'time' and
        'profile_range'. Should contain complex-valued radar returns from FFT processing.
    offset : int, optional
        Number of time steps between profile pairs used for displacement calculation.
        Larger offsets provide longer baselines but reduce temporal resolution.
        Default is 1.
    bin_size : int, optional
        Number of range samples to bin together vertically before analysis. Larger
        bins improve coherence but reduce spatial resolution. Default is 20.
    min_depth_for_ezz_fit : float, optional
        Minimum depth in meters for strain rate fitting. Default is 0.0 (surface).
    max_depth_for_ezz_fit : float, optional
        Maximum depth in meters for strain rate fitting. Default is 800.0.
    lower_limit_on_fit : float, optional
        **Deprecated**: Use max_depth_for_ezz_fit instead. Maintained for backward
        compatibility.

    Returns
    -------
    xarray.Dataset
        Dataset containing displacement analysis results with the following variables:
        
        - coherence : Complex coherence between profile pairs
        - phase : Phase difference in radians  
        - displacement : Displacement in meters
        - velocity : Vertical velocity in meters/year
        - strain_rate : Vertical strain rate in 1/year
        - phase_variance : Phase uncertainty in radians²
        - disp_variance : Displacement uncertainty in meters²
        - velocity_variance : Velocity uncertainty in (meters/year)²
        - strain_rate_variance : Strain rate uncertainty in (1/year)²
        - surface_intercept : Surface velocity from linear fit in meters/year
        - r_squared : R² value for strain rate fit

    Examples
    --------
    Basic displacement analysis with default parameters:
    
    >>> profiles = dataset.profile.isel(chirp_num=0, attenuator_setting_pair=0)
    >>> displacement_data = profiles.displacement_timeseries()
    >>> print(displacement_data.displacement.dims)
    
    Custom analysis with larger offset and different bin size:
    
    >>> displacement_data = profiles.displacement_timeseries(
    ...     offset=3, 
    ...     bin_size=30,
    ...     max_depth_for_ezz_fit=500
    ... )
    >>> displacement_data.strain_rate.plot()
    
    High-resolution analysis with small bins:
    
    >>> displacement_data = profiles.displacement_timeseries(
    ...     bin_size=10,
    ...     min_depth_for_ezz_fit=50,
    ...     max_depth_for_ezz_fit=300
    ... )

    Notes
    -----
    The method implements radar interferometry following Rosen et al. (2000), using
    the Cramer-Rao bound to estimate phase and displacement uncertainties. The time
    interval between output profiles is offset × dt, where dt is the measurement
    interval.

    Strain rates are computed using weighted least squares fitting of the velocity
    profile, with weights derived from displacement uncertainties. The fitting range
    can be controlled via the depth limit parameters.

    References
    ----------
    P. A. Rosen et al., "Synthetic aperture radar interferometry," Proceedings of 
    the IEEE, vol. 88, no. 3, pp. 333-382, March 2000, doi: 10.1109/5.838084.
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
    Compute displacement and strain analysis from pairs of complex ApRES profiles.

    This function performs radar interferometry analysis between two sets of complex 
    ApRES profiles to determine ice displacement, vertical velocities, coherence, and 
    strain rates. It implements standard radar interferometry techniques with uncertainty
    estimation based on the Cramer-Rao bound.

    The analysis bins the profiles vertically for improved coherence estimation, computes
    phase differences and converts them to displacement measurements, then performs
    weighted least squares fitting to estimate strain rates over a specified depth range.

    Parameters
    ----------
    profile1_unaligned : xarray.DataArray
        Complex ApRES profiles from the first measurement epoch. Should have dimensions
        including 'profile_range' and may include 'time' for multiple measurements.
        Values should be complex-valued radar returns from profile computation.
    profile2_unaligned : xarray.DataArray  
        Complex ApRES profiles from the second measurement epoch, with the same
        structure as profile1_unaligned. Time difference between epochs determines
        the displacement sensitivity.
    bin_size : int, optional
        Number of range samples to average vertically before coherence analysis.
        Larger bins improve coherence estimates but reduce spatial resolution.
        Default is 20.
    min_depth_for_ezz_fit : float, optional
        Minimum depth in meters for strain rate fitting. Default is 0.0 (surface).
    max_depth_for_ezz_fit : float, optional
        Maximum depth in meters for strain rate fitting. Default is 800.0.
    lower_limit_on_fit : float, optional
        **Deprecated**: Use max_depth_for_ezz_fit instead. Maintained for backward
        compatibility only.

    Returns
    -------
    xarray.Dataset
        Comprehensive dataset containing displacement analysis results:
        
        **Primary measurements:**
        - coherence : Complex coherence between profile pairs (dimensionless)
        - phase : Phase difference in radians
        - displacement : Displacement in meters
        - velocity : Vertical velocity in meters/year
        
        **Uncertainty estimates:**
        - phase_variance : Phase uncertainty in radians²
        - disp_variance : Displacement uncertainty in meters²  
        - velocity_variance : Velocity uncertainty in (meters/year)²
        
        **Strain analysis:**
        - strain_rate : Vertical strain rate in 1/year
        - strain_rate_variance : Strain rate uncertainty in (1/year)²
        - surface_intercept : Surface velocity from linear fit in meters/year
        - surface_intercept_variance : Surface velocity uncertainty in (meters/year)²
        - r_squared : Coefficient of determination for strain fit
        - sum_squared_residuals : Weighted residuals from strain fit
        
        **Coordinates:**
        - bin_depth : Depth to center of each bin in meters
        - time : Midpoint time between measurements (for time series)
        - profile_time : Time of each individual profile

    Examples
    --------
    Basic displacement analysis between two profile sets:
    
    >>> p1 = dataset.profile.isel(time=slice(0, 5))
    >>> p2 = dataset.profile.isel(time=slice(1, 6))  
    >>> displacement_data = compute_displacement(p1, p2)
    >>> displacement_data.displacement.plot()
    
    Analysis with custom binning and depth range:
    
    >>> displacement_data = compute_displacement(
    ...     p1, p2,
    ...     bin_size=30,
    ...     min_depth_for_ezz_fit=100,
    ...     max_depth_for_ezz_fit=500
    ... )
    >>> print(f"Strain rate: {displacement_data.strain_rate.values:.2e} /year")
    
    Single profile pair analysis:
    
    >>> p1 = dataset.profile.isel(time=0)
    >>> p2 = dataset.profile.isel(time=5)
    >>> displacement_data = compute_displacement(p1, p2, bin_size=25)

    Notes
    -----
    The function performs several key processing steps:
    
    1. **Profile alignment**: Handles temporal alignment and coordinate matching
    2. **Vertical binning**: Averages profiles in depth bins for coherence analysis
    3. **Coherence computation**: Calculates complex coherence between binned profiles
    4. **Phase analysis**: Extracts phase differences and estimates uncertainties
    5. **Displacement conversion**: Converts phase to displacement using radar wavelength
    6. **Velocity calculation**: Derives velocities from displacement and time intervals  
    7. **Strain fitting**: Performs weighted least squares to estimate strain rates

    Uncertainty estimation follows the Cramer-Rao bound theory, providing realistic
    error bounds on displacement and derived quantities. The weighting in strain rate
    fitting uses inverse variance weighting for optimal parameter estimation.

    The depth range for strain fitting should be chosen to represent a region of
    approximately linear velocity variation. Ice near the surface or bed may show
    non-linear behavior that violates the linear strain assumption.

    Raises
    ------
    TypeError
        If input profiles are not xarray DataArrays
    ValueError
        If profiles have incompatible dimensions or coordinates
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
    Convert radar interferometric phase differences to displacement ranges.

    This function converts phase differences measured between radar profiles into
    physical displacement distances. It implements the first-order approximation
    commonly used in radar interferometry, which is accurate for most ApRES
    applications.

    The conversion is based on the fundamental relationship between phase change
    and path length difference in radar interferometry, accounting for the two-way
    travel path and the wavelength of the radar signal.

    Parameters
    ----------
    phi : array-like
        Phase differences in radians computed between radar profile pairs.
        Can be xarray.DataArray, numpy array, or scalar values.
    lambdac : float, optional
        Wavelength at the center frequency in meters. For ApRES systems operating
        at 200-400 MHz, the default value of 0.5608 m corresponds to the center
        frequency of 300 MHz in ice. Default is 0.5608.

    Returns
    -------
    array-like
        Displacement corresponding to the input phase differences, in meters.
        Positive values indicate motion away from the radar (typically ice
        thickening or surface lowering).

    Examples
    --------
    Convert phase differences to displacement:
    
    >>> phase_diff = np.array([0.1, 0.5, 1.0])  # radians
    >>> displacement = phase2range(phase_diff)
    >>> print(f"Displacements: {displacement} m")
    
    Use with xarray data:
    
    >>> phase = coherence_data.phase
    >>> displacement = phase2range(phase, lambdac=0.56)
    >>> displacement.plot()
    
    Custom wavelength for different frequency:
    
    >>> # For 250 MHz center frequency
    >>> lambda_250MHz = 3e8 / (250e6 * np.sqrt(3.18))  # ~0.67 m in ice
    >>> displacement = phase2range(phase_diff, lambdac=lambda_250MHz)

    Notes
    -----
    The function implements the first-order approximation from Brennan et al. (2014):
    
    .. math:: r = \\frac{\\lambda_c \\phi}{4\\pi}
    
    where:
    - r is the displacement (range change)
    - λc is the wavelength at center frequency  
    - φ is the phase difference
    - The factor of 4π accounts for two-way travel and the 2π phase cycle
    
    This first-order method is accurate for typical ApRES applications. A more
    precise method exists that accounts for additional terms in the phase equation,
    but the improvement is negligible for most use cases (see referenced notebook).
    
    The default wavelength assumes:
    - Center frequency of 300 MHz
    - Ice permittivity of 3.18
    - λc = c / (f_c * √εr) where c is speed of light in vacuum
    
    For different radar systems or ice properties, adjust lambdac accordingly.

    References
    ----------
    Brennan, P. V., et al. "Phase-sensitive FMCW radar system for high-precision
    Antarctic ice shelf profile monitoring." IET Radar, Sonar & Navigation 8.7
    (2014): 776-786.
    
    See Also
    --------
    compute_displacement : Full displacement analysis including uncertainty
    displacement_timeseries : Time series displacement analysis
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
    Convert complex profile data to decibel representation.
    
    This function converts the amplitude of complex radar profiles to a logarithmic
    decibel scale, which is commonly used for displaying radar data. The conversion
    makes it easier to visualize the dynamic range of backscatter intensities.
    
    This function is automatically bound as a method to xarray DataArrays, so it
    can be called directly on profile data as `dataarray.dB()`.
    
    Parameters
    ----------
    self : xarray.DataArray
        Complex-valued radar profile data. The function uses the absolute value
        before logarithmic conversion.
    
    Returns
    -------
    xarray.DataArray
        Profile data converted to decibels (dB). Values represent 20*log10(|amplitude|)
        following standard radar conventions.
        
    Examples
    --------
    Convert profile data to dB for visualization:
    
    >>> profile_dB = dataset.profile.dB()
    >>> profile_dB.plot()
    
    Compare linear and logarithmic representations:
    
    >>> fig, (ax1, ax2) = plt.subplots(1, 2)
    >>> dataset.profile.abs().plot(ax=ax1, y='profile_range')
    >>> dataset.profile.dB().plot(ax=ax2, y='profile_range')
    
    Notes
    -----
    The decibel conversion uses the standard formula for power quantities:
    dB = 20 * log10(|amplitude|)
    
    This is appropriate for radar amplitude data where the power is proportional
    to the square of the amplitude. The factor of 20 (rather than 10) accounts
    for the amplitude-to-power relationship.
    
    See Also
    --------
    numpy.log10 : Underlying logarithm function
    numpy.abs : Absolute value function for complex data
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
    """
    Add computed radar profiles to an existing ApRES dataset.

    This function computes radar profiles from the chirp data in a dataset and merges
    them with the existing data variables. It provides a convenient way to add profile
    data to datasets that were loaded without profile computation or to recompute
    profiles with different parameters.

    The function automatically uses radar constants stored in the dataset attributes
    if available, ensuring consistency with the original data processing parameters.

    Parameters
    ----------
    self : xarray.Dataset
        ApRES dataset containing chirp data with 'chirp' data variable. Should have
        appropriate dimensions including 'chirp_time'.
    **kwargs
        Additional keyword arguments passed to the computeProfile function. Can include
        parameters like 'pad_factor', 'demean', 'max_range', 'drop_noisy_chirps', etc.
        See computeProfile documentation for full parameter list.

    Returns
    -------
    xarray.Dataset
        Dataset with profile data added. If profiles already exist, the 'profile_range'
        dimension is dropped first and then new profiles are merged. The original
        dataset structure and attributes are preserved.

    Examples
    --------
    Add profiles to a dataset loaded without them:
    
    >>> # Load data without profile computation
    >>> fd = from_dats()
    >>> ds = fd.load_all(directory='data/', computeProfiles=False)
    >>> # Add profiles later
    >>> ds_with_profiles = ds.addProfileToDs()
    
    Recompute profiles with different parameters:
    
    >>> ds_new = ds.addProfileToDs(
    ...     pad_factor=4,
    ...     demean=True,
    ...     max_range=500,
    ...     drop_noisy_chirps=True
    ... )
    
    Add profiles with custom processing:
    
    >>> ds_processed = ds.addProfileToDs(
    ...     crop_chirp_start=0.1,
    ...     crop_chirp_end=0.9,
    ...     detrend=True,
    ...     stack=True
    ... )

    Notes
    -----
    - If the dataset contains radar constants in its attributes, these are automatically
      passed to the profile computation function for consistency
    - If 'profile' data already exists, it is replaced with newly computed profiles
    - All original data variables, coordinates, and attributes are preserved
    - The function uses the 'override' combine_attrs strategy to handle attribute conflicts

    See Also
    --------
    computeProfile : Underlying function for profile computation
    from_dats.load_all : Load data with profile computation included
    """

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
    Compute complex radar profiles from raw chirp data using FFT processing.

    This function converts time-domain chirp signals from ApRES into frequency-domain
    radar profiles through a series of signal processing steps. The processing pipeline
    includes windowing, padding, FFT transformation, and range referencing to produce
    complex profiles suitable for interferometric analysis.

    The function implements the standard processing chain for FMCW (Frequency Modulated
    Continuous Wave) radar data, converting beat frequencies in the de-ramped chirps
    into range-resolved reflectivity profiles.

    Parameters
    ----------
    self : xarray.DataArray
        Input chirp data array with dimensions including 'chirp_time'. Should contain
        de-ramped chirp signals with voltage values from the radar receiver.
    pad_factor : int, optional
        Zero-padding factor applied before FFT to improve range resolution.
        Higher values increase computational cost but improve interpolation.
        Default is 2.
    drop_noisy_chirps : bool, optional
        Whether to exclude chirps that exceed the clipping threshold. Helps remove
        corrupted or saturated measurements. Default is False.
    clip_threshold : float, optional
        Voltage threshold for identifying noisy/clipped chirps. Only used when
        drop_noisy_chirps is True. Default is 1.2 volts.
    min_chirps : int, optional
        Minimum number of valid chirps required to retain a burst after noise
        filtering. Default is 0 (no minimum).
    demean : bool, optional
        Whether to remove the DC component (mean) from each chirp before processing.
        Recommended for most applications. Default is True.
    detrend : bool, optional
        Whether to remove linear trends from chirp data. Can help with drift issues
        but may affect signal content. Default is False.
    stack : bool, optional
        Whether to average all chirps within a burst before FFT processing.
        Improves SNR but removes information about individual chirps. Default is False.
    scale_for_window : bool, optional
        Whether to apply amplitude scaling to compensate for windowing effects.
        Default is True. Set to False only for comparison with legacy methods.
    crop_chirp_start : float, optional
        Start time in seconds for cropping chirps before processing. Default is 0
        (use full chirp from beginning).
    crop_chirp_end : float, optional
        End time in seconds for cropping chirps. Default is 1 (use full chirp).
        Values < 1 crop the chirp duration.
    max_range : float, optional
        Maximum range in meters for the output profile. If None, uses half the
        full range to avoid aliasing. Default is None.
    constants : dict, optional
        Dictionary of radar system constants to override defaults. Can include
        'c' (speed of light), 'ep' (ice permittivity), 'K' (chirp rate), etc.
        Missing constants are filled with default values.

    Returns
    -------
    xarray.DataArray
        Complex radar profile with 'profile_range' coordinate in meters.
        The profile contains complex reflectivity values where amplitude represents
        backscatter strength and phase contains interferometric information.
        
        Attributes include:
        - long_name : 'profile'
        - units : '-' (dimensionless)
        - description : Processing description
        - constants : Dictionary of radar constants used

    Examples
    --------
    Basic profile computation with default settings:
    
    >>> chirps = dataset.chirp.isel(time=0, attenuator_setting_pair=0)
    >>> profile = chirps.computeProfile()
    >>> profile.abs().plot()  # Plot amplitude
    
    High-resolution processing with more padding:
    
    >>> profile = chirps.computeProfile(
    ...     pad_factor=4,
    ...     demean=True,
    ...     max_range=800
    ... )
    
    Quality-controlled processing:
    
    >>> profile = chirps.computeProfile(
    ...     drop_noisy_chirps=True,
    ...     clip_threshold=1.0,
    ...     min_chirps=15,
    ...     stack=True
    ... )
    
    Custom radar constants:
    
    >>> constants = {'c': 3e8, 'ep': 3.2, 'f_c': 300e6}
    >>> profile = chirps.computeProfile(constants=constants)

    Notes
    -----
    The processing follows these main steps:
    1. Chirp cropping and quality filtering
    2. DC removal and/or detrending if requested
    3. Windowing with Blackman window
    4. Zero-padding to improve resolution
    5. FFT transformation to frequency domain
    6. Range referencing and phase correction
    7. Conversion to range coordinates
    8. Cropping to specified maximum range

    The range calculation assumes a linear frequency sweep and uses the radar
    equation for FMCW systems. The complex output preserves both amplitude and
    phase information needed for interferometric analysis.

    For large datasets, consider using dask arrays to manage memory usage during
    processing. The function automatically handles dask arrays and loads them
    before FFT operations.
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
    """
    Return default radar system constants for ApRES profile computation.

    This function provides standard parameter values for ApRES radar systems,
    which are used in profile computation when specific constants are not
    provided. These defaults correspond to typical ApRES configurations and
    can be overridden by supplying custom values.

    Returns
    -------
    dict
        Dictionary containing radar system constants with the following keys:
        
        **Frequency parameters:**
        - 'f_1' : Starting frequency in Hz (200 MHz)
        - 'f_2' : Ending frequency in Hz (400 MHz)  
        - 'f_c' : Center frequency in Hz (300 MHz)
        - 'B' : Bandwidth in Hz (200 MHz)
        
        **Timing parameters:**
        - 'T' : Chirp duration in seconds (1 s)
        - 'K' : Frequency sweep rate in Hz/s (2×10⁸ Hz/s)
        - 'dt' : Sampling time step in seconds (25 μs)
        
        **Physical constants:**
        - 'c' : Speed of light in vacuum in m/s (3×10⁸ m/s)
        - 'ep' : Relative permittivity of ice (3.18)

    Examples
    --------
    Get default constants for profile computation:
    
    >>> constants = default_constants()
    >>> print(f"Center frequency: {constants['f_c']/1e6:.0f} MHz")
    >>> print(f"Bandwidth: {constants['B']/1e6:.0f} MHz")
    
    Override specific values:
    
    >>> constants = default_constants()
    >>> constants['ep'] = 3.2  # Different ice permittivity
    >>> constants['c'] = 3e8   # Explicit scientific notation
    >>> profile = chirps.computeProfile(constants=constants)
    
    Use in profile computation:
    
    >>> # These are equivalent:
    >>> profile1 = chirps.computeProfile()  # Uses defaults
    >>> profile2 = chirps.computeProfile(constants=default_constants())

    Notes
    -----
    The default values are based on standard ApRES system specifications:
    
    - **Frequency range**: 200-400 MHz is typical for ApRES systems
    - **Chirp duration**: 1 second provides good range resolution
    - **Sampling rate**: 40 kHz (dt = 25 μs) is standard for ApRES
    - **Ice permittivity**: 3.18 is a commonly used value for glacier ice
    
    These values can be modified for:
    - Different radar systems or configurations
    - Varying ice conditions or compositions
    - Custom processing requirements
    - Comparison with historical data using different parameters
    
    The sweep rate K is calculated as K = B/T, representing the linear frequency
    change rate during the chirp. All frequency and timing relationships assume
    a linear frequency modulation (LFM) chirp waveform.

    See Also
    --------
    computeProfile : Function that uses these constants
    phase2range : Uses wavelength derived from these constants
    """
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
    """
    Bind xapres utility functions as methods to xarray DataArray and Dataset classes.

    This function extends xarray objects by adding custom methods from the xapres
    package, allowing for convenient access to ApRES-specific functionality directly
    on data objects. This enables method chaining and a more intuitive API.

    The function is called automatically when xapres is imported, so users don't
    need to call it explicitly. After import, the bound methods become available
    on all xarray objects in the session.

    Methods Added
    -------------
    **DataArray methods:**
    - dB() : Convert profile amplitudes to decibel scale
    - sonify() : Convert chirp data to audio playback  
    - displacement_timeseries() : Compute displacement analysis from profiles
    - computeProfile() : Generate radar profiles from chirp data
    - computeStrainRates() : Calculate strain rates from displacement data
    
    **Dataset methods:**
    - addProfileToDs() : Add computed profiles to existing datasets

    Examples
    --------
    After importing xapres, these methods become available:
    
    >>> import xapres
    >>> # Load some data
    >>> ds = xapres.load.generate_xarray('data/')
    >>> 
    >>> # Use bound methods directly
    >>> profile_db = ds.profile.dB()
    >>> displacement_data = ds.profile.displacement_timeseries()
    >>> ds_with_new_profiles = ds.addProfileToDs(pad_factor=4)
    
    Method chaining example:
    
    >>> processed_data = (ds.chirp
    ...                   .computeProfile(demean=True)
    ...                   .displacement_timeseries(bin_size=30))

    Notes
    -----
    This approach follows the xarray accessor pattern but uses direct method binding
    for simplicity. The bound methods maintain their original signatures and
    documentation, with the 'self' parameter automatically passed as the first
    argument.

    The binding is performed using Python's setattr() function on the class objects,
    making the methods available to all instances of DataArray and Dataset created
    after the xapres import.

    Potential conflicts with other packages that bind methods with the same names
    are resolved on a last-import-wins basis. If you encounter method conflicts,
    you can always access the original functions directly from the xapres.utils
    module.

    See Also
    --------
    Individual function documentation for detailed parameter descriptions
    xarray.accessor : Alternative approach for extending xarray functionality
    """
    
    da_methods = [dB, sonify, displacement_timeseries, computeProfile, computeStrainRates]
    for method in da_methods:
        setattr(xr.DataArray, method.__name__, method)

    ds_methods = [addProfileToDs]
    for method in ds_methods:
        setattr(xr.Dataset, method.__name__, method)

def contains_dask_array(dataarray):
    return isinstance(dataarray.data, da.Array)