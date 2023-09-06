# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 20:06:36 2020 by K. Nicholls
-- Modified on Nov 10th 2022 by J. Kingslake to include xapres class


Class definitions for ApRES processing code

xapres
==============
Instantiated with 2 optional keyword arguments, loglevel and max_range

Argument:
loglevel --- allows the user to select the level of logging messages are displayed. 
The default loglevel is warning, which means that no messages are displayed. 
If you want to see detailed log messages, use loglevel = 'debug'

max_range --- the depth the computed profiles are clipped to. 

Methods:
load_single --- load a single chirp from a single burst from a single dat file
load_dat_file --- load a dat file as a DataFileObject
list_files --- recursively find  all the files in a directory or a google bucket
load_all --- load all the files found in a directory or google bucket into an xarray

load_all is the most important method. Call it, for example, as follows:

import ApRESDefs
xa = ApRESDefs.xapres(loglevel='debug', max_range=1400)
xa.load_all(directory='gs://ldeo-glaciology/GL_apres_2022', 
            remote_load = True,
            file_numbers_to_process = [0, 1], 
            bursts_to_process=[0, 1]
           )

the resulting xarray will be saved in xa.data.

DataFileObject
==============
Instantiated with one required string, giving the name of the data file
eg fileDescriptor = DataFileObject('DATAFILENAME.DAT')

Methods:
    ExtractBurst(BurstNumber (integer))
        Output is an instance of a BurstObject
eg Burst = fileDescriptor.ExtractBurst(3)

Instance variables:
    Filename           : Name of data file
    BurstLocationList  : Python list of byte offset of each burst in file
    NoBurstsInFile     : Number of bursts in the file (len(BurstLocationList))
    
BurstObject.
============
Typically instantiated with a call to the ExtractBurst method on a DataFileObject
eg Burst = fileDescriptor.ExtractBurst(3)

Methods:
    ExtractChirp(ChirpList (Python list))
        Output is an instance of a ChirpObject, in which all chirps in the ChirpList
        have been averaged
    PlotBurst()
        Plots the full raw burst as a time series

Instance variables:
    v         : Array containing burst data
    Filename  : Name of data file
    Header    : Burst header (Python dictionary), with additional entries:
    BurstNo   : Burst number in data file

ChirpObject
===========
Typically instantiated with a call to the ExtractChirp method on a BurstBbject
eg Chirp = Burst.ExtractChirp([1,3])

Methods:
    FormProfile(StartFreq, StopFreq, padfactor, ref)
        StartFreq, StopFreq: start and end frequencies to use (eg 2e8 and 4e8)
        padfactor:           zero padding for the fft (eg. 2)
        ref:                 whether or not to apply Paul Brennan's reference
                             phase (1 or 0, for yes or no)
        Returns and instance of a ProfileObject
    PlotChirp()
        Plots the chirp as function of time

Instance variables:
    vdat:       Array containing chirp data
    t:          Array containing time for chirp samples
    ChirpList:  List of chirps averaged to make vdat
    Filename:   Name of data file
    BurstNo:    Number of burst within data file
    Header:     Burst header, as created by ExtractBurst method on FileDataObject

ProfileObject.
==============
Typically instantiated with a call to the FormProfile method on a ChirpObject
eg Profile = Chirp.FormProfile(StartFreq, StopFreq, padfactor, ref)

Methods:
    PlotProfile(MaxDepth (double))
        MaxDepth:  Maximum depth (in metres) to which to plot profile
        
Instance variables:
    Range:     Array with depth in metres each profile depth bin
    Profile:   Array containing profile (complex double)
    F0:        Start frequency used to form profile
    F1:        End frequency used to form profile
    pad:       Pad factor used when zeropadding
    ChirpList: List of chirps averaged to form profile
    Filename:  Name of original data file 
    BurstNo:   Number of burst in data file
    Header:    Burst header, as produced using ExtractBurst method on DataFileObject 
    rad2m:     radians to metres of range conversion factor
    bin2m:     bin to metres of range conversion factor
"""
import gcsfs
import numpy as np
import matplotlib.pyplot as plt
import math
import warnings
import copy
import numpy as np
import sys
import pandas as pd
import xarray as xr
from tqdm import tqdm
import glob
import os
import logging
from tqdm.notebook import trange, tqdm
sys.path.append(os.path.join(os.path.dirname(__file__), "lib"))
from utils import *


def load_zarr(site = "A101", directory = "gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/"):
    """Load ApRES data stored in a zarr directory as an xarray and add functionality"""
    
    import numpy as np
    import xarray as xr
        
    ds = xr.open_dataset(directory + site,
        engine = 'zarr', 
        chunks = {}) 
    
    # add db function as new bound method of DataArrays
    xr.DataArray.db = lambda self : 20*np.log10(np.abs(self))
    
    
    # add sonify function (only if soundfile and sounddevice are installed)
    def sonify(self, play = True, save = False, wav_filename = "chirp"):

    
        # make sure the input is just one chirp    
        if self.size != self.chirp_time.size: 
            raise BaseException('sonify only works for single chirps.')    

        #cut out the start and end to record popping
        chirp = self.isel(chirp_time =slice(5000,-500))


        if play:
            samplerate = chirp.chirp_time.size / (   (chirp.chirp_time[-1] - chirp.chirp_time[0]) /np.timedelta64(1, 's'))
            sd.play(chirp, samplerate = samplerate )

        if save:
            sf.write(F"{wav_filename} .wav", chirp, samplerate = samplerate)

 
    try:     
        import soundfile as sf
        import sounddevice as sd

        xr.DataArray.sonify = sonify 
    except ImportError:
        print("sounddevice and soundfile are required to sonify the chirps. pip install them if you need this feature") 
    
    
    return ds




class xapres():
    def __init__(self, loglevel='warning', max_range = None):
        self._setup_logging(loglevel)
        self.max_range = max_range
        
    def load_single(self, dat_filename, remote_load=False, burst_number=0, chirp_num=0):
        """Load a single chirp, from a single burst, from a single data file."""
        
        self.files_to_be_processed = dat_filename
        self.logger.debug(f"Load dat file {dat_filename} with remote_load = {remote_load}")
        self.single_dat = self.load_dat_file(dat_filename,remote_load)
        self.logger.debug(f"Extract burst number {burst_number}")
        self.single_burst = self.single_dat.ExtractBurst(burst_number)
        self.logger.debug(f"Extract chirp number {chirp_num}")
        self.single_chirp = self.single_burst.ExtractChirp([chirp_num]) 
        self.logger.debug(f"Form profile for chirp number {chirp_num}")
        self.single_profile = self.single_chirp.FormProfile()
        self.logger.debug(f"Finish call to load_single.")
    
    def load_dat_file(self, dat_filename, remote_load=False):
        """Return a DataFileObject, given a filename."""
        return DataFileObject(dat_filename,remote_load)
    
    
    def list_files(self, directory=None, remote_load=False):    
        """Recursively list all the .DAT files in a given location dir. 
        
        Arguments:
        dir -- the directory that will be looked in recursivly to find .DAT files.
        remote_load -- True/False, indicating if the data are stored locally or rmeotely. Default is False.
        """

        self.logger.debug(f"Find all the dat files in the directory {directory} with remote_load = {remote_load}")

        
        if directory is None:
            directory = os.getcwd()
                        
        if remote_load:
            fs = gcsfs.GCSFileSystem()
            dat_filenames = fs.glob(directory + '/**/*.[dD][aA][tT]', recursive = True)
        else:
            dat_filenames = glob.glob(directory + '/**/*.[dD][aA][tT]',recursive = True)
        
        self.dat_filenames = dat_filenames
        
        self.logger.debug(f"Finish call to list_files. Found {len(dat_filenames)} files")

        
        return dat_filenames
       
    
    def load_all(self, directory=None, 
                 remote_load=False, 
                 file_numbers_to_process = None, 
                 file_names_to_process = None, 
                 bursts_to_process = "All"):
        """Put all the data from all the .DAT files found recursively in 'directory', in one xarray."""
        
        
        self.file_numbers_to_process = file_numbers_to_process
        self.file_names_to_process = file_names_to_process
        #self.bursts_to_process = bursts_to_process
       
        ###### List files ######
        self.list_files(directory, remote_load)    # adds self.dat_filenames
    
    
        ###### Subset files ######
        if file_numbers_to_process is not None and file_names_to_process is not None:
            self.logger.debug("Throwing a ValueError because file_numbers_to_process and file_names_to_process cannot both be supplied to load_all")
            raise ValueError("file_numbers_to_process and file_names_to_process cannot both be supplied to load_all. You need to supply just one (or neither) of these.") 
        
        elif file_numbers_to_process is not None:
            if file_numbers_to_process == "All":
                self.logger.debug("Selecting all dats file because file_numbers_to_process == \"all\"")
                self.dat_filenames_to_process = self.dat_filenames
            else:
                self.logger.debug(f"Subset files to {file_numbers_to_process}")
                self.dat_filenames_to_process = [self.dat_filenames[i] for i in file_numbers_to_process]
        
        elif file_names_to_process is not None:
            if file_names_to_process == "All":
                self.logger.debug("Selecting all dats file because file_names_to_process == \"all\"")
                self.dat_filenames_to_process = self.dat_filenames
            else:                 
                self.logger.debug("Subset files to list of files supplied in file_names_to_process")
                self.dat_filenames_to_process = file_names_to_process
                              
        elif file_numbers_to_process is None and file_names_to_process is None:      # default is all the dat files    
            self.logger.debug("Selecting all dats file because neither file_numbers_to_process nor file_names_to_process were supplied")
            self.dat_filenames_to_process = self.dat_filenames          
        
        ## loop through the dat files, putting individual xarrays in a list
        self.logger.debug("Starting loop over dat files")
        list_of_multiBurstxarrays = []   
        for dat_filename in self.dat_filenames_to_process:
            self.logger.debug(f"Load dat file {dat_filename}")


            dat = self.load_dat_file(dat_filename,remote_load)
            
            multiBurstxarray = self._all_bursts_in_dat_to_xarray(dat,bursts_to_process)
        
            list_of_multiBurstxarrays.append(multiBurstxarray)
            self.logger.debug(f"Finished processing file {dat_filename}")
        
        self.logger.debug(f"Concatenating all the multi-burst xarrays to create xapres.data")
        # concatenate all the xarrays in the list along the time dimension
        self.data = xr.concat(list_of_multiBurstxarrays,dim='time')     
        
        self._add_attrs()
        
        self.logger.debug(f"Finish call to load_all. Call xapres.data to see the xarray this produced.")

    def _all_bursts_in_dat_to_xarray(self,dat,bursts_selected):
        """Take data from all the bursts in .DAT file and put it in an xarray.
        
        Arguments:
        dat -- a DataFileObject  
        """
        self.logger.debug(f"This dat file has {dat.NoBurstsInFile} bursts.")
        self.logger.debug(f"bursts_to_process = {bursts_selected} at the start of _all_bursts_in_dat_to_xarray.")
        # choose which bursts to process (default is all of the burst in each dat file)
        if bursts_selected == "All":
            bursts_to_process = range(dat.NoBurstsInFile)
            self.logger.debug("bursts_to_process set to \"All\"")
        else:
            bursts_to_process = bursts_selected
        
        if any(np.array(bursts_to_process) > (dat.NoBurstsInFile-1)):
            self.logger.debug(f"The burst numbers requested in bursts_to_process ({bursts_to_process}) is greater than the number of bursts in\
                the dat file ({dat.NoBurstsInFile}), so we will just process all the bursts.")
            bursts_to_process = range(dat.NoBurstsInFile)

        self.logger.debug(f"bursts_to_process = {list(bursts_to_process)} after initial parse in _all_bursts_in_dat_to_xarray.")

        self.logger.debug(f"Start loop over burst numbers {list(bursts_to_process)} in dat file {dat.Filename}")     
        
        list_of_singleBurst_xarrays = []

        for burst_number in bursts_to_process:#tqdm(bursts_to_process):
            self.logger.debug(f"Extract burst number {burst_number}")
            burst = dat.ExtractBurst(burst_number)

            singleBurst_xarray = self._burst_to_xarray(burst)

            list_of_singleBurst_xarrays.append(singleBurst_xarray)
        self.logger.debug(f"Concatenating all the single-burst xarrays from dat file {dat.Filename}")
        return xr.concat(list_of_singleBurst_xarrays,dim='time') 
    





    def _burst_to_xarray(self,burst):
        """Return an xarray containing all data from one burst with appropriate coordinates"""
        self.logger.debug(f"Put all chirps and profiles from burst number {burst.BurstNo} in 3D arrays")
        chirps , profiles = self._burst_to_3d_arrays(burst)
        chirp_time, profile_range = self._coords_from_burst(burst)
        time = self._timestamp_from_burst(burst)

        chirps = chirps[None,:,:,:]
        profiles = profiles[None,:,:,:]

        xarray_out = xr.Dataset(
            data_vars=dict(
                chirp           = (["time","chirp_time", "chirp_num", "attenuator_setting_pair"], chirps),
                profile         = (["time", "profile_range", "chirp_num", "attenuator_setting_pair"], profiles),
                latitude        = (["time"], [burst.Header['Latitude']]),
                longitude       = (["time"], [burst.Header['Longitude']]),  
                battery_voltage = (["time"], [burst.Header['BatteryVoltage']]), 
                temperature_1   = (["time"], [burst.Header['Temp1']]),
                temperature_2   = (["time"], [burst.Header['Temp2']])
            ),
            coords=dict(
                time                  = [time],
                chirp_time            = chirp_time,
                profile_range         = profile_range, 
                chirp_num             = np.arange(burst.Header['NSubBursts']),
                filename              = (["time"], [burst.Filename]),
                burst_number          = (["time"], [burst.BurstNo]),
                AFGain                = (["attenuator_setting_pair"], burst.Header['AFGain'][0:burst.Header['nAttenuators']]),
                attenuator            = (["attenuator_setting_pair"], burst.Header['Attenuator1'][0:burst.Header['nAttenuators']])
            ),
        )
        return xarray_out
    
    def _burst_to_3d_arrays(self, burst):
        """Put all chirps and their corresponding profiles in 3D numpy arrays.
        
        First extract every chirp in a burst and compute a profile from each one. 
        
        Then put the results in two 3D numpy arrays. The 3D array for the chirp data (chirp_3d) has the following dimension lengths:
            - [1] burst.Header['N_ADC_SAMPLES'] (=40001, by default) --> number of samples in each chirp  
            - [2] burst.Header['NSubBursts'] --> the number of chirps per attenuator settingp pair 
            - [3] burst.Header['nAttenuators'] --> the number of attenuator settings. 
            
        The 3D array for the profile data (profile_3d) has the same lengths in the dimensions [2] and [3], but a differnt
        length of dimension [1], equal to the length of the profile obtained from the fft processing.
        
        Keyword arguments:
        burst -- a `BurstObject` produced by DataFileObject()
        max_range -- the range to use to crop profile_3d, float or int
        
        Returns: 
        chirp_3d -- 3D numpy array containing all the chirps in the supplied burst
        cropped_profile_3d -- 3D numpy array containing all the profiles from all the chirps in this burs        
        """
        
        self.logger.debug(f"Set max range from _burst_to_3d_arrays")
        self._set_max_range(burst)
        
        # pre-allocate 3D numpy arrays
        chirp_3d = np.zeros((burst.Header['N_ADC_SAMPLES'],burst.Header['NSubBursts'],burst.Header['nAttenuators']))
        test_profile = burst.ExtractChirp([0]).FormProfile()
        profile_3d = np.zeros((len(test_profile.Profile),burst.Header['NSubBursts'],burst.Header['nAttenuators']), dtype=complex)
    
        # initialize counters for second and third dimensions of our 3D arrays respectively 
        chirp_counter = 0
        setting_counter = 0
        
        # loop over all chirps in this burst
        for i in np.arange(burst.Header['NChirps']):
            chirp_new = burst.ExtractChirp([i]) 
            profile_new = chirp_new.FormProfile()

            chirp_3d[:,chirp_counter,setting_counter] = chirp_new.vdat
            profile_3d[:,chirp_counter,setting_counter] = profile_new.Profile  

            # keep track of which pair of settings we are using
            setting_counter += 1
            if setting_counter >= burst.Header['nAttenuators']: # if the counter reaches nAttenuators, reset it to zero 
                setting_counter = 0
                chirp_counter += 1

        # crop the profile to a reasonable depth to save space 
        n = np.argmin(profile_new.Range<=self.max_range)
        cropped_profile_3d = profile_3d[:n,:,:]    

        return chirp_3d, cropped_profile_3d

    def _coords_from_burst(self, burst):
        """Return the time vector and depth vector from the first chirp in a burst. They should be the same for the whole burst."""
        
        self.logger.debug(f"Set max range from _coords_from_burst")
        self._set_max_range(burst)
        
        chirp = burst.ExtractChirp([0])

        profile = chirp.FormProfile()
        n = np.argmin(profile.Range<=self.max_range)
        cropped_range = profile.Range[:n] 
        return  chirp.t, cropped_range

    def _timestamp_from_burst(self,burst):
        """Return the time stamp of a burst"""  
        return pd.to_datetime(burst.Header["Time stamp"])  
    
    def _set_max_range(self,burst):
        
        # if not supplied with a max_range, use what is defined in the header
        if self.max_range is None:
            self.max_range = burst.Header['maxDepthToGraph']
        else: 
            self.max_range = self.max_range 
        self.logger.debug(f"Max_range has been set to {self.max_range}")
 
    def _add_attrs(self):
        """Add attributes to the xarray self.data"""
        self.logger.debug("Adding attributes to the xapres.data")
        self.data.time.attrs['long_name'] = 'time of burst'
        self.data.chirp_time.attrs['long_name'] = 'time of samples during chirps'
        self.data.chirp_time.attrs['name'] = 'time of samples during chirps'
        self.data.chirp_time.attrs['units'] = 'seconds'

        self.data.profile_range.attrs['long_name'] = 'depth'
        self.data.profile_range.attrs['units'] = 'meters'

        self.data.chirp_num.attrs['long_name'] = 'chirp number'
        self.data.chirp_num.attrs['description'] = 'the number of each chirp within each burst'


        self.data.AFGain.attrs['long_name'] = 'audio-frequency gain control setting'
        self.data.AFGain.attrs['units'] = 'decibels'

        self.data.attenuator.attrs['long_name'] = 'radio-frequency attenuator setting'
        self.data.attenuator.attrs['units'] = 'decibels'

        self.data.chirp.attrs['long_name'] = 'de-ramped chirp'
        self.data.chirp.attrs['units'] = 'volts'
        self.data.chirp.attrs['description'] = 'voltage from the analog-to-digital converter after the received signal has been mixed with the transmitted signal and the result has been filtered to leave only the low frequency compponent corresponding to the differences in the frequencies of the Tx and Rx signals'

        self.data.profile.attrs['long_name'] = 'profile' 
        self.data.profile.attrs['units'] = '-'
        self.data.profile.attrs['description'] = 'complex profile computed from the fourier transform of the de-ramped chirp'

        self.data.latitude.attrs['units'] = 'degrees'
        self.data.latitude.attrs['long_name'] = 'latitude of burst'
        
        
        self.data.longitude.attrs['units'] = 'degrees'
        self.data.longitude.attrs['long_name'] = 'longitude of burst'
        
        self.data.battery_voltage.attrs['units'] = 'volts'
        self.data.battery_voltage.attrs['long_name'] = 'battery voltage'
        
        self.data.temperature_1.attrs['units'] = 'celsius'
        self.data.temperature_1.attrs['long_name'] = 'temperature measured inside the ApRES unit in one location'
        self.data.temperature_2.attrs['units'] = 'celsius'
        self.data.temperature_2.attrs['long_name'] = 'temperature measured inside the ApRES unit in a second location'
                
        self.data.filename.attrs['description'] = 'the name of the file that contains each burst'
        self.data.burst_number.attrs['description'] = 'the number of each burst within each file'
        
        #self.data.attrs['time created'] = pd.Timestamp.now()  

    def dB(self, da):
        '''Returns decibels from the DataArray, da, which needs be ApRES complex profile (or collection of them.'''
        
        decibels = 20*np.log10(np.abs(da))
        
        
        #try:
        #    decibels = 20*np.log10(np.abs(self.data.profile))
        #except AttributeError:
        #   raise AttributeError("The xarray xapres.data does not yet exist, run xapres.load_all() to create it.")
        
        return decibels
    
    def _setup_logging(self,loglevel):
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {loglevel}")


        self.logger = logging.getLogger("default_logger")
        # Clear old logging handlers to avoid a build up of multiple similar handlers with repeat calls
        self.logger.handlers.clear()
        
        # Set stream logging level to loglevel
        self.logger.setLevel(level=numeric_level)
        

        logStreamFormatter = logging.Formatter(
            fmt=f"%(levelname)-8s %(asctime)s \t %(filename)s @function %(funcName)s line %(lineno)s - %(message)s", 
            datefmt="%H:%M:%S"
        )
        consoleHandler = logging.StreamHandler(stream=sys.stdout)
        consoleHandler.setFormatter(logStreamFormatter)
        consoleHandler.setLevel(level=numeric_level)
               
        
        self.logger.addHandler(consoleHandler)
        
        self.logger.debug(f"Stream logging level set to {loglevel.upper()}")
        self.logger.debug('Add console handler to logger')
     
           
        logFileFormatter = logging.Formatter(
            fmt=f"%(levelname)s %(asctime)s (%(relativeCreated)d) \t %(pathname)s F%(funcName)s L%(lineno)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        fileHandler = logging.FileHandler(filename='test6.log')
        fileHandler.setFormatter(logFileFormatter)
        fileHandler.setLevel(level=numeric_level)

        self.logger.addHandler(fileHandler)
        
        self.logger.debug(f"File logging level set to {loglevel.upper()}")
   
        
    def _try_logging(self):
        """Simple test of the diffrent logging levels. Not for use by users"""
        self.logger.debug("debugging something")
        self.logger.info("some message")
        self.logger.error("something went wrong")
    
    def coherence(self, s1, s2):
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

    def generate_range_diff(self, data1, data2, win_cor, step, range_ext=None, win_wrap=10, thresh=0.9, uncertainty='CR'):
        """
        Input data:
        data1, data2: xarray.DataArrays describing the "profile" variable
        range
        """

        # Get times of burst:
        t1 = data1.time.data
        t2 = data2.time.data
        dt = (t2-t1)/ np.timedelta64(1, 's')
        self.logger.info(f"Time between bursts : {dt}s")
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
    


class DataFileObject:
    
    def __init__(self, Filename, remote_load=False):
        self.BurstLocationList = []
        self.Filename = Filename
        self.remote_load = remote_load
        
        # try to looad locally, then resort to loading remotely (would remove need for "remote" arg if propagated to rest of classes)
        try:
            datafile = open(self.Filename, "rb")
        except TypeError:
            print("no local file found, trying remote load..")
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode='rb')
            if datafile is not None:
                print("remote load successful")
            else:
                print("remote load failed")
        except:
            print("file could not be found, see error message")
            

        inbuff = datafile.read()
        datafile.close()     
        
        a = "*** Burst Header ***"
        b = a.encode()
        locn = inbuff.find(b)
        while locn != -1:
            self.BurstLocationList.append(locn)
            locn = inbuff.find(b, locn + len(b))
       
        self.NoBurstsInFile = len(self.BurstLocationList)

    def ExtractBurst(self, BurstNo):
        Burst = BurstObject()
        Burst.Filename = self.Filename
        Burst.BurstNo = BurstNo
        Burst.Header = {"BurstNo":BurstNo}
        
        if self.remote_load:
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode = 'rb')
        else: 
            datafile = open(self.Filename, "rb")
 
        
        datafile.seek(self.BurstLocationList[BurstNo])
        inbuff = datafile.read(2000)
        locn = 0
        locn1 = inbuff.find(b'\x0D\x0A', locn, locn + 80)
        inline = inbuff[locn:locn1].decode()
        while inline.count("End Header") == 0:
            tmp = inline.split("=")
            if len(tmp) == 2:
                if tmp[0].lower() == "rxant" or \
                   tmp[0].lower() == "txant" or \
                   tmp[0].lower() == "afgain":
                   Burst.Header[tmp[0]] = \
                        [int(x) for x in tmp[1].split(',') if x]
                        
                elif tmp[0].lower() == "triples" or \
                     tmp[0].lower() == "attenuator1" or \
                     tmp[0].lower() == "batterycheck":
                     Burst.Header[tmp[0]] = \
                        [float(x) for x in tmp[1].split(',') if x]
                                      
                elif tmp[0].lower() == "latitude" or \
                     tmp[0].lower() == "longitude" or \
                     tmp[0].lower() == "temp1" or \
                     tmp[0].lower() == "temp2" or \
                     tmp[0].lower() == "batteryvoltage" or \
                     tmp[0].lower() == "tstepup" or \
                     tmp[0].lower() == "tstepdn" or \
                     tmp[0].lower() == "fsc" or \
                     tmp[0].lower() == "sw_issue" or \
                     tmp[0].lower() == "er_ice" or \
                     tmp[0].lower() == "position_depth_conversion" or \
                     tmp[0].lower() == "maxdepthtograph":
                    Burst.Header[tmp[0]] = float(tmp[1])
                    
                elif tmp[0].lower() == "rmb_issue" or \
                     tmp[0].lower() == "vab_issue" or \
                     tmp[0].lower() == "reg00" or \
                     tmp[0].lower() == "reg01" or \
                     tmp[0].lower() == "reg02" or \
                     tmp[0].lower() == "reg03" or \
                     tmp[0].lower() == "reg0b" or \
                     tmp[0].lower() == "reg0c" or \
                     tmp[0].lower() == "reg0d" or \
                     tmp[0].lower() == "reg0e":
                    Burst.Header[tmp[0]] = tmp[1]
                elif tmp[0].lower() == "time stamp":
                    Burst.Header[tmp[0]] = tmp[1]
                else:
                    Burst.Header[tmp[0]] = int(tmp[1])
                    
            locn = locn1 + 2
            locn1 = inbuff.find(b'\x0D\x0A', locn,locn + 80)
            inline = inbuff[locn:locn1].decode()

        # Re-open the file for binary read to get burst data
        
        if self.remote_load:
            fs = gcsfs.GCSFileSystem()
            datafile = fs.open(self.Filename, mode = 'rb')
        else: 
            datafile = open(self.Filename, "rb")      
        
        datafile.seek(self.BurstLocationList[BurstNo] + locn1 + 2)
        NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"] * Burst.Header["NSubBursts"]
        if Burst.Header["Average"] == 1:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"]
            inbuff = datafile.read(NsamplesInBurst * 4)
            Burst.v = np.frombuffer(inbuff, dtype=np.float32)/2**16*2.5-1.25
        elif Burst.Header["Average"] == 2:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"]
            inbuff = datafile.read(NsamplesInBurst * 4)
            Burst.v = np.frombuffer(inbuff, dtype=np.uint32)
        else:
            NsamplesInBurst = Burst.Header["N_ADC_SAMPLES"] * Burst.Header["NSubBursts"] *  Burst.Header["nAttenuators"]  # because each subburst conatins a certain number of chirps given by Burst.Header["nAttenuators"]
            inbuff = datafile.read(NsamplesInBurst * 2)
            Burst.v = (np.frombuffer(inbuff, dtype=np.uint16))/2**16*2.5-1.25
        datafile.close()
        
        if "SamplingFreqMode" in Burst.Header:
            if Burst.Header["SamplingFreqMode"] == 0:
                Burst.Header["dt"] = 1/40000
            else:
                Burst.Header["dt"] = 1/80000
        else:
            Burst.Header["dt"] = 1/40000
            
        Burst.Header["NChirps"] = Burst.Header["NSubBursts"] * Burst.Header["nAttenuators"] *\
            Burst.Header["TxAnt"].count(1) * Burst.Header["RxAnt"].count(1)
        
        if "FreqStepUp" in Burst.Header:
            Burst.Header["K"] = Burst.Header["FreqStepUp"] / Burst.Header["TStepUp"]
        else:
            Burst.Header["K"] = 2e8
            Burst.Header["StartFreq"] = 2e8
            Burst.Header["StopFreq"] = 4e8
            
        Burst.Header["c0"] = 3e8
        if  not ("ER_ICE" in Burst.Header):
            Burst.Header["ER_ICE"] = 3.18
        Burst.Header["CentreFreq"] = (Burst.Header["StartFreq"] + Burst.Header["StopFreq"])/2
        Burst.Header["B"] = (Burst.Header["StopFreq"] - Burst.Header["StartFreq"])
        
        # deal out attenuator settings to the chirps
        setting_counter = 0
        Burst.Header['Attenuator1_allChirps'] = []
        Burst.Header['AFGain_allChirps'] =[]
        for chirp in range(Burst.Header['NChirps']):
            Burst.Header['Attenuator1_allChirps'].append(Burst.Header['Attenuator1'][setting_counter])
            Burst.Header['AFGain_allChirps'].append(Burst.Header['AFGain'][setting_counter])
            
            # keep track of which setting to use
            setting_counter += 1
            if setting_counter >= Burst.Header['nAttenuators']: # if the counter reaches nAttenuators, reset it to zero 
                setting_counter = 0


        return(Burst)
        
class BurstObject:
    
    def __init__(self):
        self.v = 0
        self.Filename = ''
        self.Header = 0
        self.BurstNo = 0
        
    def ExtractChirp(self, ChirpList):
        Chirp = ChirpObject()
        Chirp.t = np.array(range(self.Header["N_ADC_SAMPLES"])) * self.Header["dt"]
        Chirp.Header = copy.deepcopy(self.Header)  # we dont want chirp-specific new entries in the chirp.header updating the burst.header
        Chirp.ChirpList = ChirpList
        
        Chirp.Header['Attenuator1_thisChirp'] = [self.Header['Attenuator1_allChirps'][i] for i in ChirpList]
        Chirp.Header['AFGain_thisChirp'] = [self.Header['AFGain_allChirps'][i] for i in ChirpList]

        if any(i != Chirp.Header['Attenuator1_thisChirp'][0] for i in Chirp.Header['Attenuator1_thisChirp'])\
            or any(i != Chirp.Header['AFGain_thisChirp'][0] for i in Chirp.Header['AFGain_thisChirp']):
            warnings.warn('This is stacking over chirps with different attenuator settings.')


        if self.Header["Average"] == 1:
            Chirp.vdat = self.v          
        elif self.Header["Average"] == 2:
            Chirp.vdat = self.v          
        else:
            Chirp.vdat = np.zeros(self.Header["N_ADC_SAMPLES"])
            no = 0
            for ind in ChirpList:
                if ind < self.Header["NChirps"]:
                    no += 1
                    chirpoffset = ind * self.Header["N_ADC_SAMPLES"]
                    Chirp.vdat = Chirp.vdat + self.v[chirpoffset:chirpoffset + self.Header["N_ADC_SAMPLES"]]
                else:
                    print('chirp index > number of chirps.')
            Chirp.vdat = Chirp.vdat/no

        return(Chirp)
        
    def PlotBurst(self):
        t = np.array(range(len(self.v))) * self.Header["dt"]
        plt.plot(t, self.v)
        plt.axis([0,t[-1],-1.25,1.25])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (V)")
        plt.grid()
        return(0)
    
class ChirpObject:
    def __init__(self):
        self.vdat = 0
        self.t = 0
        self.ChirpList = 0
        self.Filename = ''
        self.BurstNo = 0
        self.Header = 0
        
    def PlotChirp(self):
        plt.plot(self.t, self.vdat)
        plt.axis([0,self.t[-1],-1.25,1.25])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (V)")
        plt.grid("on")
        return(0)

    def FormProfile(self, F0=200000000, F1=400000000, pad=2, ref=1):
        Profile = ProfileObject()
        if self.Header["StartFreq"] > F0:
            F0 = self.Header["StartFreq"]
        if self.Header["StopFreq"] < F1:
            F1 = self.Header["StopFreq"]
        T0 = (F0-self.Header["StartFreq"])/self.Header["K"]
        T1 = (F1-self.Header["StartFreq"])/self.Header["K"]
        chirp = self.vdat[math.ceil(T0/self.Header["dt"]):\
                          math.floor(T1/self.Header["dt"])]

        Nt = len(chirp)
        Nt = math.floor(Nt/2) * 2
        winchirp = np.multiply(chirp[0:Nt],np.blackman(Nt))
        Nfft = math.floor(Nt*pad)

        padchirp = np.zeros(Nfft)
        padchirp[0:math.floor(Nt/2)-1] = winchirp[math.floor(Nt/2):-1]
        padchirp[-math.floor(Nt/2):-1] = winchirp[0:math.floor(Nt/2)-1]
        Profile.Profile = np.fft.fft(padchirp)/Nfft * math.sqrt(2*pad) 
        Profile.bin2m = self.Header["c0"]/(2.*(T1-T0)*pad*math.sqrt(self.Header["ER_ICE"])*self.Header["K"])
        Profile.Range = np.asarray([i for i in range(Nfft)]) * Profile.bin2m       
        Profile.Profile = Profile.Profile[0:math.floor(Nfft/2)-1]
        Profile.Range = Profile.Range[0:math.floor(Nfft/2)-1]
        if ref == 1:
            m = np.asarray([i for i in range(len(Profile.Profile))])/pad
            phiref = 2*math.pi*self.Header["CentreFreq"]*m/self.Header["B"] -\
             m * m * 2*math.pi * self.Header["K"]/2/self.Header["B"]**2
            Profile.Profile = Profile.Profile * np.exp(phiref*(-1j))
        
        Profile.BurstNo = self.BurstNo
        Profile.Header = self.Header
        Profile.Filename = self.Filename
        Profile.ChirpList = self.ChirpList
        Profile.pad = pad
        Profile.rad2m = self.Header["CentreFreq"]*math.sqrt(self.Header["c0"])/ \
            (4.*math.pi*self.Header["ER_ICE"])
        return(Profile)  
       
class ProfileObject:
    
    def __init__(self):
        self.Range = 0
        self.Profile = 0
        self.F0 = 2e8
        self.F1 = 4e8
        self.pad = 2
        self.ChirpList = 0
        self.Filename = ''
        self.BurstNo = 0
        self.Header = 0
        self.rad2m = 0
        self.bin2m = 0

    def PlotProfile(self,dmax):
        plt.plot(self.Range, 20*np.log10(np.abs(self.Profile)))
        plt.xlim(0, dmax)
        plt.xlabel("Range (m)")
        plt.ylabel("Amplitude (dB)")
        plt.grid("on")

        return(0)
        
    
