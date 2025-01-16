"""
For loading ApRES data from .dat files or zarr directories into an xarray dataset.
"""

import os
import glob
import sys
import logging

import gcsfs
import numpy as np
import xarray as xr
import pandas as pd
from tqdm import tqdm
import datetime

#print(f"once we run load we are in  {os.getcwd()}")


#sys.path.append("../bas-apres")
import apres as ap

def load_zarr(directory = "gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/A101"):
    """Load ApRES data stored in a zarr directory. """
    
    return xr.open_dataset(directory,
            engine = 'zarr', 
            chunks = {}) 

def generate_xarray(directory=None, 
                 file_numbers_to_process=None, 
                 file_names_to_process=None, 
                 bursts_to_process="All",
                 attended=False, 
                 polarmetric=False,
                 max_range = None,
                 computeProfiles = True,
                 addProfileToDs_kwargs = {},
                 loglevel = 'warning'
                 ):
    """ Load data from multiple .dat files into an xarray dataset.

    This is simple wrapper for from_dats.load_all. This slightly simplifies the process of loading ApRES data into an xarray because it avoids having to initialize the from_dats object.
    """

    fd = from_dats(loglevel=loglevel)
    
    fd.load_all(directory=directory, 
                 file_numbers_to_process=file_numbers_to_process, 
                 file_names_to_process=file_names_to_process, 
                 bursts_to_process=bursts_to_process,
                 attended=attended, 
                 disable_progress_bar = True, 
                 polarmetric=polarmetric,
                 max_range = max_range,
                 computeProfiles = computeProfiles,
                 addProfileToDs_kwargs = addProfileToDs_kwargs,
                )

    return fd.data

class from_dats():
    """
    An object containing ApRES data loaded from a dat file or many dat files, along with information about the data. 
    
        Can be instantiated with 2 optional keyword arguments, loglevel and max_range
    
        Argument:
            loglevel --- allows the user to select the level of logging messages are displayed. 
            The default loglevel is warning, which means that no messages are displayed. 
            If you want to see detailed log messages, use loglevel = 'debug'
        
        Methods:
            load --- load a single dat file into an xarray
            list_files --- recursively find  all the files in a directory or a google bucket
            load_all --- load all the files found in a directory or google bucket into an xarray
    
        load_all is the most important method. Call it, for example, as follows:
    
            import ApRESDefs
            fd_unattended = ApRESDefs.load.from_dats(loglevel='debug')
            fd_unattended.load_all(directory='gs://ldeo-glaciology/GL_apres_2022', 
                        file_numbers_to_process = [0, 1], 
                        bursts_to_process=[0, 1]
           )
    
        the resulting xarray will be saved in xa.data.
    
    """
    def __init__(self, loglevel='warning'):
        self._setup_logging(loglevel)
        
    def load(self,
            dat_filename,
            bursts_to_process="All",
            attended=False,
            polarmetric=False,
            max_range = None,
            computeProfiles = True,
            addProfileToDs_kwargs = {}):
        
        self.max_range = max_range
        self.attended = attended
        self.polarmetric = polarmetric

        if attended is True:
            self.data = self.all_bursts_at_waypoint_to_xarray(dat_filename=dat_filename)
        else:    
            self.data = self.all_bursts_in_dat_to_xarray(dat_filename, bursts_to_process) 
        
        self.correct_temperature()

        self.correct_chirp_dtype()

        self.add_attrs()

        if computeProfiles is True:
            self.data = self.data.addProfileToDs(**addProfileToDs_kwargs)

        return self.data

    def list_files(self, 
                   directory=None, 
                   search_suffix=""
                   ):    
        """Recursively list all the .DAT files in a given location dir. 
        
        Arguments:
        directory -- the directory that will be looked in recursivly to find .DAT files.
        search_suffix -- a string that can be used to search for files with a specific suffix.
        """

        self.logger.debug(f"Find all the dat files in the directory {directory}")

        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        self.is_this_a_remote_load()

        if self.remote_load:   
            fs = gcsfs.GCSFileSystem()
            dat_filenames_without_gs_prefix = fs.glob(self.directory + '/**/*' + search_suffix + '.[dD][aA][tT]', recursive = True)
            dat_filenames = ['gs://' + x for x in dat_filenames_without_gs_prefix]
        else:
            dat_filenames = glob.glob(self.directory + '/**/*' + search_suffix  +'.[dD][aA][tT]', recursive = True)
        
        self.dat_filenames = dat_filenames
        
        self.logger.debug(f"Finish call to list_files. Found {len(dat_filenames)} files")

        return dat_filenames
          
    def load_all(self,
                 directory=None, 
                 file_numbers_to_process=None, 
                 file_names_to_process=None, 
                 bursts_to_process="All",
                 disable_progress_bar = True, 
                 attended=False, 
                 polarmetric=False,
                 max_range = None,
                 computeProfiles = True,
                 addProfileToDs_kwargs = {}
                 ):
        
        """
        Method to recursively load ApRES .dat files and put them in an xarray dataset. It also 
        computes profiles from the chirp data and includes them in the output. 
        
        This method has two modes. One for unattended ApRES data and one for attended data. 
        
        For unattended data (i.e. attended=False, the default), it puts all the data from all 
        the .DAT files found recursively in 'directory', into one xarray. The most important 
        dimension of this xarray is 'time', which is the time of each burst. 

        In attended mode, the method locates the dat files corresponding to each waypoint. 
        It does this based on a user-supplied list of directories 'directory'. The method groups the 
        data by waypoint (and optionally antenna orientation).

        Parameters
        ----------
        directory : str or list, optional
            Directory or list of directories containing .DAT files. 
            If attended is False, this should be a single directory which will be search recusrivley for dat fies.  
            If attended is True, this should be a list of directories, one for each waypoint. Default is None.
        file_numbers_to_process : list, optional
            List of file numbers to process. If None, all files will be processed. Default is None.
        file_names_to_process : list, optional
            List of file names to process. If None, all files will be processed. Default is None. 
            Note that you can set either file_numbers_to_process or file_names_to_process, but not both. 
        bursts_to_process : str or list, optional
            Bursts to process from within each dat file. Default is "All".
        attended : bool, optional
            If True, load data in attended mode. Default is False.
        polarmetric : bool, optional
            If True, load data in polarmetric mode - the xarray dataset outputted will have an antenna-orientation dimension corrosponding to HH, HV, VH, and VV. 
            It designates dat files to each orientation based on the files names containing HH, HV, VH, or VV.
            Default is False.
         max_range : float, optional
            A range value used to crop the profiles. Only used in the legacy fft processing. Default is None.
        computeProfiles : bool, optional
            If True, compute profiles from the chirp data. Default is True.
        addProfileToDs_kwargs : dict, optional
            Additional keyword arguments for addProfileToDs method. 
            The following can be set:
                pad_factor = 2
                drop_noisy_chirps = False,
                clip_threshold = 1.2,
                min_chirps = 20,
                demean = False,
                detrend = False,
                stack = False,
                crop_chirp_start = None,
                crop_chirp_end = None,
                max_range = None

        Returns
        -------
        xarray.Dataset
            The loaded data as an xarray dataset.

        Raises
        ------
        ValueError
            If attended mode is True and directory_list is None.
            If both file_numbers_to_process and file_names_to_process are supplied.
        """   

        self.file_numbers_to_process = file_numbers_to_process
        self.file_names_to_process = file_names_to_process
        self.attended = attended
        self.burst_load_counter = 0   # this will increment each time a burst is loaded by _burst_to_xarray
        self.polarmetric = polarmetric
        self.max_range = max_range
        self.computeProfiles = computeProfiles
        #self.bursts_to_process = bursts_to_process
        
        if directory is None:
            self.directory = os.getcwd()
        else:
            self.directory = directory

        if isinstance(bursts_to_process, int):
            bursts_to_process = [bursts_to_process]


        self.is_this_a_remote_load()
        
        self.logger.debug(f"Start call to load_all with remote_load = {self.remote_load}, directory = {directory}, file_numbers_to_process = {file_numbers_to_process}, file_names_to_process = {file_names_to_process}, bursts_to_process = {bursts_to_process}, attended = {attended}")
      
        if attended is True and isinstance(directory, str):
            directory_list = [directory]
        elif attended is True and isinstance(directory, list):
            directory_list = directory

        if attended is False:
            self.list_files(directory)    # adds self.dat_filenames
            if not self.dat_filenames:
                self.logger.debug("No dat files found. Returning None.")
                return None
            self.subset_files()   # adds self.dat_filenames_to_process
        
            # Loop through the dat files, putting individual xarrays in a list.
            self.logger.debug("Attended is False, so starting loop over dat files")
                        
            list_of_multiBurstxarrays = [
                self.all_bursts_in_dat_to_xarray(dat_filename, bursts_to_process) 
                for dat_filename in tqdm(self.dat_filenames_to_process, disable=disable_progress_bar)
            ]

            self.logger.debug(f"Attended is False, so concatenating all the multi-burst xarrays along the time dimension, to create xapres.data")
            self.data = xr.concat(list_of_multiBurstxarrays, dim='time') 
        
        elif attended is True:
                    
            # loop over the waypoints, as defined by the directories in the user-supplied list, directory_names
            self.logger.debug("Attended is True, so starting loop over directories (each corresponding to a waypoint)")
            list_of_singlewaypoint_xarrays = [
                self.all_bursts_at_waypoint_to_xarray(directory=directory, waypoint_number=waypoint_number)
                for waypoint_number, directory in enumerate(directory_list, start=1)
            ]   
            
            self.logger.debug(f"Attended is True, so concatenating all the single-waypoint xarrays along the waypoint dimension, to create xapres.data")
            
            lengths_of_ant_dims = [len(x.attenuator_setting_pair) for x in list_of_singlewaypoint_xarrays]
            min_length = min(lengths_of_ant_dims)   
            cropped_datasets = [x.isel(attenuator_setting_pair=(min_length-1)) for x in list_of_singlewaypoint_xarrays]

            try: 
                self.data = xr.concat(cropped_datasets, dim='waypoint')
            except ValueError as e:
                return list_of_singlewaypoint_xarrays
            
        self.correct_temperature()

        self.correct_chirp_dtype()

        self.add_attrs()

        if  self.computeProfiles is True:
            self.logger.debug(f"Call addProfileToDs to add the profiles to the xarray")
            self.data = self.data.addProfileToDs(**addProfileToDs_kwargs)

        self.logger.debug(f"Finish call to load_all. Call xapres.data to see the xarray this produced.")

        return self.data

    def subset_files(self):
        """Subset files based on either file_numbers_to_process or file_names_to_process. Throws an error if both are supplied. This only gets used for unattended data."""
        if self.file_numbers_to_process is not None and self.file_names_to_process is not None:
            self.logger.debug("Throwing a ValueError because file_numbers_to_process and file_names_to_process cannot both be supplied to load_all")
            raise ValueError("file_numbers_to_process and file_names_to_process cannot both be supplied to load_all. You need to supply just one (or neither) of these.") 
        
        elif self.file_numbers_to_process is not None:
            if self.file_numbers_to_process == "All":
                self.logger.debug("Selecting all dats file because file_numbers_to_process == \"all\"")
                self.dat_filenames_to_process = self.dat_filenames
            else:
                self.logger.debug(f"Subset files to {self.file_numbers_to_process}")
                self.dat_filenames_to_process = [self.dat_filenames[i] for i in self.file_numbers_to_process]
        
        elif self.file_names_to_process is not None:
            if self.file_names_to_process == "All":
                self.logger.debug("Selecting all dats file because file_names_to_process == \"all\"")
                self.dat_filenames_to_process = self.dat_filenames
            else:                 
                self.logger.debug("Subset files to list of files supplied in file_names_to_process")
                self.dat_filenames_to_process = self.file_names_to_process
                            
        elif self.file_numbers_to_process is None and self.file_names_to_process is None:      # default is all the dat files    
            self.logger.debug("Selecting all dats file because neither file_numbers_to_process nor file_names_to_process were supplied")
            self.dat_filenames_to_process = self.dat_filenames          
        
    def all_bursts_in_dat_to_xarray(self, 
                                     dat_filename, 
                                     bursts_selected,
                                     ):
        """Take data from all the bursts in one .DAT file and put it in an xarray, concatenated by time. 
        Only gets used for unattended data because attended data should only have one burst per dat file.
        
        Arguments:
        dat_filename -- the filename of the dat file
        bursts_selected -- a list of the burst numbers to process, or "All" to process all the bursts in the dat file.
        """
        self.bursts_selected = bursts_selected
        
        self.logger.debug("Attended is False. Generating xarray for unattended data")   
        self.logger.debug(f"bursts_to_process = {bursts_selected} at the start of _all_bursts_in_dat_to_xarray.")

        with ap.ApRESFile(dat_filename) as self.f: 
            self.f.read()
        
        self.header_cleaning()
        
        self.subset_bursts_to_process()

        filename = os.path.basename(self.f.path)
        folder_name = os.path.basename(os.path.dirname(self.f.path))
        bursts = [self.f.bursts[burst_number] for burst_number in self.bursts_to_process]  

        list_of_singleBurst_xarrays = [
            xr.Dataset(
                        data_vars=dict(
                            chirp           = (["time", "chirp_num", "chirp_time", "attenuator_setting_pair"], self.burst_data(burst)/2**16*2.5-1.25),
                            latitude        = (["time"], [burst.header['Latitude']]),
                            longitude       = (["time"], [burst.header['Longitude']]),  
                            battery_voltage = (["time"], [burst.header['BatteryVoltage']]), 
                            temperature_1   = (["time"], [burst.header['Temp1']]),
                            temperature_2   = (["time"], [burst.header['Temp2']])
                        ),
                        coords=dict(    
                            time                  = [self._timestamp_from_burst(burst)],
                            chirp_time            = self.chirptime_from_burst(burst),
                            chirp_num             = np.arange(burst.header['NSubBursts']),
                            filename              = (["time"], [filename]),
                            folder_name           = (["time"], [folder_name]),                          
                            burst_number          = (["time"], [burstNo]),
                            AFGain                = (["attenuator_setting_pair"], burst.header['AFGain'][0:int(burst.header['nAttenuators'])]),
                            attenuator            = (["attenuator_setting_pair"], burst.header['Attenuator1'][0:int(burst.header['nAttenuators'])]),                   
                            orientation           = (["time"], [self._get_orientation(filename)]),
                        ),
                    )
                    for burstNo, burst in enumerate(bursts)
        ]
        
        self.current_burst = bursts[-1]
        
        return xr.concat(list_of_singleBurst_xarrays, dim='time') 

    def all_bursts_at_waypoint_to_xarray(self, 
                                          directory=None,
                                          waypoint_number=1,
                                          dat_filename=None):   
        """This is the attended equivalent to _all_bursts_in_dat_to_xarray"""
        
        if self.polarmetric is True:
            orientations = ['HH', 'HV', 'VH', 'VV']
        else:
            orientations = [""]

        list_of_singleorientation_attended_xarrays = []

        for orientation in orientations:
            self.logger.debug(f"Looking for files with orientation {orientation} in directory {directory}")
            
            if dat_filename is not None:   # this is the caes when this function is supplied with a dat_filename instead of a directory (e.g., when called by from_dats.load())
                files = [dat_filename]
            else:
                files = self.list_files(directory=directory, search_suffix=orientation)

            if len(files) != 1:
                raise Exception(f'There should by one dat file for each orientation in each directory when loading in attended mode. We found {len(files)} files.')
            
            with ap.ApRESFile(files[0]) as self.f: 
                self.f.read()

            self.header_cleaning()
            
            if len(self.f.bursts) != 1:
                raise Exception(f'There should only be one burst in each dat file when loading in attended mode. We found {len(self.f.bursts)} bursts in this dat file: {files}.')
        
            burst = self.f.bursts[0]
            filename = os.path.basename(self.f.path)
            folder_name = os.path.basename(os.path.dirname(self.f.path))

            singleorientation_attended_xarray = xr.Dataset(
                data_vars=dict(
                    chirp           = (["orientation", "waypoint", "chirp_num", "chirp_time", "attenuator_setting_pair"], self.burst_data(burst)/2**16*2.5-1.25),
                    latitude        = (["orientation", "waypoint"], np.array(burst.header['Latitude'], ndmin = 2)),
                    longitude       = (["orientation", "waypoint"], np.array(burst.header['Latitude'], ndmin = 2)),
                    battery_voltage = (["orientation", "waypoint"], np.array(burst.header['BatteryVoltage'], ndmin = 2)),
                    temperature_1   = (["orientation", "waypoint"], np.array(burst.header['Temp1'], ndmin = 2)),
                    temperature_2   = (["orientation", "waypoint"], np.array(burst.header['Temp2'], ndmin = 2))
                ),
                coords=dict(
                    time                  = (["orientation", "waypoint"], np.array(self._timestamp_from_burst(burst), ndmin = 2)),
                    chirp_time            = self.chirptime_from_burst(burst),
                    chirp_num             = np.arange(burst.header['NSubBursts']),
                    filename              = (["orientation", "waypoint"], np.array(filename, ndmin = 2)), 
                    folder_name           = (["orientation", "waypoint"], np.array(folder_name, ndmin = 2)),
                    AFGain                = (["attenuator_setting_pair"], burst.header['AFGain'][0:burst.header['nAttenuators']]),
                    attenuator            = (["attenuator_setting_pair"], burst.header['Attenuator1'][0:burst.header['nAttenuators']]),                   
                    orientation           = [orientation],
                    waypoint              = [waypoint_number]
                ),
            )

            self.current_burst = burst

            list_of_singleorientation_attended_xarrays.append(singleorientation_attended_xarray)
        
        return xr.concat(list_of_singleorientation_attended_xarrays, dim='orientation')
    
    def header_cleaning(self):
        """Clean the header of the dat file. This is done line by line, with different cleaning procedures for different header values."""
        def line_by_line_cleaning(header):
        
            to_int_list = ["rxant", "txant", "afgain"]
            to_float_list = ["triples", "attenuator1", "batterycheck"]
            to_float = ["latitude", "longitude", "temp1", "temp2", "batteryvoltage", "tstepup", "tstepdn", "fsc", "sw_issue", "er_ice", "position_depth_conversion", "maxdepthtograph"]
            do_nothing = ["time stamp", "rmb_issue", "vab_issue", "reg00", "reg01", "reg02", "reg03", "reg0b", "reg0c", "reg0d", "reg0e"] 

            for key, value in header.items():
                kl = key.lower()
                if kl in to_int_list:
                    header[key] = [int(x) for x in value.split(',') if x]
                    continue
                if kl in to_float_list:
                    header[key] = [float(x) for x in value.split(',') if x]
                    continue
                if kl in to_float:
                    header[key] = float(value)
                    continue
                if kl in do_nothing:
                    continue

                header[key] = int(value)
  
            if "FreqStepUp" in header:
                header["K"] = header["FreqStepUp"] / header["TStepUp"]
            else:
                header["K"] = 2e8
            header["StartFreq"] = 2e8
            header["StopFreq"] = 4e8
                
            header["c0"] = 3e8
            if  not ("ER_ICE" in header):
                header["ER_ICE"] = 3.18
            
            header["CentreFreq"] = (header["StartFreq"] + header["StopFreq"])/2
            header["B"] = (header["StopFreq"] - header["StartFreq"])
                
            return header

        for i, b in enumerate(self.f.bursts):
            self.f.bursts[i].header = line_by_line_cleaning(b.header)

    def subset_bursts_to_process(self):
        """Subset the bursts to process based on the user-supplied list of burst numbers.
        The default is all of the burst in each dat file."""

        if self.bursts_selected == "All":
            self.bursts_to_process = range(len(self.f.bursts))
            self.logger.debug("bursts_to_process set to \"All\"")
        else:
            self.bursts_to_process = self.bursts_selected
        
        if any(np.array(self.bursts_to_process) > (len(self.f.bursts)-1)):
            self.logger.debug(f"The burst numbers requested in bursts_to_process ({self.bursts_to_process}) is greater than the number of bursts in\
                the dat file ({len(self.f.bursts)}), so we will just process all the bursts.")
            self.bursts_to_process = range(len(self.f.bursts))
        self.logger.debug(f"After the initial parse in _all_bursts_in_dat_to_xarray, bursts_to_process = {list(self.bursts_to_process)}.")

    def _timestamp_from_burst(self, burst):
        """Return the time stamp of a burst"""  
        return pd.to_datetime(burst.header["Time stamp"])  

    def chirptime_from_burst(self, burst):
        """Return the time vector for the chirps in a burst"""
        return np.arange(burst.header["N_ADC_SAMPLES"]) * burst.header['TStepUp']

    def _get_orientation(self, filename):
        """Get the orientation of the antenna from the filename"""
        orientation = filename[-6:-4]
        #add function to make sure orientation is capitalized

        #if no orientation at end of filename, mark as unknown
        if orientation not in ['HH','HV','VH','VV']:
            orientation = 'unknown'

        return orientation
     
    def burst_data(self, burst):
        """Reshape the burst data to have the correct dimensions for the xarray.
        If only one attenuator setting is used, add an extra dimension to the end.
        If attended mode is True, add two extra dimensions to the start.
        If attended mode is False add only one extra dimension at the start.
        
        returns the burst data (numpy array)"""
        
        data = burst.data

        if burst.header['nAttenuators'] == 1:
            data = np.expand_dims(data, axis=-1)
        
        if self.attended:
            data = np.expand_dims(data, axis=(0, 1))
        else:
            data = np.expand_dims(data, axis=0)

        return data
    
    def add_attrs(self):
        """Add attributes to the xarray self.data"""
        self.logger.debug("Adding attributes to the xapres.data")
        self.data.time.attrs['long_name'] = 'time of burst'
        self.data.chirp_time.attrs['long_name'] = 'time of samples during chirps'
        self.data.chirp_time.attrs['name'] = 'time of samples during chirps'
        self.data.chirp_time.attrs['units'] = 'seconds'

        self.data.chirp_num.attrs['long_name'] = 'chirp number'
        self.data.chirp_num.attrs['description'] = 'the number of each chirp within each burst'

        self.data.AFGain.attrs['long_name'] = 'audio-frequency gain control setting'
        self.data.AFGain.attrs['units'] = 'decibels'

        self.data.attenuator.attrs['long_name'] = 'radio-frequency attenuator setting'
        self.data.attenuator.attrs['units'] = 'decibels'

        self.data.chirp.attrs['long_name'] = 'de-ramped chirp'
        self.data.chirp.attrs['units'] = 'volts'
        self.data.chirp.attrs['description'] = 'voltage from the analog-to-digital converter after the received signal has been mixed with the transmitted signal and the result has been filtered to leave only the low frequency compponent corresponding to the differences in the frequencies of the Tx and Rx signals'

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
        
        if 'burst_number' in self.data:
            self.data.burst_number.attrs['description'] = 'the number of each burst within each file'
        
        self.data.attrs['constants'] = {'c': self.current_burst.header['c0'],
                                         'K': self.current_burst.header['K'],
                                         'f_1': self.current_burst.header['StartFreq'],
                                         'f_2': self.current_burst.header['StopFreq'],
                                         'dt': self.current_burst.header['TStepUp'], 
                                         'ep': self.current_burst.header['ER_ICE'],
                                         'B': self.current_burst.header['B'],
                                         'f_c': self.current_burst.header['CentreFreq']}

        #self.data.attrs['header from last burst'] = self.current_burst.header

        self.data.orientation.attrs['description'] = 'HH, HV, VH, or VV antenna orientation as described in Ersahadi et al 2022 doi:10.5194/tc-16-1719-2022'
        
        self.data.attrs["processing"] = f"Created on {datetime.datetime.now() }"

        if self.attended:
            self.data.waypoint.attrs['description'] = 'the number of the waypoint where the data was collected'
       
    def correct_temperature(self, threshold=300, correction = -512):
        """Correct temperature values that are above a certain threshold. This appears to result from the temperature data being in the wrong format."""
        self.logger.debug(f"Correct temperature values above {threshold} by adding {correction}")
        T1 = self.data['temperature_1']
        T2 = self.data['temperature_2']

        self.data['temperature_1'] = xr.where(T1>threshold, T1+correction, T1, keep_attrs=True)
        self.data['temperature_2'] = xr.where(T2>threshold, T2+correction, T2, keep_attrs=True) 

    def correct_chirp_dtype(self):
        """ Convert the chirp data type to float64 if it is not already. This is necessary for plotting."""
        if not np.issubdtype(self.data.chirp.chirp_time.dtype, 'float64'):
            self.data.chirp['chirp_time'] = self.data.chirp.chirp_time.values.astype('float64')/1e9

    def _setup_logging(self, loglevel):
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
        fileHandler = logging.FileHandler(filename='generating_xarrays_from_dats.log')
        fileHandler.setFormatter(logFileFormatter)
        fileHandler.setLevel(level=numeric_level)

        self.logger.addHandler(fileHandler)
        
        self.logger.debug(f"File logging level set to {loglevel.upper()}")
        
    def is_this_a_remote_load(self, filename = None):
        """Detect if we are trying to load from Google Cloud Storage based on the filename being supplied"""
        
        name_to_check_for_gs = filename if filename is not None else self.directory

        self.remote_load = "gs://" in name_to_check_for_gs
        self.logger.debug(f"remote_load set to {self.remote_load}")
