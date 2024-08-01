# xapres

A package for processing data from the Autonomous phase-sensitive Radio-Echo Sounder (ApRES) using xarray. The core ApRES processing code is adapted from code written by Keith Nicholls, British Antarctic Survey, UK. This package uses Keith's code to process the raw ApRES data from ApRES surveys. It then restructures the data into a convenient format using xarray. This simplifies the challenge of dealing with multiple attenuator and gain settings, makes stacking a straight-forward xarray operation, and allows us to store the large datasets in efficient zarr format in cloud storage.

The structure of the resulting xarray depends on if the ApRES data were collected in attended or unattended mode.

## Installation

```
pip install xapres
```

## Usage

See the notebooks/guides directory for examples of how to use both the core processing code and how to restructure the resulting profiles and chirps into an xarray.

The most useful guide is notebooks/guides/UsingXaPRES.ipynb.

### Quick start
A common thing that you may want to do with XApRES is to gather multiple ApRES measurements, which are stored in .dat files, into one xarray. 

The fastest way to do this is:
```
# install the package
!pip install xapres
# import the package
import xapres as xa
# load the chirps, perform an fft, and put them all in an xarray
directory = 'data/sample/multi-burst-dat-file/'
data = xa.load.generate_xarray(directory=directory)
# stack the chirps in each burst, select one of the attenuator pairs, compute the decibels and plot
data.profile.mean(dim='chirp_num').isel(attenuator_setting_pair=0).dB().plot(x='time', yincrease=False)
```

You just need to change `directory` to the location of your .DAT files and the code will search recursively through the directory and its sub-directories to find and process all the .DAT they contain. 


### Writing and loading from zarr

The package includes the capability to write data to zarr stores, which can be accessed efficiently without immediately loading all the data to disk. This is particularly useful when performing analysis in the cloud, but can be useful when inspecting the data locally too. 

For example, to lazily access (meaning that no data is downloaded immediately) some recent ApRES data from Greenland, use:

```
import xarray as xr
def reload(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/{site}'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        chunks={}) 
    return ds
A101 = reload("A101")
A103 = reload("A103")
A104 = reload("A104")
```

Alternatively, you can use a function built-in to the package which loads these data, and also adds some functionality to the xarray it returns: 

```
import xapres as xa
ds = xa.load.load_zarr(site="A101")
ds
```

The dataset `ds` containing various variables, including the complex depth profile. 

### Functionality added to xarrays by XApRES

You can compute decibels from complex depth profiles using 

```
ds.profile.dB().compute()
```

You can also sonify the chirp. Because the frequencies contained in the chirps are mostly in the audible range, you can play them through the computer's speakers and hear what the ApRES data sound like:

```
ds.chirp.isel(time=300, chirp_num=0, attenuator_setting_pair=0).sonify()
```

This plays the chirp and if you set `save = True` as an input it will save the audio file as a .wav file. 

Contributions to this project are very welcome! Please feel free to contact us through github issues. 


--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
