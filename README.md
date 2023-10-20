# xapres_package

A package for processing data from the Autonomous phase-sensitive Radio-Echo Sounder (ApRES) using xarray. The core ApRES processing code is adapted from code written by Keith Nicholls, British Antarctic Survey, UK. This package uses Keith's code to process the raw ApRES data from ApRES surveys. It then restructures the data into a convenient format using xarray. This simplifies the challenge of dealing with multiple attenuator and gain settings, makes stacking a straight-forward xarray operation, and allows us to store the large datasets in efficient zarr format in cloud storage.

The structure of the resulting xarray depends on if the ApRES data were collected in attended or unattended mode.

## Installation

```
pip install xapres-package
```

## Usage
See the notebooks/guides directory for examples of how to use both the core processing code and how to restructure the resulting profiles and chirps into an xarray.

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
import xapres_package as xa
ds = xa.load_zarr("A101")
ds
```

Now you can compute decibels from a complex profile simply using 

```
ds.profile.isel(time=300, chirp_num=0, attenuator_setting_pair=0).dB()
```

You can also sonify the chirp. Because the frequencies contained in the chirps largely in the audible range, you can play them through the computer's speakers and hear what the ApRES data sound like:

```
ds.chirp.isel(time=300, chirp_num=0, attenuator_setting_pair=0).sonify()
```

This plays the chirp and if you set `save = True` as an input it will save the audio file as a .wav file. 

Contributions to this project are very welcome! Please feel free to contact us through github issues. 




--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
