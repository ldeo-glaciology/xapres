# xapres_package

A package for processing data from the Autonomous phase-sensitive Radio-Echo Sounder (ApRES) using xarray. The core ApRES processing code is adapted from code written by Keith Nicholls, British Antarctic Survey, UK. This package uses Keith's code to process the raw ApRES data from autonomous ApRES surveys. It then restructures the data into a convenient format using xarray. This simplifies the challenge of dealing with multiple attenuator and gain settings, makes stacking a straight-forward xarray operation, and allows us to store the large datasets in efficient zarr format in cloud storage. 

## Installation

```
pip install xapres-package
```

## Usage
See the notebooks/guides directory for examples of how to use both the core processing code and how to restructure the resulting profiles and chirps into an xarray.

The package includes the capability to write data to zarr stores, which can be accessed efficiently without immediately loading all the data to disk. This is particularly useful when performing analysis in the cloud, but can be useful when inspecting the data locally too. 

For example, to lazily access (meaning that no data is downloaded immediately) some recent data from Greenland  use:

```
def reload(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs/{site}'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        chunks={}) 
    return ds
A101 = reload("A101')
A103 = reload("A103')
A104 = reload("A104')
```

--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
