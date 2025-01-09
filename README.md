<!-- SPHINX-START-proj-desc -->
# xApRES

A package for loading and processing data collected using the Autonomous phase-sensitive Radio-Echo Sounder (ApRES). 

ApRES is an ice-penetrating radar used for measuring ice-shelf basal melt rates, englacial deformation, and other englacial changes like water content. ApRES was developed by University College London and The British Antarctic Survey. 

The xApRES package uses [xarray](https://xarray.pydata.org/en/stable/), a package that simplifies indexing, grouping, aggregating, and plotting multi-dimensional data (like ApRES data). Using xarray in xapres simplifies the challenge of dealing with multiple measurements taken using different attenuator and gain settings, using different antenna orientations, and collected in different locations. xApRES also simplifies many common processing steps, for example averaging chirps (i.e. stacking) stacking or resampling in time. Finally, using xarray as the basis for storage and processing of ApRES data helps when dealing with very large ApRES datasets because it facilitates cloud computing and storage. 

A key goal is to allow the loading, processing and plotting of full ApRES datasets, collected in either unattended or attended mode, with only a few lines of code. 

For example, loading raw ApRES data (stored in .dat files), performing an fft to get depth profiles and plotting can be achieved in one command as follows:


<!-- SPHINX-END-proj-desc -->

```
xapres.load.generate_xarray(directory='data/thwaites/')\
    .profile.mean(dim='chirp_num')\
    .isel(attenuator_setting_pair=0)\
    .dB()\
    .plot.line(hue="time", xlim = (0, 2500));
```

![Example plot](docs/src/images/plot_for_intro.png) 

## Installation

```
pip install xapres
```

## Documentation

[https://ldeo-glaciology.github.io/xapres](https://ldeo-glaciology.github.io/xapres) 

## Quick start


A common use of xApRES is gathering multiple ApRES measurements stored in .dat files into one xarray dataset. 

The fastest way to do this is:
```
import xapres as xa
directory = 'data/sample/thwaites/'
data = xa.load.generate_xarray(directory=directory)
```

A typical next step is to stack the data (i.e. average the chirps in each burst) and plot the resulting profile. 

```
profiles = data.profile.mean(dim='chirp_num')
profiles.isel(attenuator_setting_pair=0).dB().plot(x='time', yincrease=False)
```
    
Finally, you can compute velocities and strain rates from the displacement timeseries and plot them. 

```
w = profiles.displacement_timeseries()
w.velocity.plot(y = 'bin_depth', yincrease=False, xlim = (-2, 7), ylim = (1200,0))
```

You just need to change `directory` to the location of your .dat files and the code will search recursively through the directory and its sub-directories to find and process all the .dat files they contain. 

## Contributing
Any and all contributions to this project are very welcome! Please feel free to contact us through github issues. 


