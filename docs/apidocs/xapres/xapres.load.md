# {py:mod}`xapres.load`

```{py:module} xapres.load
```

```{autodoc2-docstring} xapres.load
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`from_dats <xapres.load.from_dats>`
  - ```{autodoc2-docstring} xapres.load.from_dats
    :summary:
    ```
* - {py:obj}`DataFileObject <xapres.load.DataFileObject>`
  - ```{autodoc2-docstring} xapres.load.DataFileObject
    :summary:
    ```
* - {py:obj}`BurstObject <xapres.load.BurstObject>`
  - ```{autodoc2-docstring} xapres.load.BurstObject
    :summary:
    ```
* - {py:obj}`ChirpObject <xapres.load.ChirpObject>`
  - ```{autodoc2-docstring} xapres.load.ChirpObject
    :summary:
    ```
* - {py:obj}`ProfileObject <xapres.load.ProfileObject>`
  - ```{autodoc2-docstring} xapres.load.ProfileObject
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`load_zarr <xapres.load.load_zarr>`
  - ```{autodoc2-docstring} xapres.load.load_zarr
    :summary:
    ```
* - {py:obj}`generate_xarray <xapres.load.generate_xarray>`
  - ```{autodoc2-docstring} xapres.load.generate_xarray
    :summary:
    ```
````

### API

````{py:function} load_zarr(directory='gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/A101')
:canonical: xapres.load.load_zarr

```{autodoc2-docstring} xapres.load.load_zarr
```
````

````{py:function} generate_xarray(directory=None, remote_load=False, file_numbers_to_process=None, file_names_to_process=None, bursts_to_process='All', attended=False, polarmetric=False, legacy_fft=False, corrected_pad=False, max_range=None, computeProfiles=True, addProfileToDs_kwargs={}, loglevel='warning')
:canonical: xapres.load.generate_xarray

```{autodoc2-docstring} xapres.load.generate_xarray
```
````

`````{py:class} from_dats(loglevel='warning')
:canonical: xapres.load.from_dats

```{autodoc2-docstring} xapres.load.from_dats
```

```{rubric} Initialization
```

```{autodoc2-docstring} xapres.load.from_dats.__init__
```

````{py:method} load_single(dat_filename, burst_number=0, chirp_num=0)
:canonical: xapres.load.from_dats.load_single

```{autodoc2-docstring} xapres.load.from_dats.load_single
```

````

````{py:method} load_dat_file(dat_filename)
:canonical: xapres.load.from_dats.load_dat_file

```{autodoc2-docstring} xapres.load.from_dats.load_dat_file
```

````

````{py:method} list_files(directory=None, search_suffix='')
:canonical: xapres.load.from_dats.list_files

```{autodoc2-docstring} xapres.load.from_dats.list_files
```

````

````{py:method} load_all(directory=None, remote_load=None, file_numbers_to_process=None, file_names_to_process=None, bursts_to_process='All', attended=False, polarmetric=False, legacy_fft=False, corrected_pad=False, max_range=None, computeProfiles=True, addProfileToDs_kwargs={})
:canonical: xapres.load.from_dats.load_all

```{autodoc2-docstring} xapres.load.from_dats.load_all
```

````

````{py:method} subset_files()
:canonical: xapres.load.from_dats.subset_files

```{autodoc2-docstring} xapres.load.from_dats.subset_files
```

````

````{py:method} _all_bursts_in_dat_to_xarray(dat, bursts_selected)
:canonical: xapres.load.from_dats._all_bursts_in_dat_to_xarray

```{autodoc2-docstring} xapres.load.from_dats._all_bursts_in_dat_to_xarray
```

````

````{py:method} _all_bursts_at_waypoint_to_xarray(directory, waypoint_number)
:canonical: xapres.load.from_dats._all_bursts_at_waypoint_to_xarray

```{autodoc2-docstring} xapres.load.from_dats._all_bursts_at_waypoint_to_xarray
```

````

````{py:method} _burst_to_xarray_unattended(burst: xarray.Dataset)
:canonical: xapres.load.from_dats._burst_to_xarray_unattended

```{autodoc2-docstring} xapres.load.from_dats._burst_to_xarray_unattended
```

````

````{py:method} _burst_to_xarray_attended(burst, waypoint_number)
:canonical: xapres.load.from_dats._burst_to_xarray_attended

```{autodoc2-docstring} xapres.load.from_dats._burst_to_xarray_attended
```

````

````{py:method} _burst_to_3d_arrays(burst)
:canonical: xapres.load.from_dats._burst_to_3d_arrays

```{autodoc2-docstring} xapres.load.from_dats._burst_to_3d_arrays
```

````

````{py:method} _coords_from_burst(burst)
:canonical: xapres.load.from_dats._coords_from_burst

```{autodoc2-docstring} xapres.load.from_dats._coords_from_burst
```

````

````{py:method} _timestamp_from_burst(burst)
:canonical: xapres.load.from_dats._timestamp_from_burst

```{autodoc2-docstring} xapres.load.from_dats._timestamp_from_burst
```

````

````{py:method} _get_orientation(filename)
:canonical: xapres.load.from_dats._get_orientation

```{autodoc2-docstring} xapres.load.from_dats._get_orientation
```

````

````{py:method} _set_max_range(burst)
:canonical: xapres.load.from_dats._set_max_range

```{autodoc2-docstring} xapres.load.from_dats._set_max_range
```

````

````{py:method} _add_attrs()
:canonical: xapres.load.from_dats._add_attrs

```{autodoc2-docstring} xapres.load.from_dats._add_attrs
```

````

````{py:method} correct_temperature(threshold=300, correction=-512)
:canonical: xapres.load.from_dats.correct_temperature

```{autodoc2-docstring} xapres.load.from_dats.correct_temperature
```

````

````{py:method} _setup_logging(loglevel)
:canonical: xapres.load.from_dats._setup_logging

```{autodoc2-docstring} xapres.load.from_dats._setup_logging
```

````

````{py:method} is_this_a_remote_load(filename=None)
:canonical: xapres.load.from_dats.is_this_a_remote_load

```{autodoc2-docstring} xapres.load.from_dats.is_this_a_remote_load
```

````

````{py:method} _try_logging()
:canonical: xapres.load.from_dats._try_logging

```{autodoc2-docstring} xapres.load.from_dats._try_logging
```

````

`````

`````{py:class} DataFileObject(Filename, remote_load=False)
:canonical: xapres.load.DataFileObject

```{autodoc2-docstring} xapres.load.DataFileObject
```

```{rubric} Initialization
```

```{autodoc2-docstring} xapres.load.DataFileObject.__init__
```

````{py:method} ExtractBurst(BurstNo)
:canonical: xapres.load.DataFileObject.ExtractBurst

```{autodoc2-docstring} xapres.load.DataFileObject.ExtractBurst
```

````

`````

`````{py:class} BurstObject()
:canonical: xapres.load.BurstObject

```{autodoc2-docstring} xapres.load.BurstObject
```

```{rubric} Initialization
```

```{autodoc2-docstring} xapres.load.BurstObject.__init__
```

````{py:method} ExtractChirp(ChirpList)
:canonical: xapres.load.BurstObject.ExtractChirp

```{autodoc2-docstring} xapres.load.BurstObject.ExtractChirp
```

````

````{py:method} PlotBurst()
:canonical: xapres.load.BurstObject.PlotBurst

```{autodoc2-docstring} xapres.load.BurstObject.PlotBurst
```

````

`````

`````{py:class} ChirpObject()
:canonical: xapres.load.ChirpObject

```{autodoc2-docstring} xapres.load.ChirpObject
```

```{rubric} Initialization
```

```{autodoc2-docstring} xapres.load.ChirpObject.__init__
```

````{py:method} PlotChirp()
:canonical: xapres.load.ChirpObject.PlotChirp

```{autodoc2-docstring} xapres.load.ChirpObject.PlotChirp
```

````

````{py:method} FormProfile(F0=200000000, F1=400000000, pad=2, ref=1, corrected_pad=False)
:canonical: xapres.load.ChirpObject.FormProfile

```{autodoc2-docstring} xapres.load.ChirpObject.FormProfile
```

````

`````

`````{py:class} ProfileObject()
:canonical: xapres.load.ProfileObject

```{autodoc2-docstring} xapres.load.ProfileObject
```

```{rubric} Initialization
```

```{autodoc2-docstring} xapres.load.ProfileObject.__init__
```

````{py:method} PlotProfile(dmax)
:canonical: xapres.load.ProfileObject.PlotProfile

```{autodoc2-docstring} xapres.load.ProfileObject.PlotProfile
```

````

`````
