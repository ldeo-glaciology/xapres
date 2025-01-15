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

````{py:function} generate_xarray(directory=None, file_numbers_to_process=None, file_names_to_process=None, bursts_to_process='All', attended=False, polarmetric=False, corrected_pad=False, max_range=None, computeProfiles=True, addProfileToDs_kwargs={}, loglevel='warning')
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

````{py:method} load(dat_filename, bursts_to_process='All', corrected_pad=False, attended=False, polarmetric=False, max_range=None, computeProfiles=True, addProfileToDs_kwargs={})
:canonical: xapres.load.from_dats.load

```{autodoc2-docstring} xapres.load.from_dats.load
```

````

````{py:method} list_files(directory=None, search_suffix='')
:canonical: xapres.load.from_dats.list_files

```{autodoc2-docstring} xapres.load.from_dats.list_files
```

````

````{py:method} load_all(directory=None, file_numbers_to_process=None, file_names_to_process=None, bursts_to_process='All', disable_progress_bar=True, attended=False, polarmetric=False, corrected_pad=False, max_range=None, computeProfiles=True, addProfileToDs_kwargs={})
:canonical: xapres.load.from_dats.load_all

```{autodoc2-docstring} xapres.load.from_dats.load_all
```

````

````{py:method} subset_files()
:canonical: xapres.load.from_dats.subset_files

```{autodoc2-docstring} xapres.load.from_dats.subset_files
```

````

````{py:method} all_bursts_in_dat_to_xarray(dat_filename, bursts_selected)
:canonical: xapres.load.from_dats.all_bursts_in_dat_to_xarray

```{autodoc2-docstring} xapres.load.from_dats.all_bursts_in_dat_to_xarray
```

````

````{py:method} all_bursts_at_waypoint_to_xarray(directory, waypoint_number)
:canonical: xapres.load.from_dats.all_bursts_at_waypoint_to_xarray

```{autodoc2-docstring} xapres.load.from_dats.all_bursts_at_waypoint_to_xarray
```

````

````{py:method} header_cleaning()
:canonical: xapres.load.from_dats.header_cleaning

```{autodoc2-docstring} xapres.load.from_dats.header_cleaning
```

````

````{py:method} subset_bursts_to_process()
:canonical: xapres.load.from_dats.subset_bursts_to_process

```{autodoc2-docstring} xapres.load.from_dats.subset_bursts_to_process
```

````

````{py:method} _timestamp_from_burst(burst)
:canonical: xapres.load.from_dats._timestamp_from_burst

```{autodoc2-docstring} xapres.load.from_dats._timestamp_from_burst
```

````

````{py:method} chirptime_from_burst(burst)
:canonical: xapres.load.from_dats.chirptime_from_burst

```{autodoc2-docstring} xapres.load.from_dats.chirptime_from_burst
```

````

````{py:method} _get_orientation(filename)
:canonical: xapres.load.from_dats._get_orientation

```{autodoc2-docstring} xapres.load.from_dats._get_orientation
```

````

````{py:method} burst_data(burst)
:canonical: xapres.load.from_dats.burst_data

```{autodoc2-docstring} xapres.load.from_dats.burst_data
```

````

````{py:method} add_attrs()
:canonical: xapres.load.from_dats.add_attrs

```{autodoc2-docstring} xapres.load.from_dats.add_attrs
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

`````
