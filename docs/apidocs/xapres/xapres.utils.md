# {py:mod}`xapres.utils`

```{py:module} xapres.utils
```

```{autodoc2-docstring} xapres.utils
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`displacement_timeseries <xapres.utils.displacement_timeseries>`
  - ```{autodoc2-docstring} xapres.utils.displacement_timeseries
    :summary:
    ```
* - {py:obj}`compute_displacement <xapres.utils.compute_displacement>`
  - ```{autodoc2-docstring} xapres.utils.compute_displacement
    :summary:
    ```
* - {py:obj}`combine_profiles <xapres.utils.combine_profiles>`
  - ```{autodoc2-docstring} xapres.utils.combine_profiles
    :summary:
    ```
* - {py:obj}`bin_profiles <xapres.utils.bin_profiles>`
  - ```{autodoc2-docstring} xapres.utils.bin_profiles
    :summary:
    ```
* - {py:obj}`compute_coherence <xapres.utils.compute_coherence>`
  - ```{autodoc2-docstring} xapres.utils.compute_coherence
    :summary:
    ```
* - {py:obj}`phase2range <xapres.utils.phase2range>`
  - ```{autodoc2-docstring} xapres.utils.phase2range
    :summary:
    ```
* - {py:obj}`computeStrainRates <xapres.utils.computeStrainRates>`
  - ```{autodoc2-docstring} xapres.utils.computeStrainRates
    :summary:
    ```
* - {py:obj}`dB <xapres.utils.dB>`
  - ```{autodoc2-docstring} xapres.utils.dB
    :summary:
    ```
* - {py:obj}`sonify <xapres.utils.sonify>`
  - ```{autodoc2-docstring} xapres.utils.sonify
    :summary:
    ```
* - {py:obj}`addProfileToDs <xapres.utils.addProfileToDs>`
  - ```{autodoc2-docstring} xapres.utils.addProfileToDs
    :summary:
    ```
* - {py:obj}`computeProfile <xapres.utils.computeProfile>`
  - ```{autodoc2-docstring} xapres.utils.computeProfile
    :summary:
    ```
* - {py:obj}`default_constants <xapres.utils.default_constants>`
  - ```{autodoc2-docstring} xapres.utils.default_constants
    :summary:
    ```
* - {py:obj}`add_methods_to_xarrays <xapres.utils.add_methods_to_xarrays>`
  - ```{autodoc2-docstring} xapres.utils.add_methods_to_xarrays
    :summary:
    ```
* - {py:obj}`contains_dask_array <xapres.utils.contains_dask_array>`
  - ```{autodoc2-docstring} xapres.utils.contains_dask_array
    :summary:
    ```
````

### API

````{py:function} displacement_timeseries(self: xarray.DataArray, offset: int = 1, bin_size: int = 20, lower_limit_on_fit: float = 800.0)
:canonical: xapres.utils.displacement_timeseries

```{autodoc2-docstring} xapres.utils.displacement_timeseries
```
````

````{py:function} compute_displacement(profile1_unaligned: xarray.DataArray, profile2_unaligned: xarray.DataArray, bin_size: int = 20, lower_limit_on_fit: float = 800.0)
:canonical: xapres.utils.compute_displacement

```{autodoc2-docstring} xapres.utils.compute_displacement
```
````

````{py:function} combine_profiles(profile1_unaligned, profile2_unaligned)
:canonical: xapres.utils.combine_profiles

```{autodoc2-docstring} xapres.utils.combine_profiles
```
````

````{py:function} bin_profiles(profiles, bin_size)
:canonical: xapres.utils.bin_profiles

```{autodoc2-docstring} xapres.utils.bin_profiles
```
````

````{py:function} compute_coherence(b1_binned, b2_binned)
:canonical: xapres.utils.compute_coherence

```{autodoc2-docstring} xapres.utils.compute_coherence
```
````

````{py:function} phase2range(phi, lambdac=0.5608, rc=None, K=200000000.0, ci=168230000.0)
:canonical: xapres.utils.phase2range

```{autodoc2-docstring} xapres.utils.phase2range
```
````

````{py:function} computeStrainRates(self, lower_limit_on_fit=800)
:canonical: xapres.utils.computeStrainRates

```{autodoc2-docstring} xapres.utils.computeStrainRates
```
````

````{py:function} dB(self)
:canonical: xapres.utils.dB

```{autodoc2-docstring} xapres.utils.dB
```
````

````{py:function} sonify(self, play=True, save=False, wav_filename='chirp')
:canonical: xapres.utils.sonify

```{autodoc2-docstring} xapres.utils.sonify
```
````

````{py:function} addProfileToDs(self: xarray.Dataset, **kwargs)
:canonical: xapres.utils.addProfileToDs

```{autodoc2-docstring} xapres.utils.addProfileToDs
```
````

````{py:function} computeProfile(self: xarray.DataArray, pad_factor=2, drop_noisy_chirps=False, clip_threshold=1.2, min_chirps=0, demean=True, detrend=False, stack=False, scale_for_window=True, crop_chirp_start=0, crop_chirp_end=1, max_range=None, constants={})
:canonical: xapres.utils.computeProfile

```{autodoc2-docstring} xapres.utils.computeProfile
```
````

````{py:function} default_constants()
:canonical: xapres.utils.default_constants

```{autodoc2-docstring} xapres.utils.default_constants
```
````

````{py:function} add_methods_to_xarrays()
:canonical: xapres.utils.add_methods_to_xarrays

```{autodoc2-docstring} xapres.utils.add_methods_to_xarrays
```
````

````{py:function} contains_dask_array(dataarray)
:canonical: xapres.utils.contains_dask_array

```{autodoc2-docstring} xapres.utils.contains_dask_array
```
````
