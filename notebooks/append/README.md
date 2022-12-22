# Appending to zarr stores

A collection of notebooks used ot examine the approach of appending each ApRES data file (in the form of an xarray created by `ApRESDefs.xapres`) to a zarr store, one-by-one. 

This approach runs into problems when loading from the zarr store that you have just appended to, does not appear to yield an xarray of the expected length. This problem is discussed [here](https://github.com/pydata/xarray/issues/5878). 

zarr_append_not_working.ipynb shows this the clearest. 

The issues described in these notebooks and hte issue linked above has not been solved, but the final solution to getting the full zarr stores does not use append (see the ../to_zarr directory). 