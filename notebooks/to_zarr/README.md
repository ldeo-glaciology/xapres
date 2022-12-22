# Data handling

A collection of notebooks for using `ApRESDefs.xapres` to write ApRES data to a zarr store.

The approach I finally settled on is described in to_individual_zarr.ipynb and write_big_zarrs.ipynb: We load each dat file in turn and compute the profile using `ApRESDefs.xapres`, then rechunk to a chunk size of 1 in the time dimension and -1 in the other dimensions, then write each one to an individual zarr store. 

Then we load all the individual zarr stores, calculate the stacked profiles (by taking the average along the chirp_num dimension), and rechunk to a large chunk size along the time dimension. Rechunking in this way does not cause too much difficulty memory-wise because it is gathering together many small chunks and combining them. Also we avoid this being an issue by making sure that all the chunks are the same size throughout the individual zarr files. 

Finally, we write the resulting xarray (~200GB) to a single zarr store. 

To load the store you can run:
```
def reload(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs/{site}'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        chunks={}) 
    return ds

ds_101 = reload("A101")
ds_103 = reload("A103")
ds_104 = reload("A104")

```