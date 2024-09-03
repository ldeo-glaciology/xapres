import pytest
from xapres import load 

# test the displacement calculation
def test_displacement_calculation():
    from_zarr = load.load_zarr()  # lazily load a large APRES dataset from Greenland
    p1 = from_zarr.isel(time=2000).profile_stacked # select a profile 
    p2 = from_zarr.isel(time=2100).profile_stacked # select a different profile 
       
    utils.compute_displacement(p1, p2)    # calculate the displacement between the two profiles

    t = from_zarr.sel(time='2022-07-17').profile_stacked   # select all the profiles on a specfic date
    results = t.displacement_timeseries(bin_size = 30, offset = 3) # compute a time series of displacement from these data. Use non-default values for offset and bin_size 

    assert (abs(results)>0).all().load()