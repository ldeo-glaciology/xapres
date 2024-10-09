import pytest
import xapres
import xarray as xr
def test_bound_methods_are_added_correctly():
    assert xr.DataArray.dB
    assert xr.DataArray.sonify
    assert xr.DataArray.displacement_timeseries