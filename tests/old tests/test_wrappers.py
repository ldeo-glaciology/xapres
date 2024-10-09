import pytest
from xapres import load

#  Test `generate_xarray` and `load_zarr` wrappers
def test_wrappers():
    
    from_DAT_unattended = load.generate_xarray(directory='data/sample/single_dat_file/', 
                file_numbers_to_process = [0], 
                bursts_to_process=[0],
                )

    from_zarr = load.load_zarr()
