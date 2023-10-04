import pytest
from xapres_package import load

## test polarmetric local loading by loading the same waypoint twice as if it is two different ones and chaeking 
# that you get the same thing twice.

# tests the attended option and the polarmetric option from locally stored dat files (as opposed to cloud stored dat files)

def test_polarmetric_load():
    
    xa = load.load_from_dat()
    xa.load_all(attended=True, 
                directory=["data/sample/polarmetric", "data/sample/polarmetric"], 
                polarmetric=True)
    
    assert len(xa.data.waypoint) == 2
    assert all(xa.data.isel(waypoint=0).filename.values == xa.data.isel(waypoint=1).filename.values)
