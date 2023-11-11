import pytest
from xapres import load

## test polarmetric local loading by loading the same waypoint twice as if it is two different ones and chaeking 
# that you get the same thing twice.

# tests the attended option and the polarmetric option from locally stored dat files (as opposed to cloud stored dat files)

def test_polarmetric_load():
    
    fs = load.from_dats()
    fs.load_all(attended=True, 
                directory=["data/sample/polarmetric", "data/sample/polarmetric"], 
                polarmetric=True)
    
    assert len(fs.data.waypoint) == 2
    assert all(fs.data.isel(waypoint=0).filename.values == fs.data.isel(waypoint=1).filename.values)
