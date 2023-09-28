import pytest
from xapres_package import ApRESDefs

## test polarmetric local loading by loading the same waypoint twice as if it is two different ones and chacking that you get the same thing twice.
def test_file_selection_methods():
    
    xa = ApRESDefs.xapres()
    xa.load_all(attended=True, 
                directory=["data/sample/polarmetric", "data/sample/polarmetric"], 
                polarmetric=True)
    
    assert len(xa.data.waypoint) == 2
    assert all(xa.data.isel(waypoint=0).filename.values == xa.data.isel(waypoint=1).filename.values)
