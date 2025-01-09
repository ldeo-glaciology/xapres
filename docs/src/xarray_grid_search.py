from typing import Callable
import pandas as pd
import xarray as xr
from tqdm.autonotebook import tqdm
import xarray as xr

def gridSearch(function: Callable, **kwargs) -> xr.core.dataset.Dataset:
    """
    Perform a grid search by iterating over all combinations of input parameters and running a given function.

    Parameters:
    function (callable): The function to be executed for each combination of input parameters. This function should return either an xarray dataset, an xarray datarray, or a numpy array.
    **kwargs: Keyword arguments representing the input parameters and their corresponding values.

    Returns:
    xr_unstacked (xarray.core.dataset.Dataset): The concatenated and unstacked xarray dataset containing the results of the grid search.

    Example:
    #### Define a function to be executed for each combination of input parameters
    def my_function(param1, param2):
        ##### Perform some computation using the input parameters
        result = param1 + param2
        return result

    #### Perform a grid search by iterating over all combinations of input parameters
    results = gridSearch(my_function, param1=[1, 2, 3], param2=[4, 5])
    
    """

    # extract the names of the parameters
    p_names = [x for x in kwargs] 

    # extract the values of the parameters
    p_values_list = [x for x in kwargs.values()]
    
    # create a multiIndex from the parameter names and values
    multiIndex = pd.MultiIndex.from_product(p_values_list, names=p_names)


    #loop over every conbimation of parameters stored in multiIndex
    xr_out_list = []
    for mi in tqdm(multiIndex):

        # create a dictionary of inputs for the function from the values stored in multiIndex
        inputs = {p_names[x]: mi[x] for x in range(len(p_names))}

        # run the function with this combination of inputs
        single_iteration_result = function(**inputs)

        # add coordinates to the result and store as as either a DataSet or a dataArray
        if isinstance(single_iteration_result, xr.core.dataset.Dataset):
            xr_out_new = single_iteration_result.assign_coords(inputs)
        else:
            xr_out_new = xr.DataArray(single_iteration_result, coords=inputs)    # use this line if the function returns a data array, or a numpy array
        
        # append the result to a list
        xr_out_list.append(xr_out_new)

    # concatenate the list of results into a single xarray
    xr_stacked = xr.concat(xr_out_list, dim='stacked_dim')

    # add the multiIndex to the xarray
    mindex_coords = xr.Coordinates.from_pandas_multiindex(multiIndex, 'stacked_dim')
    xr_stacked = xr_stacked.assign_coords(mindex_coords)

    # unstack the xarray - i.e. separate the multiIndex into separate dimensions
    xr_unstacked = xr_stacked.unstack()

    # convert to a dataset if the result is a data array
    if isinstance(xr_unstacked, xr.DataArray):
        xr_unstacked.name = 'result'
        xr_unstacked = xr_unstacked.to_dataset()

    return xr_unstacked

### Define two functions to test the gridSearch function
def returns_an_xrdataset(a: float = 1111.0, b: float = 5555.0, c: float = 2222.0) -> xr.core.dataset.Dataset:
    # load an example dataset
    ds = xr.tutorial.open_dataset("air_temperature")
    # perform some computation on the dataset using the model parameters
    ds['air'] = a * ds['air'] + b * ds['air'] + c
    return ds

def returns_a_float(a: float = 1111.0, b: float = 5555.0, c: float = 2222.0) -> float:
    return a + b + c

# Test the function as follows. This should return a dataset with the outputs of the functions for every combination of the supplied parameter values.
# All other parameters (the ones not specified in the function call) are set to their default values.
#ds = gridSearch(returns_an_xrdataset, a=[0.1, 0.2])
#ds = gridSearch(returns_an_xrdataset, a=[0.1, 0.2], b=[0.3, 0.4])
#ds = gridSearch(returns_an_xrdataset, a=[0.1, 0.2], b=[0.3, 0.4], c=[0.5, 0.6])

#ds = gridSearch(returns_a_float, a=[0.1, 0.2])
#ds = gridSearch(returns_a_float, a=[0.1, 0.2], b=[0.3, 0.4])
#ds = gridSearch(returns_a_float, a=[0.1, 0.2], b=[0.3, 0.4], c=[0.5, 0.6])
