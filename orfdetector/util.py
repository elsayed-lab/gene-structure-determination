"""
Miscellaneous helper functions
"""
import numpy as np

def max_filter(arr, width=10):
    """
    Moves a sliding window along a 1d array and sets all values except for the
    max in each window to 0.

    This effectively filters out everything except the local maxima along an
    array.

    For example, given: 
    
        0 0 0 1 2 5 1 0 0 0 0 0 4 15 1 0 0 0 0 0 0 84 9 0 0 0 0 0

    max_filter(x, 5) returns:
    
        0 0 0 0 0 5 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 84 0 0 0 0 0 0

    Arguments
    ---------
    arr: numpy.array
        1-dimensional range of numbers to filter on.
    width: int
        Size of filtering window to use.

    Returns
    -------
    arr: np.array
        Filtered array of same length as input.
    """
    arr = arr.copy()
    strides = rolling_window(arr, width)

    # return np.apply_along_axis(lambda x: np.where(x != max(x), 0, x), 1, strides)

    # there is probably a better way to do this..
    for i in range(strides.shape[0]):
        strides[i,][strides[i,] != max(strides[i,])] = 0

    return arr

def rolling_window(x, width):
    """
    Returns an array of windows for a specified width across a given array.
    
    Source: http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html

    Arguments
    ---------
    x: np.array
        1-dimensional array
    width: int
        Size of rolling window.
    """
    shape = x.shape[:-1] + (x.shape[-1] - width + 1, width)
    strides = x.strides + (x.strides[-1],)
    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)

