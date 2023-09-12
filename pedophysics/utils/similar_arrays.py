import numpy as np

def arrays_are_similar(a, b):
    # Check if the two arrays are of the same shape
    if a.shape != b.shape:
        return False
    
    # Check if the non-NaN elements are close to each other
    non_nan_match = np.isclose(a[~np.isnan(a)], b[~np.isnan(b)]).all()

    # Check if the NaN locations are the same in both arrays
    nan_match = np.isnan(a) == np.isnan(b)

    return non_nan_match and nan_match.all()