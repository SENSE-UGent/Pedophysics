import numpy as np

def arrays_are_similar(array1, array2, tol=1e-5):
    """
    Check if two numpy arrays are similar with a given tolerance.

        Parameters
    ----------
    array1 (numpy.ndarray): The first array to compare.
    array2 (numpy.ndarray): The second array to compare.
    tol (float): The tolerance for the comparison. Default is 1e-5.

    Returns
    -------
    bool: 
        Returns True if the arrays are considered similar, otherwise False.

    Example
    -------
    >>> a = np.array([1.0, 2.0, np.nan, 4.0])
    >>> b = np.array([0.999, 2.001, np.nan, 4.0])

    >>> similar = arrays_are_similar(a, b)

    """

    # Check if both arrays are close within the tolerance specified
    return np.allclose(array1, array2, atol=tol, equal_nan=True)