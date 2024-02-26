import numpy as np
from pedophysics import instruments

def FrequencyPerm(soil): 
    """
    Set missing values of soil.df.frequency_perm and return 

    If any value of the frequency_perm attribute is missing (NaN), 
    it will be computed using the `Inst2FreqP` function from the `instruments` module. 

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - frequency_perm : array-like
            Frequency of dielectric permittivity measurement [Hz]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state. Includes: frequency_perm

    Returns
    -------
    np.ndarray
        soil.df.frequency_perm.values: an array of updated frequency of dielectric permittivity measurement values

    Notes
    -----
    This function modifies the soil object in-place, updating the `df` dataframe 
    if necessary.

    External functions
    --------
    Inst2FreqP : Set missing values of soil.df.frequency_perm and return

    Example
    -------
    >>> sample = Soil(instrument = 'GPR')
    >>> sample.df.frequency_perm
    0   NaN
    Name: frequency_perm, dtype: float64
    >>> FrequencyPerm(sample)
    >>> sample.df.frequency_perm
    0    1e9
    Name: frequency_perm, dtype: float64
    """

    # Check if any value of frequency_perm is missing
    if (np.isnan(soil.df.frequency_perm)).any(): 
        instruments.Inst2FreqP(soil)

    return soil.df.frequency_perm.values