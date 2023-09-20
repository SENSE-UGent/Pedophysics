import numpy as np

def AirPerm(soil): 
    """
    Return and set missing values of the soil.df.air_perm attribute.

    If any value of the air_perm attribute is missing (NaN), 
    it will be set to a default value of 1.2. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - air_perm : array-like
            Soil air real relative dielectric permittivity phase [-]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    numpy.ndarray
        An array of updated soil air real relative dielectric permittivity phase values.

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.air_perm
    0   NaN
    Name: air_perm, dtype: float64
    >>> AirPerm(sample)
    >>> sample.df.air_perm
    0    1.2
    Name: air_perm, dtype: float64
    """

    # Check if any value of air_perm is missing
    if (np.isnan(soil.df.air_perm)).any:

        soil.info['air_perm'] = ["Set as 1.2 by default" if np.isnan(soil.df.air_perm[x]) or soil.info.air_perm[x] == "Set as 1.2 by default"
                                     else soil.info.air_perm[x] for x in range(soil.n_states)]

        soil.df['air_perm'] = [1.2 if np.isnan(soil.df.air_perm[x]) else soil.df.air_perm[x] for x in range(soil.n_states)]

    return soil.df.air_perm.values