import numpy as np

def SolidPerm(soil):
    """
    Set missing values of soil.df.solid_perm and return

    If any value of the solid_perm attribute is missing (NaN), 
    it will be set to a default value of 4. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - solid_perm : array-like
            Soil solid real relative dielectric permittivity phase [-]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing the qualitative information about all array-like soil attributes for each state
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.solid_perm.values: an array of updated soil solid real relative dielectric permittivity phase values

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.solid_perm
    0   NaN
    Name: solid_perm, dtype: float64
    >>> SolidPerm(sample)
    >>> sample.df.solid_perm
    0    4
    Name: solid_perm, dtype: float64
    """

    # Check if any value of solid_perm is missing
    if (np.isnan(soil.df.solid_perm)).any:  

        soil.info['solid_perm'] = ["Set as 4 by default" if np.isnan(soil.df.solid_perm[x]) or soil.info.solid_perm[x] == "Set as 4 by default"
                                     else soil.info.solid_perm[x] for x in range(soil.n_states)]
        
        soil.df.loc[np.isnan(soil.df['solid_perm']), ['solid_perm']] = 4

    return soil.df.solid_perm.values