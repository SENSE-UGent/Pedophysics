import numpy as np

def BulkPermInf(soil):
    """
    Set missing values of soil.df.bulk_perm_inf and return.

    If any value of the bulk_perm_inf attribute is missing (NaN), 
    it will be set to a default value of 5. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - bulk_perm_inf : array-like
            Soil bulk real relative permittivity at infinite frequency [-]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.bulk_perm_inf.values: an array of updated soil bulk real relative permittivity at infinite frequency values

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.bulk_perm_inf
    0   NaN
    Name: bulk_perm_inf, dtype: float64
    >>> BulkPermInf(sample)
    >>> sample.df.bulk_perm_inf
    0    5
    Name: bulk_perm_inf, dtype: float64
    """
    
    # Check if any value of bulk_perm_inf is missing
    if (np.isnan(soil.df.bulk_perm_inf)).any():

        soil.info['bulk_perm_inf'] = ["Set as 5 by default" if np.isnan(soil.df.bulk_perm_inf[x]) or soil.info.bulk_perm_inf[x] == "Set as 5 by default"
                                     else soil.info.bulk_perm_inf[x] for x in range(soil.n_states)]

        soil.df['bulk_perm_inf'] = [5 if np.isnan(soil.df.bulk_perm_inf[x]) else soil.df.bulk_perm_inf[x] for x in range(soil.n_states)]

    return soil.df.bulk_perm_inf.values