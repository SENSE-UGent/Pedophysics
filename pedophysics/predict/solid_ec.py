import numpy as np

def SolidEC(soil):
    """
    Set missing values of soil.df.solid_ec and return

    If any value of the solid_ec attribute is missing (NaN), 
    it will be set to a default value of 0. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - solid_ec : array-like
            Soil solid real electrical conductivity [S/m]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.solid_ec.values: an array of updated soil solid real electrical conductivity values

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.solid_ec
    0   NaN
    Name: solid_ec, dtype: float64
    >>> SolidEC(sample)
    >>> sample.df.solid_ec
    0    0
    Name: solid_ec, dtype: float64
    """

    # Check if any value of solid_ec is missing
    if (np.isnan(soil.df.solid_ec)).any():

        soil.info['solid_ec'] = ["Set as zero by default" if np.isnan(soil.df.solid_ec[x]) or soil.info.solid_ec[x] == "Set as zero by default"
                                 else soil.info.solid_ec[x] for x in range(soil.n_states)]

        soil.df['solid_ec'] = [0 if np.isnan(soil.df.solid_ec[x]) else soil.df.solid_ec[x] for x in range(soil.n_states)]

    return soil.df.solid_ec.values