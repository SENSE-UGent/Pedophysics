import numpy as np

def Temperature(soil):
    """
    Return and set missing values of the soil.df.temperature attribute.

    If any value of the temperature attribute is missing (NaN), 
    it will be set to a default value of 298.15. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - temperature : array-like
            Soil bulk temperature [K]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        An array of updated soil bulk temperature values

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.temperature
    0   NaN
    Name: temperature, dtype: float64
    >>> Temperature(sample)
    >>> sample.df.temperature
    0    298.15
    Name: temperature, dtype: float64
    """

    # Check if any value of solid_ec is missing
    if (np.isnan(soil.df.temperature)).any():

        soil.info['temperature'] = ["Set as 298.15 K by default" if np.isnan(soil.df.temperature[x]) or soil.info.temperature[x] == "Set as 298.15 K by default"
                                     else soil.info.temperature[x] for x in range(soil.n_states)]
        
        soil.df.loc[(np.isnan(soil.df['temperature'])), ['temperature']] = 298.15

    return soil.df.temperature.values