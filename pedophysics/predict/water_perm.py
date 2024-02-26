import numpy as np
from pedophysics.pedophysical_models.water_perm import *

def WaterPerm(soil):
    """
    Calculate or set missing values of soil.df.water_perm and return

    Determines the soil water phase real dielectric permittivity using either the MalmbergMaryott 
    function or the Olhoeft function based on the soil's salinity and
    frequency permittivity. 
    If the permittivity cannot be determined by either method, it defaults 
    to a value of 80.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - water_perm : array-like
            Soil water phase real dielectric permittivity [-]
        - temperature : array-like
            Soil bulk temperature [K]
        - salinity : array-like
            Soil salinity (NaCl) of the bulk pore fluid [mol/L]
        - frequency_perm : array-like
            Frequency of dielectric permittivity measurement [Hz]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.water_perm.values: array containing the updated soil water phase real dielectric permittivity values.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.

    External functions
    --------
    MalmbergMaryott : Calculate soil water phase real dielectric permittivity using the Malmberg & Maryott model and return
    Olhoeft : Calculate soil water phase real dielectric permittivity using the Olhoeft (1986) model and return

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.water_perm
    0   NaN
    Name: water_perm, dtype: float64
    >>> WaterPerm(sample)
    >>> sample.df.water_perm
    0    80
    Name: water_perm, dtype: float64
    """
    if (np.isnan(soil.df.water_perm)).any(): # Go over if any value is missing 

        soil.info['water_perm'] = ["Calculated using MalmbergMaryott function (RMSE = 0.0046)" if np.isnan(soil.df.water_perm[x]) & ((soil.df.salinity[x] == 0) or np.isnan(soil.df.salinity[x])) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) 
                                    or soil.info.water_perm[x] == "Calculated using MalmbergMaryott function (RMSE = 0.0046)"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
        
        soil.df['water_perm'] = [MalmbergMaryott(soil.df.temperature.values[x]) if np.isnan(soil.df.water_perm[x]) & ((soil.df.salinity[x] == 0) or np.isnan(soil.df.salinity[x])) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) else soil.df.water_perm[x] for x in range(soil.n_states)]
        
        soil.info['water_perm'] = ["Calculated using Olhoeft function" if np.isnan(soil.df.water_perm[x]) & ~np.isnan(soil.df.salinity[x]) & (soil.df.frequency_perm[x] < 100e6)
                                    or soil.info.water_perm[x] == "Calculated using Olhoeft function"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
                
        soil.df['water_perm'] = [Olhoeft(soil.df.temperature.values[x], soil.df.salinity[x]) if np.isnan(soil.df.water_perm[x]) & ~np.isnan(soil.df.salinity[x]) & (soil.df.frequency_perm[x] < 100e6) else soil.df.water_perm[x] for x in range(soil.n_states)]
                
        soil.info['water_perm'] = ["Set as 80 by default" if np.isnan(soil.df.water_perm[x])
                                    or soil.info.water_perm[x] == "Set as 80 by default"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
               
        soil.df['water_perm'] = [80 if np.isnan(soil.df.water_perm[x]) else soil.df.water_perm[x] for x in range(soil.n_states)]

    return soil.df.water_perm.values