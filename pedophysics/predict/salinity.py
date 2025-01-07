import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.water_ec import SenGoode
from .temperature import *
from .water_ec import *

def Salinity(soil):
    """
    Calculate missing values of soil.df.salinity and return

    If any value of the salinity attribute is missing (NaN), it will first compute 
    the missing values by optimizing the SenGoode function based on the soil's water 
    electrical conductivity and temperature.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - temperature : array-like
            Soil bulk temperature [K]
        - salinity : array-like
            Soil salinity (NaCl) of the bulk pore fluid [mol/L]
        - water_ec : array-like
            Soil water real electrical conductivity [S/m]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.salinity.values: an array of soil salinity (NaCl) of the bulk pore fluid values

    Notes
    -----
    This function modifies the soil object in-place, updating the `df` dataframe and `info`
    dataframe if necessary.

    External functions
    --------
    WaterEC : Compute missing values of soil.df.water_ec and return  
    Temperature : Set missing values of soil.df.temperature and return 
    SenGoode : Calculate soil water real electrical conductivity using the Sen and Goode model and return

    Example
    -------
    >>> sample = Soil(water_ec = 0.1)
    >>> sample.df.salinity
    0   NaN
    Name: salinity, dtype: float64
    >>> Salinity(sample)
    >>> sample.df.salinity
    0    0.00846
    Name: salinity, dtype: float64
    """

    if any(np.isnan(soil.df.salinity[x])for x in range(soil.n_states)):  # Go over if any value is missing 

        WaterEC(soil)
        Temperature(soil)
        sal = []

        def objective_salinity(salinity, water_ec, temperature):
            return (SenGoode(temperature, salinity) - water_ec)**2

        for x in range(soil.n_states):
            result = minimize(objective_salinity, 0.01, args=(soil.df.water_ec[x], soil.df.temperature[x]), bounds=[(0, 1)])
            sal.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))

        missing_salinity_before = soil.df['salinity'].isna()

        soil.df['salinity'] = [sal[x] if np.isnan(soil.df.salinity[x]) 
                               else soil.df.salinity[x] for x in range(soil.n_states)]
        
        missing_salinity_after = soil.df['salinity'].isna()

        soil.info['salinity'] = [str(soil.info.salinity[x]) + (
                "--> Calculated using SenGood function in predict.Salinity"
                if missing_salinity_before[x] and not missing_salinity_after[x]
                else "--> Provide salinity; otherwise, water_ec"
                if missing_salinity_before[x] and missing_salinity_after[x]
                else "")
            if missing_salinity_before[x]
            else soil.info.salinity[x]
            for x in range(soil.n_states)]
        

    return soil.df.salinity.values 
