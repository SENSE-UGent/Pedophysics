import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.water_ec import SenGoode
from .temperature import *
from .water_ec import *

def Salinity(soil):
    """
    Return or compute missing values of the soil.salinity attribute.

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
            Data Frame containing the qualitative information about all array-like soil attributes for each state
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        An array of soil salinity (NaCl) of the bulk pore fluid values

    Notes
    -----
    This function modifies the soil object in-place, updating the `df` dataframe and `info`
    dataframe if necessary.

    See Also
    --------
    WaterEC : Function to return or compute the water's electrical conductivity.
    Temperature : Function to return or compute the soil's temperature.
    SenGoode : Function to calculte water_ec based on soil.salinity and soil.temperature

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

        soil.info['salinity'] = ["Calculated using SenGood function in predict.Salinity" if np.isnan(soil.df.salinity[x]) or soil.info.salinity[x] == "Calculated using SenGood function in predict.Salinity"
                                 else soil.info.salinity[x] for x in range(soil.n_states)]

        soil.df['salinity'] = [sal[x] if np.isnan(soil.df.salinity[x]) else soil.df.salinity[x] for x in range(soil.n_states)]

    return soil.df.salinity.values 

