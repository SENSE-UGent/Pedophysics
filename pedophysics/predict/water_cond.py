import numpy as np
from pedophysics.pedophysical_models.water_cond import SenGoode
from .temperature import *

def WaterCond(soil):
    '''
    
    '''
    if (np.isnan(soil.water_cond)).any():  # Go over if any value is missing 
        Temperature(soil)
        soil.df['water_cond'] = [SenGoode(soil.df.temperature.values[x], soil.df.salinity.values[x]) if np.isnan(soil.df.water_cond[x]) else soil.df.water_cond[x] for x in range(soil.n_states)]

    return soil.df.water_cond.values