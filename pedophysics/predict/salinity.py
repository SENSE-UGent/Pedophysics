import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.water_ec import SenGoode
from .temperature import *
from .water_ec import *

def Salinity(soil):
    '''
    
    '''
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

