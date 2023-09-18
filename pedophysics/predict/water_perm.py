import numpy as np
from pedophysics.pedophysical_models.water_perm import *

def WaterPerm(soil):
    '''
        Soil water phase real dielectric permittivity prediction
    '''
    if (np.isnan(soil.df.water_perm)).any(): # Go over if any value is missing 

        soil.info['water_perm'] = ["Calculated using MalmbergMaryott function (RMSE = 0.0046)" if np.isnan(soil.df.water_perm[x]) & (soil.df.salinity[x] == 0) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) 
                                    or soil.info.water_perm[x] == "Calculated using MalmbergMaryott function (RMSE = 0.0046)"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
        
        soil.df['water_perm'] = [MalmbergMaryott(soil.df.temperature.values[x]) if np.isnan(soil.df.water_perm[x]) & (soil.df.salinity[x] == 0) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) else soil.df.water_perm[x] for x in range(soil.n_states)]
        
        soil.info['water_perm'] = ["Calculated using Olhoeft function" if np.isnan(soil.df.water_perm[x]) & ~np.isnan(soil.df.salinity[x]) & (soil.df.frequency_perm[x] < 100e6)
                                    or soil.info.water_perm[x] == "Calculated using Olhoeft function"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
                
        soil.df['water_perm'] = [Olhoeft(soil.df.temperature.values[x], soil.df.salinity[x]) if np.isnan(soil.df.water_perm[x]) & ~np.isnan(soil.df.salinity[x]) & (soil.df.frequency_perm[x] < 100e6) else soil.df.water_perm[x] for x in range(soil.n_states)]
                
        soil.info['water_perm'] = ["Set as 80 by default" if np.isnan(soil.df.water_perm[x])
                                    or soil.info.water_perm[x] == "Set as 80 by default"
                                    else soil.info.water_perm[x] for x in range(soil.n_states)]
               
        soil.df['water_perm'] = [80 if np.isnan(soil.df.water_perm[x]) else soil.df.water_perm[x] for x in range(soil.n_states)]

    return soil.df.water_perm.values