import numpy as np
from pedophysics.pedophysical_models.water_perm import *

def WaterPerm(soil):
    '''
        Soil water phase real dielectric permittivity prediction
    '''
    if (np.isnan(soil.water_perm) == True).any: # Go over if any value is missing 

        soil.df['water_perm'] = [MalmbergMaryott(soil.df.temperature.values[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) & (np.isnan(soil.salinity[x]) == True)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df['water_perm'] = [Olhoeft(soil.df.temperature.values[x], soil.salinity[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (np.isnan(soil.salinity[x]) == False)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df['water_perm'] = [Stogryn(soil.df.temperture.values[x], soil.salinity[x], soil.df.frequency_perm.values[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (np.isnan(soil.salinity[x]) == False) & (soil.df.frequency_perm[x] >= 100e6)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df.loc[(np.isnan(soil.df['water_perm']) == True), ['water_perm']] = 80

    return soil.df.water_perm.values