import numpy as np
import pedophysical_models as pm
 
def water_perm(soil):
    '''
        Soil water phase real dielectric permittivity prediction
    '''
    if (np.isnan(soil.water_perm) == True).any: # Go over if any value is missing 

        soil.df['water_perm'] = [pm.Malmberg_Maryott56(soil.df.temperature.values[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (soil.df.frequency_perm[x]  <= 100e6) & (soil.df.frequency_perm[x] >= 1e5) & (np.isnan(soil.salinity[x]) == True)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df['water_perm'] = [pm.olhoeft(soil.df.temperature.values[x], soil.salinity[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (np.isnan(soil.salinity[x]) == False)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df['water_perm'] = [pm.stogryn(soil.df.temperture.values[x], soil.salinity[x], soil.df.frequency_perm.values[x]) if ((np.isnan(soil.df.water_perm[x]) == True) & (np.isnan(soil.salinity[x]) == False) & (soil.df.frequency_perm[x] >= 100e6)) else soil.df.water_perm[x] for x in range(soil.n_states)]
        soil.df.loc[(np.isnan(soil.df['water_perm']) == True), ['water_perm']] = 80

    return soil.df.water_perm