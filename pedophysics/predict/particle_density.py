import numpy as np
from pedophysics.pedotransfer_functions.particle_density import Schjonnen
from .texture import Texture

def ParticleDensity(soil):
    '''
    
    '''
    if (np.isnan(soil.df.particle_density)).any():  # Go over if any value is missing 
        Texture(soil)

        soil.info['particle_density'] = ["Calculated using Schjonnen function (RMSE = 0.011 g/cm3)" if np.isnan(soil.df.particle_density[x]) 
                                         or soil.info.particle_density[x] == "Calculated using Schjonnen function (RMSE = 0.011 g/cm3)"
                                         else soil.info.particle_density[x] for x in range(soil.n_states)]
        
        soil.df['particle_density'] = [Schjonnen(soil.df.clay.values[x], soil.df.orgm.values[x]) if np.isnan(soil.df.particle_density[x])  
                                       else soil.df.particle_density[x] for x in range(soil.n_states)]
        
        soil.info['particle_density'] = ["Set as 2.65 by default" if np.isnan(soil.df.particle_density[x]) or soil.info.particle_density[x] == "Set as 2.65 by default"
                                     else soil.info.particle_density[x] for x in range(soil.n_states)]
        
        soil.df['particle_density'] = [2.65 if np.isnan(soil.df.particle_density[x]) else soil.df.particle_density[x] for x in range(soil.n_states)]

    return soil.df.particle_density.values
