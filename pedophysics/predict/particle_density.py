import numpy as np
from pedophysics import pedotransfer_functions as ptf
from .texture import Texture

def ParticleDensity(soil):
    '''
    
    '''
    if (np.isnan(soil.particle_density)).any():  # Go over if any value is missing 
        Texture(soil)
        soil.df['particle_density'] = [ptf.Schjonnen(soil.df.clay.values[x], soil.df.orgm.values[x]) if (np.isnan(soil.df.particle_density[x]) == True) else soil.df.particle_density[x] for x in range(soil.n_states)]
        soil.df.loc[(np.isnan(soil.df['particle_density'])), ['particle_density']] = 2.65

    return soil.df.particle_density.values