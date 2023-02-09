from predict import texture
import numpy as np
import pedotransfer_functions as ptf


def particle_density(soil):
    '''
    
    '''
    texture.texture(soil)
    if (np.isnan(soil.particle_density) == True).any:  # Go over if any value is missing 
        soil.df['particle_density'] = [ptf.schjonnpd(soil.df.clay.values[x], soil.df.orgm.values[x]) if (np.isnan(soil.df.particle_density[x]) == True) else soil.df.particle_density[x] for x in range(soil.n_states)]
        soil.df.loc[(np.isnan(soil.df['particle_density']) == True), ['particle_density']] = 2.65

    return soil.df.particle_density