import numpy as np

def AirPerm(soil): 
    '''
    
    '''
    if (np.isnan(soil.df.air_perm)).any:

        soil.info['air_perm'] = ["Set as 1.2 by default" if np.isnan(soil.df.air_perm[x]) or soil.info.air_perm[x] == "Set as 1.2 by default"
                                     else soil.info.air_perm[x] for x in range(soil.n_states)]

        soil.df['air_perm'] = [1.2 if np.isnan(soil.df.air_perm[x]) else soil.df.air_perm[x] for x in range(soil.n_states)]

    return soil.df.air_perm.values