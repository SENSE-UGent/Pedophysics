import numpy as np

def SolidPerm(soil):
    '''
    
    '''
    if (np.isnan(soil.df.solid_perm)).any:  # Go over if any value is missing

        soil.info['solid_perm'] = ["Set as 4 by default" if np.isnan(soil.df.solid_perm[x]) or soil.info.solid_perm[x] == "Set as 4 by default"
                                     else soil.info.solid_perm[x] for x in range(soil.n_states)]
        
        soil.df.loc[np.isnan(soil.df['solid_perm']), ['solid_perm']] = 4

    return soil.df.solid_perm.values