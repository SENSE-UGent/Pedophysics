import numpy as np

def solid_perm(soil):
    '''
    
    '''
    if (np.isnan(soil.solid_perm) == True).any:  # Go over if any value is missing
        soil.df.loc[(np.isnan(soil.df['solid_perm']) == True), ['solid_perm']] = 4

    return soil.df.solid_perm