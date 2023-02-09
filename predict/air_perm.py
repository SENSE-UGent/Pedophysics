import numpy as np

def air_perm(soil):    #### This is currently a non-state variable
    '''
    
    '''
    if (np.isnan(soil.air_perm) == True).any:  # Go over if any value is missing
        soil.df.loc[(np.isnan(soil.df['air_perm']) == True), ['air_perm']] = 1.2

    return soil.df.air_perm