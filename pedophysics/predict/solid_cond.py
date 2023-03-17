import numpy as np

def SolidCond(soil):
    '''
    
    '''
    if (np.isnan(soil.solid_cond)).any():  # Go over if any value is missing 
        soil.df['solid_cond'] = [0 if np.isnan(soil.df.solid_cond[x]) else soil.df.solid_cond[x] for x in range(soil.n_states)]

    return soil.df.solid_cond.values