import numpy as np
import additional_functions as adf

def frequency_perm(soil): 
    '''

    '''
    if (np.isnan(soil.frequency_perm) == True).any: # Go over if any value is missing 
        adf.inst_to_freq(soil)

    return soil.df.frequency_perm