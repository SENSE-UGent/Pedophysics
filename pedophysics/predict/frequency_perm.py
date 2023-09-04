import numpy as np
from pedophysics import instruments

def FrequencyPerm(soil): 
    '''

    '''
    if (np.isnan(soil.df.frequency_perm)).any(): # Go over if any value is missing 
        instruments.Inst2FreqP(soil)

    return soil.df.frequency_perm.values