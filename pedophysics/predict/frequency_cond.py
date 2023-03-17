import numpy as np
from pedophysics import instruments

def FrequencyCond(soil): 
    '''

    '''
    if (np.isnan(soil.frequency_cond)).any(): # Go over if any value is missing 
        instruments.Inst2FreqC(soil)
        soil.df['frequency_cond'] = [90 if np.isnan(soil.df.frequency_cond[x]) else soil.df.frequency_cond[x] for x in range(soil.n_states)]

    return soil.df.frequency_cond.values