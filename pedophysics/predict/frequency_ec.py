import numpy as np
from pedophysics import instruments

def FrequencyEC(soil): 
    '''

    '''
    if (np.isnan(soil.df.frequency_ec)).any(): # Go over if any value is missing 
        instruments.Inst2FreqC(soil)

        soil.info['frequency_ec'] = ["Set as 0Hz (direct current) by default" if np.isnan(soil.df.frequency_ec[x]) or soil.info.frequency_ec[x] == "Set as 0Hz (direct current) by default"
                                     else soil.info.frequency_ec[x] for x in range(soil.n_states)]
        
        soil.df['frequency_ec'] = [0 if np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    return soil.df.frequency_ec.values