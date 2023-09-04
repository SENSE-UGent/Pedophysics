import numpy as np

def SolidEC(soil):
    '''
    
    '''
    if (np.isnan(soil.df.solid_ec)).any():  # Go over if any value is missing 

        soil.info['solid_ec'] = ["Set as zero by default" if np.isnan(soil.df.solid_ec[x]) or soil.info.solid_ec[x] == "Set as zero by default"
                                 else soil.info.solid_ec[x] for x in range(soil.n_states)]

        soil.df['solid_ec'] = [0 if np.isnan(soil.df.solid_ec[x]) else soil.df.solid_ec[x] for x in range(soil.n_states)]

    return soil.df.solid_ec.values