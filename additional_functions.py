import numpy as np


def inst_to_freq(soil):
    """ 


    """
    if soil.instrument == 'GPR':
        soil.df.loc[(np.isnan(soil.df['frequency_perm']) == True), ['frequency_perm']] = 1e9

    if soil.instrument == 'TDR':
        soil.df.loc[ (np.isnan(soil.df['frequency_perm']) == True), ['frequency_perm']] = 200e6

    if soil.instrument == 'HydraProbe':
        soil.df.loc[(np.isnan(soil.df['frequency_perm']) == True), ['frequency_perm']] = 50e6  
#    frequency_perm = [1e9 if (instrument[x] == 'GPR' and np.isnan(frequency_perm[x]) == True) else frequency_perm[x] for x in range(n_states)]
#    frequency_perm = [200e6 if (instrument[x] == 'TDR' and np.isnan(frequency_perm[x]) == True) else frequency_perm[x] for x in range(n_states)]
#    frequency_perm = [50e6 if (instrument[x] == 'HydraProbe' and np.isnan(frequency_perm[x]) == True) else frequency_perm[x] for x in range(n_states)]
# anote anotee
    return soil.df.frequency_perm
