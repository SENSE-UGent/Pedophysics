import numpy as np

def BulkPermInf(soil):
    '''
    
    '''
    if (np.isnan(soil.df.bulk_perm_inf)).any():  # Go over if any value is missing

        soil.info['bulk_perm_inf'] = ["Set as 5 by default" if np.isnan(soil.df.bulk_perm_inf[x]) or soil.info.bulk_perm_inf[x] == "Set as 5 by default"
                                     else soil.info.bulk_perm.inf[x] for x in range(soil.n_states)]

        soil.df['bulk_perm_inf'] = [5 if np.isnan(soil.df.bulk_perm_inf[x]) else soil.df.bulk_perm_inf[x] for x in range(soil.n_states)]

    return soil.df.bulk_perm_inf.values