import numpy as np

def BulkPermInf(soil):
    '''
    
    '''
    if (np.isnan(soil.bulk_perm_inf) == True).any:  # Go over if any value is missing
        soil.df.loc[(np.isnan(soil.df['bulk_perm_inf']) == True), ['bulk_perm_inf']] = 5

    return soil.df.bulk_perm_inf.values