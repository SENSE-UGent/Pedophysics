import numpy as np

def Inst2FreqP(soil):
    """ 

    """
    if soil.instrument == 'GPR':
        soil.df.loc[np.isnan(soil.df['frequency_perm']), ['frequency_perm']] = 1e9

    elif soil.instrument == 'TDR':
        soil.df.loc[np.isnan(soil.df['frequency_perm']), ['frequency_perm']] = 200e6

    elif soil.instrument == 'HydraProbe':
        soil.df.loc[np.isnan(soil.df['frequency_perm']), ['frequency_perm']] = 50e6  

    return soil.df.frequency_perm


def Inst2FreqC(soil):
    """ 

    """
    if soil.instrument == 'EMI Dualem':
        soil.df.loc[np.isnan(soil.df['frequency_cond']), ['frequency_cond']] = 9e3

    elif soil.instrument == 'EMI EM38-DD':
        soil.df.loc[np.isnan(soil.df['frequency_cond']), ['frequency_cond']] = 16e3

#    elif soil.instrument == 'Miller 400D':
#        soil.df.loc[np.isnan(soil.df['frequency_cond']), ['frequency_cond']] = 82  

#    elif soil.instrument == 'Miller 400D':
#        soil.df.loc[np.isnan(soil.df['frequency_cond']), ['frequency_cond']] = 82  

    return soil.df.frequency_cond