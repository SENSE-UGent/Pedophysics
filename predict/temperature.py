import numpy as np

def temperature(soil):
    '''

    '''
    if (np.isnan(soil.temperature) == True).any: # Go over if any value is missing 
        soil.df.loc[(np.isnan(soil.df['temperature']) == True), ['temperature']] = 298.15

    return soil.df.temperature