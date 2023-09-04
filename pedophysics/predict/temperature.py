import numpy as np

def Temperature(soil):
    '''

    '''
    if (np.isnan(soil.df.temperature)).any():

        soil.info['temperature'] = ["Set as 298.15 K by default" if np.isnan(soil.df.temperature[x]) or soil.info.temperature[x] == "Set as 298.15 K by default"
                                     else soil.info.temperature[x] for x in range(soil.n_states)]
        
        soil.df.loc[(np.isnan(soil.df['temperature'])), ['temperature']] = 298.15

    return soil.df.temperature.values