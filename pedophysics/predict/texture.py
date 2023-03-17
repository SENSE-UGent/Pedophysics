import numpy as np

def Texture(soil):
    '''
    
    '''
    if (np.isnan(soil.sand) == True).any() or (np.isnan(soil.silt) == True).any() or (np.isnan(soil.clay) == True).any() :  # Go over if any value is missing 

        # Complete a third fraction if just two are given
        soil.df['sand'] = [100 - soil.df.clay[x] - soil.df.silt[x] if ((np.isnan(soil.df.sand[x]) == True) & (np.isnan(soil.df.silt[x]) == False) & (np.isnan(soil.df.clay[x]) == False)) else soil.df.sand[x] for x in range(soil.n_states)]
        soil.df['silt'] = [100 - soil.df.clay[x] - soil.df.sand[x] if ((np.isnan(soil.df.silt[x]) == True) & (np.isnan(soil.df.sand[x]) == False) & (np.isnan(soil.df.clay[x]) == False)) else soil.df.silt[x] for x in range(soil.n_states)]
        soil.df['clay'] = [100 - soil.df.sand[x] - soil.df.silt[x] if ((np.isnan(soil.df.clay[x]) == True) & (np.isnan(soil.df.silt[x]) == False) & (np.isnan(soil.df.sand[x]) == False)) else soil.df.clay[x] for x in range(soil.n_states)]

        # TODO warn if sabd+silt+clay is diferent than 100

        # Complete a all fractions if texture is given
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == 'Sand'), ['sand', 'silt', 'clay']] = 95, 3, 2
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Loamy sand"), ['sand', 'silt', 'clay']] = 82, 12, 6
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Sandy loam"), ['sand', 'silt', 'clay']] = 65, 25, 10
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Loam"), ['sand', 'silt', 'clay']] = 40, 40, 20
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Silt loam"), ['sand', 'silt', 'clay']] = 20, 65, 15
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Silt"), ['sand', 'silt', 'clay']] = 8, 86, 6
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Sandy clay loam"), ['sand', 'silt', 'clay']] = 60, 25, 15
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Clay loam"), ['sand', 'silt', 'clay']] = 30, 35, 35
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Silty clay loam"), ['sand', 'silt', 'clay']] = 10, 55, 35
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Sandy clay"), ['sand', 'silt', 'clay']] = 50, 10, 40
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Clay"), ['sand', 'silt', 'clay']] = 15, 20, 65
        soil.df.loc[(np.isnan(soil.df['sand']) == True) & (np.isnan(soil.df['silt']) == True) & (np.isnan(soil.df['clay']) == True) & (soil.texture == "Silty clay"), ['sand', 'silt', 'clay']] = 7, 48, 45
