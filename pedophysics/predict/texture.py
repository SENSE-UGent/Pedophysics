import numpy as np
import warnings

def Texture(soil):
    '''
    
    '''
    if (np.isnan(soil.df.sand)).any() or (np.isnan(soil.df.silt)).any() or (np.isnan(soil.df.clay)).any() :  # Go over if any value is missing 
        
    # Warn texture fractions that does not sum 100
        for x in range(soil.n_states):
                    if ~np.isnan(soil.df.sand[x]) & ~np.isnan(soil.df.silt[x]) & ~np.isnan(soil.df.clay[x]):
                        total_percent = soil.df.sand[x] + soil.df.silt[x] + soil.df.clay[x]
                        if total_percent != 100:
                            warnings.warn(f"Total percentage of texture fractions in state: {x} is equal to {total_percent}")

        # Complete a third fraction if just two are given
        soil.info['sand'] = ["Fraction calculated using: 100 - clay - silt" if (np.isnan(soil.df.sand[x]) & ~np.isnan(soil.df.silt[x]) & ~np.isnan(soil.df.clay[x])) 
                                    or (soil.info.sand[x] == "Fraction completed using: 100 - clay - silt") else soil.info.sand[x] for x in range(soil.n_states)]
        soil.df['sand'] = [100 - soil.df.clay[x] - soil.df.silt[x] if ((np.isnan(soil.df.sand[x]) ) & (~np.isnan(soil.df.silt[x])  ) & (~np.isnan(soil.df.clay[x])  )) else soil.df.sand[x] for x in range(soil.n_states)]
        
        soil.info['silt'] = ["Fraction calculated using: 100 - clay - sand" if (np.isnan(soil.df.silt[x]) & ~np.isnan(soil.df.sand[x]) & ~np.isnan(soil.df.clay[x]))
                             or (soil.info.silt[x] == "Fraction completed using: 100 - clay - sand") else soil.info.silt[x] for x in range(soil.n_states)]
        soil.df['silt'] = [100 - soil.df.clay[x] - soil.df.sand[x] if ((np.isnan(soil.df.silt[x]) ) & (~np.isnan(soil.df.sand[x])  ) & (~np.isnan(soil.df.clay[x])  )) else soil.df.silt[x] for x in range(soil.n_states)]
        
        soil.info['clay'] = ["Fraction calculated using: 100 - sand - silt" if (np.isnan(soil.df.clay[x]) & ~np.isnan(soil.df.silt[x]) & ~np.isnan(soil.df.sand[x]))
                             or (soil.info.clay[x] == "Fraction completed using: 100 - sand - silt") else soil.info.clay[x] for x in range(soil.n_states)]
        soil.df['clay'] = [100 - soil.df.sand[x] - soil.df.silt[x] if ((np.isnan(soil.df.clay[x]) ) & (~np.isnan(soil.df.silt[x])  ) & (~np.isnan(soil.df.sand[x])  )) else soil.df.clay[x] for x in range(soil.n_states)]

    # Create a dictionary mapping soil textures to their corresponding fractions
    texture_to_fractions = {
        "Sand": (95, 3, 2),
        "Loamy sand": (82, 12, 6),
        "Sandy loam": (65, 25, 10),
        "Loam": (40, 40, 20),
        "Silt loam": (20, 65, 15),
        "Silt": (8, 86, 6),
        "Sandy clay loam": (60, 25, 15),
        "Clay loam": (30, 35, 35),
        "Silty clay loam": (10, 55, 35),
        "Sandy clay": (50, 10, 40),
        "Clay": (15, 20, 65),
        "Silty clay": (7, 48, 45)
    }

    # Go over each texture and assign the corresponding fractions where needed
    for texture, fractions in texture_to_fractions.items():
        soil.info.loc[(np.isnan(soil.df['sand'])) & (np.isnan(soil.df['silt'])) & (np.isnan(soil.df['clay'])) & (soil.texture == texture), ['sand', 'silt', 'clay']] = ('Fraction calculated using soil.texture', 'Fraction calculated using soil.texture', 'Fraction calculated using soil.texture')
        soil.df.loc[(np.isnan(soil.df['sand'])) & (np.isnan(soil.df['silt'])) & (np.isnan(soil.df['clay'])) & (soil.texture == texture), ['sand', 'silt', 'clay']] = fractions