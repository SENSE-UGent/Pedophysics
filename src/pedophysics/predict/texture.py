import warnings
import numpy as np

def Texture(soil):
    """
    Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return

    If any value of the sand, silt, or clay attribute is missing, this function will:
    1. Warn if the sum of texture fractions does not equate to 100%.
    2. Calculate missing texture fraction if only two out of three are given.
    3. Assign default texture fractions based on the `texture` attribute of the soil object.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - clay : array-like
            Soil clay content [g/g]*100
        - sand : array-like
            Soil sand content [g/g]*100
        - silt : array-like
            Soil silt content [g/g]*100
        - texture : str
            Soil texture according to USDA convention: "Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Silty clay loam", "Sandy clay", "Clay", "Silty clay"
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Notes
    -----
    This function modifies the soil object in-place, updating the `df` and `info` dataframes.

    External functions
    --------
    warnings.warn : Issue a warning, or maybe ignore it or raise an exception.

    Example
    -------
    >>> sample = Soil(clay = 20,
                    silt = 40)
    >>> sample.df.sand
    0   NaN
    Name: sand, dtype: float64
    >>> Texture(sample)
    >>> sample.df.sand
    0    40
    Name: sand , dtype: float64
    """

    # Check if any value of sand, silt, or clay is missing. 
    if (np.isnan(soil.df.sand)).any() or (np.isnan(soil.df.silt)).any() or (np.isnan(soil.df.clay)).any() : 

        # Warn texture fractions that does not sum 100
        if any(not np.isnan(soil.df.sand[x]) and not np.isnan(soil.df.silt[x]) and not np.isnan(soil.df.clay[x]) and soil.df.sand[x] + soil.df.silt[x] + soil.df.clay[x] != 100 for x in range(soil.n_states)):
            states_warns = []
            total_percents = []
            for x in range(soil.n_states):
                if ~np.isnan(soil.df.sand[x]) & ~np.isnan(soil.df.silt[x]) & ~np.isnan(soil.df.clay[x]):
                    total_percent = soil.df.sand[x] + soil.df.silt[x] + soil.df.clay[x]
                    if total_percent != 100:
                        total_percents.append(total_percent)
                        states_warns.append(x)

            warnings.warn(f"Total percentage of texture fractions in states: {states_warns} are equal to {total_percents}")

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