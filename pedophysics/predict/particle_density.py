import numpy as np
from pedophysics.pedotransfer_functions.particle_density import Schjonnen
from .texture import Texture

def ParticleDensity(soil):
    """
    Return or compute missing values of the soil.df.particle_density attribute.

    If any value of the particle_density attribute is missing (NaN), it will first
    be computed using the Schjonnen function based on the clay and organic matter values.
    If it remains missing, it's set to a default value of 2.65. 

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - particle_density : array-like
            Soil particle density [kg/m**3]
        - clay : array-like
            Soil clay content [g/g]*100
        - orgm : array-like
            Soil organic matter [g/g]*100
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        An array of updated soil particle density values

    Notes
    -----
    This function modifies the soil object in-place, updating the `df` dataframe and `info`
    dataframe if necessary.

    External functions
    --------
    Texture : Function to calculate missing clay, silt or sand attributes based on soil.texture
    Schjonnen : Function to calculate particle_density based on soil.df.clay and soil.df.orgm

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.particle_density
    0   NaN
    Name: particle_density, dtype: float64
    >>> ParticleDensity(sample)
    >>> sample.df.particle_density
    0    2.65
    Name: particle_density, dtype: float64
    """

    # Check if any value of particle_density is missing
    if (np.isnan(soil.df.particle_density)).any(): 
        Texture(soil)

        soil.info['particle_density'] = ["Calculated using Schjonnen function (RMSE = 0.011 g/cm3)" if np.isnan(soil.df.particle_density[x]) 
                                         or soil.info.particle_density[x] == "Calculated using Schjonnen function (RMSE = 0.011 g/cm3)"
                                         else soil.info.particle_density[x] for x in range(soil.n_states)]
        
        soil.df['particle_density'] = [Schjonnen(soil.df.clay.values[x], soil.df.orgm.values[x]) if np.isnan(soil.df.particle_density[x])  
                                       else soil.df.particle_density[x] for x in range(soil.n_states)]
        
        soil.info['particle_density'] = ["Set as 2.65 by default" if np.isnan(soil.df.particle_density[x]) or soil.info.particle_density[x] == "Set as 2.65 by default"
                                     else soil.info.particle_density[x] for x in range(soil.n_states)]
        
        soil.df['particle_density'] = [2.65 if np.isnan(soil.df.particle_density[x]) else soil.df.particle_density[x] for x in range(soil.n_states)]

    return soil.df.particle_density.values
