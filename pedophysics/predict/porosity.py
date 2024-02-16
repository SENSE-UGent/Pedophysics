import numpy as np
from pedophysics.predict.particle_density import ParticleDensity

def Porosity(soil):
    """

    """
    ParticleDensity(soil)
    # Check if any value of porosity is missing
    if any(np.isnan(soil.df.porosity[x]) and not np.isnan(soil.df.bulk_density[x]) for x in range(soil.n_states)): 
        
        soil.info['porosity'] = ["Calculated based on bulk density" if np.isnan(soil.df.porosity[x]) or soil.info.porosity[x] == "Calculated based on bulk density"
                                     else soil.info.porosity[x] for x in range(soil.n_states)]
        
        soil.df['porosity'] = [round(1 - soil.df.bulk_density[x]/soil.df.particle_density[x], soil.roundn) if np.isnan(soil.df.porosity[x]) else soil.df.porosity[x] for x in range(soil.n_states)]

    return soil.df.porosity.values