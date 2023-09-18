import numpy as np

from .water_from_ec import WaterFromEC
from .water_from_perm import WaterFromPerm
from ..frequency_perm import FrequencyPerm

def Water(soil):
    """
        Soil volumetric water content prediction
    """
    # Condition to obtain water from bulk_perm 
    FrequencyPerm(soil)

    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) for x in range(soil.n_states)) and (np.isnan(soil.df.frequency_perm)).all():
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Unmodified value. Please provide soil.frequency_perm" if True else soil.info.bulk_perm[x] for x in range(soil.n_states)]
    
    elif any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) for x in range(soil.n_states)) and not (np.isnan(soil.df.frequency_perm).all()):
        WaterFromPerm(soil) 

    # Condition to obtain water from bulk_ec
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)):
        WaterFromEC(soil)        

    return soil.df.water.values